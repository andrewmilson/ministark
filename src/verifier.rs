use crate::air::AirConfig;
use crate::challenges::Challenges;
use crate::channel::VerifierChannelArtifacts;
use crate::composer::DeepCompositionCoeffs;
use crate::constraints::AlgebraicItem;
use crate::constraints::CompositionItem;
use crate::fri;
use crate::fri::FriVerifier;
use crate::hints::Hints;
use crate::merkle::MatrixMerkleTree;
use crate::random::draw_multiple;
use crate::random::PublicCoin;
use crate::stark::Stark;
use crate::utils::horner_evaluate;
use crate::utils::FieldVariant;
use crate::Air;
use crate::Proof;
use alloc::collections::BTreeMap;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use ministark_gpu::utils::bit_reverse_index;
use snafu::Snafu;

#[allow(clippy::too_many_lines)]
pub fn default_verify<S: Stark>(
    this: &S,
    proof: Proof<S>,
    required_security_bits: u32,
) -> Result<VerifierChannelArtifacts<S::Fq>, VerificationError> {
    use VerificationError::*;

    if proof.security_level_bits() < required_security_bits {
        return Err(InvalidProofSecurity);
    }

    let Proof {
        options,
        base_trace_commitment,
        extension_trace_commitment,
        composition_trace_commitment,
        execution_trace_ood_evals,
        composition_trace_ood_evals,
        trace_queries,
        trace_len,
        fri_proof,
        pow_nonce,
        ..
    } = proof;

    let air = Air::new(trace_len, this.get_public_inputs(), options);
    let mut public_coin = this.gen_public_coin(&air);

    public_coin.reseed_with_digest(&base_trace_commitment);
    let num_challenges = air.num_challenges();
    let air_challenges = Challenges::new(draw_multiple(&mut public_coin, num_challenges));
    let air_hints = air.gen_hints(&air_challenges);

    let extension_trace_commitment = extension_trace_commitment.map(|commitment| {
        public_coin.reseed_with_digest(&commitment);
        commitment
    });

    let num_composition_coeffs = air.num_composition_constraint_coeffs();
    let composition_coeffs = draw_multiple(&mut public_coin, num_composition_coeffs);
    public_coin.reseed_with_digest(&composition_trace_commitment);

    let z = public_coin.draw();
    let ood_evals = [
        execution_trace_ood_evals.clone(),
        composition_trace_ood_evals.clone(),
    ]
    .concat();
    public_coin.reseed_with_field_elements(&ood_evals);
    // execution trace ood evaluation map
    let trace_ood_eval_map = air
        .trace_arguments()
        .into_iter()
        .zip(execution_trace_ood_evals)
        .collect::<BTreeMap<(usize, isize), S::Fq>>();
    let calculated_ood_constraint_evaluation = ood_constraint_evaluation::<S::AirConfig>(
        &composition_coeffs,
        &air_challenges,
        &air_hints,
        &trace_ood_eval_map,
        &air,
        z,
    );

    let provided_ood_constraint_evaluation = horner_evaluate(&composition_trace_ood_evals, &z);

    if calculated_ood_constraint_evaluation != provided_ood_constraint_evaluation {
        return Err(InconsistentOodConstraintEvaluations);
    }

    let deep_coeffs = this.gen_deep_coeffs(&mut public_coin, &air);
    let fri_verifier = FriVerifier::<S::Fq, S::Digest, S::MerkleTree>::new(
        &mut public_coin,
        options.into_fri_options(),
        fri_proof,
        trace_len - 1,
    )?;

    if options.grinding_factor != 0 {
        if !public_coin.verify_proof_of_work(options.grinding_factor, pow_nonce) {
            return Err(FriProofOfWork);
        }
        public_coin.reseed_with_int(pow_nonce);
    }

    let lde_domain_size = air.trace_len() * air.lde_blowup_factor();
    let query_positions =
        Vec::from_iter(public_coin.draw_queries(options.num_queries.into(), lde_domain_size));

    let base_trace_rows = trace_queries
        .base_trace_values
        .chunks(S::AirConfig::NUM_BASE_COLUMNS)
        .collect::<Vec<_>>();
    let extension_trace_rows = if S::AirConfig::NUM_EXTENSION_COLUMNS == 0 {
        Vec::new()
    } else {
        trace_queries
            .extension_trace_values
            .chunks(S::AirConfig::NUM_EXTENSION_COLUMNS)
            .collect::<Vec<_>>()
    };

    let composition_trace_rows = trace_queries
        .composition_trace_values
        .chunks(air.ce_blowup_factor())
        .collect::<Vec<&[S::Fq]>>();

    // base trace positions
    S::MerkleTree::verify_rows(
        &base_trace_commitment,
        &query_positions,
        &base_trace_rows,
        trace_queries.base_trace_proof,
    )
    .map_err(|_| BaseTraceQueryDoesNotMatchCommitment)?;

    if let Some(extension_trace_commitment) = extension_trace_commitment {
        S::MerkleTree::verify_rows(
            &extension_trace_commitment,
            &query_positions,
            &extension_trace_rows,
            trace_queries.extension_trace_proof.unwrap(),
        )
        .map_err(|_| ExtensionTraceQueryDoesNotMatchCommitment)?;
    }

    // composition trace positions
    S::MerkleTree::verify_rows(
        &composition_trace_commitment,
        &query_positions,
        &composition_trace_rows,
        trace_queries.composition_trace_proof,
    )
    .map_err(|_| CompositionTraceQueryDoesNotMatchCommitment)?;

    let deep_evaluations = deep_composition_evaluations(
        &air,
        &query_positions,
        &deep_coeffs,
        &base_trace_rows,
        &extension_trace_rows,
        &composition_trace_rows,
        &trace_ood_eval_map,
        &composition_trace_ood_evals,
        z,
    );

    let fri_alphas = fri_verifier.layer_alphas.clone();
    fri_verifier.verify(&query_positions, &deep_evaluations)?;

    Ok(VerifierChannelArtifacts {
        air_challenges,
        air_hints,
        fri_alphas,
        query_positions,
    })
}

/// Errors that are returned during verification of a STARK proof
#[derive(Debug, Snafu)]
pub enum VerificationError {
    #[snafu(display("proof params do not satisfy security requirements"))]
    InvalidProofSecurity,
    #[snafu(display("constraint evaluations at the out-of-domain point are inconsistent"))]
    InconsistentOodConstraintEvaluations,
    #[snafu(context(false))]
    #[snafu(display("fri verification failed: {source}"))]
    FriVerification { source: fri::VerificationError },
    #[snafu(display("query does not resolve to the base trace commitment"))]
    BaseTraceQueryDoesNotMatchCommitment,
    #[snafu(display("query does not resolve to the extension trace commitment"))]
    ExtensionTraceQueryDoesNotMatchCommitment,
    #[snafu(display("query does not resolve to the composition trace commitment"))]
    CompositionTraceQueryDoesNotMatchCommitment,
    #[snafu(display("insufficient proof of work on fri commitments"))]
    FriProofOfWork,
}

pub fn ood_constraint_evaluation<A: AirConfig>(
    composition_coefficients: &[A::Fq],
    challenges: &Challenges<A::Fq>,
    hints: &Hints<A::Fq>,
    trace_ood_eval_map: &BTreeMap<(usize, isize), A::Fq>,
    air: &Air<A>,
    x: A::Fq,
) -> A::Fq {
    use AlgebraicItem::*;
    use CompositionItem::*;
    air.composition_constraint()
        .graph_eval(&mut |leaf| match leaf {
            Item(X) => FieldVariant::Fq(x),
            &Item(Constant(v)) => v,
            &Item(Challenge(i)) => FieldVariant::Fq(challenges[i]),
            &Item(Hint(i)) => FieldVariant::Fq(hints[i]),
            &Item(Periodic(col)) => {
                let trace_len = air.trace_len();
                let point = x.pow([(trace_len / col.interval_size()) as u64]);
                let coeffs = col
                    .coeffs()
                    .iter()
                    .map(FieldVariant::as_fq)
                    .collect::<Vec<_>>();
                FieldVariant::Fq(horner_evaluate(&coeffs, &point))
            }
            &Item(Trace(i, j)) => FieldVariant::Fq(trace_ood_eval_map[&(i, j)]),
            &CompositionCoeff(i) => FieldVariant::Fq(composition_coefficients[i]),
        })
        .as_fq()
}

#[allow(clippy::too_many_arguments)]
pub fn deep_composition_evaluations<A: AirConfig>(
    air: &Air<A>,
    query_positions: &[usize],
    composition_coeffs: &DeepCompositionCoeffs<A::Fq>,
    base_trace_rows: &[&[A::Fp]],
    extension_trace_rows: &[&[A::Fq]],
    composition_trace_rows: &[&[A::Fq]],
    execution_trace_ood_evals_map: &BTreeMap<(usize, isize), A::Fq>,
    composition_trace_ood_evals: &[A::Fq],
    z: A::Fq,
) -> Vec<A::Fq> {
    let trace_domain = air.trace_domain();
    let g = trace_domain.group_gen();
    let g_inv = trace_domain.group_gen_inv();
    let z_n = z.pow([air.ce_blowup_factor() as u64]);
    let lde_domain = air.lde_domain();
    let lde_domain_size = lde_domain.size();
    let xs = query_positions
        .iter()
        .map(|pos| lde_domain.element(bit_reverse_index(lde_domain_size, *pos)))
        .collect::<Vec<A::Fp>>();

    let mut evals = vec![A::Fq::zero(); query_positions.len()];

    let num_columns = A::NUM_BASE_COLUMNS + A::NUM_EXTENSION_COLUMNS;
    let base_column_range = 0..A::NUM_BASE_COLUMNS;
    let extension_column_range = A::NUM_BASE_COLUMNS..num_columns;

    for (i, (&x, eval)) in xs.iter().zip(&mut evals).enumerate() {
        // execution trace
        for (j, ((column, offset), ood_eval)) in execution_trace_ood_evals_map.iter().enumerate() {
            let trace_value = if base_column_range.contains(column) {
                A::Fq::from(base_trace_rows[i][*column])
            } else if extension_column_range.contains(column) {
                extension_trace_rows[i][column - A::NUM_BASE_COLUMNS]
            } else {
                panic!("column {column} does not exist");
            };

            let alpha = composition_coeffs.execution_trace[j];
            let shift = if *offset >= 0 { g } else { g_inv }.pow([offset.unsigned_abs() as u64]);
            *eval += alpha * (trace_value - ood_eval) / (A::Fq::from(x) - z * shift);
        }

        // composition trace
        for (j, value) in composition_trace_rows[i].iter().enumerate() {
            let alpha = composition_coeffs.composition_trace[j];
            let ood_eval = composition_trace_ood_evals[j];
            *eval += alpha * (*value - ood_eval) / (A::Fq::from(x) - z_n);
        }
    }

    // adjust degree
    let (alpha, beta) = composition_coeffs.degree;
    for (x, eval) in xs.iter().zip(&mut evals) {
        *eval *= alpha + beta * x;
    }

    evals
}
