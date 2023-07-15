use crate::air::AirConfig;
use crate::challenges::Challenges;
use crate::composer::DeepCompositionCoeffs;
use crate::constraints::AlgebraicItem;
use crate::constraints::CompositionItem;
use crate::fri;
use crate::fri::FriVerifier;
use crate::hints::Hints;
use crate::merkle;
use crate::merkle::MerkleTree;
use crate::random::draw_multiple;
use crate::random::PublicCoin;
use crate::utils::horner_evaluate;
use crate::utils::FieldVariant;
use crate::Air;
use crate::Proof;
use crate::StarkExtensionOf;
use alloc::collections::BTreeMap;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use ministark_gpu::GpuFftField;
use snafu::Snafu;

pub trait Verifiable {
    type Fp: GpuFftField + FftField;
    type Fq: StarkExtensionOf<Self::Fp>;
    type AirConfig: AirConfig<Fp = Self::Fp, Fq = Self::Fq>;
    type Digest: Digest;
    type PublicCoin: PublicCoin<Digest = Self::Digest, Field = Self::Fq>;

    fn get_public_inputs(&self) -> <Self::AirConfig as AirConfig>::PublicInputs;

    fn gen_public_coin(&self, air: &Air<Self::AirConfig>) -> Self::PublicCoin {
        let mut seed = Vec::new();
        air.public_inputs().serialize_compressed(&mut seed).unwrap();
        air.trace_len().serialize_compressed(&mut seed).unwrap();
        air.options().serialize_compressed(&mut seed).unwrap();
        PublicCoin::new(Self::Digest::digest(&seed))
    }

    fn gen_deep_coeffs(
        &self,
        public_coin: &mut Self::PublicCoin,
        air: &Air<Self::AirConfig>,
    ) -> DeepCompositionCoeffs<Self::Fq> {
        let num_execution_trace = air.trace_arguments().len();
        let num_composition_trace = air.ce_blowup_factor();
        DeepCompositionCoeffs {
            execution_trace: draw_multiple(public_coin, num_execution_trace),
            composition_trace: draw_multiple(public_coin, num_composition_trace),
            degree: (public_coin.draw(), public_coin.draw()),
        }
    }

    #[allow(clippy::too_many_lines)]
    fn verify(&self, proof: Proof<Self::Fp, Self::Fq>) -> Result<(), VerificationError> {
        use VerificationError::*;

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

        let air = Air::new(trace_len, self.get_public_inputs(), options);
        let mut public_coin = self.gen_public_coin(&air);

        let base_trace_commitment = Output::<Self::Digest>::from_iter(base_trace_commitment);
        public_coin.reseed_with_hash(&base_trace_commitment);
        let num_challenges = air.num_challenges();
        let challenges = Challenges::new(draw_multiple(&mut public_coin, num_challenges));
        let hints = air.gen_hints(&challenges);

        let extension_trace_commitment = extension_trace_commitment.map(|commitment| {
            let commitment = Output::<Self::Digest>::from_iter(commitment);
            public_coin.reseed_with_hash(&commitment);
            commitment
        });

        let num_composition_coeffs = air.num_composition_constraint_coeffs();
        let composition_coeffs = draw_multiple(&mut public_coin, num_composition_coeffs);
        let composition_trace_commitment =
            Output::<Self::Digest>::from_iter(composition_trace_commitment);
        public_coin.reseed_with_hash(&composition_trace_commitment);

        let z = public_coin.draw();
        for eval in &execution_trace_ood_evals {
            public_coin.reseed_with_field_element(eval);
        }
        // execution trace ood evaluation map
        let trace_ood_eval_map = air
            .trace_arguments()
            .into_iter()
            .zip(execution_trace_ood_evals)
            .collect::<BTreeMap<(usize, isize), Self::Fq>>();
        let calculated_ood_constraint_evaluation = ood_constraint_evaluation::<Self::AirConfig>(
            &composition_coeffs,
            &challenges,
            &hints,
            &trace_ood_eval_map,
            &air,
            z,
        );

        for eval in &composition_trace_ood_evals {
            public_coin.reseed_with_field_element(eval);
        }
        let provided_ood_constraint_evaluation = horner_evaluate(&composition_trace_ood_evals, &z);

        if calculated_ood_constraint_evaluation != provided_ood_constraint_evaluation {
            return Err(InconsistentOodConstraintEvaluations);
        }

        let deep_coeffs = self.gen_deep_coeffs(&mut public_coin, &air);
        let fri_verifier = FriVerifier::<Self::Fq, Self::Digest>::new(
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
            .chunks(Self::AirConfig::NUM_BASE_COLUMNS)
            .collect::<Vec<_>>();
        let extension_trace_rows = if Self::AirConfig::NUM_EXTENSION_COLUMNS == 0 {
            Vec::new()
        } else {
            trace_queries
                .extension_trace_values
                .chunks(Self::AirConfig::NUM_EXTENSION_COLUMNS)
                .collect::<Vec<_>>()
        };

        let composition_trace_rows = trace_queries
            .composition_trace_values
            .chunks(air.ce_blowup_factor())
            .collect::<Vec<&[Self::Fq]>>();

        // base trace positions
        verify_positions::<Self::Digest>(
            &base_trace_commitment,
            &query_positions,
            &base_trace_rows,
            trace_queries.base_trace_proofs,
        )
        .map_err(|_| BaseTraceQueryDoesNotMatchCommitment)?;

        if let Some(extension_trace_commitment) = extension_trace_commitment {
            // extension trace positions
            verify_positions::<Self::Digest>(
                &extension_trace_commitment,
                &query_positions,
                &extension_trace_rows,
                trace_queries.extension_trace_proofs,
            )
            .map_err(|_| ExtensionTraceQueryDoesNotMatchCommitment)?;
        }

        // composition trace positions
        verify_positions::<Self::Digest>(
            &composition_trace_commitment,
            &query_positions,
            &composition_trace_rows,
            trace_queries.composition_trace_proofs,
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

        Ok(fri_verifier.verify(&query_positions, &deep_evaluations)?)
    }
}

/// Errors that are returned during verification of a STARK proof
#[derive(Debug, Snafu)]
pub enum VerificationError {
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
            &Item(Trace(i, j)) => FieldVariant::Fq(trace_ood_eval_map[&(i, j)]),
            &CompositionCoeff(i) => FieldVariant::Fq(composition_coefficients[i]),
        })
        .as_fq()
}

pub fn verify_positions<D: Digest>(
    commitment: &Output<D>,
    positions: &[usize],
    rows: &[&[impl CanonicalSerialize]],
    proofs: Vec<merkle::Proof>,
) -> Result<(), merkle::Error> {
    for ((position, proof), row) in positions.iter().zip(proofs).zip(rows) {
        let proof = proof.parse::<D>();
        let expected_leaf = &proof[0];
        let mut row_bytes = Vec::with_capacity(row.compressed_size());
        row.serialize_compressed(&mut row_bytes).unwrap();
        let actual_leaf = D::digest(&row_bytes);

        if *expected_leaf != actual_leaf {
            println!("WTF");
            return Err(merkle::Error::InvalidProof);
        }

        MerkleTree::<D>::verify(commitment, &proof, *position)?;
    }

    Ok(())
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
    let lde_domain = air.lde_domain();
    let xs = query_positions
        .iter()
        .map(|pos| lde_domain.element(*pos))
        .collect::<Vec<A::Fp>>();

    let mut evals = vec![A::Fq::zero(); query_positions.len()];

    let num_columns = A::NUM_BASE_COLUMNS + A::NUM_EXTENSION_COLUMNS;
    let base_column_range = 0..A::NUM_BASE_COLUMNS;
    let extension_column_range = A::NUM_BASE_COLUMNS..num_columns;

    // add execution trace
    for (i, (&x, eval)) in xs.iter().zip(&mut evals).enumerate() {
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
    }

    // add composition trace
    let z_n = z.pow([air.ce_blowup_factor() as u64]);
    for ((&x, row), eval) in xs.iter().zip(composition_trace_rows).zip(&mut evals) {
        for (i, value) in row.iter().enumerate() {
            let alpha = composition_coeffs.composition_trace[i];
            let ood_eval = composition_trace_ood_evals[i];
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
