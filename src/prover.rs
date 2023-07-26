use crate::air::AirConfig;
use crate::challenges::Challenges;
use crate::channel::ProverChannel;
use crate::composer::DeepPolyComposer;
use crate::fri::FriProver;
use crate::merkle::MatrixMerkleTree;
use crate::merkle::MerkleTree;
use crate::random::draw_multiple;
use crate::stark::Stark;
use crate::trace::Queries;
use crate::utils::GpuAllocator;
use crate::utils::GpuVec;
use crate::Air;
use crate::Matrix;
use crate::Proof;
use crate::ProofOptions;
use crate::Trace;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ministark_gpu::utils::bit_reverse;
use std::time::Instant;

#[allow(clippy::too_many_lines)]
pub fn default_prove<S: Stark>(
    this: &S,
    options: ProofOptions,
    witness: S::Witness,
) -> Result<Proof<S::Fp, S::Fq, S::Digest, S::MerkleTree>, ProvingError> {
    let now = Instant::now();
    let trace = this.generate_trace(witness);
    println!(
        "Generated execution trace (cols={}, rows={}) in {:.0?}",
        trace.base_columns().num_cols(),
        trace.base_columns().num_rows(),
        now.elapsed(),
    );

    let now = Instant::now();
    let air = Air::new(trace.len(), this.get_public_inputs(), options);
    let public_coin = this.gen_public_coin(&air);
    let mut channel = ProverChannel::<S>::new(&air, public_coin);
    println!("Init air: {:?}", now.elapsed());

    let lde_blowup_factor = air.lde_blowup_factor();
    let ce_blowup_factor = air.ce_blowup_factor();

    // reverse engineering TODOs:
    // 1. don't create leaf hash if row is single item
    // 2. commit to out of order rows
    let now = Instant::now();
    let trace_xs = air.trace_domain();
    let lde_xs = air.lde_domain();
    let base_trace = trace.base_columns();
    assert_eq!(S::AirConfig::NUM_BASE_COLUMNS, base_trace.num_cols());
    let base_trace_polys = base_trace.interpolate(trace_xs);
    let mut base_trace_lde = base_trace_polys.bit_reversed_evaluate(lde_xs);
    let base_trace_ce_lde =
        reduduce_blowup_factor(&base_trace_lde, lde_blowup_factor, ce_blowup_factor);
    println!("made it herre");
    let base_trace_tree = S::MerkleTree::from_matrix(&base_trace_lde);
    channel.commit_base_trace(base_trace_tree.root());
    let num_challenges = air.num_challenges();
    let challenges = Challenges::new(draw_multiple(&mut channel.public_coin, num_challenges));
    let hints = air.gen_hints(&challenges);
    println!("Base trace: {:?}", now.elapsed());

    let now = Instant::now();
    let extension_trace = trace.build_extension_columns(&challenges);
    let num_extension_cols = extension_trace.as_ref().map_or(0, Matrix::num_cols);
    assert_eq!(S::AirConfig::NUM_EXTENSION_COLUMNS, num_extension_cols);
    let extension_trace_polys = extension_trace.as_ref().map(|t| t.interpolate(trace_xs));
    let mut extension_trace_lde = extension_trace_polys
        .as_ref()
        .map(|p| p.bit_reversed_evaluate(lde_xs));
    let extension_trace_ce_lde = extension_trace_lde
        .as_ref()
        .map(|t| reduduce_blowup_factor(t, lde_blowup_factor, ce_blowup_factor));
    let extension_trace_tree = extension_trace_lde.as_ref().map(S::MerkleTree::from_matrix);
    if let Some(t) = extension_trace_tree.as_ref() {
        channel.commit_extension_trace(t.root());
    }
    println!("Extension trace: {:?}", now.elapsed());

    #[cfg(debug_assertions)]
    this.validate_constraints(&challenges, &hints, base_trace, extension_trace.as_ref());
    drop((trace, extension_trace));

    let now = Instant::now();
    let num_composition_coeffs = air.num_composition_constraint_coeffs();
    let composition_coeffs = draw_multiple(&mut channel.public_coin, num_composition_coeffs);
    let ce_lde_xs = air.ce_domain();
    let x_lde = ce_lde_xs.elements().collect::<Vec<_>>();
    println!("X lde: {:?}", now.elapsed());
    let now = Instant::now();
    // TODO: improve
    let composition_evals = S::AirConfig::eval_constraint(
        air.composition_constraint(),
        &challenges,
        &hints,
        &composition_coeffs,
        air.ce_blowup_factor(),
        x_lde.to_vec_in(GpuAllocator),
        &base_trace_ce_lde,
        extension_trace_ce_lde.as_ref(),
    );
    println!("Constraint eval: {:?}", now.elapsed());
    let now = Instant::now();
    let composition_poly = composition_evals.into_polynomials(air.ce_domain());
    let composition_trace_cols = air.ce_blowup_factor();
    let composition_trace_polys = Matrix::from_rows(
        GpuVec::try_from(composition_poly)
            .unwrap()
            .chunks(composition_trace_cols)
            .map(<[S::Fq]>::to_vec)
            .collect(),
    );
    let composition_trace_lde = composition_trace_polys.bit_reversed_evaluate(air.lde_domain());
    let composition_trace_tree = S::MerkleTree::from_matrix(&composition_trace_lde);

    channel.commit_composition_trace(composition_trace_tree.root());
    println!("Constraint composition polys: {:?}", now.elapsed());

    let now = Instant::now();
    let z = channel.get_ood_point();
    let mut deep_poly_composer = DeepPolyComposer::new(
        &air,
        z,
        &base_trace_polys,
        extension_trace_polys.as_ref(),
        composition_trace_polys,
    );
    let (execution_trace_oods, composition_trace_oods) = deep_poly_composer.get_ood_evals();
    channel.send_execution_trace_ood_evals(execution_trace_oods);
    channel.send_composition_trace_ood_evals(composition_trace_oods);

    let deep_coeffs = this.gen_deep_coeffs(&mut channel.public_coin, &air);
    let deep_composition_poly = deep_poly_composer.into_deep_poly(deep_coeffs);
    // let deep_xs = Radix2EvaluationDomain::new(lde_xs.size());
    let deep_composition_lde = deep_composition_poly.into_bit_reversed_evaluations(lde_xs);
    println!("Deep composition: {:?}", now.elapsed());

    let now = Instant::now();
    let fri_options = options.into_fri_options();
    let mut fri_prover = FriProver::<S::Fq, S::Digest, S::MerkleTree>::new(fri_options);
    fri_prover.build_layers(&mut channel, deep_composition_lde.try_into().unwrap());

    channel.grind_fri_commitments();

    let query_positions = Vec::from_iter(channel.get_fri_query_positions());
    let fri_proof = fri_prover.into_proof(&query_positions);
    println!("FRI: {:?}", now.elapsed());

    let queries = Queries::new(
        &base_trace_lde,
        extension_trace_lde.as_ref(),
        &composition_trace_lde,
        &base_trace_tree,
        extension_trace_tree.as_ref(),
        &composition_trace_tree,
        &query_positions,
    );
    Ok(channel.build_proof(queries, fri_proof))
}

/// Errors that can occur during the proving stage
#[derive(Debug)]
pub enum ProvingError {
    Fail,
    // TODO
}

// TODO: this is a bit hacky fix
fn reduduce_blowup_factor<F: Field>(
    bit_rev_matrix: &Matrix<F>,
    from_blowup: usize,
    to_blowup: usize,
) -> Matrix<F> {
    let degree_bound = bit_rev_matrix.num_rows() / from_blowup;
    let n = degree_bound * to_blowup;
    let mut columns = Vec::new();
    for column in &bit_rev_matrix.0 {
        let mut column = column[0..n].to_vec_in(GpuAllocator);
        bit_reverse(&mut column);
        columns.push(column);
    }
    Matrix::new(columns)
}
