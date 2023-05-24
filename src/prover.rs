use crate::air::AirConfig;
use crate::channel::ProverChannel;
use crate::composer::DeepPolyComposer;
use crate::fri::FriProver;
// use crate::merkle::ElementCRHScheme;
use crate::trace::Queries;
use crate::utils::GpuAllocator;
use crate::utils::GpuVec;
use crate::Air;
use crate::Matrix;
use crate::Proof;
use crate::ProofOptions;
use crate::StarkExtensionOf;
use crate::Trace;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ministark_gpu::GpuFftField;
use sha2::Sha256;

/// Errors that can occur during the proving stage
#[derive(Debug)]
pub enum ProvingError {
    Fail,
    // TODO
}

pub trait Prover {
    type Fp: GpuFftField + FftField;
    type Fq: StarkExtensionOf<Self::Fp>;
    type AirConfig: AirConfig<Fp = Self::Fp, Fq = Self::Fq>;
    type Trace: Trace<Fp = Self::Fp, Fq = Self::Fq>;
    // type HashScheme: ElementCRHScheme<Self::Fp>;

    fn new(options: ProofOptions) -> Self;

    fn get_pub_inputs(&self, trace: &Self::Trace) -> <Self::AirConfig as AirConfig>::PublicInputs;

    fn options(&self) -> ProofOptions;

    async fn generate_proof(
        &self,
        trace: Self::Trace,
    ) -> Result<Proof<Self::AirConfig>, ProvingError> {
        let options = self.options();
        let trace_info = trace.info();
        let pub_inputs = self.get_pub_inputs(&trace);
        let air = Air::new(trace_info.trace_len, pub_inputs, options);
        let mut channel = ProverChannel::<Self::AirConfig, Sha256>::new(&air);

        let trace_xs = air.trace_domain();
        let lde_xs = air.lde_domain();
        let base_trace = trace.base_columns();
        let base_trace_polys = base_trace.interpolate(trace_xs);
        assert_eq!(Self::Trace::NUM_BASE_COLUMNS, base_trace_polys.num_cols());
        let base_trace_lde = base_trace_polys.evaluate(lde_xs);
        let base_trace_lde_tree = base_trace_lde.commit_to_rows::<Sha256>();
        channel.commit_base_trace(base_trace_lde_tree.root());
        let challenges = air.gen_challenges(&mut channel.public_coin);
        let hints = air.gen_hints(&challenges);

        let extension_trace = trace.build_extension_columns(&challenges);
        let num_extension_columns = extension_trace.as_ref().map_or(0, Matrix::num_cols);
        assert_eq!(Self::Trace::NUM_EXTENSION_COLUMNS, num_extension_columns);
        let extension_trace_polys = extension_trace.as_ref().map(|t| t.interpolate(trace_xs));
        let extension_trace_lde = extension_trace_polys.as_ref().map(|p| p.evaluate(lde_xs));
        let extension_trace_tree = extension_trace_lde.as_ref().map(Matrix::commit_to_rows);
        if let Some(t) = extension_trace_tree.as_ref() {
            channel.commit_extension_trace(t.root());
        }

        // #[cfg(all(feature = "std", debug_assertions))]
        // air.validate_constraints(&challenges, &hints, base_trace,
        // extension_trace.as_ref());
        drop((base_trace, extension_trace));

        let composition_constraint_coeffs =
            air.gen_composition_constraint_coeffs(&mut channel.public_coin);
        let x_lde = lde_xs.elements().collect::<Vec<_>>();
        let composition_evals = Self::AirConfig::eval_constraint(
            air.composition_constraint(),
            &challenges,
            &hints,
            &composition_constraint_coeffs,
            air.lde_blowup_factor(),
            x_lde.to_vec_in(GpuAllocator),
            &base_trace_lde,
            extension_trace_lde.as_ref(),
        );
        let composition_poly = composition_evals.into_polynomials(air.lde_domain());
        let composition_trace_cols = air.ce_blowup_factor();
        let composition_trace_polys = Matrix::from_rows(
            GpuVec::try_from(composition_poly)
                .unwrap()
                .chunks(composition_trace_cols)
                .map(<[Self::Fq]>::to_vec)
                .collect(),
        );
        let composition_trace_lde = composition_trace_polys.evaluate(air.lde_domain());
        let composition_trace_lde_tree = composition_trace_lde.commit_to_rows();
        channel.commit_composition_trace(composition_trace_lde_tree.root());

        let mut deep_poly_composer = DeepPolyComposer::new(
            &air,
            channel.get_ood_point(),
            &base_trace_polys,
            extension_trace_polys.as_ref(),
            composition_trace_polys,
        );
        let (execution_trace_oods, composition_trace_oods) = deep_poly_composer.get_ood_evals();
        channel.send_execution_trace_ood_evals(execution_trace_oods);
        channel.send_composition_trace_ood_evals(composition_trace_oods);
        let deep_coeffs = air.gen_deep_composition_coeffs(&mut channel.public_coin);
        let deep_composition_poly = deep_poly_composer.into_deep_poly(deep_coeffs);
        let deep_composition_lde = deep_composition_poly.into_evaluations(lde_xs);

        let mut fri_prover = FriProver::<Self::Fq, Sha256>::new(air.options().into_fri_options());
        #[cfg(feature = "std")]
        let now = std::time::Instant::now();
        fri_prover.build_layers(&mut channel, deep_composition_lde.try_into().unwrap());
        #[cfg(feature = "std")]
        println!("yo {:?}", now.elapsed());

        channel.grind_fri_commitments();

        let query_positions = channel.get_fri_query_positions();
        let fri_proof = fri_prover.into_proof(&query_positions);

        let queries = Queries::new(
            &base_trace_lde,
            extension_trace_lde.as_ref(),
            &composition_trace_lde,
            &base_trace_lde_tree,
            extension_trace_tree.as_ref(),
            &composition_trace_lde_tree,
            &query_positions,
        );
        Ok(channel.build_proof(queries, fri_proof))
    }
}
