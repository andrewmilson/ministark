use crate::channel::ProverChannel;
use crate::composer::ConstraintComposer;
use crate::composer::DeepPolyComposer;
use crate::fri::FriProver;
use crate::matrix::GroupItem;
use crate::matrix::MatrixGroup;
use crate::trace::Queries;
use crate::utils::Timer;
use crate::Air;
use crate::Proof;
use crate::ProofOptions;
use crate::Trace;
use ark_ff::FftField;
use ark_ff::Field;
use gpu_poly::GpuField;
use sha2::Sha256;

/// Errors that can occur during the proving stage
#[derive(Debug)]
pub enum ProvingError {
    Fail,
    // TODO
}

pub trait Prover {
    type Fp: GpuField<FftField = Self::Fp> + FftField;
    type Fq: GpuField<FftField = Self::Fp>;
    type Air: Air<Fp = Self::Fp, Fq = Self::Fq>;
    type Trace: Trace<Fp = Self::Fp, Fq = Self::Fq>;

    fn new(options: ProofOptions) -> Self;

    fn get_pub_inputs(&self, trace: &Self::Trace) -> <Self::Air as Air>::PublicInputs;

    fn options(&self) -> ProofOptions;

    fn generate_proof(&self, trace: Self::Trace) -> Result<Proof<Self::Air>, ProvingError> {
        let _timer = Timer::new("proof generation");

        let options = self.options();
        let trace_info = trace.info();
        let pub_inputs = self.get_pub_inputs(&trace);
        let air = Self::Air::new(trace_info, pub_inputs, options);
        air.validate();
        let mut channel = ProverChannel::<Self::Air, Sha256>::new(&air);

        let trace_xs = air.trace_domain();
        let lde_xs = air.lde_domain();
        let base_trace = trace.base_columns();
        let base_trace_polys = base_trace.interpolate(trace_xs);
        assert_eq!(Self::Trace::NUM_BASE_COLUMNS, base_trace_polys.num_cols());
        let base_trace_lde = base_trace_polys.evaluate(lde_xs);
        let base_trace_lde_tree = base_trace_lde.commit_to_rows();
        channel.commit_base_trace(base_trace_lde_tree.root());
        let challenges = air.get_challenges(&mut channel.public_coin);

        let extension_trace = trace.build_extension_columns(&challenges);
        let num_extension_columns = extension_trace.as_ref().map_or(0, |t| t.num_cols());
        assert_eq!(Self::Trace::NUM_EXTENSION_COLUMNS, num_extension_columns);
        let extension_trace_polys = extension_trace.as_ref().map(|t| t.interpolate(trace_xs));
        let extension_trace_lde = extension_trace_polys.as_ref().map(|p| p.evaluate(lde_xs));
        let extension_trace_tree = extension_trace_lde.as_ref().map(|lde| lde.commit_to_rows());
        if let Some(t) = extension_trace_tree.as_ref() {
            channel.commit_extension_trace(t.root())
        }

        #[cfg(debug_assertions)]
        air.validate_constraints(&challenges, base_trace, extension_trace.as_ref());
        #[cfg(not(debug_assertions))]
        drop((base_trace, extension_trace));

        let composition_coeffs = air.get_constraint_composition_coeffs(&mut channel.public_coin);
        let constraint_coposer = ConstraintComposer::new(&air, composition_coeffs);
        // TODO: move commitment here
        let (composition_trace_lde, composition_trace_polys, composition_trace_lde_tree) =
            constraint_coposer.build_commitment(
                &challenges,
                &base_trace_lde,
                extension_trace_lde.as_ref(),
            );
        channel.commit_composition_trace(composition_trace_lde_tree.root());

        let g = trace_xs.group_gen;
        let z = channel.get_ood_point();
        let mut execution_trace_polys = MatrixGroup::new(vec![GroupItem::Fp(&base_trace_polys)]);
        if let Some(extension_trace_polys) = extension_trace_polys.as_ref() {
            execution_trace_polys.append(GroupItem::Fp(extension_trace_polys))
        }
        let ood_execution_trace_evals = execution_trace_polys.evaluate_at(z);
        let ood_execution_trace_evals_next = execution_trace_polys.evaluate_at(z * g);
        channel.send_ood_trace_states(&ood_execution_trace_evals, &ood_execution_trace_evals_next);
        let z_n = z.pow([composition_trace_polys.num_cols() as u64]);
        let ood_composition_trace_evals = composition_trace_polys.evaluate_at(z_n);
        channel.send_ood_constraint_evaluations(&ood_composition_trace_evals);

        let deep_coeffs = air.get_deep_composition_coeffs(&mut channel.public_coin);
        let mut deep_poly_composer = DeepPolyComposer::new(&air, deep_coeffs, z);
        deep_poly_composer.add_execution_trace_polys(
            base_trace_polys,
            extension_trace_polys,
            ood_execution_trace_evals,
            ood_execution_trace_evals_next,
        );
        deep_poly_composer
            .add_composition_trace_polys(composition_trace_polys, ood_composition_trace_evals);
        let deep_composition_poly = deep_poly_composer.into_deep_poly();
        let deep_composition_lde = deep_composition_poly.into_evaluations(lde_xs);

        let mut fri_prover = FriProver::<Self::Fp, Sha256>::new(air.options().into_fri_options());
        fri_prover.build_layers(&mut channel, deep_composition_lde.try_into().unwrap());

        channel.grind_fri_commitments();

        let query_positions = channel.get_fri_query_positions();
        let fri_proof = fri_prover.into_proof(&query_positions);

        let queries = Queries::new(
            &base_trace_lde,
            extension_trace_lde.as_ref(),
            &composition_trace_lde,
            base_trace_lde_tree,
            extension_trace_tree,
            composition_trace_lde_tree,
            &query_positions,
        );

        Ok(channel.build_proof(queries, fri_proof))
    }
}
