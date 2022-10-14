use crate::channel::ProverChannel;
use crate::composer::ConstraintComposer;
use crate::composer::DeepPolyComposer;
use crate::merkle::MerkleTree;
use crate::utils::Timer;
use crate::Air;
use crate::Matrix;
use crate::Trace;
use crate::TraceInfo;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::UniformRand;
use ark_ff::Zero;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use fast_poly::GpuField;
use sha2::Sha256;

// TODO: include ability to specify:
// - base field
// - extension field
// - hashing function
// - determine if grinding factor is appropriate
// - fri folding factor
// - fri max remainder size
#[derive(Debug, Clone, Copy, CanonicalSerialize, CanonicalDeserialize)]
pub struct ProofOptions {
    pub num_queries: u8,
    // would be nice to make this clear as LDE blowup factor vs constraint blowup factor
    pub blowup_factor: u8,
}

impl ProofOptions {
    pub fn new(num_queries: u8, blowup_factor: u8) -> Self {
        ProofOptions {
            num_queries,
            blowup_factor,
        }
    }
}

/// A proof generated by a mini-stark prover
#[derive(Debug, Clone)]
pub struct Proof {
    options: ProofOptions,
    trace_info: TraceInfo,
    commitments: Vec<u64>,
}

/// Errors that can occur during the proving stage
#[derive(Debug)]
pub enum ProvingError {
    Fail,
    // /// This error occurs when a transition constraint evaluated over a specific execution
    // trace /// does not evaluate to zero at any of the steps.
    // UnsatisfiedTransitionConstraintError(usize),
    // /// This error occurs when polynomials built from the columns of a constraint evaluation
    // /// table do not all have the same degree.
    // MismatchedConstraintPolynomialDegree(usize, usize),
}

pub trait Prover {
    type Fp: GpuField;
    type Air: Air<Fp = Self::Fp>;
    type Trace: Trace<Fp = Self::Fp>;

    fn new(options: ProofOptions) -> Self;

    fn get_pub_inputs(&self, trace: &Self::Trace) -> <Self::Air as Air>::PublicInputs;

    fn options(&self) -> ProofOptions;

    /// Return value is of the form `(lde, polys, merkle_tree)`
    fn build_trace_commitment(
        &self,
        trace: &Matrix<Self::Fp>,
        trace_domain: Radix2EvaluationDomain<Self::Fp>,
        lde_domain: Radix2EvaluationDomain<Self::Fp>,
    ) -> (Matrix<Self::Fp>, Matrix<Self::Fp>, MerkleTree<Sha256>) {
        let trace_polys = {
            let _timer = Timer::new("trace interpolation");
            trace.interpolate_columns(trace_domain)
        };
        let trace_lde = {
            let _timer = Timer::new("trace low degree extension");
            trace_polys.evaluate(lde_domain)
        };
        let merkle_tree = {
            let _timer = Timer::new("trace commitment");
            trace_lde.commit_to_rows()
        };
        (trace_lde, trace_polys, merkle_tree)
    }

    fn generate_proof(&self, trace: Self::Trace) -> Result<Proof, ProvingError> {
        let _timer = Timer::new("proof generation");

        let options = self.options();
        let trace_info = trace.info();
        let pub_inputs = self.get_pub_inputs(&trace);
        let air = Self::Air::new(trace_info.clone(), pub_inputs, options);
        let mut channel = ProverChannel::<Self::Air, Sha256>::new(&air);

        {
            // TODO: move into validation section
            let ce_blowup_factor = air.ce_blowup_factor();
            let lde_blowup_factor = air.lde_blowup_factor();
            assert!(ce_blowup_factor <= lde_blowup_factor, "constraint evaluation blowup factor {ce_blowup_factor} is larger than the lde blowup factor {lde_blowup_factor}");
        }

        let trace_domain = air.trace_domain();
        let lde_domain = air.lde_domain();
        let (base_trace_lde, base_trace_polys, base_trace_lde_tree) =
            self.build_trace_commitment(trace.base_columns(), trace_domain, lde_domain);

        channel.commit_base_trace(base_trace_lde_tree.root());
        let num_challenges = air.num_challenges();
        let challenges = channel.get_challenges::<Self::Fp>(num_challenges);

        #[cfg(debug_assertions)]
        let mut execution_trace = trace.base_columns().clone();
        let mut execution_trace_lde = base_trace_lde;
        let mut execution_trace_polys = base_trace_polys;
        let mut extension_trace_tree = None;

        if let Some(extension_trace) = trace.build_extension_columns(&challenges) {
            let (extension_trace_lde, extension_trace_polys, extension_trace_lde_tree) =
                self.build_trace_commitment(&extension_trace, trace_domain, lde_domain);
            channel.commit_extension_trace(extension_trace_lde_tree.root());
            #[cfg(debug_assertions)]
            execution_trace.append(extension_trace);
            execution_trace_lde.append(extension_trace_lde);
            execution_trace_polys.append(extension_trace_polys);
            extension_trace_tree = Some(extension_trace_lde_tree);
        }

        #[cfg(debug_assertions)]
        air.validate_constraints(&challenges, &execution_trace);

        let composition_coeffs = channel.get_constraint_composition_coeffs();
        let constraint_coposer = ConstraintComposer::new(&air, composition_coeffs);
        // TODO: move commitment here
        let (composition_trace_lde, composition_trace_polys, composition_trace_lde_tree) =
            constraint_coposer.build_commitment(&challenges, &execution_trace_lde);
        channel.commit_composition_trace(composition_trace_lde_tree.root());

        let g = trace_domain.group_gen;
        let z = channel.get_ood_point();
        let ood_execution_trace_evals = execution_trace_polys.evaluate_at(z);
        let ood_execution_trace_evals_next = execution_trace_polys.evaluate_at(z * g);
        channel.send_ood_trace_states(ood_execution_trace_evals, ood_execution_trace_evals_next);
        let z_n = z.pow([execution_trace_polys.num_cols() as u64]);
        let ood_composition_trace_evals = composition_trace_polys.evaluate_at(z_n);
        channel.send_ood_constraint_evaluations(ood_composition_trace_evals);

        let deep_coeffs = channel.get_deep_composition_coeffs();
        let deep_poly_composer = DeepPolyComposer::new(&air, deep_coeffs, z);
        deep_poly_composer.add_execution_trace_polys(
            execution_trace_polys,
            ood_execution_trace_evals,
            ood_execution_trace_evals_next,
        );
        deep_poly_composer
            .add_composition_trace_polys(composition_trace_polys, ood_composition_trace_evals);

        Ok(Proof {
            options,
            trace_info,
            commitments: Vec::new(),
        })
    }
}
