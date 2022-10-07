use crate::Constraint;
use crate::ProofOptions;
use crate::TraceInfo;
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use fast_poly::GpuField;

pub trait Air {
    type Fp: GpuField;
    type PublicInputs: CanonicalSerialize;

    // TODO: could make this borrow info and options if so inclined
    fn new(info: TraceInfo, inputs: Self::PublicInputs, options: ProofOptions) -> Self;

    fn pub_inputs(&self) -> &Self::PublicInputs;

    fn trace_info(&self) -> &TraceInfo;

    fn options(&self) -> &ProofOptions;

    fn domain_offset(&self) -> Self::Fp {
        Self::Fp::GENERATOR
    }

    fn trace_len(&self) -> usize {
        self.trace_info().trace_len
    }

    fn ce_blowup_factor(&self) -> usize {
        // the blowup factor is the maximum of one or the degree of the highest
        // degree transition constraint.
        let highest_degree = self
            .transition_constraints()
            .into_iter()
            .map(|constraint| constraint.degree())
            .max()
            .unwrap_or(0);
        std::cmp::max(highest_degree, 1)
    }

    fn lde_blowup_factor(&self) -> usize {
        self.options().blowup_factor as usize
    }

    fn trace_domain(&self) -> Radix2EvaluationDomain<Self::Fp> {
        let trace_len = self.trace_len();
        Radix2EvaluationDomain::new(trace_len).unwrap()
    }

    /// Constraint evaluation domain
    fn ce_domain(&self) -> Radix2EvaluationDomain<Self::Fp> {
        let offset = self.domain_offset();
        let trace_len = self.trace_len();
        let ce_blowup_factor = self.ce_blowup_factor();
        Radix2EvaluationDomain::new_coset(trace_len * ce_blowup_factor, offset).unwrap()
    }

    /// Low degree extension domain
    fn lde_domain(&self) -> Radix2EvaluationDomain<Self::Fp> {
        let offset = self.domain_offset();
        let trace_len = self.trace_len();
        let lde_blowup_factor = self.lde_blowup_factor();
        Radix2EvaluationDomain::new_coset(trace_len * lde_blowup_factor, offset).unwrap()
    }

    fn boundary_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        Vec::new()
    }

    fn transition_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        Vec::new()
    }

    fn terminal_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        Vec::new()
    }

    fn num_challenges(&self) -> usize {
        let all_constraints = vec![
            self.boundary_constraints(),
            self.transition_constraints(),
            self.terminal_constraints(),
        ]
        .concat();
        // TODO: change get_challenge_indices to a constraint iterator and extract the
        // constraint with the highest index
        all_constraints
            .iter()
            .flat_map(|constraint| constraint.get_challenge_indices())
            .max()
            .unwrap_or(0)
    }
}
