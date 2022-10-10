use crate::Constraint;
use crate::ProofOptions;
use crate::TraceInfo;
use ark_ff::batch_inversion;
use ark_ff::FftField;
use ark_ff::One;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use fast_poly::allocator::PageAlignedAllocator;
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

    fn boundary_constraints(&self) -> &[Constraint<Self::Fp>] {
        &[]
    }

    fn transition_constraints(&self) -> &[Constraint<Self::Fp>] {
        &[]
    }

    fn terminal_constraints(&self) -> &[Constraint<Self::Fp>] {
        &[]
    }

    fn transition_constraint_divisor(&self) -> Vec<Self::Fp, PageAlignedAllocator> {
        let trace_domain = self.trace_domain();
        let lde_domain = self.lde_domain();
        let n = lde_domain.size();
        // TODO: make this easier to understand
        // Transition constraints apply to all rows of execution trace except the last
        // row. We need to change the inverse zerofier from being the
        // evaluations of the polynomial `1/((x - o^0)...(x - o^(n-1)))` to
        // `1/((x - o^0)...(x - o^(n-2)))`. This is achieved by performing the
        // dot product of the zerofier evaluations and the evaluations of the polynomial
        // (x - o^(n-1)). Note that o^(n-1) is the inverse of `o`.
        let last_trace_x = trace_domain.group_gen_inv;
        let mut divisor = Vec::with_capacity_in(n, PageAlignedAllocator);
        divisor.resize(n, Self::Fp::zero());
        // TODO: make parallel
        for (i, lde_x) in lde_domain.elements().enumerate() {
            divisor[i] = trace_domain.evaluate_vanishing_polynomial(lde_x);
        }
        batch_inversion(&mut divisor);
        // TODO: make parallel
        for (i, lde_x) in lde_domain.elements().enumerate() {
            divisor[i] *= lde_x - last_trace_x;
        }
        divisor
    }

    fn boundary_constraint_divisor(&self) -> Vec<Self::Fp, PageAlignedAllocator> {
        let lde_domain = self.lde_domain();
        let n = lde_domain.size();
        // TODO: make this easier to understand
        // Evaluations of the polynomial `1/(x - o^0)` over the FRI domain
        // Context: boundary constraints have to be 0 in the first row.
        let first_trace_x = Self::Fp::one();
        let mut divisor = Vec::with_capacity_in(n, PageAlignedAllocator);
        divisor.resize(n, Self::Fp::zero());
        // TODO: make parallel
        for (i, lde_x) in lde_domain.elements().enumerate() {
            divisor[i] = lde_x - first_trace_x;
        }
        batch_inversion(&mut divisor);
        divisor
    }

    fn terminal_constraint_divisor(&self) -> Vec<Self::Fp, PageAlignedAllocator> {
        let trace_domain = self.trace_domain();
        let lde_domain = self.lde_domain();
        let n = lde_domain.size();
        // TODO: make this easier to understand
        // Evaluations of the polynomial `1/(x - o^(n-1))` over the FRI domain
        // Context: terminal constraints have to be 0 in the last row.
        let last_trace_x = trace_domain.group_gen_inv;
        let mut divisor = Vec::with_capacity_in(n, PageAlignedAllocator);
        divisor.resize(n, Self::Fp::zero());
        // TODO: make parallel
        for (i, lde_x) in lde_domain.elements().enumerate() {
            divisor[i] = lde_x - last_trace_x;
        }
        batch_inversion(&mut divisor);
        divisor
    }

    fn num_challenges(&self) -> usize {
        // TODO: change get_challenge_indices to a constraint iterator and extract the
        // constraint with the highest index
        [
            self.boundary_constraints(),
            self.transition_constraints(),
            self.terminal_constraints(),
        ]
        .into_iter()
        .flatten()
        .flat_map(|constraint| constraint.get_challenge_indices())
        .max()
        .map(|max| max + 1)
        .unwrap_or(0)
    }
}
