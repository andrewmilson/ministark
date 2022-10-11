use crate::challenges::Challenges;
use crate::constraint::Element;
use crate::Constraint;
use crate::Matrix;
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

    // Constraint evaluation blowup factor
    fn ce_blowup_factor(&self) -> usize {
        // the blowup factor is the maximum of one or the degree of the highest
        // degree transition constraint.
        let highest_degree = self
            .transition_constraints()
            .iter()
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
        self.all_constraint_elements()
            .iter()
            .filter_map(|element| match element {
                Element::Challenge(index) => Some(index + 1),
                _ => None,
            })
            .max()
            .unwrap_or(0)
    }

    fn all_constraint_elements(&self) -> Vec<Element> {
        // TODO: change get_challenge_indices to a constraint iterator and extract the
        // constraint with the highest index
        let mut indicies = [
            self.boundary_constraints(),
            self.transition_constraints(),
            self.terminal_constraints(),
        ]
        .into_iter()
        .flatten()
        .flat_map(|constraint| constraint.get_elements())
        .collect::<Vec<Element>>();
        indicies.sort();
        indicies.dedup();
        indicies
    }

    fn validate(&self, challenges: &Challenges<Self::Fp>, full_trace: &Matrix<Self::Fp>) {
        let mut col_indicies = vec![false; full_trace.num_cols()];
        let mut challenge_indicies = vec![false; challenges.len()];
        for element in self.all_constraint_elements() {
            match element {
                Element::Curr(i) | Element::Next(i) => col_indicies[i] = true,
                Element::Challenge(i) => challenge_indicies[i] = true,
            }
        }
        for (index, exists) in col_indicies.into_iter().enumerate() {
            if !exists {
                // TODO: make assertion
                println!("WARN: no constraints for column {index}");
            }
        }
        for (index, exists) in challenge_indicies.into_iter().enumerate() {
            if !exists {
                // TODO: make assertion
                println!("WARN: challenge at index {index} never used");
            }
        }

        let trace_rows = full_trace.rows();
        let first_row = trace_rows.first().unwrap();
        let last_row = trace_rows.last().unwrap();

        // check boundary constraints
        for (i, constraint) in self.boundary_constraints().iter().enumerate() {
            let eval = constraint.evaluate(challenges, first_row, &[]);
            assert!(eval.is_zero(), "boundary {i} mismatch");
        }

        // check terminal constraints
        for (i, constraint) in self.terminal_constraints().iter().enumerate() {
            let eval = constraint.evaluate(challenges, last_row, &[]);
            assert!(eval.is_zero(), "terminal {i} mismatch");
        }

        // check transition constraints
        for (i, window) in trace_rows.windows(2).enumerate() {
            let curr_row = &window[0];
            let next_row = &window[1];

            if i == 0 {
                println!("Curr:");
                print_row(curr_row);
                println!("Next:");
                print_row(next_row);
            }

            for (j, constraint) in self.terminal_constraints().iter().enumerate() {
                let eval = constraint.evaluate(challenges, curr_row, next_row);
                assert!(eval.is_zero(), "transition {j} mismatch at row {i}");
            }
        }
    }
}

fn print_row<F: GpuField>(row: &[F]) {
    for val in row {
        print!("{val}, ");
    }
    println!()
}
