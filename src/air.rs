use crate::challenges::Challenges;
use crate::composer::DeepCompositionCoeffs;
use crate::constraint::Element;
use crate::hints::Hints;
use crate::random::PublicCoin;
use crate::utils;
use crate::utils::fill_vanishing_polynomial;
use crate::utils::Timer;
use crate::Constraint;
use crate::ProofOptions;
use crate::StarkExtensionOf;
use crate::TraceInfo;
use ark_ff::batch_inversion;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::UniformRand;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use gpu_poly::prelude::*;
use gpu_poly::GpuFftField;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Deref;

pub trait Air {
    type Fp: GpuFftField;
    type Fq: StarkExtensionOf<Self::Fp>;
    // TODO: consider removing clone requirement
    type PublicInputs: CanonicalSerialize + CanonicalDeserialize + Clone;

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

    /// Constraint evaluation blowup factor
    /// Must be a power of two.
    fn ce_blowup_factor(&self) -> usize {
        let max_boundary_constraint_degree = self
            .boundary_constraints()
            .iter()
            .map(|constraint| constraint.degree())
            .max()
            .unwrap_or(0);
        let boundary_ce_blowup_factor = utils::ceil_power_of_two(max_boundary_constraint_degree);

        let max_terminal_constraint_degree = self
            .terminal_constraints()
            .iter()
            .map(|constraint| constraint.degree())
            .max()
            .unwrap_or(0);
        let terminal_ce_blowup_factor = utils::ceil_power_of_two(max_terminal_constraint_degree);

        let max_transition_constraint_degree = self
            .transition_constraints()
            .iter()
            .map(|constraint| constraint.degree())
            .max()
            .unwrap_or(0);
        // TODO: improve explanation of why we negate these constraint degrees by 1
        // Transition constraints must evaluate to zero in all execution trace rows
        // except the last. These rows are divided out from the transition constraint
        // evaluations which has the effect of reducing the overall degree of the
        // transition constraint evaluations by `trace_len - 1`. Therefore the
        // total constraint evaluation degree is `constraint_degree * (trace_len - 1) -
        // (trace_len - 1) = (constraint_degree - 1) * (trace_len - 1)`
        let transition_ce_blowup_factor =
            utils::ceil_power_of_two(max_transition_constraint_degree.saturating_sub(1));

        [
            transition_ce_blowup_factor,
            terminal_ce_blowup_factor,
            boundary_ce_blowup_factor,
        ]
        .into_iter()
        .max()
        .unwrap()
    }

    /// Returns a degree that all constraints polynomials must be normalized to.
    fn composition_degree(&self) -> usize {
        let trace_len = self.trace_len();
        let ce_domain_size = trace_len * self.ce_blowup_factor();
        ce_domain_size - 1
    }

    fn lde_blowup_factor(&self) -> usize {
        self.options().lde_blowup_factor as usize
    }

    /// Validate properties of this air
    fn validate(&self) {
        let ce_blowup_factor = self.ce_blowup_factor();
        let lde_blowup_factor = self.lde_blowup_factor();
        assert!(
            ce_blowup_factor <= lde_blowup_factor,
            "constraint evaluation blowup factor {ce_blowup_factor} is 
            larger than the lde blowup factor {lde_blowup_factor}"
        );
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

    fn boundary_constraints(&self) -> &[Constraint<Self::Fq>] {
        &[]
    }

    fn transition_constraints(&self) -> &[Constraint<Self::Fq>] {
        &[]
    }

    fn terminal_constraints(&self) -> &[Constraint<Self::Fq>] {
        &[]
    }

    fn transition_constraint_divisor(&self) -> Divisor<Self::Fp> {
        let _timer = Timer::new("===TRANSITION DIVISOR===");
        let trace_domain = self.trace_domain();
        let last_trace_x = trace_domain.group_gen_inv;
        let degree = trace_domain.size() - 1;
        let lde_domain = self.lde_domain();
        let n = lde_domain.size();

        let mut lde = Vec::with_capacity_in(n, PageAlignedAllocator);
        lde.resize(n, Self::Fp::zero());

        // evaluates `(x - t_0)(x - t_1)...(x - t_n-1)` over the lde domain
        fill_vanishing_polynomial(&mut lde, &trace_domain, &lde_domain);

        // invert the vanishing polynomial evaluations
        // i.e. evaluations of `1 / (x - t_0)(x - t_1)...(x - t_n-1)`
        batch_inversion(&mut lde);

        // transition constraints apply to all rows except the last
        // multiplies out the last term of the vanishing polynomial
        // i.e. evaluations of `1 / (x - t_0)(x - t_1)...(x - t_n-2)`
        // Note: `t^(n-1) = t^(-1)`
        #[cfg(feature = "parallel")]
        let chunk_size = std::cmp::max(n / rayon::current_num_threads(), 1024);
        #[cfg(not(feature = "parallel"))]
        let chunk_size = n;
        ark_std::cfg_chunks_mut!(lde, chunk_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                let mut lde_x = lde_domain.element(i * chunk_size);
                chunk.iter_mut().for_each(|coeff| {
                    *coeff *= lde_x - last_trace_x;
                    lde_x *= &lde_domain.group_gen
                })
            });

        Divisor { lde, degree }
    }

    fn boundary_constraint_divisor(&self) -> Divisor<Self::Fp> {
        let _timer = Timer::new("===BOUNDARY DIVISOR===");

        let first_trace_x = Self::Fp::one();
        let lde_domain = self.lde_domain();
        let n = lde_domain.size();
        let mut lde = Vec::with_capacity_in(n, PageAlignedAllocator);
        lde.resize(n, Self::Fp::zero());

        #[cfg(feature = "parallel")]
        let chunk_size = std::cmp::max(n / rayon::current_num_threads(), 1024);
        #[cfg(not(feature = "parallel"))]
        let chunk_size = n;

        // evaluates `(x - t_0)` over the lde domain
        ark_std::cfg_chunks_mut!(lde, chunk_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                let mut lde_x = lde_domain.group_gen.pow([(i * chunk_size) as u64]);
                chunk.iter_mut().for_each(|coeff| {
                    *coeff = lde_domain.offset * lde_x - first_trace_x;
                    lde_x *= &lde_domain.group_gen
                })
            });

        // invert the evaluations
        // i.e. evaluations of `1 / (x - t_0)`
        batch_inversion(&mut lde);

        Divisor { lde, degree: 1 }
    }

    fn terminal_constraint_divisor(&self) -> Divisor<Self::Fp> {
        let _timer = Timer::new("===TERMINAL DIVISOR===");

        let last_trace_x = self.trace_domain().group_gen_inv();
        let lde_domain = self.lde_domain();
        let n = lde_domain.size();
        let mut lde = Vec::with_capacity_in(n, PageAlignedAllocator);
        lde.resize(n, lde_domain.offset);

        #[cfg(feature = "parallel")]
        let chunk_size = std::cmp::max(n / rayon::current_num_threads(), 1024);
        #[cfg(not(feature = "parallel"))]
        let chunk_size = n;

        // evaluates `(x - t_n-1)` over the lde domain
        // Note: `t^(n-1) = t^(-1)`
        ark_std::cfg_chunks_mut!(lde, chunk_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                let mut lde_x = lde_domain.group_gen.pow([(i * chunk_size) as u64]);
                chunk.iter_mut().for_each(|coeff| {
                    *coeff = *coeff * lde_x - last_trace_x;
                    lde_x *= &lde_domain.group_gen
                })
            });

        // invert the evaluations
        // i.e. evaluations of `1 / (x - t_n-1)`
        batch_inversion(&mut lde);

        Divisor { lde, degree: 1 }
    }

    fn get_challenges(&self, public_coin: &mut PublicCoin<impl Digest>) -> Challenges<Self::Fq> {
        // TODO: change get_challenge_indices to a constraint iterator and extract the
        // constraint with the highest index
        let num_challenges = self
            .all_constraint_elements()
            .iter()
            .filter_map(|element| match element {
                Element::Challenge(index) => Some(index + 1),
                _ => None,
            })
            .max()
            .unwrap_or(0);

        if num_challenges == 0 {
            Challenges::default()
        } else {
            let mut rng = public_coin.draw_rng();
            Challenges::new(&mut rng, num_challenges)
        }
    }

    fn get_hints(&self, _challenges: &Challenges<Self::Fq>) -> Hints<Self::Fq> {
        Hints::default()
    }

    // TODO: make this generic
    fn get_constraint_composition_coeffs(
        &self,
        public_coin: &mut PublicCoin<impl Digest>,
    ) -> Vec<(Self::Fq, Self::Fq)> {
        let mut rng = public_coin.draw_rng();
        (0..self.num_constraints())
            .map(|_| (Self::Fq::rand(&mut rng), Self::Fq::rand(&mut rng)))
            .collect()
    }

    // TODO: make this generic
    /// Output is of the form `(trace_coeffs, composition_coeffs,
    /// degree_adjustment_coeffs)`
    fn get_deep_composition_coeffs(
        &self,
        public_coin: &mut PublicCoin<impl Digest>,
    ) -> DeepCompositionCoeffs<Self::Fq> {
        let mut rng = public_coin.draw_rng();

        // execution trace coeffs
        let trace_info = self.trace_info();
        let mut base_trace_coeffs = Vec::new();
        for _ in 0..trace_info.num_base_columns {
            base_trace_coeffs.push((
                Self::Fq::rand(&mut rng),
                Self::Fq::rand(&mut rng),
                Self::Fq::rand(&mut rng),
            ));
        }

        let mut extension_trace_coeffs = Vec::new();
        for _ in 0..trace_info.num_extension_columns {
            extension_trace_coeffs.push((
                Self::Fq::rand(&mut rng),
                Self::Fq::rand(&mut rng),
                Self::Fq::rand(&mut rng),
            ));
        }

        // composition trace coeffs
        let num_composition_trace_cols = self.ce_blowup_factor();
        let mut composition_trace_coeffs = Vec::new();
        for _ in 0..num_composition_trace_cols {
            composition_trace_coeffs.push(Self::Fq::rand(&mut rng));
        }

        DeepCompositionCoeffs {
            base_trace: base_trace_coeffs,
            extension_trace: extension_trace_coeffs,
            constraints: composition_trace_coeffs,
            degree: (Self::Fq::rand(&mut rng), Self::Fq::rand(&mut rng)),
        }
    }

    fn all_constraint_elements(&self) -> Vec<Element> {
        // TODO: change get_challenge_indices to a constraint iterator and extract the
        // constraint with the highest index
        let mut indicies: Vec<Element> = [
            self.boundary_constraints(),
            self.transition_constraints(),
            self.terminal_constraints(),
        ]
        .into_iter()
        .flatten()
        .flat_map(|constraint| constraint.get_elements())
        .collect();
        indicies.sort();
        indicies.dedup();
        indicies
    }

    #[cfg(debug_assertions)]
    fn validate_constraints(
        &self,
        challenges: &Challenges<Self::Fq>,
        hints: &Hints<Self::Fq>,
        base_trace: &crate::Matrix<Self::Fp>,
        extension_trace: Option<&crate::Matrix<Self::Fq>>,
    ) {
        use crate::matrix::GroupItem;
        use crate::matrix::MatrixGroup;

        let mut execution_trace = MatrixGroup::new(vec![GroupItem::Fp(base_trace)]);
        if let Some(extension_trace) = extension_trace.as_ref() {
            execution_trace.append(GroupItem::Fq(extension_trace))
        }

        let mut col_indicies = vec![false; execution_trace.num_cols()];
        let mut challenge_indicies = vec![false; challenges.len()];
        let mut hint_indicies = vec![false; hints.len()];

        for element in self.all_constraint_elements() {
            match element {
                Element::Curr(i) | Element::Next(i) => col_indicies[i] = true,
                Element::Challenge(i) => challenge_indicies[i] = true,
                Element::Hint(i) => hint_indicies[i] = true,
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

        for (index, exists) in hint_indicies.into_iter().enumerate() {
            if !exists {
                // TODO: make assertion
                println!("WARN: hint at index {index} never used");
            }
        }

        let trace_rows = execution_trace.rows();
        let first_row = trace_rows.first().unwrap();
        let last_row = trace_rows.last().unwrap();

        // check boundary constraints
        for (i, constraint) in self.boundary_constraints().iter().enumerate() {
            let eval = constraint.evaluate(challenges, hints, first_row, &[]);
            assert!(eval.is_zero(), "boundary {i} mismatch");
        }

        // check terminal constraints
        for (i, constraint) in self.terminal_constraints().iter().enumerate() {
            let eval = constraint.evaluate(challenges, hints, last_row, &[]);
            assert!(eval.is_zero(), "terminal {i} mismatch");
        }

        // check transition constraints
        for (i, [curr, next]) in trace_rows.array_windows::<2>().enumerate() {
            for (j, constraint) in self.transition_constraints().iter().enumerate() {
                let eval = constraint.evaluate(challenges, hints, curr, next);
                assert!(eval.is_zero(), "transition {j} mismatch at row {i}");
            }
        }
    }

    fn num_constraints(&self) -> usize {
        //Vec<(Self::Fp, Self::Fp)> {
        self.boundary_constraints().len()
            + self.transition_constraints().len()
            + self.terminal_constraints().len()
    }
}

pub struct Divisor<F> {
    pub lde: GpuVec<F>,
    pub degree: usize,
}

impl<F: GpuField> Deref for Divisor<F> {
    type Target = GpuVec<F>;

    fn deref(&self) -> &Self::Target {
        &self.lde
    }
}
