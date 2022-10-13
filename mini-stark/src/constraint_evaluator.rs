use crate::air::Divisor;
use crate::challenges::Challenges;
use crate::Air;
use crate::Constraint;
use crate::Matrix;
use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::PLANNER;
use fast_poly::stage::MulPowStage;
use fast_poly::utils::buffer_no_copy;
use fast_poly::GpuVec;
use std::collections::BTreeMap;
use std::iter::zip;
use std::ops::Deref;

pub struct ConstraintEvaluator<'a, A: Air> {
    air: &'a A,
    domain: Radix2EvaluationDomain<A::Fp>,
    composition_coefficients: Vec<(A::Fp, A::Fp)>,
    domain_elements: Vec<A::Fp, PageAlignedAllocator>,
}

impl<'a, A: Air> ConstraintEvaluator<'a, A> {
    pub fn new(air: &'a A, composition_coefficients: Vec<(A::Fp, A::Fp)>) -> Self {
        let domain = air.lde_domain();
        let mut domain_elements = Vec::with_capacity_in(domain.size(), PageAlignedAllocator);
        // TODO: parallelise
        for element in domain.elements() {
            domain_elements.push(element);
        }
        ConstraintEvaluator {
            air,
            domain,
            composition_coefficients,
            domain_elements,
        }
    }

    fn generate_quotients(
        &self,
        challenges: &Challenges<A::Fp>,
        constraints: &[Constraint<A::Fp>],
        trace_lde: &Matrix<A::Fp>,
        divisor: &Divisor<A::Fp>,
    ) -> Matrix<A::Fp> {
        let trace_step = self.air.lde_blowup_factor() as usize;
        let mut all_evaluations = Matrix::join(
            constraints
                .iter()
                .map(|constraint| constraint.evaluate_symbolic(challenges, trace_step, trace_lde))
                .collect(),
        );

        // Generate quotients
        let library = &PLANNER.library;
        let command_queue = &PLANNER.command_queue;
        let command_buffer = command_queue.new_command_buffer();
        let multiplier = MulPowStage::<A::Fp>::new(library, divisor.len(), 0);
        let divisor_buffer = buffer_no_copy(command_queue.device(), divisor);
        // TODO: let's move GPU stuff out of here and make it readable in here.
        for evaluations in &mut all_evaluations.0 {
            let mut evaluations_buffer = buffer_no_copy(command_queue.device(), evaluations);
            multiplier.encode(command_buffer, &mut evaluations_buffer, &divisor_buffer, 0);
        }
        command_buffer.commit();
        command_buffer.wait_until_completed();

        all_evaluations
    }

    // output of the form `(boundary, transition, terminal)`
    pub fn evaluate(
        self,
        challenges: &Challenges<A::Fp>,
        trace_lde: &Matrix<A::Fp>,
    ) -> Matrix<A::Fp> {
        let air = self.air;

        let boundary_divisor = air.boundary_constraint_divisor();
        let transition_divisor = air.transition_constraint_divisor();
        let terminal_divisor = air.terminal_constraint_divisor();

        let boundary_constraints = air.boundary_constraints();
        let transition_constraints = air.transition_constraints();
        let terminal_constraints = air.terminal_constraints();

        let boundary_quotients = self.generate_quotients(
            challenges,
            boundary_constraints,
            trace_lde,
            &boundary_divisor,
        );
        let transition_quotients = self.generate_quotients(
            challenges,
            transition_constraints,
            trace_lde,
            &transition_divisor,
        );
        let terminal_quotients = self.generate_quotients(
            challenges,
            terminal_constraints,
            trace_lde,
            &terminal_divisor,
        );

        let boundary_iter =
            zip(boundary_constraints, boundary_quotients.0).map(|(c, q)| (c, q, &boundary_divisor));
        let transition_iter = zip(transition_constraints, transition_quotients.0)
            .map(|(c, q)| (c, q, &transition_divisor));
        let terminal_iter = zip(terminal_constraints, terminal_quotients.0)
            .map(|(c, q)| (c, q, &transition_divisor));

        let mut coefficients = self.composition_coefficients.into_iter();
        let mut groups = BTreeMap::new();
        for (constraint, quotient, divisor) in
            boundary_iter.chain(transition_iter).chain(terminal_iter)
        {
            let trace_degree = air.trace_len() - 1;
            let evaluation_degree = constraint.degree() * trace_degree - divisor.degree;
            let target_degree = air.composition_degree();
            println!("target: {}", target_degree);
            let degree_adjustment = target_degree - evaluation_degree;
            let group = groups
                .entry(degree_adjustment)
                .or_insert_with(|| DegreeAdjustmentGroup {
                    degree_adjustment,
                    columns: Vec::new(),
                    coefficients: Vec::new(),
                });
            group.columns.push(quotient);
            group.coefficients.push(coefficients.next().unwrap());
        }

        // TODO: GPU
        let lde_domain = air.lde_domain();
        let mut accumulator = Matrix::new(vec![]);
        for (_, group) in groups.into_iter() {
            let DegreeAdjustmentGroup {
                degree_adjustment,
                columns,
                coefficients,
            } = group;

            println!("Adjustment: {}", degree_adjustment);
            let mut alpha_cols = Vec::new();
            let mut beta_cols = Vec::new();
            for (column, (alpha, beta)) in zip(columns, coefficients) {
                let mut alpha_col = Vec::with_capacity_in(column.len(), PageAlignedAllocator);
                let mut beta_col = Vec::with_capacity_in(column.len(), PageAlignedAllocator);
                for v in column {
                    alpha_col.push(alpha * v);
                    beta_col.push(beta * v);
                }
                alpha_cols.push(alpha_col);
                beta_cols.push(beta_col);
            }

            let beta_coeffs = Matrix::new(
                alpha_cols
                    .iter()
                    .map(|col| col.to_vec_in(PageAlignedAllocator))
                    .collect(),
            )
            .interpolate_columns(lde_domain);
            let alpha_coeffs = Matrix::new(
                alpha_cols
                    .iter()
                    .map(|col| col.to_vec_in(PageAlignedAllocator))
                    .collect(),
            )
            .interpolate_columns(lde_domain);

            for (alpha_col, beta_col) in zip(&beta_coeffs.0, &alpha_coeffs.0) {
                let alpha_poly = DensePolynomial::from_coefficients_slice(alpha_col);
                let beta_poly = DensePolynomial::from_coefficients_slice(beta_col);
                println!("alpha degree is: {}", alpha_poly.degree());
                println!("beta degree is: {}", beta_poly.degree());
            }

            let alpha_col = Matrix::new(alpha_cols).sum_columns();
            let mut beta_col = Matrix::new(beta_cols).sum_columns();
            for (v, x) in zip(&mut beta_col.0[0], lde_domain.elements()) {
                *v *= x.pow([degree_adjustment as u64])
            }

            let beta_coeffs = Matrix::new(
                beta_col
                    .iter()
                    .map(|col| col.to_vec_in(PageAlignedAllocator))
                    .collect(),
            )
            .interpolate_columns(lde_domain);

            for beta_col in &beta_coeffs.0 {
                let beta_poly = DensePolynomial::from_coefficients_slice(beta_col);
                println!("BETA degree is: {}", beta_poly.degree());
            }

            accumulator.append(alpha_col);
            accumulator.append(beta_col);
        }

        accumulator.sum_columns()
    }
}

struct DegreeAdjustmentGroup<F> {
    degree_adjustment: usize,
    columns: Vec<GpuVec<F>>,
    coefficients: Vec<(F, F)>,
}
