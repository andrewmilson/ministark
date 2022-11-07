use crate::air::Divisor;
use crate::challenges::Challenges;
use crate::matrix::GroupItem;
use crate::matrix::MatrixGroup;
use crate::merkle::MerkleTree;
use crate::utils::synthetic_divide;
use crate::utils::Timer;
use crate::Air;
use crate::Column;
use crate::Constraint;
use crate::Matrix;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use gpu_poly::prelude::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use sha2::Sha256;
use std::collections::BTreeMap;

pub struct ConstraintComposer<'a, A: Air> {
    air: &'a A,
    composition_coeffs: Vec<(A::Fp, A::Fp)>,
}

impl<'a, A: Air> ConstraintComposer<'a, A> {
    pub fn new(air: &'a A, composition_coeffs: Vec<(A::Fp, A::Fp)>) -> Self {
        ConstraintComposer {
            air,
            composition_coeffs,
        }
    }

    #[cfg(feature = "gpu")]
    fn generate_quotients_gpu(
        &self,
        mut constraint_evaluations: Matrix<A::Fp>,
        divisor: &Divisor<A::Fp>,
    ) -> Matrix<A::Fp> {
        // Generate quotients
        let _timer = Timer::new("=QQ====QUOTIENT EVALUATIOOOONS===QQ=");
        let library = &PLANNER.library;
        let command_queue = &PLANNER.command_queue;
        let command_buffer = command_queue.new_command_buffer();
        let multiplier = MulPowStage::<A::Fp>::new(library, divisor.len());
        let divisor_buffer = buffer_no_copy(command_queue.device(), divisor);
        // TODO: let's move GPU stuff out of here and make it readable in here.
        for evaluations in &mut constraint_evaluations.0 {
            let mut evaluations_buffer = buffer_no_copy(command_queue.device(), evaluations);
            multiplier.encode(
                command_buffer,
                &mut evaluations_buffer,
                &divisor_buffer,
                1,
                0,
            );
        }
        command_buffer.commit();
        command_buffer.wait_until_completed();

        drop(_timer);

        constraint_evaluations
    }

    #[cfg(not(feature = "gpu"))]
    fn generate_quotients_cpu(
        &self,
        mut constraint_evaluations: Matrix<A::Fp>,
        divisor: &Divisor<A::Fp>,
    ) -> Matrix<A::Fp> {
        let n = constraint_evaluations.num_rows();
        #[cfg(not(feature = "parallel"))]
        let chunk_size = n;
        #[cfg(feature = "parallel")]
        let chunk_size = std::cmp::max(n / rayon::current_num_threads().next_power_of_two(), 1024);

        for column in &mut constraint_evaluations.0 {
            ark_std::cfg_chunks_mut!(column, chunk_size)
                .enumerate()
                .for_each(|(chunk_offset, chunk)| {
                    let offset = chunk_offset * chunk_size;
                    for (i, quotient) in chunk.iter_mut().enumerate() {
                        *quotient *= divisor[offset + i]
                    }
                });
        }

        constraint_evaluations
    }

    fn generate_quotients(
        &self,
        challenges: &Challenges<A::Fp>,
        constraints: &[Constraint<A::Fp>],
        trace_lde: &Matrix<A::Fp>,
        divisor: &Divisor<A::Fp>,
    ) -> Matrix<A::Fp> {
        let trace_step = self.air.lde_blowup_factor();
        let _timer = Timer::new("=======CONSTRAINT EVALUATIOOOONS======");
        let all_evaluations =
            Constraint::evaluate_symbolic(constraints, challenges, trace_step, trace_lde);
        drop(_timer);

        #[cfg(feature = "gpu")]
        return self.generate_quotients_gpu(all_evaluations, divisor);
        #[cfg(not(feature = "gpu"))]
        return self.generate_quotients_cpu(all_evaluations, divisor);
    }

    pub fn evaluate(
        &mut self,
        challenges: &Challenges<A::Fp>,
        base_trace_lde: &Matrix<A::Fp>,
        extension_trace_lde: Option<&Matrix<A::Fp>>,
    ) -> Matrix<A::Fp> {
        // create a matrix group with all the LDEs we need for composition
        let mut lde_columns = MatrixGroup::default();

        // add execution trace LDE
        lde_columns.append(GroupItem::Fp(base_trace_lde));
        extension_trace_lde.map(|t| lde_columns.append(GroupItem::Fq(t)));

        let boundary_constraints = self.air.boundary_constraints();
        let boundary_divisor_idx = lde_columns.num_cols();
        let boundary_divisor = self.air.boundary_constraint_divisor();
        let _boundary_divisor_matrix = Matrix::new(vec![boundary_divisor.lde]);
        // add boundary constraint divisor LDE
        lde_columns.append(GroupItem::Fp(&_boundary_divisor_matrix));
        let boundary_iter = boundary_constraints
            .iter()
            .map(|c| (c, boundary_divisor_idx.curr(), boundary_divisor.degree));

        let transition_constraints = self.air.transition_constraints();
        let transition_divisor_idx = lde_columns.num_cols();
        let transition_divisor = self.air.transition_constraint_divisor();
        let _transition_divisor_matrix = Matrix::new(vec![transition_divisor.lde]);
        // add transition constraint divisor LDE
        lde_columns.append(GroupItem::Fp(&_transition_divisor_matrix));
        let transition_iter = transition_constraints
            .iter()
            .map(|c| (c, transition_divisor_idx.curr(), transition_divisor.degree));

        let terminal_constraints = self.air.terminal_constraints();
        let terminal_divisor_idx = lde_columns.num_cols();
        let terminal_divisor = self.air.terminal_constraint_divisor();
        let _terminal_divisor_matrix = Matrix::new(vec![terminal_divisor.lde]);
        // add terminal constraint divisor LDE
        lde_columns.append(GroupItem::Fp(&_terminal_divisor_matrix));
        let terminal_iter = terminal_constraints
            .iter()
            .map(|c| (c, terminal_divisor_idx.curr(), terminal_divisor.degree));

        // add degree adjustment LDEs
        let trace_degree = self.air.trace_len() - 1;
        let composition_degree = self.air.composition_degree();
        let lde_domain = self.air.lde_domain();
        let mut degree_adjustment_matricies = Vec::new();
        let mut degree_adjustment_map = BTreeMap::<usize, Constraint<A::Fp>>::new();
        for (constraint, _, divisor_degree) in boundary_iter
            .clone()
            .chain(transition_iter.clone())
            .chain(terminal_iter.clone())
        {
            let evaluation_degree = constraint.degree() * trace_degree - divisor_degree;
            assert!(evaluation_degree <= composition_degree);
            let degree_adjustment = composition_degree - evaluation_degree;

            degree_adjustment_map
                .entry(degree_adjustment)
                .or_insert_with(|| {
                    if degree_adjustment == 0 {
                        Constraint::from(A::Fp::one())
                    } else {
                        let col_idx = lde_columns.num_cols() + degree_adjustment_matricies.len();
                        let mut domain = lde_domain;
                        // TODO: this is hacky. fix
                        domain.offset = domain.offset.pow([degree_adjustment as u64]);
                        domain.group_gen = domain.group_gen.pow([degree_adjustment as u64]);
                        let column = domain.elements().collect::<Vec<_>>();
                        let matrix = Matrix::new(vec![column.to_vec_in(PageAlignedAllocator)]);
                        degree_adjustment_matricies.push(matrix);
                        col_idx.curr()
                    }
                });
        }

        // add degree adjustment LDEs
        for degree_adjustment_matrix in &degree_adjustment_matricies {
            lde_columns.append(GroupItem::Fp(&degree_adjustment_matrix));
        }

        let mut composition_constraint = Constraint::zero();
        for (constraint, divisor, divisor_degree) in
            boundary_iter.chain(transition_iter).chain(terminal_iter)
        {
            let evaluation_degree = constraint.degree() * trace_degree - divisor_degree;
            assert!(evaluation_degree <= composition_degree);
            let degree_adjustment = composition_degree - evaluation_degree;
            let degree_adjustor = degree_adjustment_map.get(&degree_adjustment).unwrap();

            let (alpha, beta) = self.composition_coeffs.pop().unwrap();
            composition_constraint += constraint * divisor * (degree_adjustor * alpha + beta);
        }

        let lde_step = self.air.lde_blowup_factor();
        lde_columns.evaluate_symbolic(&[composition_constraint], challenges, lde_step)

        // let boundary_divisor_index =

        // let boundary_constraints = air.boundary_constraints();
        // let transition_constraints = air.transition_constraints();
        // let terminal_constraints = air.terminal_constraints();

        // // let _timer = Timer::new("generating quotients");
        // let boundary_quotients = self.generate_quotients(
        //     challenges,
        //     boundary_constraints,
        //     trace_lde,
        //     &boundary_divisor,
        // );
        // let transition_quotients = self.generate_quotients(
        //     challenges,
        //     transition_constraints,
        //     trace_lde,
        //     &transition_divisor,
        // );
        // let terminal_quotients = self.generate_quotients(
        //     challenges,
        //     terminal_constraints,
        //     trace_lde,
        //     &terminal_divisor,
        // );
        // // drop(_timer);

        // let boundary_iter =
        //     zip(boundary_constraints, boundary_quotients.0).map(|(c, q)| (c,
        // q, &boundary_divisor)); let transition_iter =
        // zip(transition_constraints, transition_quotients.0)
        //     .map(|(c, q)| (c, q, &transition_divisor));
        // let terminal_iter =
        //     zip(terminal_constraints, terminal_quotients.0).map(|(c, q)| (c,
        // q, &terminal_divisor));

        // let _timer = Timer::new("asjusting degree");
        // let trace_degree = air.trace_len() - 1;
        // let composition_degree = air.composition_degree();
        // let mut groups = BTreeMap::new();
        // for (constraint, quotient, divisor) in
        //     boundary_iter.chain(transition_iter).chain(terminal_iter)
        // {
        //     // TODO: handle case when degree is 0?
        //     let evaluation_degree = constraint.degree() * trace_degree -
        // divisor.degree;     #[cfg(debug_assertions)]
        //     self.validate_quotient_degree(&quotient, evaluation_degree);
        //     assert!(evaluation_degree <= composition_degree);
        //     let degree_adjustment = composition_degree - evaluation_degree;

        //     let group = groups
        //         .entry(degree_adjustment)
        //         .or_insert_with(|| DegreeAdjustmentGroup {
        //             degree_adjustment,
        //             columns: Vec::new(),
        //             coeffs: Vec::new(),
        //         });
        //     group.columns.push(quotient);
        //     // TODO: don't use pop. use index
        //     group.coeffs.push(self.composition_coeffs.pop().unwrap());
        // }

        // // TODO: GPU
        // let lde_domain = air.lde_domain();
        // let mut accumulator = Matrix::new(vec![]);
        // for (_, group) in groups.into_iter() {
        //     let DegreeAdjustmentGroup {
        //         degree_adjustment,
        //         columns,
        //         coeffs,
        //     } = group;

        //     // TODO this step can be skipped
        //     let (alpha_cols, beta_cols) = ark_std::cfg_iter!(columns)
        //         .zip(coeffs)
        //         .map(|(column, (alpha, beta))| {
        //             let n = column.len();
        //             let mut alpha_col = Vec::with_capacity_in(n,
        // PageAlignedAllocator);             let mut beta_col =
        // Vec::with_capacity_in(n, PageAlignedAllocator);

        //             for v in column {
        //                 alpha_col.push(alpha * v);
        //                 beta_col.push(beta * v);
        //             }

        //             (alpha_col, beta_col)
        //         })
        //         .unzip();

        //     let mut alpha_col = Matrix::new(alpha_cols).sum_columns();

        //     let _timer = Timer::new("===DEG ADJUSTOR===");

        //     if degree_adjustment != 0 {
        //         // TODO: make parallel. also this is hacky and needs to go
        //         // modify domain to go from x to x^degree_adjustment
        //         let mut adjust_domain = lde_domain;
        //         adjust_domain.offset =
        // adjust_domain.offset.pow([degree_adjustment as u64]);
        //         adjust_domain.group_gen =
        // adjust_domain.group_gen.pow([degree_adjustment as u64]);
        //         for (v, x) in zip(&mut alpha_col.0[0],
        // adjust_domain.elements()) {             *v *= x
        //         }
        //     }

        //     // TODO: multithreaded but not faster. hmmm..
        //     // let adjust_offset = lde_domain.offset.pow([degree_adjustment
        // as u64]);     // let adjust_group_gen =
        // lde_domain.group_gen.pow([degree_adjustment as u64]);     //
        // Radix2EvaluationDomain::distribute_powers_and_mul_by_const(
        //     //     &mut beta_col.0[0],
        //     //     adjust_group_gen,
        //     //     adjust_offset,
        //     // );

        //     drop(_timer);

        //     accumulator.append(alpha_col);
        //     accumulator.append(Matrix::new(beta_cols));
        // }
        // drop(_timer);

        // accumulator.sum_columns()
    }

    fn trace_polys(&self, composed_evaluations: Matrix<A::Fp>) -> Matrix<A::Fp> {
        assert_eq!(composed_evaluations.num_cols(), 1);
        let mut composition_poly = composed_evaluations.into_polynomials(self.air.lde_domain());

        let composition_poly_degree = composition_poly.column_degrees()[0];
        assert_eq!(composition_poly_degree, self.air.composition_degree());
        assert_eq!(composition_poly_degree, self.air.composition_degree());
        composition_poly.0[0].truncate(composition_poly_degree + 1);

        let num_composition_trace_cols = self.air.ce_blowup_factor();
        assert_eq!(
            composition_poly.num_rows() / self.air.trace_len(),
            num_composition_trace_cols
        );
        let composition_trace_poly = if num_composition_trace_cols == 1 {
            composition_poly
        } else {
            Matrix::from_rows(
                composition_poly.0[0]
                    .chunks(num_composition_trace_cols)
                    .map(|chunk| chunk.to_vec())
                    .collect(),
            )
        };

        composition_trace_poly
    }

    /// builds a commitment to the composed trace polynomial.
    /// Output is of the form `(lde, poly, lde_merkle_tree)`
    pub fn build_commitment(
        mut self,
        challenges: &Challenges<A::Fp>,
        base_trace_lde: &Matrix<A::Fp>,
        extension_trace_lde: Option<&Matrix<A::Fp>>,
    ) -> (Matrix<A::Fp>, Matrix<A::Fp>, MerkleTree<Sha256>) {
        let _timer = Timer::new("constraint evaluation");
        let composed_evaluations = self.evaluate(challenges, base_trace_lde, extension_trace_lde);
        drop(_timer);
        let composition_trace_polys = self.trace_polys(composed_evaluations);
        let composition_trace_lde = composition_trace_polys.evaluate(self.air.lde_domain());
        let merkle_tree = composition_trace_lde.commit_to_rows();
        (composition_trace_lde, composition_trace_polys, merkle_tree)
    }

    /// Ensures the post-division degree of the polynomial matches the expected
    #[cfg(debug_assertions)]
    fn validate_quotient_degree(&self, quotient: &[A::Fp], expected_degree: usize) {
        use ark_poly::univariate::DensePolynomial;
        use ark_poly::DenseUVPolynomial;
        use ark_poly::Polynomial;

        let coeffs = self.air.lde_domain().ifft(quotient);
        let poly = DensePolynomial::from_coefficients_vec(coeffs);

        // TODO: set up errors
        assert_eq!(expected_degree, poly.degree());
    }
}

pub struct DeepPolyComposer<'a, A: Air> {
    air: &'a A,
    composition_coeffs: DeepCompositionCoeffs<A::Fp>,
    z: A::Fp,
    poly: GpuVec<A::Fp>,
}

impl<'a, A: Air> DeepPolyComposer<'a, A> {
    pub fn new(air: &'a A, composition_coeffs: DeepCompositionCoeffs<A::Fp>, z: A::Fp) -> Self {
        let poly = Vec::with_capacity_in(air.trace_len(), PageAlignedAllocator);
        DeepPolyComposer {
            air,
            composition_coeffs,
            z,
            poly,
        }
    }

    pub fn add_execution_trace_polys(
        &mut self,
        base_trace_polys: Matrix<A::Fp>,
        extension_trace_polys: Option<Matrix<A::Fp>>,
        ood_evals: Vec<A::Fp>,
        ood_evals_next: Vec<A::Fp>,
    ) {
        assert!(self.poly.is_empty());

        let trace_domain = self.air.trace_domain();
        let next_z = self.z * trace_domain.group_gen();
        let n = trace_domain.size();

        let mut t1_composition = Vec::with_capacity_in(n, PageAlignedAllocator);
        t1_composition.resize(n, A::Fp::zero());
        let mut t2_composition = Vec::with_capacity_in(n, PageAlignedAllocator);
        t2_composition.resize(n, A::Fp::zero());

        // TODO: clean up code
        for (i, poly) in base_trace_polys.iter().enumerate() {
            let (alpha, beta, _) = self.composition_coeffs.base_trace[i];
            let ood_eval = ood_evals[i];
            let ood_eval_next = ood_evals_next[i];

            acc_trace_poly::<A::Fp>(&mut t1_composition, poly, ood_eval, alpha);
            acc_trace_poly::<A::Fp>(&mut t2_composition, poly, ood_eval_next, beta);
        }

        if let Some(extension_trace_polys) = extension_trace_polys {
            let num_base_columns = self.air.trace_info().num_base_columns;
            for (i, poly) in extension_trace_polys.iter().enumerate() {
                let (alpha, beta, _) = self.composition_coeffs.extension_trace[i];
                let ood_eval = ood_evals[num_base_columns + i];
                let ood_eval_next = ood_evals_next[num_base_columns + i];

                acc_trace_poly::<A::Fp>(&mut t1_composition, poly, ood_eval, alpha);
                acc_trace_poly::<A::Fp>(&mut t2_composition, poly, ood_eval_next, beta);
            }
        }

        // TODO: multithread
        synthetic_divide(&mut t1_composition, 1, self.z);
        synthetic_divide(&mut t2_composition, 1, next_z);

        for (t1, t2) in t1_composition.into_iter().zip(t2_composition) {
            self.poly.push(t1 + t2)
        }

        // TODO:
        // Check that the degree has reduced by 1 as a result of the divisions
        assert!(self.poly.last().unwrap().is_zero());
    }

    pub fn add_composition_trace_polys(&mut self, mut polys: Matrix<A::Fp>, ood_evals: Vec<A::Fp>) {
        assert!(!self.poly.is_empty());

        let z_n = self.z.pow([polys.num_cols() as u64]);

        ark_std::cfg_iter_mut!(polys.0)
            .zip(ood_evals)
            .for_each(|(poly, ood_eval)| {
                poly[0] -= ood_eval;
                synthetic_divide(poly, 1, z_n);
            });

        for (i, poly) in polys.0.into_iter().enumerate() {
            let alpha = self.composition_coeffs.constraints[i];
            for (lhs, rhs) in self.poly.iter_mut().zip(poly) {
                *lhs += rhs * alpha;
            }
        }

        // Check that the degree has reduced by 1 as a result of the divisions
        assert!(self.poly.last().unwrap().is_zero());
    }

    pub fn into_deep_poly(mut self) -> Matrix<A::Fp> {
        let (alpha, beta) = self.composition_coeffs.degree;

        // TODO: consider making multithreaded
        // Adjust the degree
        // P(x) * (alpha + x * beta)
        let mut last = A::Fp::zero();
        for coeff in &mut self.poly {
            let tmp = *coeff;
            *coeff *= alpha;
            *coeff += last * beta;
            last = tmp;
        }

        Matrix::new(vec![self.poly])
    }
}

pub struct DeepCompositionCoeffs<F> {
    /// Base trace polynomial composition coefficients
    pub base_trace: Vec<(F, F, F)>,
    /// Base trace polynomial composition coefficients
    pub extension_trace: Vec<(F, F, F)>,
    /// Composition poly trace column composition coefficients
    pub constraints: Vec<F>,
    /// Degree adjustment composition coefficients
    pub degree: (F, F),
}

struct DegreeAdjustmentGroup<F> {
    degree_adjustment: usize,
    columns: Vec<GpuVec<F>>,
    coeffs: Vec<(F, F)>,
}

/// Computes (P(x) - value) * k
/// Source https://github.com/novifinancial/winterfell
fn acc_trace_poly<F: GpuField>(acc: &mut [F], poly: &[F], value: F, k: F) {
    // P(x) * k
    ark_std::cfg_iter_mut!(acc)
        .zip(poly)
        .for_each(|(v, coeff)| *v += k * coeff);
    // -value * k
    acc[0] -= value * k;
}
