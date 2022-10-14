use crate::air::Divisor;
use crate::challenges::Challenges;
use crate::merkle::MerkleTree;
use crate::Air;
use crate::Constraint;
use crate::Matrix;
use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::PLANNER;
use fast_poly::stage::MulPowStage;
use fast_poly::utils::buffer_no_copy;
use fast_poly::GpuField;
use fast_poly::GpuVec;
use sha2::Sha256;
use std::collections::BTreeMap;
use std::iter::zip;

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
        &mut self,
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

        let composition_degree = air.composition_degree();
        let mut groups = BTreeMap::new();
        for (constraint, quotient, divisor) in
            boundary_iter.chain(transition_iter).chain(terminal_iter)
        {
            // TODO: handle case when degree is 0?
            let trace_degree = air.trace_len() - 1;
            let evaluation_degree = constraint.degree() * trace_degree - divisor.degree;
            let degree_adjustment = composition_degree - evaluation_degree;
            let group = groups
                .entry(degree_adjustment)
                .or_insert_with(|| DegreeAdjustmentGroup {
                    degree_adjustment,
                    columns: Vec::new(),
                    coeffs: Vec::new(),
                });
            group.columns.push(quotient);
            group.coeffs.push(self.composition_coeffs.pop().unwrap());
        }

        // TODO: GPU
        let lde_domain = air.lde_domain();
        let mut accumulator = Matrix::new(vec![]);
        for (_, group) in groups.into_iter() {
            let DegreeAdjustmentGroup {
                degree_adjustment,
                columns,
                coeffs,
            } = group;

            let mut alpha_cols = Vec::new();
            let mut beta_cols = Vec::new();
            for (column, (alpha, beta)) in zip(columns, coeffs) {
                let mut alpha_col = Vec::with_capacity_in(column.len(), PageAlignedAllocator);
                let mut beta_col = Vec::with_capacity_in(column.len(), PageAlignedAllocator);
                for v in column {
                    alpha_col.push(alpha * v);
                    beta_col.push(beta * v);
                }
                alpha_cols.push(alpha_col);
                beta_cols.push(beta_col);
            }

            let alpha_col = Matrix::new(alpha_cols).sum_columns();
            let mut beta_col = Matrix::new(beta_cols).sum_columns();
            for (v, x) in zip(&mut beta_col.0[0], lde_domain.elements()) {
                *v *= x.pow([degree_adjustment as u64])
            }

            accumulator.append(alpha_col);
            accumulator.append(beta_col);
        }

        accumulator.sum_columns()
    }

    fn trace_polys(&self, composed_evaluations: Matrix<A::Fp>) -> Matrix<A::Fp> {
        assert!(composed_evaluations.num_cols() == 1);
        let mut composition_poly = composed_evaluations.interpolate_columns(self.air.lde_domain());

        let composition_poly_degree = composition_poly.column_degrees()[0];
        assert_eq!(composition_poly_degree, self.air.composition_degree());
        composition_poly.0[0].truncate(composition_poly_degree + 1);

        let num_composition_trace_cols = self.air.ce_blowup_factor();
        assert_eq!(
            composition_poly.num_rows() / self.air.trace_len(),
            num_composition_trace_cols
        );
        let composition_trace_poly = Matrix::from_rows(
            composition_poly.0[0]
                .chunks(num_composition_trace_cols)
                .map(|chunk| chunk.to_vec())
                .collect(),
        );

        composition_trace_poly
    }

    /// builds a commitment to the composed trace polynomial.
    /// Output is of the form `(lde, poly, lde_merkle_tree)`
    pub fn build_commitment(
        mut self,
        challenges: &Challenges<A::Fp>,
        execution_trace_lde: &Matrix<A::Fp>,
    ) -> (Matrix<A::Fp>, Matrix<A::Fp>, MerkleTree<Sha256>) {
        let composed_evaluations = self.evaluate(challenges, execution_trace_lde);
        let composition_trace_polys = self.trace_polys(composed_evaluations);
        let composition_trace_lde = composition_trace_polys.evaluate(self.air.lde_domain());
        let merkle_tree = composition_trace_lde.commit_to_rows();
        (composition_trace_lde, composition_trace_polys, merkle_tree)
    }
}

pub struct DeepPolyComposer<'a, A: Air> {
    air: &'a A,
    coeffs: DeepCompositionCoeffs<A::Fp>,
    z: A::Fp,
}

impl<'a, A: Air> DeepPolyComposer<'a, A> {
    pub fn new(air: &'a A, coeffs: DeepCompositionCoeffs<A::Fp>, z: A::Fp) -> Self {
        DeepPolyComposer { air, coeffs, z }
    }

    pub fn add_execution_trace_polys(
        &mut self,
        polys: Matrix<A::Fp>,
        ood_evals: Vec<A::Fp>,
        ood_evals_next: Vec<A::Fp>,
    ) {
        let trace_domain = self.air.trace_domain();
        let next_z = self.z * trace_domain.group_gen();
        let n = trace_domain.size();

        let mut t1_composition = Vec::with_capacity_in(n, PageAlignedAllocator);
        t1_composition.resize(n, A::Fp::zero());
        let mut t2_composition = Vec::with_capacity_in(n, PageAlignedAllocator);
        t2_composition.resize(n, A::Fp::zero());

        for (i, poly) in polys.iter().enumerate() {
            let (alpha, beta, _) = self.coeffs.trace.pop().unwrap();
            let ood_eval = ood_evals[i];
            let ood_eval_next = ood_evals_next[i];

            // compute T' += (T - T(z)) * alpha
            acc_trace_poly::<A::Fp>(&mut t1_composition, poly, ood_eval, alpha);
            
            // compute T'' += (T - T(z * g)) * beta
            acc_trace_poly::<A::Fp>(&mut t2_composition, poly, ood_eval_next, beta);
        }

        
    }


    pub fn add_composition_trace_polys(&mut self, polys: Matrix<A::Fp>, ood_evals: Vec<A::Fp>) {}
}

pub struct DeepCompositionCoeffs<F> {
    /// Execution trace polynomial composition coefficients
    pub trace: Vec<(F, F, F)>,
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
        .for_each(|(v, coeff)| *v += coeff * k)
    // -value * k
    acc[0] -= value * k;
}
