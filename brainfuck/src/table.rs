use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_poly::EvaluationDomain;
use ark_poly::Evaluations;
use ark_poly::GeneralEvaluationDomain;
use ark_poly::Polynomial;
use legacy_algebra::Multivariate;

pub trait Table<F>
where
    F: FftField,
    F::BasePrimeField: FftField,
{
    /// The width of the table before extension
    const BASE_WIDTH: usize;

    /// The width of the table after extension
    const EXTENSION_WIDTH: usize;

    /// Returns the (non power of two) length of the execution trace
    fn len(&self) -> usize;

    /// Returns the power of two height of the execution trace
    fn height(&self) -> usize;

    // TODO: remove this is not the place. Should have a proof struct with trace
    // length
    fn set_height(&mut self, height: usize);

    /// Returns true if the table has no rows in it, otherwise false
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Pads the execution trace to n rows in length
    fn pad(&mut self, n: usize);

    /// TODO
    fn base_boundary_constraints() -> Vec<Multivariate<F>>;

    /// TODO
    fn base_transition_constraints() -> Vec<Multivariate<F>>;

    /// TODO
    fn extension_boundary_constraints(challenges: &[F]) -> Vec<Multivariate<F>>;

    /// TODO
    fn extension_transition_constraints(challenges: &[F]) -> Vec<Multivariate<F>>;

    /// TODO
    fn extension_terminal_constraints(
        &self,
        challenges: &[F],
        terminals: &[F],
    ) -> Vec<Multivariate<F>>;

    // TODO
    fn interpolant_degree(&self) -> usize;

    // TODO
    fn max_degree(&self) -> usize {
        // TODO: This NEEDS to be improved...
        let transition_constraints = Self::extension_transition_constraints(&[F::one(); 30]);
        let mut max_degree = 1;
        for air in transition_constraints {
            let degree_bounds = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH * 2];
            let degree =
                air.symbolic_degree_bound(&degree_bounds) - (self.height().saturating_sub(1));
            max_degree = max_degree.max(degree);
        }
        max_degree
    }

    fn boundary_quotients<D: EvaluationDomain<F>>(
        &self,
        codeword_len: usize,
        codewords: &[Evaluations<F, D>],
        challenges: &[F],
    ) -> Vec<Evaluations<F, D>> {
        println!("boundary_quotient");
        if codewords.is_empty() {
            return Vec::new();
        }

        let eval_domain = codewords[0].domain();
        // Evaluations of the polynomial `1/(x - o^0)` over the FRI domain
        let mut boundary_zerofier_inv = eval_domain
            .elements()
            .map(|x| x - F::one())
            .collect::<Vec<F>>();
        ark_ff::batch_inversion(&mut boundary_zerofier_inv);
        let mut quotient_codewords = Vec::new();

        let boundary_constraints = Self::extension_boundary_constraints(challenges);
        for constraint in boundary_constraints {
            println!("constraint");
            let mut quotient_codeword = Vec::new();
            for i in 0..codeword_len {
                let point = codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect::<Vec<F>>();
                quotient_codeword.push(constraint.evaluate(&point) * boundary_zerofier_inv[i]);
            }
            quotient_codewords.push(Evaluations::from_vec_and_domain(
                quotient_codeword,
                eval_domain,
            ));
        }

        quotient_codewords
    }

    // fn transition_quotients<D: EvaluationDomain<F>>(
    //     &self,
    //     codeword_len: usize,
    //     codewords: &[Evaluations<F, D>],
    //     challenges: &[F],
    // ) -> Vec<Evaluations<F, D>> {
    //     println!("trans_quotient");
    //     let mut quotient_codewords = Vec::new();
    //     // Evaluations of the polynomial (x - o^0)...(x - o^(n-1)) over the FRI
    // domain     // (x - o^0)...(x - o^(n-1)) = x^n - 1
    //     // Note: codeword_len = n * expansion_factor
    //     let omega = F::BasePrimeField::get_root_of_unity(codeword_len as
    // u64).unwrap();     let subgroup_zerofier = (0..codeword_len)
    //         .map(|i| {
    //             (F::BasePrimeField::GENERATOR * omega.pow([i as
    // u64])).pow([self.height() as u64])
    //                 - F::BasePrimeField::one()
    //         })
    //         .collect::<Vec<F::BasePrimeField>>();
    //     let mut subgroup_zerofier_inv = subgroup_zerofier;
    //     ark_ff::batch_inversion(&mut subgroup_zerofier_inv);

    //     // dot product of the inverse zerofier codeword and the codeword defined
    // by     // the evaluations of the polynomial (x - o^(n-1)). Note that
    // o^(n-1) is the     // inverse of `o`.
    //     let last_omicron = F::BasePrimeField::get_root_of_unity(self.height() as
    // u64)         .unwrap()
    //         .inverse()
    //         .unwrap();
    //     let zerofier_inv = (0..codeword_len)
    //         .map(|i| {
    //             subgroup_zerofier_inv[i]
    //                 * (F::BasePrimeField::GENERATOR * omega.pow([i as u64]) -
    //                   last_omicron)
    //         })
    //         .map(F::from_base_prime_field)
    //         .collect::<Vec<F>>();

    //     let row_step = codeword_len / self.height();
    //     println!("Row step is: {}", row_step);
    //     let transition_constraints =
    // Self::extension_transition_constraints(challenges);     for constraint in
    // transition_constraints {         println!("constraint");
    //         let mut quotient_codeword = Vec::new();
    //         // let combination_codeword = Vec::new();
    //         for i in 0..codeword_len {
    //             // if i % 1024 == 0 {
    //             //     println!("pos:{i} {codeword_len}");
    //             // }
    //             let point_lhs = codewords
    //                 .iter()
    //                 .map(|codeword| codeword[i])
    //                 .collect::<Vec<F>>();
    //             // TODO: HELP: why do we just wrap around here. Don't we need the
    // codeword to be             // extended so we have the right point?
    //             // Right. We are dealing with roots of unity so the evaluation is
    // o^(n-1)*o=o^0             let point_rhs = codewords
    //                 .iter()
    //                 .map(|codeword| codeword[(i + row_step) % codeword_len])
    //                 .collect::<Vec<F>>();
    //             let point = vec![point_lhs, point_rhs].concat();
    //             let evaluation = constraint.evaluate(&point);
    //             // combination_codeword.push(evaluation);
    //             quotient_codeword.push(evaluation * zerofier_inv[i]);
    //         }
    //         quotient_codewords.push(Evaluations::from_vec_and_domain(
    //             quotient_codeword,
    //             codewords[0].domain(),
    //         ));
    //     }

    //     quotient_codewords
    // }

    fn transition_quotients<D: EvaluationDomain<F>>(
        &self,
        codeword_len: usize,
        codewords: &[Evaluations<F, D>],
        challenges: &[F],
    ) -> Vec<Evaluations<F, D>> {
        println!("trans_quotient");
        if codewords.is_empty() {
            return Vec::new();
        }

        // TODO: move to params.
        let eval_domain = codewords[0].domain();
        let interp_domain = GeneralEvaluationDomain::new_subgroup(self.height()).unwrap();
        // Evaluations of the polynomial 1/((x - o^0)...(x - o^(n-1))) over the FRI
        // domain (x - o^0)...(x - o^(n-1)) = x^n - 1
        let mut subgroup_zerofier_inv = eval_domain
            .elements()
            .map(|e| interp_domain.evaluate_vanishing_polynomial(e))
            .collect::<Vec<F>>();
        ark_ff::batch_inversion(&mut subgroup_zerofier_inv);

        // Transition constraints apply to all rows of execution trace except the last
        // row. We need to change the inverse zerofier from being the
        // evaluations of the polynomial `1/((x - o^0)...(x - o^(n-1)))` to
        // `1/((x - o^0)...(x - o^(n-2)))`. This is achieved by performing the
        // dot product of the zerofier evaluations and the evaluations of the polynomial
        // (x - o^(n-1)). Note that o^(n-1) is the inverse of `o`.
        let last_x = interp_domain.element(self.height() - 1);
        let transition_zerofier_inv = subgroup_zerofier_inv
            .into_iter()
            .enumerate()
            .map(|(i, z)| z * (eval_domain.element(i) - last_x))
            .collect::<Vec<F>>();

        let row_step = eval_domain.size() / self.height();
        let row_step_old = codeword_len / self.height();
        println!("Row step is: {} {}", row_step, row_step_old);
        let transition_constraints = Self::extension_transition_constraints(challenges);
        let mut quotient_codewords = (0..transition_constraints.len())
            .map(|_| Evaluations::from_vec_and_domain(vec![F::zero(); codeword_len], eval_domain))
            .collect::<Vec<Evaluations<F, D>>>();

        for i in 0..codeword_len {
            // if i % 1024 == 0 {
            //     println!("pos:{i} {codeword_len}");
            // }
            let (point_lhs, point_rhs): (Vec<F>, Vec<F>) = codewords
                .iter()
                .map(|codeword| (codeword[i], codeword[(i + row_step) % codeword_len]))
                .unzip();
            // TODO: HELP: why do we just wrap around here. Don't we need the codeword to be
            // extended so we have the right point?
            // Right. We are dealing with roots of unity so the evaluation is o^(n-1)*o=o^0
            let point = vec![point_lhs, point_rhs].concat();
            for (j, constraint) in transition_constraints.iter().enumerate() {
                let evaluation = constraint.evaluate(&point);
                // combination_codeword.push(evaluation);
                quotient_codewords[j].evals[i] = evaluation * transition_zerofier_inv[i];
            }
        }

        quotient_codewords
    }

    fn terminal_quotients<D: EvaluationDomain<F>>(
        &self,
        codeword_len: usize,
        codewords: &[Evaluations<F, D>],
        challenges: &[F],
        terminals: &[F],
    ) -> Vec<Evaluations<F, D>> {
        println!("term_quotient");
        if codewords.is_empty() {
            return Vec::new();
        }

        let eval_domain = codewords[0].domain();
        let interp_domain = GeneralEvaluationDomain::<F>::new_subgroup(self.height()).unwrap();
        let last_interp_x = interp_domain.element(self.height() - 1);
        // evaluations of the polynomial (x - o^(n-1)). Note that o^(n-1) is the
        // inverse of `o`.
        let mut terminal_zerofier_inv = eval_domain
            .elements()
            .map(|x| x - last_interp_x)
            .collect::<Vec<F>>();
        ark_ff::batch_inversion(&mut terminal_zerofier_inv);

        let mut quotient_codewords = Vec::new();
        let terminal_constraints = self.extension_terminal_constraints(challenges, terminals);
        for constraint in terminal_constraints {
            println!("constraint");
            let mut quotient_codeword = Vec::new();
            for i in 0..codeword_len {
                let point = codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect::<Vec<F>>();
                // for i in Self::BASE_WIDTH {
                //     i * 30 % BASE_WIDTH
                // }
                quotient_codeword.push(constraint.evaluate(&point) * terminal_zerofier_inv[i]);
            }
            quotient_codewords.push(Evaluations::from_vec_and_domain(
                quotient_codeword,
                eval_domain,
            ));
        }

        quotient_codewords
    }

    fn all_quotients<D: EvaluationDomain<F>>(
        &self,
        codeword_len: usize,
        codewords: &[Evaluations<F, D>], // &[Vec<F>],
        challenges: &[F],
        terminals: &[F],
    ) -> Vec<Evaluations<F, D>> {
        // println!("pos:{codeword_len}");
        let boundary_quotients = self.boundary_quotients(codeword_len, codewords, challenges);
        println!("BOUNDARY");
        // for codeword in &boundary_quotients {
        //     determine_codeword_degree(codeword);
        // }
        let transition_quotients = self.transition_quotients(codeword_len, codewords, challenges);
        println!("TRANSITION {}-{}", Self::BASE_WIDTH, Self::EXTENSION_WIDTH);
        // for codeword in &transition_quotients {
        //     determine_codeword_degree(codeword);
        // }
        let terminal_quotients =
            self.terminal_quotients(codeword_len, codewords, challenges, terminals);
        println!("TERMINALS");
        // for codeword in &transition_quotients {
        //     determine_codeword_degree(codeword);
        // }
        vec![boundary_quotients, transition_quotients, terminal_quotients].concat()
    }

    fn boundary_quotient_degree_bounds(&self, challenges: &[F]) -> Vec<usize> {
        let max_degrees = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH];
        Self::extension_boundary_constraints(challenges)
            .into_iter()
            // TODO: improve this comment. It's late. Can't think
            // -1 at the end since the boundary is divided out
            .map(|constraint| constraint.symbolic_degree_bound(&max_degrees) - 1)
            .collect()
    }

    fn transition_quotient_degree_bounds(&self, challenges: &[F]) -> Vec<usize> {
        let max_degrees = vec![self.interpolant_degree(); 2 * Self::EXTENSION_WIDTH];
        Self::extension_transition_constraints(challenges)
            .into_iter()
            // TODO: improve this comment. It's late. Can't think
            // TODO: can cause overflow
            // divide out all 0 roots. +1 at the end since the last point is not checked
            .map(|constraint| constraint.symbolic_degree_bound(&max_degrees) - self.height() + 1)
            .collect()
    }

    fn terminal_quotient_degree_bounds(&self, challenges: &[F], terminals: &[F]) -> Vec<usize> {
        let max_degrees = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH];
        self.extension_terminal_constraints(challenges, terminals)
            .into_iter()
            // TODO: improve this comment. It's late. Can't think
            // -1 at the end since the terminal is divided out
            .map(|constraint| constraint.symbolic_degree_bound(&max_degrees) - 1)
            .collect()
    }

    fn all_quotient_degree_bounds(&self, challenges: &[F], terminals: &[F]) -> Vec<usize> {
        let boundary_degree_bounds = self.boundary_quotient_degree_bounds(challenges);
        let transition_degree_bounds = self.transition_quotient_degree_bounds(challenges);
        let terminal_degree_bounds = self.terminal_quotient_degree_bounds(challenges, terminals);
        vec![
            boundary_degree_bounds,
            transition_degree_bounds,
            terminal_degree_bounds,
        ]
        .concat()
    }

    // //
    // fn get_base_columns(&self) -> [Vec<Fx>; Self::BASE_WIDTH];
    // //
    // fn get_extension_columns(&self) -> [Vec<Fx>; Self::EXTENSION_WIDTH];

    // TODO
    fn set_matrix(&mut self, matrix: Vec<[F::BasePrimeField; Self::BASE_WIDTH]>);

    fn extend(&mut self, challenges: &[F], initials: &[F]);

    /// Computes the low degree extension of the base columns
    fn base_lde<D: EvaluationDomain<F::BasePrimeField>>(
        &mut self,
        eval_domain: D,
    ) -> Vec<Evaluations<F::BasePrimeField, D>>;

    /// Computes the low degree extension of all columns
    fn extension_lde<D: EvaluationDomain<F>>(&mut self, eval_domain: D) -> Vec<Evaluations<F, D>>;
}

fn determine_codeword_degree<F: FftField, D: EvaluationDomain<F>>(evaluations: &Evaluations<F, D>) {
    // let poly = DensePolynomial::from_coefficients_vec(
    //     legacy_algebra::number_theory_transform::inverse_number_theory_transform(evaluations),
    // );
    println!("Degree is: {}", evaluations.clone().interpolate().degree());
}

// fn determine_codeword_degree<F>(evaluations: &[F])
// where
//     F: Field,
//     F::BasePrimeField: FftField,
// {
//     let mut coeffs =
//         legacy_algebra::number_theory_transform::inverse_number_theory_transform(evaluations);
//     while coeffs.last().map_or(false, |v| v.is_zero()) {
//         coeffs.pop();
//     }

//     let mut degree = if coeffs.is_empty() {
//         0
//     } else {
//         coeffs.len() - 1
//     };

//     println!("Degree is {}", degree);
// }
