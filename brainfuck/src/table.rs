use crate::OpCode;
use legacy_algebra::batch_inverse;
use legacy_algebra::ExtensionOf;
use legacy_algebra::Felt;
use legacy_algebra::Multivariate;
use legacy_algebra::PrimeFelt;
use legacy_algebra::StarkFelt;
use std::marker::PhantomData;

pub trait Table<F, E = F>
where
    F: StarkFelt,
    E: Felt<BaseFelt = F> + ExtensionOf<F>,
{
    /// The width of the table before extension
    const BASE_WIDTH: usize;

    /// The width of the table after extension
    const EXTENSION_WIDTH: usize;

    /// Returns the (non power of two) length of the execution trace
    fn len(&self) -> usize;

    /// Returns the power of two height of the execution trace
    fn height(&self) -> usize;

    /// Returns true if the table has no rows in it, otherwise false
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Pads the execution trace to n rows in length
    fn pad(&mut self, n: usize);

    /// TODO
    fn base_boundary_constraints() -> Vec<Multivariate<E>>;

    /// TODO
    fn base_transition_constraints() -> Vec<Multivariate<E>>;

    /// TODO
    fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>>;

    /// TODO
    fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>>;

    /// TODO
    fn extension_terminal_constraints(
        &self,
        challenges: &[E],
        terminals: &[E],
    ) -> Vec<Multivariate<E>>;

    // TODO
    fn interpolant_degree(&self) -> usize;

    // TODO
    fn max_degree(&self) -> usize {
        // TODO: This NEEDS to be improved...
        let transition_constraints = Self::extension_transition_constraints(&[E::one(); 30]);
        let mut max_degree = 1;
        println!("{}, {}", Self::BASE_WIDTH, Self::EXTENSION_WIDTH);
        for air in transition_constraints {
            let degree_bounds = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH * 2];
            let degree = air.symbolic_degree_bound(&degree_bounds) - (self.len().saturating_sub(1));
            max_degree = max_degree.max(degree);
        }
        max_degree
    }

    fn boundary_quotients(
        &self,
        codeword_len: usize,
        codewords: &[Vec<E>],
        challenges: &[E],
    ) -> Vec<Vec<E>> {
        println!("boundary_quotient");
        // TODO: HELP: trying to understand zerofier here
        let mut quotient_codewords = Vec::new();
        let omega = F::get_root_of_unity(codeword_len.ilog2());
        // Evaluations of the polynomial (x - o^0) over the FRI domain
        let zerofier = (0..codeword_len)
            .map(|i| F::GENERATOR * omega.pow(&[i as u64]) - F::one())
            .collect::<Vec<F>>();
        let zerofier_inv = batch_inverse(&zerofier)
            .into_iter()
            .map(|z| z.unwrap().into())
            .collect::<Vec<E>>();

        let boundary_constraints = Self::extension_boundary_constraints(challenges);
        for constraint in boundary_constraints {
            println!("constraint");
            let mut quotient_codeword = Vec::new();
            for i in 0..codeword_len {
                let point = codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect::<Vec<E>>();
                quotient_codeword.push(constraint.evaluate(&point) * zerofier_inv[i]);
            }
            quotient_codewords.push(quotient_codeword);
        }

        quotient_codewords
    }

    fn transition_quotients(
        &self,
        codeword_len: usize,
        codewords: &[Vec<E>],
        challenges: &[E],
    ) -> Vec<Vec<E>> {
        println!("trans_quotient");
        let mut quotient_codewords = Vec::new();
        // Evaluations of the polynomial (x - o^0)...(x - o^(n-1)) over the FRI domain
        // (x - o^0)...(x - o^(n-1)) = x^n - 1
        // Note: codeword_len = n * expansion_factor
        let omega = F::get_root_of_unity(codeword_len.ilog2());
        let subgroup_zerofier = (0..codeword_len)
            .map(|i| {
                (F::GENERATOR * omega.pow(&[i as u64])).pow(&[self.height() as u64]) - F::one()
            })
            .collect::<Vec<F>>();
        let subgroup_zerofier_inv = batch_inverse(&subgroup_zerofier)
            .into_iter()
            .map(|z| z.unwrap())
            .collect::<Vec<F>>();
        // Transition constraints apply to all rows of execution trace except the last
        // row. We need to change the inverse zerofier from being the
        // evaluations of the polynomial `1/((x - o^0)...(x - o^(n-1)))` to
        // `1/((x - o^0)...(x - o^(n-2)))`. This is achieved by performing the
        // dot product of the inverse zerofier codeword and the codeword defined by
        // the evaluations of the polynomial (x - o^(n-1)). Note that o^(n-1) is the
        // inverse of `o`.
        let last_omicron = F::get_root_of_unity(self.height().ilog2())
            .inverse()
            .unwrap();
        let zerofier_inv = (0..codeword_len)
            .map(|i| {
                subgroup_zerofier_inv[i] * (F::GENERATOR * omega.pow(&[i as u64]) - last_omicron)
            })
            .map(E::from)
            .collect::<Vec<E>>();

        let row_step = codeword_len / self.height();
        let transition_constraints = Self::extension_transition_constraints(challenges);
        for constraint in transition_constraints {
            println!("constraint");
            let mut quotient_codeword = Vec::new();
            // let combination_codeword = Vec::new();
            for i in 0..codeword_len {
                if i % 1024 == 0 {
                    println!("pos:{i} {codeword_len}");
                }
                let point_lhs = codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect::<Vec<E>>();
                // TODO: HELP: why do we just wrap around here. Don't we need the codeword to be
                // extended so we have the right point?
                // Right. We are dealing with roots of unity so the evaluation is o^(n-1)*o=o^0
                let point_rhs = codewords
                    .iter()
                    .map(|codeword| codeword[(i + row_step) % codeword_len])
                    .collect::<Vec<E>>();
                let point = vec![point_lhs, point_rhs].concat();
                let evaluation = constraint.evaluate(&point);
                // combination_codeword.push(evaluation);
                quotient_codeword.push(evaluation * zerofier_inv[i]);
            }
            quotient_codewords.push(quotient_codeword);
        }

        quotient_codewords
    }

    fn terminal_quotients(
        &self,
        codeword_len: usize,
        codewords: &[Vec<E>],
        challenges: &[E],
        terminals: &[E],
    ) -> Vec<Vec<E>> {
        println!("term_quotient");
        let mut quotient_codewords = Vec::new();
        let omega = F::get_root_of_unity(codeword_len.ilog2());
        let last_omicron = F::get_root_of_unity(self.height().ilog2())
            .inverse()
            .unwrap();
        // evaluations of the polynomial (x - o^(n-1)). Note that o^(n-1) is the
        // inverse of `o`.
        let zerofier = (0..codeword_len)
            .map(|i| F::GENERATOR * omega.pow(&[i as u64]) - F::one())
            .collect::<Vec<F>>();
        let zerofier_inv = batch_inverse(&zerofier)
            .into_iter()
            .map(|z| z.unwrap().into())
            .collect::<Vec<E>>();

        let terminal_constraints = self.extension_terminal_constraints(challenges, terminals);
        for constraint in terminal_constraints {
            println!("constraint");
            let mut quotient_codeword = Vec::new();
            for i in 0..codeword_len {
                let point = codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect::<Vec<E>>();
                quotient_codeword.push(constraint.evaluate(&point) * zerofier_inv[i]);
            }
            quotient_codewords.push(quotient_codeword);
        }

        quotient_codewords
    }

    fn all_quotients(
        &self,
        codeword_len: usize,
        codewords: &[Vec<E>],
        challenges: &[E],
        terminals: &[E],
    ) -> Vec<Vec<E>> {
        println!("pos:{codeword_len}");
        let boundary_quotients = self.boundary_quotients(codeword_len, codewords, challenges);
        let transition_quotients = self.transition_quotients(codeword_len, codewords, challenges);
        let terminal_quotients =
            self.terminal_quotients(codeword_len, codewords, challenges, terminals);
        vec![boundary_quotients, transition_quotients, terminal_quotients].concat()
    }

    fn boundary_quotient_degree_bounds(&self, challenges: &[E]) -> Vec<usize> {
        let max_degrees = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH];
        Self::extension_boundary_constraints(challenges)
            .into_iter()
            // TODO: improve this comment. It's late. Can't think
            // -1 at the end since the boundary is divided out
            .map(|constraint| constraint.symbolic_degree_bound(&max_degrees) - 1)
            .collect()
    }

    fn transition_quotient_degree_bounds(&self, challenges: &[E]) -> Vec<usize> {
        let max_degrees = vec![self.interpolant_degree(); 2 * Self::EXTENSION_WIDTH];
        Self::extension_transition_constraints(challenges)
            .into_iter()
            // TODO: improve this comment. It's late. Can't think
            // divide out all 0 roots. +1 at the end since the last point is not checked
            .map(|constraint| constraint.symbolic_degree_bound(&max_degrees) - self.height() + 1)
            .collect()
    }

    fn terminal_quotient_degree_bounds(&self, challenges: &[E], terminals: &[E]) -> Vec<usize> {
        let max_degrees = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH];
        self.extension_terminal_constraints(challenges, terminals)
            .into_iter()
            // TODO: improve this comment. It's late. Can't think
            // -1 at the end since the terminal is divided out
            .map(|constraint| constraint.symbolic_degree_bound(&max_degrees) - 1)
            .collect()
    }

    fn all_quotient_degree_bounds(&self, challenges: &[E], terminals: &[E]) -> Vec<usize> {
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
    // fn get_base_columns(&self) -> [Vec<E>; Self::BASE_WIDTH];
    // //
    // fn get_extension_columns(&self) -> [Vec<E>; Self::EXTENSION_WIDTH];

    // TODO
    fn set_matrix(&mut self, matrix: Vec<[F; Self::BASE_WIDTH]>);

    fn extend(&mut self, challenges: &[E], initials: &[E]);

    /// Computes the low degree extension of the base columns
    fn base_lde(&mut self, offset: F, codeword_len: usize) -> Vec<Vec<E>>;

    /// Computes the low degree extension of all columns
    fn extension_lde(&mut self, offset: F, codeword_len: usize) -> Vec<Vec<E>>;
}
