use crate::OpCode;
use algebra::batch_inverse;
use algebra::ExtensionOf;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use algebra::StarkFelt;
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
        for air in transition_constraints {
            let degree_bounds = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH * 2];
            let degree = air.symbolic_degree_bound(&degree_bounds) - (self.len() - 1);
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
        // TODO: HELP: trying to understand zerofier here
        let mut quotient_codewords = Vec::new();
        let omega = F::get_root_of_unity(codeword_len.ilog2());
        // Evaluations of the polynomial (x - ω^0) over the FRI domain
        let zerofier = (0..codeword_len)
            .map(|i| F::GENERATOR * omega.pow(&[i as u64]) - F::one())
            .collect::<Vec<F>>();
        let zerofier_inv = batch_inverse(&zerofier)
            .into_iter()
            .map(|z| E::from(z.unwrap()))
            .collect::<Vec<E>>();

        let boundary_constraints = Self::extension_boundary_constraints(challenges);
        for constraint in boundary_constraints {
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
    ) -> Vec<Vec<E>> {
        let boundary_quotients = self.boundary_quotients(codeword_len, codewords, challenges);
        let transition_quotients = todo!();
        let terminal_quotients = todo!();
        vec![boundary_quotients, transition_quotients, terminal_quotients].concat()
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
