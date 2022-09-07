use crate::OpCode;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use std::marker::PhantomData;

pub trait Table<E: PrimeFelt> {
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
        let transition_constraints = Self::extension_transition_constraints(&[E::one(); 50]);
        let mut max_degree = 1;
        for air in transition_constraints {
            let degree_bounds = vec![self.interpolant_degree(); Self::EXTENSION_WIDTH * 2];
            let degree = air.symbolic_degree_bound(&degree_bounds) - (self.len() - 1);
            max_degree = max_degree.max(degree);
        }
        max_degree
    }

    // //
    // fn get_base_columns(&self) -> [Vec<E>; Self::BASE_WIDTH];

    // //
    // fn get_extension_columns(&self) -> [Vec<E>; Self::EXTENSION_WIDTH];

    // TODO
    fn set_matrix(&mut self, matrix: Vec<[E; Self::BASE_WIDTH]>);
}
