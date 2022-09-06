use crate::OpCode;
use algebra::Felt;
use algebra::Multivariate;
use std::marker::PhantomData;

pub trait Table<E: Felt> {
    /// The width of the table before extension
    const BASE_WIDTH: usize;

    /// The width of the table after extension
    const EXTENSION_WIDTH: usize;

    /// Returns the (non power of two) length of the execution trace
    fn len(&self) -> usize;

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
    fn extension_terminal_constraints(challenges: &[E], terminals: &[E]) -> Vec<Multivariate<E>>;

    // TODO
    fn max_degree(&self) -> usize;

    // TODO
    fn set_matrix(&mut self, matrix: Vec<[E; Self::BASE_WIDTH]>);
}
