use crate::OpCode;
use algebra::Felt;
use algebra::Multivariate;
use std::marker::PhantomData;

pub struct Table<E> {
    /// width of the matrix before extension
    pub base_width: usize,
    /// width of the matrix after extension
    pub full_width: usize,
    /// denotes the (non power of two) length of the execution trace
    pub length: usize,
    /// number of randomizers used when interpolating
    pub num_randomizers: usize,
    /// the height of the matrix (must be power of two)
    pub height: usize,
    /// two-dimensional array of field elements that represents the part of the
    /// algebraic execution trace that this table captures
    pub matrix: Vec<Vec<E>>,
}

impl<E: Felt> Table<E> {
    pub fn new(
        base_width: usize,
        full_width: usize,
        length: usize,
        num_randomizers: usize,
    ) -> Table<E> {
        Table {
            base_width,
            full_width,
            length,
            num_randomizers,
            height: length.next_power_of_two(),
            matrix: Vec::new(),
        }
    }
}

pub fn instr_zerofier<E: Felt>(curr_instr: &Multivariate<E>) -> Multivariate<E> {
    let mut accumulator = Multivariate::one();
    for opcode in OpCode::iterator() {
        let factor = curr_instr.clone() - E::from(Into::<usize>::into(opcode.clone()));
        accumulator = accumulator * factor;
    }
    accumulator
}
