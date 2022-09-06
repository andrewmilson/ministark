use crate::OpCode;
use algebra::Felt;
use algebra::Multivariate;

pub fn instr_zerofier<E: Felt>(curr_instr: &Multivariate<E>) -> Multivariate<E> {
    let mut accumulator = Multivariate::one();
    for opcode in OpCode::iterator() {
        let factor = curr_instr.clone() - E::from(Into::<usize>::into(opcode.clone()));
        accumulator = accumulator * factor;
    }
    accumulator
}
