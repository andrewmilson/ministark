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

/// returns a polynomial in X that evaluates to 0 in all instructions except
/// for one provided
pub(crate) fn if_not_instr<E: Felt>(
    instr: &OpCode,
    indeterminate: &Multivariate<E>,
) -> Multivariate<E> {
    let mut accumulator = Multivariate::one();
    for opcode in OpCode::iterator() {
        if opcode != instr {
            let factor = indeterminate.clone() - E::from(Into::<usize>::into(opcode.clone()));
            accumulator = accumulator * factor;
        }
    }
    accumulator
}

pub(crate) fn if_instr<E: Felt>(
    instr: &OpCode,
    indeterminate: &Multivariate<E>,
) -> Multivariate<E> {
    Multivariate::constant(Into::<usize>::into(instr.clone()).into()) - indeterminate.clone()
}
