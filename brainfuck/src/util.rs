use crate::OpCode;
use algebra::ExtensionOf;
use algebra::Felt;
use algebra::Multivariate;
use algebra::StarkFelt;
use mini_stark::number_theory_transform::fast_interpolate;
use mini_stark::number_theory_transform::inverse_number_theory_transform;
use mini_stark::polynomial::Polynomial;
use rand::Rng;

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

// Lifts a vector of field elements into array of extension field elements
pub(crate) fn lift<F: Felt, E: Felt + ExtensionOf<F>>(v: Vec<F>) -> Vec<E> {
    v.into_iter().map(|v| E::from(v)).collect()
}

pub(crate) fn interpolate_columns<F, const WIDTH: usize>(
    matrix: &[[F; WIDTH]],
    num_randomizers: usize,
) -> Vec<Polynomial<F>>
where
    F: Felt,
    F::BaseFelt: StarkFelt,
{
    assert!(matrix.len().is_power_of_two());
    let mut rng = rand::thread_rng();
    let n = matrix.len();
    let omicron = F::BaseFelt::get_root_of_unity(n.ilog2());
    let omega = F::BaseFelt::get_root_of_unity(n.ilog2() + 1);

    let matrix_domain = (0..n)
        .map(|i| omicron.pow(&[i as u64]).into())
        .collect::<Vec<F>>();
    // Odd indices to avoid collision with `matrix_domain`
    let randomizer_domain = (0..num_randomizers)
        .map(|i| omega.pow(&[1 + 2 * i as u64]).into())
        .collect::<Vec<F>>();
    let domain = vec![matrix_domain, randomizer_domain].concat();

    let mut polynomials = Vec::new();
    for col_idx in 0..WIDTH {
        println!("col");
        let trace_column = matrix.iter().map(|row| row[col_idx]).collect::<Vec<F>>();
        let randomizers = (0..num_randomizers)
            .map(|_| F::rand(&mut rng))
            .collect::<Vec<F>>();
        let values = vec![trace_column, randomizers].concat();
        assert_eq!(values.len(), domain.len());
        if num_randomizers == 0 {
            polynomials.push(Polynomial::new(inverse_number_theory_transform(&values)));
        } else {
            polynomials.push(fast_interpolate(&domain, &values))
        }
    }

    polynomials
}
