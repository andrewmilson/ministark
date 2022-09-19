use crate::OpCode;
use ark_ff::FftField;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use legacy_algebra::number_theory_transform::fast_interpolate;
use legacy_algebra::number_theory_transform::inverse_number_theory_transform;
use legacy_algebra::Multivariate;

pub fn instr_zerofier<F: Field>(curr_instr: &Multivariate<F>) -> Multivariate<F> {
    let mut accumulator = Multivariate::one();
    for opcode in OpCode::iterator() {
        let factor = curr_instr.clone() - F::from(Into::<usize>::into(opcode.clone()) as u64);
        accumulator = accumulator * factor;
    }
    accumulator
}

/// returns a polynomial in X that evaluates to 0 in all instructions except
/// for one provided
pub(crate) fn if_not_instr<F: Field>(
    instr: &OpCode,
    indeterminate: &Multivariate<F>,
) -> Multivariate<F> {
    let mut accumulator = Multivariate::one();
    for opcode in OpCode::iterator() {
        if opcode != instr {
            let opcode: u64 = opcode.clone().into();
            let factor = indeterminate.clone() - F::from(opcode);
            accumulator = accumulator * factor;
        }
    }
    accumulator
}

pub(crate) fn if_instr<F: Field>(
    instr: &OpCode,
    indeterminate: &Multivariate<F>,
) -> Multivariate<F> {
    Multivariate::constant(Into::<u64>::into(instr.clone()).into()) - indeterminate.clone()
}

// Lifts a vector of field elements into array of extension field elements
pub(crate) fn lift<F: Field>(v: Vec<F::BasePrimeField>) -> Vec<F> {
    v.into_iter().map(F::from_base_prime_field).collect()
}

pub(crate) fn interpolate_columns<F, const WIDTH: usize>(
    matrix: &[[F; WIDTH]],
    num_randomizers: usize,
) -> Vec<DensePolynomial<F>>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    assert!(matrix.len().is_power_of_two());
    let mut rng = rand::thread_rng();
    let n = matrix.len() as u64;
    let omicron = F::BasePrimeField::get_root_of_unity(n).unwrap();
    let omega = F::BasePrimeField::get_root_of_unity(n * 2).unwrap();

    let matrix_domain = (0..n)
        .map(|i| omicron.pow([i]))
        .map(F::from_base_prime_field)
        .collect::<Vec<F>>();
    // Odd indices to avoid collision with `matrix_domain`
    let randomizer_domain = (0..num_randomizers)
        .map(|i| omega.pow([1 + 2 * i as u64]))
        .map(F::from_base_prime_field)
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
            polynomials.push(DensePolynomial::from_coefficients_vec(
                inverse_number_theory_transform(&values),
            ));
        } else {
            polynomials.push(fast_interpolate(&domain, &values))
        }
    }

    polynomials
}
