use crate::OpCode;
use ark_ff::FftField;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Evaluations;
use ark_poly::GeneralEvaluationDomain;
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
pub fn lift<F: FftField>(
    e: Vec<Evaluations<F::BasePrimeField, impl EvaluationDomain<F::BasePrimeField>>>,
) -> Vec<Evaluations<F, GeneralEvaluationDomain<F>>> {
    e.into_iter()
        .map(|codeword| {
            let domain = codeword.domain();
            let offset = F::from_base_prime_field(domain.offset());
            Evaluations::from_vec_and_domain(
                codeword
                    .evals
                    .into_iter()
                    .map(F::from_base_prime_field)
                    .collect(),
                GeneralEvaluationDomain::new_coset(domain.size(), offset).unwrap(),
            )
        })
        .collect()
}

pub(crate) fn interpolate_columns<F: FftField, const WIDTH: usize>(
    matrix: &[[F; WIDTH]],
) -> Vec<DensePolynomial<F>> {
    assert!(matrix.len().is_power_of_two());
    let domain = GeneralEvaluationDomain::new_subgroup(matrix.len()).unwrap();
    let mut polynomials = Vec::new();
    for col_idx in 0..WIDTH {
        let column = matrix.iter().map(|row| row[col_idx]).collect::<Vec<F>>();
        let evals = Evaluations::from_vec_and_domain(column, domain);
        polynomials.push(evals.interpolate());
    }
    polynomials
}
