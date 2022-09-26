use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_poly::Evaluations;

pub fn quotient<F: FftField, D: EvaluationDomain<F>>(
    a: &Evaluations<F, D>,
    b: &Evaluations<F, D>,
) -> Evaluations<F, D> {
    let domain = a.domain();
    assert_eq!(domain, b.domain());
    let difference = a - b;
    // The polynomial (x - 1) evaluated over the FRI domain.
    // Essentially a boundary consition stating that the the permutation arguments
    // should be equal in the first itteration
    let mut boundary_zerofier_inv =
        Evaluations::from_vec_and_domain(domain.elements().map(|x| x - F::one()).collect(), domain);
    ark_ff::batch_inversion(&mut boundary_zerofier_inv.evals);
    &difference * &boundary_zerofier_inv
}
