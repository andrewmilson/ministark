use ark_ff::FftField;
use ark_ff::Field;
use num_traits::One;

// pub fn evaluation_difference<F>(a: &[F], b: &[F]) -> Vec<F>
// where
//     F: Field,
//     F::BasePrimeField: FftField,
// {

// }

pub fn quotient<F>(a: &[F], b: &[F]) -> Vec<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    assert_eq!(a.len(), b.len());
    assert!(a.len().is_power_of_two());
    let codeword_len = a.len();
    let offset = F::BasePrimeField::GENERATOR;
    let omega = F::BasePrimeField::get_root_of_unity(codeword_len as u64).unwrap();
    let difference_codeword = a.iter().zip(b).map(|(&x, &y)| x - y);
    // The polynomial (x - 1) evaluated over the FRI domain.
    // The intuition is that the permutation arguments should be equal in the first
    // itteration (1 == omicron^0)
    let zerofier = (0..codeword_len)
        .map(|i| offset * omega.pow([i as u64]) - F::BasePrimeField::one())
        .collect::<Vec<F::BasePrimeField>>();
    let mut zerofier_inv = zerofier;
    ark_ff::batch_inversion(&mut zerofier_inv);

    difference_codeword
        .zip(zerofier_inv)
        .map(|(d, z)| d * F::from_base_prime_field(z))
        .collect()
}
