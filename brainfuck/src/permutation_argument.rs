use ark_ff::FftField;
use ark_ff::Field;
use legacy_algebra::batch_inverse;
use legacy_algebra::Felt;
use legacy_algebra::StarkFelt;
use num_traits::One;

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
    let difference_codeword = b.iter().zip(b).map(|(&x, &y)| x - y);
    // The polynomial (x - 1) evaluated over the FRI domain.
    let zerofier = (0..codeword_len)
        .map(|i| offset * omega.pow(&[i as u64]) - F::BasePrimeField::one())
        .collect::<Vec<F::BasePrimeField>>();
    let mut zerofier_inv = zerofier;
    ark_ff::batch_inversion(&mut zerofier_inv);

    difference_codeword
        .zip(zerofier_inv)
        .map(|(d, z)| d * F::from_base_prime_field(z))
        .collect()
}
