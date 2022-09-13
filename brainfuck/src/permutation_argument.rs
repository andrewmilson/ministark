use algebra::batch_inverse;
use algebra::Felt;
use algebra::StarkFelt;
use num_traits::One;

pub fn quotient<E>(a: &[E], b: &[E]) -> Vec<E>
where
    E: Felt,
    E::BaseFelt: StarkFelt,
{
    assert_eq!(a.len(), b.len());
    assert!(a.len().is_power_of_two());
    let codeword_len = a.len();
    let offset = E::BaseFelt::GENERATOR;
    let omega = E::BaseFelt::get_root_of_unity(codeword_len.ilog2());
    let difference_codeword = b.iter().zip(b).map(|(&x, &y)| x - y);
    // The polynomial (x - 1) evaluated over the FRI domain.
    let zerofier = (0..codeword_len)
        .map(|i| offset * omega.pow(&[i as u64]) - E::BaseFelt::one())
        .collect::<Vec<E::BaseFelt>>();
    let zerofier_inv = batch_inverse(&zerofier);

    difference_codeword
        .zip(zerofier_inv)
        .map(|(d, z)| d * E::from(z.unwrap()))
        .collect()
}
