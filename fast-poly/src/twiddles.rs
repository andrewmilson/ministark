use crate::NttDirection;
use ark_ff::FftField;
use ark_ff::Field;

/// Fills a slice with twiddle factors
///
/// TODO: Generate of the GPU https://kieber-emmons.medium.com/9e60b974d62
///
/// [Inverse](NttDirection::Inverse) twiddles are normalized by `1 / n`.
pub fn fill_twiddles<F: FftField>(dst: &mut [F], n: usize, direction: NttDirection) {
    assert!(n.is_power_of_two(), "must be a power of two");
    let n = n as u64;
    let twiddle = F::get_root_of_unity(n).unwrap();
    let norm_factor = match direction {
        NttDirection::Forward => F::one(),
        NttDirection::Inverse => {
            F::from_base_prime_field(F::BasePrimeField::from(n).inverse().unwrap())
        }
    };
    let mut accumulator = F::one();
    for (i, v) in dst.iter_mut().enumerate() {
        // Using this method creates a data dependency over each iteration
        // preventing parallelization so is actually slower.
        *v = norm_factor * accumulator;
        accumulator *= twiddle;
    }
}
