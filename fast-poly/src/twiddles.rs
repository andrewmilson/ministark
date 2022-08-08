use crate::fields::StarkFelt;
use crate::NttDirection;

/// Fills a slice with twiddle factors
///
/// TODO: Generate of the GPU https://kieber-emmons.medium.com/9e60b974d62
///
/// [Inverse](NttDirection::Inverse) twiddles are normalized by `1 / n`.
pub fn fill_twiddles<E: StarkFelt>(dst: &mut [E], n: usize, direction: NttDirection) {
    assert!(n.is_power_of_two(), "must be a power of two");
    let twiddle = E::get_root_of_unity(n.log2());
    let norm_factor = match direction {
        NttDirection::Forward => E::one(),
        NttDirection::Inverse => E::from(n).inverse().unwrap(),
    };
    for (i, v) in dst.iter_mut().enumerate() {
        *v = norm_factor * twiddle.pow(i.try_into().ok().expect("can't calculate power"));
    }
}
