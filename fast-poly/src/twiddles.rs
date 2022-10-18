use ark_ff::FftField;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Fills a slice with twiddle factors
///
/// TODO: Generate of the GPU https://kieber-emmons.medium.com/9e60b974d62
///
/// [Inverse](FftDirection::Inverse) twiddles are normalized by `1 / n`.
pub fn fill_twiddles<F: FftField>(dst: &mut [F], root: F) {
    #[cfg(not(feature = "parallel"))]
    let chunk_size = dst.len();
    #[cfg(feature = "parallel")]
    let chunk_size = std::cmp::max(
        dst.len() / rayon::current_num_threads().next_power_of_two(),
        1024,
    );

    ark_std::cfg_chunks_mut!(dst, chunk_size)
        .enumerate()
        .for_each(|(chunk_offset, chunk)| {
            chunk[0] = root.pow([(chunk_offset * chunk_size) as u64]);
            for i in 1..chunk.len() {
                chunk[i] = chunk[i - 1] * root;
            }
        });
}
