use crate::FftDirection;
use ark_ff::FftField;
use ark_ff::Field;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

// stolen from Winterfell library
macro_rules! batch_iter_mut {
    ($e: expr, $c: expr) => {
        #[cfg(feature = "parallel")]
        {
            let batch_size = $e.len() / rayon::current_num_threads().next_power_of_two();
            if batch_size < 1 {
                $c($e, 0);
            }
            else {
                $e.par_chunks_mut(batch_size).enumerate().for_each(|(i, batch)| {
                    $c(batch, i * batch_size);
                });
            }
        }

        #[cfg(not(feature = "parallel"))]
        $c($e, 0);
    };
    ($e: expr, $min_batch_size: expr, $c: expr) => {
        #[cfg(feature = "parallel")]
        {
            println!("feature parallel");
            let batch_size = $e.len() / rayon::current_num_threads().next_power_of_two();
            if batch_size < $min_batch_size {
                $c($e, 0);
            }
            else {
                $e.par_chunks_mut(batch_size).enumerate().for_each(|(i, batch)| {
                    $c(batch, i * batch_size);
                });
            }
        }

        #[cfg(not(feature = "parallel"))]
        $c($e, 0);
    };
}

/// Fills a slice with twiddle factors
///
/// TODO: Generate of the GPU https://kieber-emmons.medium.com/9e60b974d62
///
/// [Inverse](FftDirection::Inverse) twiddles are normalized by `1 / n`.
pub fn fill_twiddles<F: FftField>(dst: &mut [F], root: F) {
    batch_iter_mut!(dst, 1024, |batch: &mut [F], batch_offset: usize| {
        batch[0] = root.pow([batch_offset as u64]);
        for i in 1..batch.len() {
            batch[i] = batch[i - 1] * root;
        }
    });
}
