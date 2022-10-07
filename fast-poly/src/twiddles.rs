use crate::FftDirection;
use ark_ff::FftField;
use ark_ff::Field;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

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
            println!("feature Serr");
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
pub fn fill_twiddles<F: FftField>(dst: &mut [F], n: usize, direction: FftDirection) {
    let n = n as u64;
    let mut root = F::get_root_of_unity(n).unwrap();
    if direction == FftDirection::Inverse {
        root.inverse_in_place();
    }

    let norm_factor = if direction == FftDirection::Inverse {
        F::from_base_prime_field(F::BasePrimeField::from(n).inverse().unwrap())
    } else {
        F::one()
    };

    batch_iter_mut!(dst, 1024, |batch: &mut [F], batch_offset: usize| {
        println!("Doint batch {batch_offset}");
        batch[0] = norm_factor * root.pow([batch_offset as u64]);
        for i in 1..batch.len() {
            batch[i] = batch[i - 1] * root;
        }
    });

    // // ark_std::cfg_iter_mut!(dst)
    // for v in dst.iter_mut() {
    //     // Using this method creates a data dependency over each iteration
    //     // preventing parallelization so is actually slower.
    //     *v = norm_factor;
    //     norm_factor *= root;
    // }
}

// /// Stolen from Winterfell
// #[cfg(not(feature = "parallel"))]
// fn prefix_scan(&mut) {
//     batch_iter_mut!(&mut result, 1024, |batch: &mut [E], batch_offset: usize|
// {         let start = b.exp((batch_offset as u64).into());
//         fill_power_series(batch, b, start);
//     });
//     result
// }

// #[cfg(feature = "parallel")]
// fn compute_powers() {}
