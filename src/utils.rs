use ark_ff::FftField;
use ark_ff::Field;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Add;
use std::time::Instant;

pub struct Timer<'a> {
    name: &'a str,
    start: Instant,
}

impl<'a> Timer<'a> {
    pub fn new(name: &'a str) -> Timer<'a> {
        let start = Instant::now();
        Timer { name, start }
    }
}

impl<'a> Drop for Timer<'a> {
    fn drop(&mut self) {
        println!("{} in {:?}", self.name, self.start.elapsed());
    }
}

pub fn interleave<T: Copy + Send + Sync + Default, const RADIX: usize>(
    source: &[T],
) -> Vec<[T; RADIX]> {
    let n = source.len() / RADIX;
    let mut res = vec![[T::default(); RADIX]; n];
    ark_std::cfg_iter_mut!(res)
        .enumerate()
        .for_each(|(i, element)| {
            for j in 0..RADIX {
                element[j] = source[i + j * n]
            }
        });
    res
}

// pub(crate) fn print_row<F: GpuField>(row: &[F]) {
//     for val in row {
//         print!("{val}, ");
//     }
//     println!()
// }

/// Rounds the input value up the the nearest power of two
pub fn ceil_power_of_two(value: usize) -> usize {
    if value.is_power_of_two() {
        value
    } else {
        value.next_power_of_two()
    }
}

// Evaluates the vanishing polynomial for `vanish_domain` over `eval_domain`
// E.g. evaluates `(x - v_0)(x - v_1)...(x - v_n-1)` over `eval_domain`
pub fn fill_vanishing_polynomial<F: FftField>(
    dst: &mut [F],
    vanish_domain: &Radix2EvaluationDomain<F>,
    eval_domain: &Radix2EvaluationDomain<F>,
) {
    let n = vanish_domain.size();
    let scaled_eval_offset = eval_domain.coset_offset().pow([n as u64]);
    let scaled_eval_generator = eval_domain.group_gen().pow([n as u64]);
    let scaled_vanish_offset = vanish_domain.coset_offset_pow_size();

    #[cfg(feature = "parallel")]
    let chunk_size = std::cmp::max(n / rayon::current_num_threads(), 1024);
    #[cfg(not(feature = "parallel"))]
    let chunk_size = n;

    ark_std::cfg_chunks_mut!(dst, chunk_size)
        .enumerate()
        .for_each(|(i, chunk)| {
            let mut acc = scaled_eval_offset * scaled_eval_generator.pow([(i * chunk_size) as u64]);
            chunk.iter_mut().for_each(|coeff| {
                *coeff = acc - scaled_vanish_offset;
                acc *= &scaled_eval_generator
            })
        });
}

// taken from arkworks-rs
/// Horner's method for polynomial evaluation
#[inline]
pub fn horner_evaluate<F: Field, T: Field>(poly_coeffs: &[F], point: &T) -> T
where
    T: for<'a> Add<&'a F, Output = T>,
{
    poly_coeffs
        .iter()
        .rfold(T::zero(), move |result, coeff| result * point + coeff)
}
