#![feature(allocator_api, int_log)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff_optimized::fp64::Fp;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::Planner;
use fast_poly::utils::bit_reverse;
use objc::rc::autoreleasepool;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use std::time::Instant;

// Number theory transform
fn ntt_control<F: FftField>(values: &[F]) -> Vec<F> {
    let domain = Radix2EvaluationDomain::<F>::new(values.len()).unwrap();
    let mut vals = values.to_owned();
    domain.fft_in_place(&mut vals);
    vals
    // if values.len() <= 1 {
    //     values.to_vec()
    // } else {
    //     let twiddle = F::get_root_of_unity(values.len() as u64).unwrap();
    //     let half = values.len() / 2;
    //     let odd_values = values
    //         .iter()
    //         .skip(1)
    //         .step_by(2)
    //         .copied()
    //         .collect::<Vec<F>>();
    //     let even_values =
    // values.iter().step_by(2).copied().collect::<Vec<F>>();     let odds =
    // ntt_control(&odd_values);     let evens = ntt_control(&even_values);
    //     (0..values.len())
    //         .map(|i| evens[i % half] + twiddle.pow([i as u64]) * odds[i %
    // half])         .collect()
    // }
}

fn gen_pcg_input<F: Field>(n: usize) -> Vec<F, PageAlignedAllocator> {
    let mut rng = Pcg64::seed_from_u64(42); //n.try_into().unwrap());
    let mut res = Vec::new_in(PageAlignedAllocator);
    res.resize(n, F::zero());
    res.iter_mut().for_each(|v| *v = F::rand(&mut rng));
    res
}

#[test]
fn ntt_2048_vals() {
    autoreleasepool(|| {
        let n = 2048;
        let mut input = gen_pcg_input::<Fp>(n);
        let now = Instant::now();
        let mut expected = ntt_control(&input);
        bit_reverse(&mut expected);
        println!("Control fft time: {:?}", now.elapsed());
        let planner = Planner::default();

        let now = Instant::now();
        let mut ntt = planner.plan_fft(n);
        ntt.process(&mut input);
        println!("Gpu fft time: {:?}", now.elapsed());

        for (i, (a, b)) in input.iter().zip(&expected).enumerate() {
            assert_eq!(a, b, "mismatch at index {i}");
        }
    });
}

#[test]
fn ntt_524288_vals() {
    autoreleasepool(|| {
        let n = 4194304;
        let mut input = gen_pcg_input::<Fp>(n);
        let now = Instant::now();
        let mut expected = ntt_control(&input);
        bit_reverse(&mut expected);
        println!("Control fft time: {:?}", now.elapsed());
        let planner = Planner::default();

        let now = Instant::now();
        let mut ntt = planner.plan_fft(n);
        ntt.process(&mut input);
        println!("Gpu fft time: {:?}", now.elapsed());

        for (i, (a, b)) in input.iter().zip(&expected).enumerate().take(10) {
            println!("{}, {}", a, b);
        }

        for (i, (a, b)) in input.iter().zip(&expected).enumerate() {
            assert_eq!(a, b, "mismatch at index {i}");
        }
    });
}
