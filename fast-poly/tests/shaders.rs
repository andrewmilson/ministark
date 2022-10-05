#![feature(allocator_api, int_log)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff_optimized::fp64::Fp;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::EvaluationDomain;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::Fft;
use fast_poly::plan::Planner;
use fast_poly::plan::PLANNER;
use fast_poly::utils::bit_reverse;
use objc::rc::autoreleasepool;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use std::time::Instant;

fn gen_pcg_input<F: Field>(n: usize) -> Vec<F, PageAlignedAllocator> {
    let mut rng = Pcg64::seed_from_u64(42); //n.try_into().unwrap());
    let mut res = Vec::new_in(PageAlignedAllocator);
    res.resize(n, F::zero());
    res.iter_mut().for_each(|v| *v = F::rand(&mut rng));
    res
}

#[test]
fn fft_2048_vals() {
    autoreleasepool(|| {
        let domain = Radix2EvaluationDomain::new(2048).unwrap();
        let mut input = gen_pcg_input::<Fp>(domain.size());
        let now = Instant::now();
        let mut expected = domain.fft(&input);
        println!("Control fft time: {:?}", now.elapsed());

        let now = Instant::now();
        let mut fft = Fft::from(domain);
        fft.process(&mut input);
        println!("Gpu fft time: {:?}", now.elapsed());

        for (i, (a, b)) in input.iter().zip(&expected).enumerate() {
            assert_eq!(a, b, "mismatch at index {i}");
        }
    });
}

#[test]
fn fft_524288_vals() {
    autoreleasepool(|| {
        let n = 33554432;
        let domain = Radix2EvaluationDomain::new(n).unwrap();
        let mut input = gen_pcg_input::<Fp>(domain.size());
        let mut expected = input.to_vec();
        let now = Instant::now();
        domain.fft_in_place(&mut expected);
        println!("Control fft time: {:?}", now.elapsed());

        let now = Instant::now();
        let mut fft = Fft::from(domain);
        fft.process(&mut input);
        println!("Gpu fft time: {:?}", now.elapsed());

        // for (i, (a, b)) in input.iter().zip(&expected).enumerate().take(10) {
        //     println!("{}, {}", a, b);
        // }

        // for (i, (a, b)) in input.iter().zip(&expected).enumerate() {
        //     assert_eq!(a, b, "mismatch at index {i}");
        // }
    });
}
