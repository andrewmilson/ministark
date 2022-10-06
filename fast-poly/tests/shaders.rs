#![feature(allocator_api, int_log)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff_optimized::fp64::Fp;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
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
fn fft() {
    autoreleasepool(|| {
        let domains = [
            Radix2EvaluationDomain::new(2048).unwrap(),
            Radix2EvaluationDomain::new(8192).unwrap(),
            Radix2EvaluationDomain::new_coset(2048, Fp::GENERATOR).unwrap(),
            Radix2EvaluationDomain::new_coset(8192, Fp::GENERATOR).unwrap(),
        ];

        for (i, domain) in domains.into_iter().enumerate() {
            let poly = DensePolynomial::<Fp>::rand(domain.size() - 1, &mut ark_std::test_rng());

            let mut evals = poly.coeffs.to_vec_in(PageAlignedAllocator);
            let mut fft = Fft::from(domain);
            fft.encode(&mut evals);
            fft.execute();

            for (j, (x, y)) in domain.elements().zip(evals).enumerate() {
                assert_eq!(poly.evaluate(&x), y, "domain ({i}) mismatch at index {j}");
            }
        }
    });
}
