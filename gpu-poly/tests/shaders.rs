#![cfg(target_arch = "aarch64")]
#![feature(allocator_api, int_log)]

use ark_ff::FftField;
use ark_ff_optimized::fp64::Fp;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use gpu_poly::prelude::*;
use objc::rc::autoreleasepool;

#[test]
fn fft() {
    autoreleasepool(|| {
        let domains = [
            Radix2EvaluationDomain::new(2048).unwrap(),
            Radix2EvaluationDomain::new(4096).unwrap(),
            Radix2EvaluationDomain::new_coset(2048, Fp::GENERATOR).unwrap(),
            Radix2EvaluationDomain::new_coset(4096, Fp::GENERATOR).unwrap(),
        ];

        for (i, domain) in domains.into_iter().enumerate() {
            let poly = DensePolynomial::rand(domain.size() - 1, &mut ark_std::test_rng());
            let mut evals = poly.coeffs.to_vec_in(PageAlignedAllocator);
            let mut fft = GpuFft::from(domain);
            fft.encode(&mut evals);
            fft.execute();

            for (j, (x, y)) in domain.elements().zip(evals).enumerate() {
                assert_eq!(poly.evaluate(&x), y, "domain ({i}) mismatch at index {j}");
            }
        }
    });
}

#[test]
fn ifft() {
    autoreleasepool(|| {
        let domains = [
            Radix2EvaluationDomain::new(2048).unwrap(),
            Radix2EvaluationDomain::new(4096).unwrap(),
            Radix2EvaluationDomain::new_coset(2048, Fp::GENERATOR).unwrap(),
            Radix2EvaluationDomain::new_coset(4096, Fp::GENERATOR).unwrap(),
        ];

        for (i, domain) in domains.into_iter().enumerate() {
            let poly = DensePolynomial::rand(domain.size() - 1, &mut ark_std::test_rng());
            let evals = poly.evaluate_over_domain_by_ref(domain).evals;

            let mut coeffs = evals.to_vec_in(PageAlignedAllocator);
            let mut ifft = GpuIfft::from(domain);
            ifft.encode(&mut coeffs);
            ifft.execute();

            for (j, (expected, actual)) in poly.coeffs.into_iter().zip(coeffs).enumerate() {
                assert_eq!(expected, actual, "domain ({i}) mismatch at index {j}");
            }
        }
    });
}
