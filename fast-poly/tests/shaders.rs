#![feature(allocator_api, int_log)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::UniformRand;
use ark_ff_optimized::fp64::Fp;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Evaluations;
use ark_poly::Polynomial;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::Fft;
use fast_poly::plan::Ifft;
use objc::rc::autoreleasepool;
use rand::SeedableRng;
use rand_pcg::Pcg64;

fn gen_pcg_input<F: Field>(n: usize) -> Vec<F, PageAlignedAllocator> {
    let mut rng = Pcg64::seed_from_u64(42); //n.try_into().unwrap());
    let mut res = Vec::new_in(PageAlignedAllocator);
    res.resize(n, F::zero());
    res.iter_mut().for_each(|v| *v = F::rand(&mut rng));
    res
}

#[test]
fn reg_fft() {
    autoreleasepool(|| {
        let domains = [
            Radix2EvaluationDomain::new(2048).unwrap(),
            Radix2EvaluationDomain::new(4096).unwrap(),
            Radix2EvaluationDomain::new_coset(2048, Fp::GENERATOR).unwrap(),
            Radix2EvaluationDomain::new_coset(4096, Fp::GENERATOR).unwrap(),
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
            // let mut rng = ark_std::test_rng();
            // let evals = (0..domain.size())
            //     .map(|_| Fp::rand(&mut rng))
            //     .collect::<Vec<Fp>>();
            // let expected_coeffs = domain.ifft(&evals);
            let poly = DensePolynomial::<Fp>::rand(domain.size() - 1, &mut ark_std::test_rng());
            let evals = poly.evaluate_over_domain_by_ref(domain).evals;

            let mut coeffs = evals.to_vec_in(PageAlignedAllocator);
            let mut ifft = Ifft::from(domain);
            ifft.encode(&mut coeffs);
            ifft.execute();

            let domain_size_inv = Fp::one(); //domain.size_inv;

            for (j, (expected, actual)) in poly.coeffs.into_iter().zip(coeffs).enumerate() {
                assert_eq!(
                    expected,
                    actual * domain_size_inv,
                    "domain ({i}) mismatch at index {j}"
                );
            }
        }
    });
}
