#![feature(array_chunks)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff_optimized::fp64::Fp;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use ark_std::UniformRand;
use mini_stark::fri::apply_drp;
use mini_stark::utils::interleave;

#[test]
fn drp() {
    // let n = 4;
    // let mut rng = ark_std::test_rng();
    // let evals = (0..n).map(|_| Fp::rand(&mut rng)).collect::<Vec<Fp>>();

    let alpha = Fp::from(123);
    let mut evals = vec![
        Fp::from(43),
        Fp::from(2),
        Fp::from(908),
        Fp::from(3242),
        Fp::from(12),
        Fp::from(44),
        Fp::from(94),
        Fp::from(94),
    ];
    let n = evals.len();
    let coset = Radix2EvaluationDomain::<Fp>::new_coset(n, Fp::GENERATOR).unwrap();
    let poly = coset.ifft(&evals);
    let (evens, odds): (Vec<Fp>, Vec<Fp>) = poly.array_chunks::<2>().map(|[a, b]| (a, b)).unzip();
    let half_coset =
        Radix2EvaluationDomain::<Fp>::new_coset(n / 2, Fp::GENERATOR.square()).unwrap();

    let evens_poly = DensePolynomial::from_coefficients_slice(&evens);
    let odds_poly = DensePolynomial::from_coefficients_slice(&odds);

    let new_coeffs = evens
        .into_iter()
        .zip(odds)
        .map(|(e, o)| e + alpha * o)
        .collect::<Vec<Fp>>();
    let new_poly = DensePolynomial::from_coefficients_slice(&new_coeffs);
    prety_print(&half_coset.fft(&new_coeffs));

    println!("HERE WE GO");
    // let domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    for x in coset.elements() {
        let x2 = x.square();
        println!(
            "{}",
            evens_poly.evaluate(&x2) + alpha * odds_poly.evaluate(&x2)
        );
        println!("dedup:{}", new_poly.evaluate(&x2))
    }

    // let new_poly = evens
    //     .into_iter()
    //     .zip(odds)
    //     .map(|(e, o)| e + alpha * o)
    //     .collect::<Vec<Fp>>();

    // prety_print(&half_coset.fft(&new_poly));

    // prety_print(
    //     &,
    // );

    // std::iter::zip(half_coset.fft(evens), half_coset.fft(odds))

    // // let interleaved_evals = interleave::<Fp, 2>(&evals);

    // // let drp1 = apply_drp(&interleaved_evals, Fp::GENERATOR, alpha);

    // // prety_print(&drp1);
    // let root = Fp::get_root_of_unity(evals.len() as u64).unwrap();

    // // let mut drp2 = evals.clone();
    // // let poly = coset.ifft(&evals);
    // let half_coset = Radix2EvaluationDomain::new_coset(4,
    // Fp::GENERATOR).unwrap(); half_coset.ifft_in_place(&mut evals);
    // let (mut evens, mut odds): (Vec<Fp>, Vec<Fp>) =
    //     evals.array_chunks::<2>().map(|[a, b]| (a, b)).unzip();
    // // half_coset
    // //     .get_coset(Fp::GENERATOR.square())
    // //     .unwrap()
    // //     .ifft_in_place(&mut odds);

    // prety_print(&evens);
    // println!("odds");
    // prety_print(&odds);

    // println!("Yest");

    // prety_print(
    //     &half_coset
    //         .fft(&evens)
    //         .into_iter()
    //         .zip(half_coset.fft(&odds))
    //         .zip(half_coset.elements())
    //         .map(|((e, o), x)| [e + x * o, e - x * o])
    //         .flatten()
    //         .collect::<Vec<_>>(),
    // );

    // let even_evals = half_coset.fft(&evens);
    // let odd_evals = half_coset.fft(&odds);

    // let combined = even_evals
    // .into_iter()
    // .zip(odd_evals)
    // .map(|(e, o)| e + alpha * o)
    // .collect::<Vec<_>>();

    // let mut combined = evens
    //     .into_iter()
    //     .zip(odds)
    //     .map(|(e, o)| e + alpha * o)
    //     .collect::<Vec<_>>();

    // let half_coset = Radix2EvaluationDomain::new_coset(2,
    // Fp::GENERATOR).unwrap(); half_coset.fft_in_place(&mut combined);

    // println!("FINITO");
    // prety_print(&combined);

    // coset.ifft_in_place(&mut drp2);
    // // let (mut even, mut odd): (Vec<Fp>, Vec<Fp>) = drp2
    // //     .array_chunks::<2>()
    // //     .copied()
    // //     .map(|[v1, v2]| (v1, v2))
    // //     .unzip();
    // let (mut even, mut odd): (Vec<Fp>, Vec<Fp>) = interleave::<Fp, 2>(&drp2)
    //     .into_iter()
    //     .map(|[v1, v2]| (v1, v2))
    //     .unzip();

    // let half_coset = Radix2EvaluationDomain::new_coset(n / 2,
    // Fp::GENERATOR).unwrap(); half_coset.fft_in_place(&mut even);
    // half_coset.fft_in_place(&mut odd);

    // println!("even:{:?}", even);
    // println!("odd:{:?}", odd);

    // let mut drp2 = even
    //     .into_iter()
    //     .zip(odd)
    //     .map(|(e, o)| e + alpha * o)
    //     .collect::<Vec<Fp>>();
    // // half_coset.fft_in_place(&mut drp2);

    // println!("drp1: {:?}", drp1);
    // println!("drp2: {:?}", drp2);

    // for (i, (v1, v2)) in drp1.into_iter().zip(drp2).enumerate() {
    //     assert_eq!(v1, v2, "mismatch at index {i}");
    // }
}

fn prety_print(vals: &[Fp]) {
    for val in vals {
        println!("{val},");
    }
}
