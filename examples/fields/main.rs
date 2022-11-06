#![feature(allocator_api)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::UniformRand;
use ark_ff_optimized::fp64::Fp;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use gpu_poly::prelude::PageAlignedAllocator;
use ministark::Matrix;

mod fq3;
use fq3::Fq3;

fn main() {
    // let mut my_elem = Fq3::from_base_prime_field_elems(&[Fp::one(), Fp::one(),
    // Fp::one()]).unwrap(); let mut extension_element = Fq3::one();
    // let generator = Fp::TWO_ADIC_ROOT_OF_UNITY;

    let n = 2048;
    let mut rng = ark_std::test_rng();
    let column = (0..n)
        .map(|_| Fq3::rand(&mut rng))
        .collect::<Vec<Fq3>>()
        .to_vec_in(PageAlignedAllocator);

    let my_matrix = Matrix::new(vec![column.to_vec_in(PageAlignedAllocator)]);

    let domain = Radix2EvaluationDomain::new(n).unwrap();
    let polys = my_matrix.interpolate(domain);
    let evals = polys.evaluate(domain);

    assert_eq!(polys.0[0], column);

    // println!("one: {}", generator);
    // println!("generator: {}", generator);
    // extension_element *= generator;
    // println!("mul: {}", extension_element.inverse().unwrap());
    // println!("mul: {}", extension_element.pow([(1 << 32) - 1]));
    // my_elem *= generator;
    // println!("test: {:?}", (my_elem * my_elem).inverse());
}
