#![feature(allocator_api)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff_optimized::fp64::Fp;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use fast_poly::allocator::PageAlignedAllocator;
use mini_stark::Matrix;

#[test]
fn degree_adjustor() {
    println!("{:?}", Fp::one());

    let n = 2048;
    let ce_domain = Radix2EvaluationDomain::new_coset(n, Fp::GENERATOR).unwrap();
    let lde_domain = Radix2EvaluationDomain::new_coset(n * 4, Fp::GENERATOR).unwrap();

    let ce_elements = ce_domain.elements().collect::<Vec<Fp>>();
    let ce_matrix = Matrix::new(vec![ce_elements.to_vec_in(PageAlignedAllocator)]);

    let ce_poly = ce_matrix.interpolate_columns(ce_domain);
    let mut lde_evals = ce_poly.evaluate(lde_domain);

    let actual_poly = DensePolynomial::from_coefficients_slice(&ce_poly.0[0]);
    println!("{:?}", actual_poly.degree());
    println!("HERE::{:?}::GONE", actual_poly);

    println!("first: {}", ce_elements[0],);
    println!("first: {}, {}", lde_evals.0[0][0], lde_domain.element(0));

    println!("second: {}", ce_elements[1]);
    println!("second: {}, {}", lde_evals.0[0][1], lde_domain.element(1));

    println!("third: {}", ce_elements[2]);
    println!("third: {}, {}", lde_evals.0[0][2], lde_domain.element(2));

    for eval in &mut lde_evals.0[0] {
        eval.square_in_place();
    }

    println!("first: {}", lde_evals.0[0][0]);
    println!("second: {}", lde_evals.0[0][1]);
    println!("third: {}", lde_evals.0[0][2]);

    let new_poly = lde_evals.interpolate_columns(lde_domain);
    let actual_poly = DensePolynomial::from_coefficients_slice(&new_poly.0[0]);
    println!("Degree:::{:?}", actual_poly);
    println!("first: {}", new_poly.0[0][0]);
    println!("second: {}", new_poly.0[0][1]);
    println!("third: {}", new_poly.0[0][2]);

    // let mut eval_col = Vec::with_capacity_in(ce_domain.size(),
    // PageAlignedAllocator); eval_col.resize(ce_domain.size(),
    // A::Fp::zero()); fill_twiddles(&mut eval_col, ce_domain.group_gen);
    // let eval_matrix = Matrix::new(vec![eval_col]);
    // let poly_matrix = eval_matrix.interpolate_columns(ce_domain);
    // poly_matrix.evaluate(self.air.lde_domain())
}
