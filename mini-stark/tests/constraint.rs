#![feature(allocator_api)]
use ark_ff::FftField;
use ark_ff::One;
use ark_ff::UniformRand;
use ark_ff::Zero;
use ark_ff_optimized::fp64::Fp;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_std::rand::Rng;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::GpuFft;
use fast_poly::GpuField;
use mini_stark::constraint::helper::are_eq;
use mini_stark::constraint::helper::is_binary;
use mini_stark::constraint::Column;
use mini_stark::Constraint;
use mini_stark::Matrix;

fn fib_matrix<F: GpuField>(n: usize) -> Matrix<F> {
    let mut columns = vec![
        Vec::with_capacity_in(n, PageAlignedAllocator),
        Vec::with_capacity_in(n, PageAlignedAllocator),
    ];
    columns[0].push(F::one());
    columns[1].push(F::one());
    for i in 1..n {
        let n0 = *columns[0].last().unwrap() + columns[1].last().unwrap();
        let n1 = n0 + columns[1].last().unwrap();
        columns[0].push(n0);
        columns[1].push(n1);
    }
    Matrix::new(columns)
}

fn binary_matrix<F: GpuField>(n: usize) -> Matrix<F> {
    let mut rng = ark_std::test_rng();
    let mut col = Vec::with_capacity_in(n, PageAlignedAllocator);
    for _ in 0..n {
        if rng.gen() {
            col.push(F::one())
        } else {
            col.push(F::zero());
        }
    }
    Matrix::new(vec![col])
}

#[test]
fn evaluate_fibonacci_constraint() {
    let n = 2048;
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = trace_domain.get_coset(Fp::GENERATOR).unwrap();
    let matrix = fib_matrix(n);
    let constraints: Vec<Constraint<Fp>> = vec![
        are_eq(0.next(), 0.curr() + 1.curr()),
        are_eq(1.next(), 0.next() + 1.curr()),
    ];

    let poly_matrix = matrix.interpolate_columns(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let mut constraint_evals = Vec::new();
    for constraint in constraints {
        constraint_evals.push(constraint.evaluate(&[], 1, &lde_matrix));
    }
    let constraint_evals_matrix = Matrix::new(constraint_evals);
    let mut constraint_matrix_poly = constraint_evals_matrix.interpolate_columns(lde_domain);
    let constraint_matrix = constraint_matrix_poly.evaluate(trace_domain);

    println!("{}, {}", &constraint_matrix[0][0], &constraint_matrix[1][0]);
    println!("{}, {}", &constraint_matrix[0][1], &constraint_matrix[1][1]);
    println!("{}, {}", &constraint_matrix[0][2], &constraint_matrix[1][2]);
    println!("{}, {}", &constraint_matrix[0][3], &constraint_matrix[1][3]);
    println!("{}, {}", &constraint_matrix[0][4], &constraint_matrix[1][4]);
    println!("{}, {}", &constraint_matrix[0][5], &constraint_matrix[1][5]);
    println!("{}, {}", &constraint_matrix[0][6], &constraint_matrix[1][6]);
    println!("{}, {}", &constraint_matrix[0][7], &constraint_matrix[1][7]);
}

#[test]
fn evaluate_binary_constraint() {
    let n = 2048;
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = trace_domain.get_coset(Fp::GENERATOR).unwrap();
    let matrix = binary_matrix(n);
    let constraints: Vec<Constraint<Fp>> = vec![is_binary(0.curr())];

    let poly_matrix = matrix.interpolate_columns(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let mut constraint_evals = Vec::new();
    for constraint in constraints {
        constraint_evals.push(constraint.evaluate(&[], 1, &lde_matrix));
    }
    let constraint_evals_matrix = Matrix::new(constraint_evals);
    let mut constraint_matrix_poly = constraint_evals_matrix.interpolate_columns(lde_domain);
    let constraint_matrix = constraint_matrix_poly.evaluate(trace_domain);

    println!("{}, {}", &constraint_matrix[0][0], &constraint_matrix[1][0]);
    println!("{}, {}", &constraint_matrix[0][1], &constraint_matrix[1][1]);
    println!("{}, {}", &constraint_matrix[0][2], &constraint_matrix[1][2]);
    println!("{}, {}", &constraint_matrix[0][3], &constraint_matrix[1][3]);
    println!("{}, {}", &constraint_matrix[0][4], &constraint_matrix[1][4]);
    println!("{}, {}", &constraint_matrix[0][5], &constraint_matrix[1][5]);
    println!("{}, {}", &constraint_matrix[0][6], &constraint_matrix[1][6]);
    println!("{}, {}", &constraint_matrix[0][7], &constraint_matrix[1][7]);
}

// #[test]
// fn evaluate_fibonacci_constraint() {
//     let n = 2048;
//     let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
//     let lde_domain = trace_domain.get_coset(Fp::GENERATOR).unwrap();
//     let mut col1 = Vec::new_in(PageAlignedAllocator);
//     // col1.resize(n, Fp::one());
//     let mut col2 = Vec::new_in(PageAlignedAllocator);
//     let three = Fp::one() + Fp::one() + Fp::one();
//     // col2.resize(n, three);
//     col1.push(Fp::one());
//     col2.push(three);
//     for _ in 1..n {
//         col1.push(*col1.last().unwrap() * three);
//         col2.push(*col2.last().unwrap() * three);
//     }
//     // println!("LES GOOOO: {}, {}", col1[1], col2[(1 + 1) % n]);
//     let matrix = Matrix::new(vec![col1, col2]);
//     let constraints: Vec<Constraint<Fp>> = vec![
//         are_eq(0.curr() * three * three, 1.next()), // are_eq(1.next(),
// 0.next() + 1.curr()),         are_eq(0.curr() * three, 1.curr()),
//     ];

//     // for constraint in &constraints {
//     //     println!("Degreeeee: {}", constraint.degree());
//     // }

//     let poly_matrix = matrix.interpolate_columns(trace_domain);
//     let lde_matrix = poly_matrix.evaluate(lde_domain);

//     // println!("v0:{}, v1:{}", matrix[0][0], matrix[0][1]);
//     let mut my_eval1 = Vec::new_in(PageAlignedAllocator);
//     let mut my_eval2 = Vec::new_in(PageAlignedAllocator);
//     for i in 0..n {
//         my_eval1.push(lde_matrix.0[0][i] * three * three - lde_matrix.0[1][(i
// + 1) % n]);         my_eval2.push(lde_matrix.0[0][i] * three -
// lde_matrix.0[1][i]);     }
//     println!("First eval:{}", my_eval1[0]);
//     println!("Second eval:{}", my_eval1[1]);
//     let my_mat = Matrix::new(vec![my_eval1, my_eval2]);
//     let my_poly = my_mat.interpolate_columns(lde_domain);
//     let my_actual_eval = my_poly.interpolate_columns(trace_domain);

//     // println!("1:{}", &my_actual_eval[0][0]);
//     // println!("2:{}", &my_actual_eval[0][1]);
//     // println!("{}", &my_actual_eval[0][2]);
//     // println!("{}", &my_actual_eval[0][3]);
//     // println!("{}", &my_actual_eval[0][4]);
//     // println!("{}", &my_actual_eval[0][5]);
//     // println!("{}", &my_actual_eval[0][6]);
//     // println!("{}", &my_actual_eval[0][7]);
//     // println!("{}", &my_actual_eval[0][2046]);
//     // println!("{}", &my_actual_eval[0][2047]);

//     let mut constraint_evals = Vec::new();
//     for constraint in constraints {
//         constraint_evals.push(constraint.evaluate(&[], 1, &lde_matrix));
//     }

//     for i in 0..n {
//         assert_eq!(constraint_evals[0][i], my_mat.0[0][i], "Mismatch at
// {i}");     }

//     println!("First eval:{}", constraint_evals[0][0]);
//     println!("Second eval:{}", constraint_evals[0][1]);
//     let constraint_evals_matrix = Matrix::new(constraint_evals);
//     let cons_mat_poly =
// constraint_evals_matrix.interpolate_columns(lde_domain);
//     let constraint_matrix = cons_mat_poly.evaluate(trace_domain);

//     println!("{}", &constraint_matrix[0][0]);
//     println!("{}", &constraint_matrix[0][1]);
//     println!("{}", &constraint_matrix[0][2]);
//     println!("{}", &constraint_matrix[0][3]);
//     println!("{}", &constraint_matrix[0][4]);
//     println!("{}", &constraint_matrix[0][5]);
//     println!("{}", &constraint_matrix[0][6]);
//     println!("{}", &constraint_matrix[0][7]);

//     // println!("{}, {}", &constraint_evals[0][0], &constraint_evals[1][0]);
//     // println!("{}, {}", &constraint_evals[0][1], &constraint_evals[1][1]);
//     // println!("{}, {}", &constraint_evals[0][2], &constraint_evals[1][2]);
//     // println!("{}, {}", &constraint_evals[0][3], &constraint_evals[1][3]);
//     // println!("{}, {}", &constraint_evals[0][4], &constraint_evals[1][4]);
//     // println!("{}, {}", &constraint_evals[0][5], &constraint_evals[1][5]);
//     // println!("{}, {}", &constraint_evals[0][6], &constraint_evals[1][6]);
//     // println!("{}, {}", &constraint_evals[0][7], &constraint_evals[1][7]);
// }
