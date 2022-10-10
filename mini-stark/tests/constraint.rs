#![feature(allocator_api)]
use ark_ff::FftField;
use ark_ff::One;
use ark_ff::Zero;
use ark_ff_optimized::fp64::Fp;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use ark_std::rand::Rng;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::GpuField;
use mini_stark::constraint::are_eq;
use mini_stark::constraint::is_binary;
use mini_stark::constraint::Challenge;
use mini_stark::constraint::Column;
use mini_stark::Constraint;
use mini_stark::Matrix;

enum TestChallenge {
    Alpha,
}

impl Challenge for TestChallenge {
    fn index(&self) -> usize {
        *self as usize
    }
}

#[test]
fn constraint_with_challenges() {
    let constraint = (TestChallenge::Alpha.get_challenge() - 0.curr()) * 1.curr();

    assert!(constraint
        .evaluate(&[Fp::one()], &[Fp::one(), Fp::from(100)], &[])
        .is_zero());
}

#[test]
fn constraint_multiplication() {
    let zero = Fp::zero();
    let one = Fp::one();
    let two = one + one;
    let three = two + one;
    let four = three + one;
    let five = four + one;
    let six = five + one;
    let seven = six + one;
    let eight = seven + one;
    let nine = eight + one;
    let ten = nine + one;
    let eleven = ten + one;
    let twelve = eleven + one;

    // checks the column values are between 0 and 10
    let is_between_0_and_10: Constraint<Fp> = (0.curr() - one)
        * (0.curr() - two)
        * (0.curr() - three)
        * (0.curr() - four)
        * (0.curr() - five)
        * (0.curr() - six)
        * (0.curr() - seven)
        * (0.curr() - eight)
        * (0.curr() - nine);

    assert!(!is_between_0_and_10.evaluate(&[], &[-two], &[]).is_zero());
    assert!(!is_between_0_and_10.evaluate(&[], &[-one], &[]).is_zero());
    assert!(!is_between_0_and_10.evaluate(&[], &[zero], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[one], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[two], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[three], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[four], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[five], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[six], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[seven], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[eight], &[]).is_zero());
    assert!(is_between_0_and_10.evaluate(&[], &[nine], &[]).is_zero());
    assert!(!is_between_0_and_10.evaluate(&[], &[ten], &[]).is_zero());
    assert!(!is_between_0_and_10.evaluate(&[], &[eleven], &[]).is_zero());
    assert!(!is_between_0_and_10.evaluate(&[], &[twelve], &[]).is_zero());
}

#[test]
fn evaluate_fibonacci_constraint() {
    let n = 2048;
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = trace_domain.get_coset(Fp::GENERATOR).unwrap();
    let matrix = gen_fib_matrix(n);
    let poly_matrix = matrix.interpolate_columns(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let constraints: Vec<Constraint<Fp>> = vec![
        are_eq(0.next(), 0.curr() + 1.curr()),
        are_eq(1.next(), 0.next() + 1.curr()),
    ];

    let constraint_evals = constraints
        .into_iter()
        .map(|constraint| constraint.evaluate_symbolic(&[], 1, &lde_matrix))
        .collect();

    let constraint_evals_matrix = Matrix::new(constraint_evals);
    let constraint_evals_poly = constraint_evals_matrix.interpolate_columns(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_evals_poly);
}

#[test]
fn evaluate_binary_constraint() {
    let n = 2048;
    let constraint: Constraint<Fp> = is_binary(0.curr());
    let blowup = constraint.degree();
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let matrix = gen_binary_valued_matrix(n);
    let poly_matrix = matrix.interpolate_columns(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);

    let constraint_eval = constraint.evaluate_symbolic(&[], 1, &lde_matrix);

    let constraint_eval_matrix = Matrix::new(vec![constraint_eval]);
    let constraint_eval_poly = constraint_eval_matrix.interpolate_columns(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_eval_poly);
}

/// Generates a matrix of fibbonacci sequence across two columns i.e.
/// ┌───────┬───────┐
/// │ Col 0 | Col 1 │
/// ├───────┼───────┤
/// │ 1     │ 1     │ #1 -> #2 ->
/// ├───────┼───────┤
/// │ 2     │ 3     │ #3 -> #4 ->
/// ├───────┼───────┤
/// │ 5     │ 8     │ #5 -> #6 ->
/// ├───────┼───────┤
/// │ ...   │ ...   │ ...
/// └───────┴───────┘
fn gen_fib_matrix<F: GpuField>(n: usize) -> Matrix<F> {
    let mut columns = vec![
        Vec::with_capacity_in(n, PageAlignedAllocator),
        Vec::with_capacity_in(n, PageAlignedAllocator),
    ];
    columns[0].push(F::one());
    columns[1].push(F::one());
    for _ in 1..n {
        let n0 = *columns[0].last().unwrap() + columns[1].last().unwrap();
        let n1 = n0 + columns[1].last().unwrap();
        columns[0].push(n0);
        columns[1].push(n1);
    }
    Matrix::new(columns)
}

/// Generates a single column matrix of zero and one values
/// ┌───────┐
/// │ Col 0 │
/// ├───────┤
/// │ 0     │
/// ├───────┤
/// │ 1     │
/// ├───────┤
/// │ 0     │
/// ├───────┤
/// │ 0     │
/// ├───────┤
/// │ 1     │
/// ├───────┤
/// │ ...   │
/// └───────┘
fn gen_binary_valued_matrix<F: GpuField>(n: usize) -> Matrix<F> {
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

fn assert_valid_over_transition_domain<F: GpuField>(
    domain: Radix2EvaluationDomain<F>,
    poly_matrix: Matrix<F>,
) {
    let mut x_values = domain.elements().collect::<Vec<F>>();
    // transition constraints apply to all rows except the last.
    x_values.pop();
    for (i, column) in poly_matrix.iter().enumerate() {
        let poly = DensePolynomial::from_coefficients_slice(column);
        for (j, x) in x_values.iter().enumerate() {
            let y = poly.evaluate(x);
            assert!(y.is_zero(), "polynomial {i} invalid at index {j}");
        }
    }
}
