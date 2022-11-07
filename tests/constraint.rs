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
use ark_std::rand::seq::SliceRandom;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use gpu_poly::allocator::PageAlignedAllocator;
use gpu_poly::GpuField;
use ministark::constraint::are_eq;
use ministark::constraint::is_binary;
use ministark::constraint::Challenge;
use ministark::constraint::Column;
use ministark::matrix::GroupItem;
use ministark::matrix::MatrixGroup;
use ministark::Constraint;
use ministark::Matrix;

#[test]
fn constraint_with_challenges() {
    let constraint = (0.get_challenge() - 0.curr()) * 1.curr();
    let challenges = [Fp::one()];
    let col_values = [Fp::one(), Fp::from(100)];

    assert!(constraint.evaluate(&challenges, &col_values, &[]).is_zero());
}

#[test]
fn symbolic_evaluation_with_challenges() {
    let n = 2048;
    let constraint = (0.curr() - 0.get_challenge()) * (0.curr() - 1.get_challenge());
    let blowup = constraint.degree();
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let alpha = Fp::from(3);
    let beta = Fp::from(7);
    let matrix = gen_binary_valued_matrix(n, alpha, beta);
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let matrix_group = MatrixGroup::new(vec![GroupItem::Fp(&lde_matrix)]);

    let constraint_eval = matrix_group.evaluate_symbolic(&[constraint], &[alpha, beta], blowup);

    let constraint_eval_poly = constraint_eval.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_eval_poly);
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
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let matrix_group = MatrixGroup::new(vec![GroupItem::Fp(&lde_matrix)]);
    let constraints: Vec<Constraint<Fp>> = vec![
        are_eq(0.next(), 0.curr() + 1.curr()),
        are_eq(1.next(), 0.next() + 1.curr()),
    ];

    let constraint_evals = matrix_group.evaluate_symbolic(&constraints, &[], 1);

    let constraint_evals_poly = constraint_evals.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_evals_poly);
}

#[test]
fn evaluate_binary_constraint() {
    let n = 2048;
    let constraint: Constraint<Fp> = is_binary(0.curr());
    let blowup = constraint.degree();
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let matrix = gen_binary_valued_matrix(n, Fp::zero(), Fp::one());
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let matrix_group = MatrixGroup::new(vec![GroupItem::Fp(&lde_matrix)]);

    let constraint_eval = matrix_group.evaluate_symbolic(&[constraint], &[], blowup);

    let constraint_eval_poly = constraint_eval.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_eval_poly);
}

#[test]
fn evaluate_permutation_constraint() {
    let n = 2048;
    let mut rng = ark_std::test_rng();
    let original_col = (0..n).map(|_| Fp::rand(&mut rng)).collect::<Vec<Fp>>();
    let mut shuffled_col = original_col.clone();
    shuffled_col.shuffle(&mut rng);
    let challenge = Fp::rand(&mut rng); // verifier challenge
    let original_product = original_col
        .iter()
        .scan(Fp::one(), |product, v| {
            let ret = *product;
            *product *= challenge - v;
            Some(ret)
        })
        .collect::<Vec<Fp>>();
    let shuffled_product = shuffled_col
        .iter()
        .scan(Fp::one(), |product, v| {
            let ret = *product;
            *product *= challenge - v;
            Some(ret)
        })
        .collect::<Vec<Fp>>();
    let matrix = Matrix::new(vec![
        original_col.to_vec_in(PageAlignedAllocator),
        shuffled_col.to_vec_in(PageAlignedAllocator),
        original_product.to_vec_in(PageAlignedAllocator),
        shuffled_product.to_vec_in(PageAlignedAllocator),
    ]);
    let alpha = 0; // first verifier challenge
    let original_col = 0;
    let shuffled_col = 1;
    let original_product = 2;
    let shuffled_product = 3;
    let constraints = vec![
        original_product.curr() * (alpha.get_challenge() - original_col.curr())
            - original_product.next(),
        shuffled_product.curr() * (alpha.get_challenge() - shuffled_col.curr())
            - shuffled_product.next(),
    ];
    let blowup = 2;
    let trace_domain = Radix2EvaluationDomain::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let matrix_group = MatrixGroup::new(vec![GroupItem::Fp(&lde_matrix)]);

    let constraint_evals = matrix_group.evaluate_symbolic(&constraints, &[challenge], blowup);

    let last_original_val = matrix.0[original_col].last().unwrap();
    let last_shuffled_val = matrix.0[shuffled_col].last().unwrap();
    let last_original_product = matrix.0[original_product].last().unwrap();
    let last_shuffled_product = matrix.0[shuffled_product].last().unwrap();
    let final_original_product = *last_original_product * (challenge - last_original_val);
    let final_shuffled_product = *last_shuffled_product * (challenge - last_shuffled_val);
    assert_eq!(final_original_product, final_shuffled_product);
    let constraint_eval_poly = constraint_evals.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_eval_poly);
}

#[test]
fn evaluate_zerofier_constraint() {
    // TODO: clean up this testcase
    let n = 2048;
    let challenges = &[Fp::from(999), Fp::from(43)];
    let curr_instr = 0;
    let permutation = 1;
    let alpha = 0;
    let a = 1;
    let instr = Fp::from(b'+');
    let constraint = curr_instr.curr()
        * (permutation.curr() * (alpha.get_challenge() - a.get_challenge() * curr_instr.curr())
            - permutation.next())
        + (curr_instr.curr() - instr) * (permutation.curr() - permutation.next());
    let blowup = 16;
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let curr_instr_column = vec![instr; n].to_vec_in(PageAlignedAllocator);
    let permutation_column = (0..n)
        .scan(Fp::one(), |acc, i| {
            let ret = *acc;
            *acc *= challenges[alpha] - challenges[a] * curr_instr_column[i];
            Some(ret)
        })
        .collect::<Vec<Fp>>()
        .to_vec_in(PageAlignedAllocator);
    let matrix = Matrix::new(vec![curr_instr_column, permutation_column]);
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    let matrix_group = MatrixGroup::new(vec![GroupItem::Fp(&lde_matrix)]);

    let constraint_eval = matrix_group.evaluate_symbolic(&[constraint], challenges, blowup);

    let constraint_eval_poly = constraint_eval.interpolate(lde_domain);
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

/// Generates a single column matrix of consisting of two values i.e.
/// ┌───────┐
/// │ Col 0 │
/// ├───────┤
/// │ 3     │
/// ├───────┤
/// │ 7     │
/// ├───────┤
/// │ 3     │
/// ├───────┤
/// │ 3     │
/// ├───────┤
/// │ 7     │
/// ├───────┤
/// │ ...   │
/// └───────┘
fn gen_binary_valued_matrix<F: GpuField>(n: usize, v1: F, v2: F) -> Matrix<F> {
    let mut rng = ark_std::test_rng();
    let mut col = Vec::with_capacity_in(n, PageAlignedAllocator);
    for _ in 0..n {
        if rng.gen() {
            col.push(v1)
        } else {
            col.push(v2);
        }
    }
    Matrix::new(vec![col])
}

fn assert_valid_over_transition_domain<F: GpuField>(
    domain: Radix2EvaluationDomain<F::FftField>,
    poly_matrix: Matrix<F>,
) {
    let mut x_values = domain.elements().map(|e| e.into()).collect::<Vec<F>>();
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
