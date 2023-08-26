#![feature(allocator_api)]
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::UniformRand;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use ark_std::rand::seq::SliceRandom;
use ministark::constraints::AlgebraicItem;
use ministark::constraints::Constraint;
use ministark::constraints::ExecutionTraceColumn;
use ministark::constraints::VerifierChallenge;
use ministark::expression::Expr;
use ministark::utils;
use ministark::utils::tests::gen_binary_valued_matrix;
use ministark::utils::tests::gen_fib_matrix;
use ministark::utils::FieldVariant;
use ministark::utils::GpuAllocator;
use ministark::Matrix;
use ministark::StarkExtensionOf;
use ministark_gpu::fields::p18446744069414584321::ark::Fp;
use ministark_gpu::GpuFftField;
use ministark_gpu::GpuField;
use num_traits::Pow;

// TODO: handle
// #[test]
// fn expressions_are_equal() {
//     use AlgebraicItem::*;
//     let mut rng = ark_std::test_rng();
//     let x = Fp::rand(&mut rng);
//     let left: Expr<AlgebraicItem<Fp>> = X.into();
//     let right: Expr<AlgebraicItem<Fp>> = X.pow(2) / X;
//     assert_eq!(left.evaluation_hash(x), right.evaluation_hash(x));
// }
// #[test]
// fn expressions_are_unequal() {
//     use AlgebraicItem::*;
//     let mut rng = ark_std::test_rng();
//     let x = Fp::rand(&mut rng);
//     let left: Expr<AlgebraicItem<Fp>> = X;
//     let right: Expr<AlgebraicItem<Fp>> = X.pow(3) / X;
//     assert_ne!(left.evaluation_hash(x), right.evaluation_hash(x));
// }

#[test]
fn constraint_degree() {
    use AlgebraicItem::*;
    let trace_degree = 2usize.pow(10) - 1;
    let expected_degree = (trace_degree + 1) - 100;
    let constraint = Constraint::<()>::new((Trace(0, 0) * X) / X.pow(100) + X);

    let (numerator, denominator) = constraint.degree(trace_degree);
    let actual_degree = numerator - denominator;

    assert_eq!(expected_degree, actual_degree);
}

#[test]
fn constraint_with_challenges() {
    // TODO: hints
    let constraint: Expr<AlgebraicItem<Fp>> = (0.challenge() - 0.curr()) * 1.curr();
    let challenges = [Fp::one()];
    let col_values = [Fp::one(), Fp::from(100)];

    use AlgebraicItem::*;
    assert!(constraint
        .eval::<FieldVariant<Fp, Fp>>(&mut |leaf| match leaf {
            X => unreachable!(),
            &Constant(v) => FieldVariant::Fp(v),
            &Hint(_) => unreachable!(),
            &Periodic(_) => todo!(),
            &Challenge(i) => FieldVariant::Fp(challenges[i]),
            &Trace(i, j) => {
                assert_eq!(0, j);
                FieldVariant::Fp(col_values[i])
            }
        })
        .is_zero());
}

#[test]
fn symbolic_evaluation_with_challenges() {
    let n = 2048;
    let constraint = Constraint::new((0.curr() - 0.challenge()) * (0.curr() - 1.challenge()));
    let (numerator_degree, denominator_degree) = constraint.degree(n);
    let blowup = utils::ceil_power_of_two((numerator_degree - denominator_degree) / n);
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let alpha = Fp::from(3);
    let beta = Fp::from(7);
    let matrix = gen_binary_valued_matrix(n, alpha, beta);
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);

    let constraint_eval = evaluate_symbolic(
        lde_domain,
        blowup,
        &[],
        &[alpha, beta],
        &constraint,
        &lde_matrix,
    );

    let constraint_eval_poly = constraint_eval.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_eval_poly);
}

#[test]
fn constraint_multiplication() {
    use AlgebraicItem::*;
    let zero = FieldVariant::Fp(Fp::zero());
    let one = FieldVariant::Fp(Fp::one());
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
    let between_0_and_10 = (0.curr() - Constant(one))
        * (0.curr() - Constant(two))
        * (0.curr() - Constant(three))
        * (0.curr() - Constant(four))
        * (0.curr() - Constant(five))
        * (0.curr() - Constant(six))
        * (0.curr() - Constant(seven))
        * (0.curr() - Constant(eight))
        * (0.curr() - Constant(nine));

    let f = |val: FieldVariant<Fp, Fp>| {
        move |leaf: &AlgebraicItem<FieldVariant<Fp, Fp>>| match leaf {
            X => one,
            Challenge(_) => unreachable!(),
            Hint(_) => unreachable!(),
            Periodic(_) => todo!(),
            &Constant(v) => v,
            &Trace(i, j) => {
                assert_eq!(0, i, "for value {val}");
                assert_eq!(0, j, "for value {val}");
                val
            }
        }
    };

    assert!(!between_0_and_10.eval(&mut f(-two)).is_zero());
    assert!(!between_0_and_10.eval(&mut f(-one)).is_zero());
    assert!(!between_0_and_10.eval(&mut f(zero)).is_zero());
    assert!(between_0_and_10.eval(&mut f(one)).is_zero());
    assert!(between_0_and_10.eval(&mut f(two)).is_zero());
    assert!(between_0_and_10.eval(&mut f(three)).is_zero());
    assert!(between_0_and_10.eval(&mut f(four)).is_zero());
    assert!(between_0_and_10.eval(&mut f(five)).is_zero());
    assert!(between_0_and_10.eval(&mut f(six)).is_zero());
    assert!(between_0_and_10.eval(&mut f(seven)).is_zero());
    assert!(between_0_and_10.eval(&mut f(eight)).is_zero());
    assert!(between_0_and_10.eval(&mut f(nine)).is_zero());
    assert!(!between_0_and_10.eval(&mut f(ten)).is_zero());
    assert!(!between_0_and_10.eval(&mut f(eleven)).is_zero());
    assert!(!between_0_and_10.eval(&mut f(twelve)).is_zero());
}

#[test]
fn evaluate_fibonacci_constraint() {
    let n = 2048;
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = trace_domain.get_coset(Fp::GENERATOR).unwrap();
    let matrix = gen_fib_matrix(n);
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);
    // TODO: don't use field variant here
    let constraints: Vec<Expr<AlgebraicItem<FieldVariant<Fp, Fp>>>> = vec![
        0.next() - (0.curr() + 1.curr()),
        1.next() - (0.next() + 1.curr()),
    ];

    let constraint_evals = Matrix::join(
        constraints
            .into_iter()
            .map(Constraint::new)
            .map(|c| evaluate_symbolic(lde_domain, 1, &[], &[], &c, &lde_matrix))
            .collect(),
    );

    let constraint_evals_poly = constraint_evals.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_evals_poly);
}

#[test]
fn evaluate_binary_constraint() {
    let n = 2048;
    // constrains column 0 values to 0 or 1
    use AlgebraicItem::*;
    let one = Constant(FieldVariant::Fp(Fp::one()));
    let constraint = Constraint::new(0.curr() * (0.curr() - one));
    let (numerator_degree, denominator_degree) = constraint.degree(n);
    let blowup = utils::ceil_power_of_two((numerator_degree - denominator_degree) / n);
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let matrix = gen_binary_valued_matrix(n, Fp::zero(), Fp::one());
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);

    let constraint_eval = evaluate_symbolic(lde_domain, blowup, &[], &[], &constraint, &lde_matrix);

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
        original_col.to_vec_in(GpuAllocator),
        shuffled_col.to_vec_in(GpuAllocator),
        original_product.to_vec_in(GpuAllocator),
        shuffled_product.to_vec_in(GpuAllocator),
    ]);
    let alpha = 0; // first verifier challenge
    let original_col = 0;
    let shuffled_col = 1;
    let original_product = 2;
    let shuffled_product = 3;
    let constraints = vec![
        original_product.curr() * (alpha.challenge() - original_col.curr())
            - original_product.next(),
        shuffled_product.curr() * (alpha.challenge() - shuffled_col.curr())
            - shuffled_product.next(),
    ];
    let blowup = 2;
    let trace_domain = Radix2EvaluationDomain::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);

    let constraint_evals = Matrix::join(
        constraints
            .into_iter()
            .map(Constraint::new)
            .map(|c| evaluate_symbolic(lde_domain, blowup, &[], &[challenge], &c, &lde_matrix))
            .collect(),
    );

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
    let constraint = Constraint::new(
        curr_instr.curr()
            * (permutation.curr() * (alpha.challenge() - a.challenge() * curr_instr.curr())
                - permutation.next())
            + (curr_instr.curr() - AlgebraicItem::Constant(FieldVariant::Fp(instr)))
                * (permutation.curr() - permutation.next()),
    );
    let blowup = 16;
    let trace_domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
    let lde_domain = Radix2EvaluationDomain::<Fp>::new_coset(n * blowup, Fp::GENERATOR).unwrap();
    let curr_instr_column = vec![instr; n].to_vec_in(GpuAllocator);
    let permutation_column = (0..n)
        .scan(Fp::one(), |acc, i| {
            let ret = *acc;
            *acc *= challenges[alpha] - challenges[a] * curr_instr_column[i];
            Some(ret)
        })
        .collect::<Vec<Fp>>()
        .to_vec_in(GpuAllocator);
    let matrix = Matrix::new(vec![curr_instr_column, permutation_column]);
    let poly_matrix = matrix.interpolate(trace_domain);
    let lde_matrix = poly_matrix.evaluate(lde_domain);

    let constraint_eval = evaluate_symbolic(
        lde_domain,
        blowup,
        &[],
        challenges,
        &constraint,
        &lde_matrix,
    );

    let constraint_eval_poly = constraint_eval.interpolate(lde_domain);
    assert_valid_over_transition_domain(trace_domain, constraint_eval_poly);
}

fn assert_valid_over_transition_domain<F: GpuField + Field>(
    domain: Radix2EvaluationDomain<F::FftField>,
    poly_matrix: Matrix<F>,
) where
    F: From<F::FftField>,
    F::FftField: FftField,
{
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

/// TODO: consider merging with ConstraintComposer::evaluate_constraint_cpu
fn evaluate_symbolic<Fp: GpuFftField + FftField, Fq: StarkExtensionOf<Fp>>(
    lde_domain: Radix2EvaluationDomain<Fp>,
    blowup_factor: usize,
    hints: &[Fq],
    challenges: &[Fq],
    constraint: &Constraint<FieldVariant<Fp, Fq>>,
    lde_matrix: &Matrix<Fq>,
) -> Matrix<Fq> {
    let blowup_factor = blowup_factor as isize;
    let xs = lde_domain.elements();
    let n = lde_domain.size();
    let mut result = Vec::with_capacity_in(n, GpuAllocator);
    result.resize(n, Fq::zero());

    for (i, (v, x)) in result.iter_mut().zip(xs).enumerate() {
        use AlgebraicItem::*;
        let eval_result = constraint.eval(&mut |leaf| match leaf {
            X => FieldVariant::Fp(x),
            &Constant(v) => v,
            &Hint(i) => FieldVariant::Fq(hints[i]),
            &Challenge(i) => FieldVariant::Fq(challenges[i]),
            &Periodic(_col) => todo!(),
            &Trace(col_idx, offset) => {
                let pos = (i as isize + blowup_factor * offset).rem_euclid(n as isize) as usize;
                let column = &lde_matrix[col_idx];
                FieldVariant::Fq(column[pos])
            }
        });

        *v = match eval_result {
            FieldVariant::Fp(v) => Fq::from(v),
            FieldVariant::Fq(v) => v,
        };
    }

    Matrix::new(vec![result])
}
