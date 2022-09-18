use legacy_algebra::fp_u128::BaseFelt;
use mini_stark::polynomial::MultivariatePolynomial;
use mini_stark::polynomial::Polynomial;
use mini_stark::StandardProofStream;
use mini_stark::Stark;
use rand::Rng;
use std::iter::once;

/// Generates a boundary constraint.
///
/// The first value corresponds to the execution cycle,
/// the second corresponds to the register and the third corresponds the value.
/// The boundary constraints are public to the prover and verifier.
fn boundary_constraint(cycle: usize, register: usize, value: u32) -> (usize, usize, BaseFelt) {
    (cycle, register, BaseFelt::from(value))
}

fn increment(counter: &mut usize) -> usize {
    *counter += 1;
    *counter
}

/// Transition constraints for read/write memory.
///
/// Alpha and beta provided by the verifier for a permutation check.
/// ref: https://hackmd.io/@bobbinth/HJr56BKKt
fn transition_constraints(
    alpha: BaseFelt,
    beta: BaseFelt,
) -> Vec<MultivariatePolynomial<BaseFelt>> {
    // constants
    let one = MultivariatePolynomial::one();
    let u16_values = MultivariatePolynomial::constant(2u128.pow(16).into());
    let alpha = MultivariatePolynomial::constant(alpha);
    let beta = MultivariatePolynomial::constant(beta);

    let mut column_counter = 1;
    // column values for current row of the trace
    let c = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let a = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let i = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    // field elements stored at a given context/address/clock cycle prior to the
    // memory operation
    let u0 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let u1 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let u2 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let u3 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    // field elements stored at a given context/address/clock cycle after the memory
    // operation
    let v0 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let v1 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let v2 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let v3 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let d0 =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let d1 =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let t = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let perm_prod = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));

    // Column values in next row of the trace
    let c_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let a_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let i_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    // field elements stored at a given context/address/clock cycle prior to the
    // memory operation
    let u0_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let u1_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let u2_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let u3_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    // field elements stored at a given context/address/clock cycle after the memory
    // operation
    let v0_prime =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let v1_prime =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let v2_prime =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let v3_prime =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let d0_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let d1_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let t_prime =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));
    let perm_prod_prime =
        MultivariatePolynomial::<BaseFelt>::lift(Polynomial::x(), increment(&mut column_counter));

    // Pseudo columns
    let n0 = (c_prime.clone() - c.clone()) * t.clone();
    let n1 = (a_prime.clone() - a.clone()) * t;

    // Constraints for `t` such that `t` guarantees:
    // - when context changes `n0 = 1`
    // - when context remains the same but address changes `(1 - n0) * n1 = 1`
    // - when context and address are unchanged `(1 - n0) * (1 - n1) = 1`
    let t_constraints = vec![
        (n0.clone() ^ 2) - n0.clone(), // n0 is always 0 or 1
        (one.clone() - n0.clone()) * (c_prime.clone() - c.clone()),
        (one.clone() - n0.clone()) * ((n1.clone() ^ 2) - n1.clone()),
        (one.clone() - n0.clone()) * (one.clone() - n1.clone()) * (a_prime.clone() - a.clone()),
    ];

    // Enforce the values of context, address and clock cycles grow monotonically:
    // 1. `n0 = 1` - when the context changes, columns `d0` and `d1` at the next row
    // should    contain the delta between the old and the new contexts.
    let context_increase_constraint = n0.clone()
        * ((c_prime.clone() - c.clone())
            - (d0_prime.clone() + d1_prime.clone() * u16_values.clone()));
    // 2. `n0 = 0`, `n1 = 1` - when the context remains the same but the address
    // changes, columns    `d0` and `d1` at the next row should contain the
    // delta between the old and the new addresses.
    let address_increase_constraint = (one.clone() - n0.clone())
        * (n1.clone())
        * ((a_prime.clone() - a.clone())
            - (d0_prime.clone() + d1_prime.clone() * u16_values.clone()));
    // 3. `n0 = 0`, `n1 = 0` - when both the context and the address remain the
    // same, columns    `d0` and `d1` at the next row should contain the delta
    // between the old and the new clock cycle.
    let clock_cycle_increase_constraint = (one.clone() - n0.clone())
        * (one.clone() - n1.clone())
        * ((i_prime.clone() - i.clone() - one.clone())
            - (d0_prime.clone() + d1_prime.clone() * u16_values.clone()));
    let _increase_constraint =
        context_increase_constraint * address_increase_constraint * clock_cycle_increase_constraint;
    // Since `n0, n1 ∈ (0, 1)` and given the fact one constraint has to hold each
    // row we can optimize the above constraint from a degree 13 polynomial to a
    // degree 5.
    let increase_constraint = n0.clone() * (c_prime.clone() - c.clone())
        + (one.clone() - n0.clone())
            * (n1.clone() * (a.clone() - a_prime.clone())
                + (one.clone() - n1.clone()) * (i_prime.clone() - i.clone() - one.clone()))
        - (d0_prime.clone() + d1_prime.clone() * u16_values.clone());
    assert_eq!(_increase_constraint.degree(), 13);
    assert_eq!(increase_constraint.degree(), 5);

    // make sure that values at a given memory address are always initialized to 0
    let init_zero_constraints = vec![
        (n0.clone() + (one.clone() - n0.clone()) * n0.clone()) * u0_prime.clone(),
        (n0.clone() + (one.clone() - n0.clone()) * n0.clone()) * u1_prime.clone(),
        (n0.clone() + (one.clone() - n0.clone()) * n0.clone()) * u2_prime.clone(),
        (n0.clone() + (one.clone() - n0.clone()) * n0.clone()) * u3_prime.clone(),
    ];

    // We need to make sure that for the same context/address combination, the `v0,
    // v1, ...` of the current row are equal to the `u0, u1, ...` of the next
    // row.
    let memory_address_consistency_constraints = vec![
        (one.clone() - n0.clone()) * (one.clone() - n1.clone()) * (u0_prime.clone() - v0.clone()),
        (one.clone() - n0.clone()) * (one.clone() - n1.clone()) * (u1_prime.clone() - v1.clone()),
        (one.clone() - n0.clone()) * (one.clone() - n1.clone()) * (u2_prime.clone() - v2.clone()),
        (one.clone() - n0.clone()) * (one.clone() - n1.clone()) * (u3_prime.clone() - v3.clone()),
    ];

    // To use the above table in permutation checks, we need to reduce each row of
    // the memory table to a single value. This can be done like so:
    let v = beta.clone()
        + alpha.clone() * c.clone()
        + (alpha.clone() ^ 2) * a.clone()
        + (alpha.clone() ^ 3) * i.clone()
        + (alpha.clone() ^ 4) * u0.clone()
        + (alpha.clone() ^ 5) * u1.clone()
        + (alpha.clone() ^ 6) * u2.clone()
        + (alpha.clone() ^ 7) * u3.clone()
        + (alpha.clone() ^ 8) * v0.clone()
        + (alpha.clone() ^ 9) * v1.clone()
        + (alpha.clone() ^ 10) * v2.clone()
        + (alpha.clone() ^ 11) * v3.clone();

    // Used for the permutation check.
    let running_product = v - perm_prod;

    t_constraints
        .into_iter()
        .chain(once(running_product))
        .chain(once(increase_constraint))
        .chain(init_zero_constraints.into_iter())
        .chain(memory_address_consistency_constraints.into_iter())
        .collect()

    // let v = MultivariatePolynomial::lift(Polynomial::x(), 1);
    // let p0 = MultivariatePolynomial::lift(Polynomial::x(), 2);
    // let v_prime = MultivariatePolynomial::lift(Polynomial::x(), 3);
    // let p0_prime = MultivariatePolynomial::lift(Polynomial::x(), 4);
    // let one = MultivariatePolynomial::one();
    // let alpha = MultivariatePolynomial::constant(alpha);

    // // Constraint makes sure the register increases by '1' each cycle.
    // Context: // - `x` symbolizes the value of the first register of the
    // current cycle. // - `y` symbolizes the value of the first register of
    // the next cycle. // Transition constraints need to evaluate to '0' to
    // represent valid transitions. vec![
    //     (v - v_prime) * (v_prime - v - one),
    //     p0_prime - p0 * (alpha + v),
    // ]
}

fn main() {
    // // Parameters for the program
    // let starting_number = 7;
    // let increment_times = 100;

    // let starting_value = 0;
    // let ending_value = 255;
    // let starting_value_constraint = boundary_constraint(0, 0, starting_value);
    // let ending_value_constraint = boundary_constraint(256, 0, ending_value);
    // let boundary_constraints = vec![starting_value_constraint,
    // ending_value_constraint];

    let mut rng = rand::thread_rng();
    let alpha = BaseFelt::new(rng.gen());
    let beta = BaseFelt::new(rng.gen());

    let transition_constraints = transition_constraints(alpha, beta);
    let transition_constraints_degree = transition_constraints
        .iter()
        .map(|transition_constraint| transition_constraint.degree() as usize)
        // .max()
        // .unwrap()
        .collect::<Vec<usize>>();
    println!(
        "Transition constraint degree: {:?}",
        transition_constraints_degree
    );

    // // let trace = (starting_value..=ending_value).collect::<Vec<u32>>();

    // // Verifier wants the:
    // // - first register of the first cycle to correspond to the program's
    // starting value. // - first register of the last cycle to correspond
    // to the program's output value.

    // // ====================================

    // // Construct a stark with a 128 bit security level
    // let security_level = 128;
    // let stark = Stark::new(
    //     4,
    //     security_level / 2,
    //     security_level,
    //     1,
    //     program_cycles,
    //     transition_constraint_degree,
    // );
    // let (transition_zerofier, transition_zerofier_codeword,
    // transition_zerofier_root) =     stark.preprocess();

    // // Prover generates the STARK proof.
    // let mut prover_proof_stream = StandardProofStream::new();
    // println!("Generating proof to prove program output
    // '{program_output}'..."); let now = Instant::now();
    // let serialized_proof = stark.prove(
    //     program_trace,
    //     vec![transition_constraint.clone()],
    //     &boundary_constraints,
    //     &transition_zerofier,
    //     &transition_zerofier_codeword,
    //     &mut prover_proof_stream,
    // );
    // println!("Proof generated in {:.2?}", now.elapsed());
    // println!("Proof size {}KB\n\n", serialized_proof.len() / 1000);

    // // Verifier receives and verifies the STARK proof.
    // let mut verifier_proof_stream = StandardProofStream::new();
    // println!("Verifying validity of proof for program output
    // {program_output}..."); let now = Instant::now();
    // let verification_result = stark.verify(
    //     &serialized_proof,
    //     vec![transition_constraint],
    //     &boundary_constraints,
    //     transition_zerofier_root,
    //     &mut verifier_proof_stream,
    // );
    // println!(
    //     "Verified validity of proof for program output '{program_output}' in
    // {:.2?}",     now.elapsed()
    // );

    // println!(
    //     "Output '{program_output}' is valid? {}",
    //     if verification_result.is_ok() {
    //         "yes"
    //     } else {
    //         "no"
    //     }
    // );
}
