use mini_stark::{
    polynomial::{MultivariatePolynomial, Polynomial},
    prime_field_u128::BaseElement as PrimeFieldElement,
    StandardProofStream, Stark,
};
use rand::Rng;

fn boundary_constraint(
    cycle: usize,
    register: usize,
    value: u32,
) -> (usize, usize, PrimeFieldElement) {
    (cycle, register, PrimeFieldElement::from(value))
}

/// Alpha provided by the verifier.
fn transition_constraints(
    alpha: PrimeFieldElement,
) -> Vec<MultivariatePolynomial<PrimeFieldElement>> {
    let mut column_counter = 1;
    let v = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let p0 = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let v_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));
    let p0_prime = MultivariatePolynomial::lift(Polynomial::x(), increment(&mut column_counter));

    let one = MultivariatePolynomial::one();
    let alpha = MultivariatePolynomial::constant(alpha);

    vec![
        (v.clone() - v_prime.clone()) * (v_prime.clone() - v.clone() - one.clone()),
        p0_prime.clone() - p0.clone() * (alpha.clone() + v.clone()),
    ]
}

fn increment(counter: &mut usize) -> usize {
    *counter += 1;
    *counter
}

fn main() {
    let starting_value = 0;
    let ending_value = 255;
    let starting_value_constraint = boundary_constraint(0, 0, starting_value);
    let ending_value_constraint = boundary_constraint(256, 0, ending_value);
    let boundary_constraints = vec![starting_value_constraint, ending_value_constraint];

    let mut rng = rand::thread_rng();
    let alpha = PrimeFieldElement::new(rng.gen());

    let transition_constraints = transition_constraints(alpha);
    let transition_constraints_degree = transition_constraints
        .iter()
        .map(|transition_constraint| transition_constraint.degree() as usize)
        .max()
        .unwrap();
    println!(
        "Transition constraint degree: {}",
        transition_constraints_degree
    );

    // let trace = (starting_value..=ending_value).collect::<Vec<u32>>();

    // // Verifier wants the:
    // // - first register of the first cycle to correspond to the program's starting value.
    // // - first register of the last cycle to correspond to the program's output value.

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
    // let (transition_zerofier, transition_zerofier_codeword, transition_zerofier_root) =
    //     stark.preprocess();

    // // Prover generates the STARK proof.
    // let mut prover_proof_stream = StandardProofStream::new();
    // println!("Generating proof to prove program output '{program_output}'...");
    // let now = Instant::now();
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
    // println!("Verifying validity of proof for program output {program_output}...");
    // let now = Instant::now();
    // let verification_result = stark.verify(
    //     &serialized_proof,
    //     vec![transition_constraint],
    //     &boundary_constraints,
    //     transition_zerofier_root,
    //     &mut verifier_proof_stream,
    // );
    // println!(
    //     "Verified validity of proof for program output '{program_output}' in {:.2?}",
    //     now.elapsed()
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
