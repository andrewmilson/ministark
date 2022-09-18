//! STARK proof and verification of a program that increases the number `7` by
//! `+1` each cycle.
#![feature(allocator_api)]

use fast_poly::allocator::PageAlignedAllocator;
use legacy_algebra::fp_u128::BaseFelt;
use mini_stark::polynomial::MultivariatePolynomial;
use mini_stark::polynomial::Polynomial;
use mini_stark::StandardProofStream;
use mini_stark::Stark;
use std::time::Instant;

/// Takes a number `x` as input and returns `x + 1`.
///
/// This is the function the stark proof will be generated over. Look at
/// [transition_constraints] to see the constraint that ensures the
/// computational integrity of this function.
fn increment_by_one(x: u32) -> u32 {
    x + 1
}

/// Generates a boundary constraint.
///
/// The first value corresponds to the execution cycle,
/// the second corresponds to the register and the third corresponds the value.
/// The boundary constraints are public to the prover and verifier.
fn boundary_constraint(cycle: usize, register: usize, value: u32) -> (usize, usize, BaseFelt) {
    (cycle, register, BaseFelt::from(value))
}

/// Generates the constraint to specify the number increases by '1' from one
/// cycle to the next.
///
/// Constraints are expressed as multivariate polynomials. Variables represent
/// register values in either the current and next execution cycle.
fn transition_constraint() -> MultivariatePolynomial<BaseFelt> {
    let x = MultivariatePolynomial::lift(Polynomial::x(), 1);
    let y = MultivariatePolynomial::lift(Polynomial::x(), 2);
    let one = MultivariatePolynomial::one();

    // Constraint makes sure the register increases by '1' each cycle. Context:
    // - `x` symbolizes the value of the first register of the current cycle.
    // - `y` symbolizes the value of the first register of the next cycle.
    // Transition constraints need to evaluate to '0' to represent valid
    // transitions.
    y - x - one
}

/// Returns program output and generates a trace of its execution.
///
/// The trace is an array representing the state of program every cycle. The
/// state of the program at each cycle is obtained by providing the values of
/// all registers. For this program only a single register is needed.
/// It's technically an "algebraic" trace due to the fact that register values
/// consist of finite field elements.
fn execute_program_with_trace(
    starting_number: u32,
    increment_times: u32,
) -> (u32, Vec<Vec<BaseFelt, PageAlignedAllocator>>) {
    let mut trace = Vec::new();

    // The starting state needs to be in the trace.
    let mut starting_register_values = Vec::new_in(PageAlignedAllocator);
    starting_register_values.push(BaseFelt::from(starting_number));
    trace.push(starting_register_values);

    let mut last_value = starting_number;

    for _ in 1..increment_times {
        last_value = increment_by_one(last_value);
        // This example only needs one register for capturing the state.
        let mut next_register_values = Vec::new_in(PageAlignedAllocator);
        next_register_values.push(BaseFelt::from(last_value));
        trace.push(next_register_values);
    }

    (last_value, trace)
}

fn main() {
    println!("IS THIS WORKING");
    // Parameters for the program
    let starting_number = 7;
    let increment_times = 150;

    // Execution and generation of the trace. Prover would provide the output along
    // with the STARK proof to the verifier. The program increments the number
    // every cycle. **Trace only needs to be known by the prover**.
    let (program_output, program_trace) =
        execute_program_with_trace(starting_number, increment_times);
    let program_cycles = program_trace.len();

    // Public inputs to determine correct execution of the program.
    // **Known by both prover and verifier**.
    // ====================================
    let transition_constraint = transition_constraint();
    let transition_constraint_degree = transition_constraint.degree() as usize;
    let transition_constraints = vec![transition_constraint];
    // Verifier wants the:
    // - first register of the first cycle to correspond to the program's starting
    //   value.
    // - first register of the last cycle to correspond to the program's output
    //   value.
    let program_input_constraint = boundary_constraint(0, 0, starting_number);
    let program_output_constraint = boundary_constraint(program_cycles - 1, 0, program_output);
    let boundary_constraints = vec![program_input_constraint, program_output_constraint];
    // ====================================

    // Construct a stark with a 128 bit security level
    let security_level = 128;
    let stark = Stark::new(
        4,
        security_level / 2,
        security_level,
        1,
        program_cycles,
        transition_constraint_degree,
    );
    let (transition_zerofier, transition_zerofier_codeword, transition_zerofier_root) =
        stark.preprocess();

    // Prover generates the STARK proof.
    let mut prover_proof_stream = StandardProofStream::new();
    println!("WHOA");
    println!("Generating proof to prove program output '{program_output}'...");
    let now = Instant::now();
    let serialized_proof = stark.prove(
        program_trace,
        transition_constraints.clone(),
        &boundary_constraints,
        &transition_zerofier,
        &transition_zerofier_codeword,
        &mut prover_proof_stream,
    );
    println!("Proof generated in {:.2?}", now.elapsed());
    println!("Proof size {}KB\n\n", serialized_proof.len() / 1000);

    // Verifier receives and verifies the STARK proof.
    let mut verifier_proof_stream = StandardProofStream::new();
    println!("Verifying validity of proof for program output {program_output}...");
    let now = Instant::now();
    let verification_result = stark.verify(
        &serialized_proof,
        transition_constraints,
        &boundary_constraints,
        transition_zerofier_root,
        &mut verifier_proof_stream,
    );
    println!(
        "Verified validity of proof for program output '{program_output}' in {:.2?}",
        now.elapsed()
    );

    println!(
        "Output '{program_output}' is valid? {}",
        if verification_result.is_ok() {
            "yes"
        } else {
            "no"
        }
    );
}
