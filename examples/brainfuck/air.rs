use crate::tables;
use crate::tables::Challenge;
use crate::tables::EvaluationArgumentHint;
use crate::vm::compile;
use crate::BrainfuckClaim;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ministark::air::AirConfig;
use ministark::challenges::Challenges;
use ministark::constraints::AlgebraicItem;
use ministark::constraints::Constraint;
use ministark::constraints::Hint;
use ministark::constraints::VerifierChallenge;
use ministark::hints::Hints;
use ministark::utils::FieldVariant;
use ministark_gpu::fields::p18446744069414584321::ark::Fp;
use ministark_gpu::fields::p18446744069414584321::ark::Fq3;
use num_traits::Pow;

pub struct BrainfuckAirConfig;

impl AirConfig for BrainfuckAirConfig {
    const NUM_BASE_COLUMNS: usize = 17;
    const NUM_EXTENSION_COLUMNS: usize = 9;

    type Fp = Fp;
    type Fq = Fq3;
    type PublicInputs = BrainfuckClaim;

    fn gen_hints(
        trace_len: usize,
        execution_info: &BrainfuckClaim,
        challenges: &Challenges<Self::Fq>,
    ) -> Hints<Self::Fq> {
        use Challenge::*;
        use EvaluationArgumentHint::*;
        let BrainfuckClaim {
            source_code,
            input,
            output,
        } = execution_info;

        let (input_eval_arg, input_eval_offset) =
            io_terminal_helper(input, challenges[Gamma.index()], trace_len);
        let (output_eval_arg, output_eval_offset) =
            io_terminal_helper(output, challenges[Delta.index()], trace_len);
        let instruction_eval_arg = compute_instruction_evaluation_argument(source_code, challenges);

        Hints::new(vec![
            (Instruction.index(), instruction_eval_arg),
            (Input.index(), input_eval_arg),
            (InputOffset.index(), input_eval_offset),
            (Output.index(), output_eval_arg),
            (OutputOffset.index(), output_eval_offset),
        ])
    }

    fn constraints(trace_len: usize) -> Vec<Constraint<FieldVariant<Self::Fp, Self::Fq>>> {
        use AlgebraicItem::*;
        let one = Constant(FieldVariant::<Fp, Fq3>::Fp(Fp::one()));
        let trace_xs = Radix2EvaluationDomain::<Fp>::new(trace_len).unwrap();
        let first_trace_x = Constant(FieldVariant::<Fp, Fq3>::Fp(trace_xs.element(0)));
        let last_trace_x = Constant(FieldVariant::Fp(trace_xs.element(trace_len - 1)));

        let transition_constraints = [
            tables::ProcessorBaseColumn::transition_constraints(),
            tables::ProcessorExtensionColumn::transition_constraints(),
            tables::MemoryBaseColumn::transition_constraints(),
            tables::MemoryExtensionColumn::transition_constraints(),
            tables::InstructionBaseColumn::transition_constraints(),
            tables::InstructionExtensionColumn::transition_constraints(),
            tables::InputExtensionColumn::transition_constraints(),
            tables::OutputExtensionColumn::transition_constraints(),
        ]
        .into_iter()
        .flatten()
        .map(|constraint| {
            // ensure constraints hold in all rows except the last
            // multiply by `(x - t_(n-1))` to remove the last term
            // NOTE: `x^trace_len - 1 = (x - t_0)(x - t_1)...(x - t_(n-1))`
            // NOTE: `t^(n-1) = t^(-1)`
            constraint * ((X - last_trace_x) / (X.pow(trace_len) - one))
        });

        let boundary_constraints = [
            tables::ProcessorBaseColumn::boundary_constraints(),
            tables::ProcessorExtensionColumn::boundary_constraints(),
            tables::MemoryBaseColumn::boundary_constraints(),
            tables::InstructionBaseColumn::boundary_constraints(),
            tables::InstructionExtensionColumn::boundary_constraints(),
            tables::InputExtensionColumn::boundary_constraints(),
            tables::OutputExtensionColumn::boundary_constraints(),
        ]
        .into_iter()
        .flatten()
        .map(|constraint| {
            // ensure constraint holds in the first row
            // symbolically divide `(x - t_0)`
            constraint / (X - first_trace_x)
        });

        let terminal_constraints = [
            tables::ProcessorExtensionColumn::terminal_constraints(),
            tables::InstructionExtensionColumn::terminal_constraints(),
            tables::InputExtensionColumn::terminal_constraints(),
            tables::OutputExtensionColumn::terminal_constraints(),
        ]
        .into_iter()
        .flatten()
        .map(|constraint| {
            // ensure constraint holds in the last row
            // symbolically divide `(x - t_(n-1))`
            // NOTE: `t^(n-1) = t^(-1)`
            constraint / (X - last_trace_x)
        });

        transition_constraints
            .chain(boundary_constraints)
            .chain(terminal_constraints)
            .map(Constraint::from)
            .collect()
    }
}

// Computes the evaluation terminal for the instruction table
fn compute_instruction_evaluation_argument(source_code: &str, challenges: &Challenges<Fq3>) -> Fq3 {
    use Challenge::Eta;
    use Challenge::A;
    use Challenge::B;
    use Challenge::C;
    let mut program = compile(source_code);
    // add padding
    program.push(0);
    // let prev_ip = None;
    let mut acc = Fq3::zero();
    for (ip, curr_instr) in program.iter().copied().enumerate() {
        let next_instr = program.get(ip + 1).copied().unwrap_or(0);
        acc = acc * challenges[Eta.index()]
            + challenges[A.index()] * Fp::from(ip as u64)
            + challenges[B.index()] * Fp::from(curr_instr as u64)
            + challenges[C.index()] * Fp::from(next_instr as u64);
    }
    acc
}

// Computes the evaluation terminal for the input and output table
// output is of the form `(evaluatoin_argument, evaluation_offset)`
fn io_terminal_helper<F: Field>(symbols: &[u8], challenge: F, trace_len: usize) -> (F, F) {
    let mut acc = F::zero();
    for symbol in symbols {
        acc = challenge * acc + F::from(*symbol as u64);
    }
    let evaluation_argument = acc;
    // from BrainSTARK
    // In every additional row, the running evaluation variable is
    // multiplied by another `challenge` factor. So we multiply by
    // `challenge^(trace_len - num_symbols)` to get the value of
    // the evaluation terminal after all 2^k trace rows.
    let offset = challenge.pow([(trace_len - symbols.len()) as u64]);
    (evaluation_argument, offset)
}
