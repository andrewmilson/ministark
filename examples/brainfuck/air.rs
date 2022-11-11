use crate::tables;
use crate::tables::Challenge;
use crate::tables::EvaluationArgumentHint;
use crate::vm::compile;
use ark_ff::Field;
use ark_ff::Zero;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use gpu_poly::fields::p18446744069414584321::Fp;
use gpu_poly::fields::p18446744069414584321::Fq3;
use ministark::challenges::Challenges;
use ministark::constraint::Hint;
use ministark::hints::Hints;
use ministark::Air;
use ministark::Constraint;
use ministark::ProofOptions;
use ministark::TraceInfo;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct ExecutionInfo {
    pub source_code: String,
    pub input: Vec<u8>,
    pub output: Vec<u8>,
}

pub struct BrainfuckAir {
    options: ProofOptions,
    trace_info: TraceInfo,
    execution_info: ExecutionInfo,
    boundary_constraints: Vec<Constraint<Fq3>>,
    transition_constraints: Vec<Constraint<Fq3>>,
    terminal_constraints: Vec<Constraint<Fq3>>,
}

impl Air for BrainfuckAir {
    type Fp = Fp;
    type Fq = Fq3;
    type PublicInputs = ExecutionInfo;

    fn new(trace_info: TraceInfo, execution_info: ExecutionInfo, options: ProofOptions) -> Self {
        BrainfuckAir {
            options,
            trace_info,
            execution_info,
            transition_constraints: vec![
                tables::ProcessorBaseColumn::transition_constraints(),
                tables::ProcessorExtensionColumn::transition_constraints(),
                tables::MemoryBaseColumn::transition_constraints(),
                tables::MemoryExtensionColumn::transition_constraints(),
                tables::InstructionBaseColumn::transition_constraints(),
                tables::InstructionExtensionColumn::transition_constraints(),
                tables::InputExtensionColumn::transition_constraints(),
                tables::OutputExtensionColumn::transition_constraints(),
            ]
            .concat(),
            boundary_constraints: vec![
                tables::ProcessorBaseColumn::boundary_constraints(),
                tables::ProcessorExtensionColumn::boundary_constraints(),
                tables::MemoryBaseColumn::boundary_constraints(),
                tables::InstructionBaseColumn::boundary_constraints(),
                tables::InstructionExtensionColumn::boundary_constraints(),
                tables::InputExtensionColumn::boundary_constraints(),
                tables::OutputExtensionColumn::boundary_constraints(),
            ]
            .concat(),
            terminal_constraints: vec![
                tables::ProcessorExtensionColumn::terminal_constraints(),
                tables::InstructionExtensionColumn::terminal_constraints(),
                tables::InputExtensionColumn::terminal_constraints(),
                tables::OutputExtensionColumn::terminal_constraints(),
            ]
            .concat(),
        }
    }

    fn get_hints(&self, challenges: &Challenges<Self::Fq>) -> Hints<Self::Fq> {
        use Challenge::*;
        use EvaluationArgumentHint::*;

        let ExecutionInfo {
            source_code,
            input,
            output,
        } = &self.execution_info;
        let trace_len = self.trace_info().trace_len;

        let (input_eval_arg, input_eval_offset) =
            io_terminal_helper(input, challenges[Gamma], trace_len);
        let (output_eval_arg, output_eval_offset) =
            io_terminal_helper(output, challenges[Delta], trace_len);
        let instruction_eval_arg = compute_instruction_evaluation_argument(source_code, challenges);

        Hints::new(vec![
            (Instruction.index(), instruction_eval_arg),
            (Input.index(), input_eval_arg),
            (InputOffset.index(), input_eval_offset),
            (Output.index(), output_eval_arg),
            (OutputOffset.index(), output_eval_offset),
        ])
    }

    fn options(&self) -> &ProofOptions {
        &self.options
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.execution_info
    }

    fn transition_constraints(&self) -> &[Constraint<Fq3>] {
        &self.transition_constraints
    }

    fn boundary_constraints(&self) -> &[Constraint<Fq3>] {
        &self.boundary_constraints
    }

    fn terminal_constraints(&self) -> &[Constraint<Fq3>] {
        &self.terminal_constraints
    }

    fn trace_info(&self) -> &TraceInfo {
        &self.trace_info
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
        acc = acc * challenges[Eta]
            + challenges[A] * Fp::from(ip as u64)
            + challenges[B] * Fp::from(curr_instr as u64)
            + challenges[C] * Fp::from(next_instr as u64);
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
