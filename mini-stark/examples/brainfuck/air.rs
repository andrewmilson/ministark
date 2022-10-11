use crate::tables;
use ark_ff::One;
use ark_ff_optimized::fp64::Fp;
use ark_serialize::CanonicalSerialize;
use mini_stark::constraint::are_eq;
use mini_stark::Air;
use mini_stark::Column;
use mini_stark::Constraint;
use mini_stark::ProofOptions;
use mini_stark::TraceInfo;

#[derive(CanonicalSerialize)]
pub struct ExecutionInfo {
    pub execution_len: usize,
    pub input: Vec<usize>,
    pub output: Vec<usize>,
}

// Quitient pre:32767 after:20470
// Quitient pre:32767 after:16376
// Quitient pre:32767 after:16376
// Quitient pre:32767 after:0
// Quitient pre:32767 after:4094
// Quitient pre:32767 after:4094
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:32767 13
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:2047
// Quitient pre:32767 after:32767
// Quitient pre:32767 after:32767
// Quitient pre:0 after:0
// Quitient pre:32767 after:0 21

pub struct BrainfuckAir {
    options: ProofOptions,
    trace_info: TraceInfo,
    execution_info: ExecutionInfo,
    // boundary_constraints: Vec<Constraint<Fp>>,
    transition_constraints: Vec<Constraint<Fp>>,
    // terminal_constraints: Vec<Constraint<Fp>>,
}

impl Air for BrainfuckAir {
    type Fp = Fp;
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
        }
    }

    fn options(&self) -> &ProofOptions {
        &self.options
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.execution_info
    }

    fn transition_constraints(&self) -> &[Constraint<Self::Fp>] {
        &self.transition_constraints
    }

    fn trace_info(&self) -> &TraceInfo {
        &self.trace_info
    }
}
