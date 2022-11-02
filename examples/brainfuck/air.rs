use crate::tables;
use ark_ff_optimized::fp64::Fp;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ministark::Air;
use ministark::Constraint;
use ministark::ProofOptions;
use ministark::TraceInfo;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct ExecutionInfo {
    pub execution_len: usize,
    pub input: Vec<usize>,
    pub output: Vec<usize>,
}

pub struct BrainfuckAir {
    options: ProofOptions,
    trace_info: TraceInfo,
    execution_info: ExecutionInfo,
    boundary_constraints: Vec<Constraint<Fp>>,
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
            // TODO
            boundary_constraints: vec![],
            //     tables::ProcessorBaseColumn::boundary_constraints(),
            //     tables::ProcessorExtensionColumn::boundary_constraints(),
            //     tables::MemoryBaseColumn::boundary_constraints(),
            //     tables::MemoryExtensionColumn::boundary_constraints(),
            //     tables::InstructionBaseColumn::boundary_constraints(),
            //     tables::InstructionExtensionColumn::boundary_constraints(),
            //     tables::InputExtensionColumn::boundary_constraints(),
            //     tables::OutputExtensionColumn::boundary_constraints(),
            // ]
            // .concat(),
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

    fn boundary_constraints(&self) -> &[Constraint<Self::Fp>] {
        &self.boundary_constraints
    }

    fn trace_info(&self) -> &TraceInfo {
        &self.trace_info
    }
}
