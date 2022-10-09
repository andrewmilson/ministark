use crate::tables;
use ark_ff::One;
use ark_ff_optimized::fp64::Fp;
use ark_serialize::CanonicalSerialize;
use mini_stark::constraint::helper::are_eq;
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

pub struct BrainfuckAir {
    options: ProofOptions,
    trace_info: TraceInfo,
    execution_info: ExecutionInfo,
}

impl Air for BrainfuckAir {
    type Fp = Fp;
    type PublicInputs = ExecutionInfo;

    fn new(trace_info: TraceInfo, execution_info: ExecutionInfo, options: ProofOptions) -> Self {
        BrainfuckAir {
            options,
            trace_info,
            execution_info,
        }
    }

    fn options(&self) -> &ProofOptions {
        &self.options
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.execution_info
    }

    // fn boundary_constraints(&self) -> Vec<Constraint<Self::Fp>> {

    // }

    fn transition_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        vec![
            tables::Processor::transition_constraints(),
            // tables::Memory::transition_constraints(),
            tables::Instruction::transition_constraints(),
            // tables::Input::transition_constraints(),
            // tables::Output::transition_constraints(),
        ]
        .concat()
    }

    // fn terminal_constraints(&self) -> Vec<Constraint<Self::Fp>> {
    //     vec![7.curr() - self.result]
    // }

    fn trace_info(&self) -> &TraceInfo {
        &self.trace_info
    }
}
