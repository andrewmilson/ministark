use crate::vm::OpCode;
use ark_ff::Zero;
use fast_poly::GpuField;
use mini_stark::constraint::Challenge;
use mini_stark::constraint::Column;
use mini_stark::Constraint;
use std::borrow::Borrow;

enum Challenges {
    A,
    B,
    C,
    D,
    E,
    F,
    Alpha,
    Beta,
    Gamma,
    Delta,
    Eta,
}

impl Challenge for Challenges {
    fn index(&self) -> usize {
        match self {
            Self::A => 0,
            Self::B => 1,
            Self::C => 2,
            Self::D => 3,
            Self::E => 4,
            Self::F => 5,
            Self::Alpha => 6,
            Self::Beta => 7,
            Self::Gamma => 8,
            Self::Delta => 9,
            Self::Eta => 10,
        }
    }
}

pub enum Processor {
    Cycle,
    Ip,
    CurrInstr,
    NextInstr,
    Mp,
    MemVal,
    MemValInv,
    InstructionPermutation,
    MemoryPermutation,
    InputEvaluation,
    OutputEvaluation,
}

impl Processor {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        // use Challenges::Alpha;
        // use Challenges::Eta;
        // use Challenges::A;
        // use Challenges::B;
        // use Challenges::C;
        use Processor::*;

        let zero = F::zero();
        let one = F::one();
        let two = one + one;
        let mem_val_is_zero = MemVal.curr() * MemValInv.curr() - one;
        let mut constraints = (Constraint::zero(), Constraint::zero(), Constraint::zero());

        use OpCode::*;
        for instr in OpCode::iterator().copied() {
            // max degree: 4
            let mut instr_constraints =
                (Constraint::zero(), Constraint::zero(), Constraint::zero());

            match instr {
                IncrementPointer => {
                    instr_constraints.0 = Ip.next() - Ip.curr() - one;
                    instr_constraints.1 = Mp.next() - Mp.curr() - one;
                }
                DecrementPointer => {
                    instr_constraints.0 = Ip.next() - Ip.curr() - one;
                    instr_constraints.1 = Mp.next() - Mp.curr() + one;
                }
                Increment => {
                    instr_constraints.0 = Ip.next() - Ip.curr() - one;
                    instr_constraints.1 = Mp.next() - Mp.curr();
                    instr_constraints.2 = MemVal.next() - MemVal.curr() - one;
                }
                Decrement => {
                    instr_constraints.0 = Ip.next() - Ip.curr() - one;
                    instr_constraints.1 = Mp.next() - Mp.curr();
                    instr_constraints.2 = MemVal.next() - MemVal.curr() + one;
                }
                Write => {
                    instr_constraints.0 = Ip.next() - Ip.curr() - one;
                    instr_constraints.1 = Mp.next() - Mp.curr();
                }
                Read => {
                    instr_constraints.0 = Ip.next() - Ip.curr() - one;
                    instr_constraints.1 = Mp.next() - Mp.curr();
                    instr_constraints.2 = MemVal.next() - MemVal.curr();
                }
                LoopBegin => {
                    instr_constraints.0 = MemVal.curr() * (Ip.next() - Ip.curr() - two)
                        + mem_val_is_zero.clone() * (Ip.next() - NextInstr.curr());
                    instr_constraints.1 = Mp.next() - Mp.curr();
                    instr_constraints.2 = MemVal.next() - MemVal.curr();
                }
                LoopEnd => {
                    instr_constraints.0 = &mem_val_is_zero * (Ip.next() - Ip.curr() - two)
                        + MemVal.curr() * (Ip.next() - NextInstr.curr());
                    instr_constraints.1 = Mp.next() - Mp.curr();
                    instr_constraints.2 = MemVal.next() - MemVal.curr();
                }
            }

            // max degree: 7
            let deselector = if_not_instr(instr, CurrInstr.curr());

            // TODO: mul assign
            // account for padding and deactivate all polynomials if curr instruction is 0
            constraints.0 = constraints.0 + &deselector * instr_constraints.0 * CurrInstr.curr();
            constraints.1 = constraints.1 + &deselector * instr_constraints.1 * CurrInstr.curr();
            constraints.2 = constraints.2 + &deselector * instr_constraints.2 * CurrInstr.curr();
        }

        vec![
            constraints.0,
            constraints.1,
            constraints.2,
            // cycle independent constraints
            Cycle.next() - Cycle.curr() - one,
            MemVal.curr() * &mem_val_is_zero,
            MemValInv.curr() * &mem_val_is_zero,
        ]
    }
}

impl Processor {
    const TRACE_START_INDEX: usize = 0;
}

impl Column for Processor {
    fn index(&self) -> usize {
        Self::TRACE_START_INDEX + *self as usize
    }
}

pub enum Memory {
    Cycle,
    Mp,
    MemVal,
    Dummy,
    Permutation,
}

impl Memory {
    const TRACE_START_INDEX: usize = 11;
}

impl Column for Memory {
    fn index(&self) -> usize {
        Self::TRACE_START_INDEX + *self as usize
    }
}

pub enum Instruction {
    Ip,
    CurrInstr,
    NextInstr,
    ProcessorPermutation,
    ProgramEvaluation,
}

impl Instruction {
    const TRACE_START_INDEX: usize = 16;

    fn boundary_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenges::A;
        use Challenges::B;
        use Challenges::C;
        use Instruction::*;
        vec![
            Ip.curr(),
            ProgramEvaluation.curr()
                - A.get_challenge() * Ip.curr()
                - B.get_challenge() * CurrInstr.curr()
                - C.get_challenge() * NextInstr.curr(),
        ]
    }

    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenges::Alpha;
        use Challenges::Eta;
        use Challenges::A;
        use Challenges::B;
        use Challenges::C;
        use Instruction::*;
        let one = F::one();
        vec![
            // instruction pointer increases by 0 or 1
            (Ip.next() - Ip.curr() - one) * (Ip.next() - Ip.curr()),
            // if address increases the next instruction in the current row must equal the current
            // instruction in the next row
            (Ip.next() - Ip.curr()) * (NextInstr.curr() - CurrInstr.next()),
            // if address is the same, then current instruction is also
            (Ip.next() - Ip.curr() - one) * (CurrInstr.next() - CurrInstr.curr()),
            // if address is the same, then next instruction is also
            (Ip.next() - Ip.curr() - one) * (NextInstr.next() - NextInstr.curr()),
            // - processor permutation changes correctly if ip changes
            // - processor permutation doesn't change if `curr_instr=0` i.e. padding
            // - processor permutation doesn't change if `ip` stays the same
            CurrInstr.curr()
                * (Ip.curr() - Ip.next() + one)
                * (ProcessorPermutation.next()
                    - ProcessorPermutation.curr()
                        * (Alpha.get_challenge()
                            - A.get_challenge() * Ip.next()
                            - B.get_challenge() * CurrInstr.next()
                            - C.get_challenge() * NextInstr.next()))
                + instr_zerofier(CurrInstr.curr())
                    * (ProcessorPermutation.next() - ProcessorPermutation.curr())
                + (Ip.curr() - Ip.next())
                    * (ProcessorPermutation.curr() - ProcessorPermutation.next()),
            // - no evaluation change if `ip` remains the same
            // - evaluation change if `ip` changes
            (Ip.next() - Ip.curr() - one) * (ProgramEvaluation.next() - ProgramEvaluation.curr())
                + (Ip.next() - Ip.curr())
                    * (ProgramEvaluation.next()
                        - ProgramEvaluation.curr() * Eta.get_challenge()
                        - A.get_challenge() * Ip.next()
                        - B.get_challenge() * CurrInstr.next()
                        - C.get_challenge() * NextInstr.next()),
        ]
    }

    fn terminal_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        todo!()
    }
}

impl Column for Instruction {
    fn index(&self) -> usize {
        Self::TRACE_START_INDEX + *self as usize
    }
}

enum Input {
    Value,
    Evaluation,
}

impl Input {
    const TRACE_START_INDEX: usize = 21;
}

impl Column for Input {
    fn index(&self) -> usize {
        Self::TRACE_START_INDEX + *self as usize
    }
}

enum Output {
    Value,
    Evaluation,
}

impl Output {
    const TRACE_START_INDEX: usize = 23;
}

impl Column for Output {
    fn index(&self) -> usize {
        Self::TRACE_START_INDEX + *self as usize
    }
}

fn instr_zerofier<F: GpuField>(instr: Constraint<F>) -> Constraint<F> {
    let mut accumulator = Constraint::from(F::one());
    for opcode in OpCode::iterator().copied() {
        let opcode: u64 = opcode.into();
        let factor = &instr - F::from(opcode);
        // TODO: mul assign
        accumulator = accumulator * factor;
    }
    accumulator
}

/// returns a polynomial in X that evaluates to 0 in all instructions except
/// for one provided
pub(crate) fn if_not_instr<F: GpuField>(
    instr: OpCode,
    indeterminate: impl Borrow<Constraint<F>>,
) -> Constraint<F> {
    let mut accumulator = Constraint::from(F::one());
    for opcode in OpCode::iterator().copied() {
        if opcode != instr {
            let opcode: u64 = opcode.into();
            let factor = indeterminate.borrow() - F::from(opcode);
            accumulator = accumulator * factor;
        }
    }
    accumulator
}
