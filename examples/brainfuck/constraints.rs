use crate::tables::Challenge;
use crate::tables::InputBaseColumn;
use crate::tables::InputExtensionColumn;
use crate::tables::InstructionBaseColumn;
use crate::tables::InstructionExtensionColumn;
use crate::tables::MemoryBaseColumn;
use crate::tables::MemoryExtensionColumn;
use crate::tables::OutputBaseColumn;
use crate::tables::OutputExtensionColumn;
use crate::tables::ProcessorBaseColumn;
use crate::tables::ProcessorExtensionColumn;
use crate::vm::OpCode;
use ark_ff::Zero;
use gpu_poly::GpuField;
use ministark::constraint::Challenge as _;
use ministark::Column;
use ministark::Constraint;
use std::borrow::Borrow;

impl ProcessorBaseColumn {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use ProcessorBaseColumn::*;
        let one = F::one();
        let two = one + one;
        let mem_val_is_zero = MemVal.curr() * MemValInv.curr() - one;
        let mut constraints = (Constraint::zero(), Constraint::zero(), Constraint::zero());

        use OpCode::*;
        for instr in OpCode::VALUES {
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
            constraints.0 += &deselector * &instr_constraints.0 * CurrInstr.curr();
            constraints.1 += &deselector * &instr_constraints.1 * CurrInstr.curr();
            constraints.2 += &deselector * &instr_constraints.2 * CurrInstr.curr();
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

impl ProcessorExtensionColumn {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenge::Alpha;
        use Challenge::Beta;
        use Challenge::Delta;
        use Challenge::Gamma;
        use Challenge::A;
        use Challenge::B;
        use Challenge::C;
        use ProcessorBaseColumn::*;
        use ProcessorExtensionColumn::*;

        vec![
            // running product for instruction table permutation
            CurrInstr.curr()
                * (InstructionPermutation.curr()
                    * (Alpha.get_challenge()
                        - A.get_challenge() * Ip.curr()
                        - B.get_challenge() * CurrInstr.curr()
                        - C.get_challenge() * NextInstr.curr())
                    - InstructionPermutation.next())
                + instr_zerofier(CurrInstr.curr())
                    * (InstructionPermutation.curr() - InstructionPermutation.next()),
            // running product for memory table permutation
            CurrInstr.curr()
                * (MemoryPermutation.curr()
                    * (Beta.get_challenge()
                        - Challenge::D.get_challenge() * Cycle.curr()
                        - Challenge::E.get_challenge() * Mp.curr()
                        - Challenge::F.get_challenge() * MemVal.curr())
                    - MemoryPermutation.next())
                * instr_zerofier(CurrInstr.curr())
                * (MemoryPermutation.curr() - MemoryPermutation.next()),
            // running evaluation for input tape
            // TODO: think can remove CurrInstr.curr() from the start here.
            CurrInstr.curr()
                * if_not_instr(OpCode::Read, CurrInstr.curr())
                * (InputEvaluation.next()
                    - Gamma.get_challenge() * InputEvaluation.curr()
                    - MemVal.next())
                + if_instr(OpCode::Read, CurrInstr.curr())
                    * (InputEvaluation.next() - InputEvaluation.curr()),
            // running evaluation for output tape
            CurrInstr.curr()
                * if_not_instr(OpCode::Write, CurrInstr.curr())
                * (OutputEvaluation.next()
                    - OutputEvaluation.curr() * Delta.get_challenge()
                    - MemVal.curr())
                + if_instr(OpCode::Write, CurrInstr.curr())
                    * (OutputEvaluation.next() - OutputEvaluation.curr()),
        ]
    }
}

impl MemoryBaseColumn {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use MemoryBaseColumn::*;
        let one = F::one();
        vec![
            // 1. memory pointer increases by one or zero
            // note: remember table is sorted by memory address
            (Mp.next() - Mp.curr() - one) * (Mp.next() - Mp.curr()),
            //
            // 2. the memory value changes only if (a.) the memory pointer does not increase or
            // (b.) the cycle count increases by one.These constraints are implied by 3.
            //
            // 3. if the memory pointer increases by one, then the memory value must be set to zero
            (Mp.next() - Mp.curr()) * MemVal.next(),
            // 4. dummy has to be zero or one
            (Dummy.next() - one) * Dummy.next(),
            // 5. if dummy is set the memory pointer can not change
            (Mp.next() - Mp.curr()) * Dummy.curr(),
            // 6. if dummy is set the memory value can not change
            (MemVal.next() - MemVal.curr()) * Dummy.curr(),
            // 7. if the memory pointer remains the same, then the cycle has to increase by one
            (Mp.next() - Mp.curr() - one) * (Cycle.next() - Cycle.curr() - one),
        ]
    }
}

impl MemoryExtensionColumn {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenge::Beta;
        use MemoryBaseColumn::*;
        use MemoryExtensionColumn::*;
        // Only progress permutation if dummy != 1
        vec![
            (Permutation.next()
                - Permutation.curr()
                    * (Beta.get_challenge()
                        - Challenge::D.get_challenge() * Cycle.curr()
                        - Challenge::E.get_challenge() * Mp.curr()
                        - Challenge::F.get_challenge() * MemVal.curr()))
                * (Dummy.curr() - F::one())
                + (Permutation.next() - Permutation.curr()) * Dummy.curr(),
        ]
    }
}

impl InstructionBaseColumn {
    pub fn boundary_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use InstructionBaseColumn::*;
        vec![Ip.curr()]
    }

    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use InstructionBaseColumn::*;
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
        ]
    }
}

impl InstructionExtensionColumn {
    pub fn boundary_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenge::A;
        use Challenge::B;
        use Challenge::C;
        use InstructionBaseColumn::*;
        use InstructionExtensionColumn::*;
        vec![
            ProgramEvaluation.curr()
                - A.get_challenge() * Ip.curr()
                - B.get_challenge() * CurrInstr.curr()
                - C.get_challenge() * NextInstr.curr(),
        ]
    }

    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenge::Alpha;
        use Challenge::Eta;
        use Challenge::A;
        use Challenge::B;
        use Challenge::C;
        use InstructionBaseColumn::*;
        use InstructionExtensionColumn::*;
        let one = F::one();
        vec![
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
}

impl InputExtensionColumn {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenge::Gamma;
        use InputBaseColumn::*;
        use InputExtensionColumn::*;
        vec![Evaluation.curr() * Gamma.get_challenge() + Value.next() - Evaluation.next()]
    }
}

impl OutputExtensionColumn {
    pub fn transition_constraints<F: GpuField>() -> Vec<Constraint<F>> {
        use Challenge::Delta;
        use OutputBaseColumn::*;
        use OutputExtensionColumn::*;
        vec![Evaluation.curr() * Delta.get_challenge() + Value.next() - Evaluation.next()]
    }
}

fn instr_zerofier<F: GpuField>(instr: Constraint<F>) -> Constraint<F> {
    let mut accumulator = Constraint::from(F::one());
    for opcode in OpCode::VALUES {
        accumulator *= &instr - F::from(opcode as u64);
    }
    accumulator
}

/// returns a polynomial in X that evaluates to 0 in all instructions except
/// for one provided
fn if_not_instr<F: GpuField>(
    instr: OpCode,
    indeterminate: impl Borrow<Constraint<F>>,
) -> Constraint<F> {
    let mut accumulator = Constraint::from(F::one());
    for opcode in OpCode::VALUES {
        if opcode != instr {
            let factor = indeterminate.borrow() - F::from(opcode as u64);
            accumulator *= factor;
        }
    }
    accumulator
}

fn if_instr<F: GpuField>(
    instr: OpCode,
    indeterminate: impl Borrow<Constraint<F>>,
) -> Constraint<F> {
    indeterminate.borrow() - F::from(instr as u64)
}
