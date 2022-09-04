use super::table::Table;
use crate::table::instr_zerofier;
use crate::OpCode;
use algebra::Felt;
use algebra::Multivariate;
use std::fs::OpenOptions;
use std::mem;
use std::ops::Mul;

pub struct ProcessorTable<E> {
    table: Table<E>,
}

impl<E: Felt> ProcessorTable<E> {
    // Column indices
    pub const CYCLE: usize = 0;
    pub const IP: usize = 1;
    pub const CURR_INSTR: usize = 2;
    pub const NEXT_INSTR: usize = 3;
    pub const MP: usize = 4;
    pub const MEM_VAL: usize = 5;
    pub const MEM_VAL_INV: usize = 6;
    // Extension columns
    pub const INSTRUCTION_PERMUTATION: usize = 7;
    pub const MEMORY_PERMUTATION: usize = 8;
    pub const INPUT_EVALUATION: usize = 9;
    pub const OUTPUT_EVALUATION: usize = 10;

    pub fn new(length: usize, num_randomizers: usize) -> ProcessorTable<E> {
        ProcessorTable {
            table: Table::new(7, 11, length, num_randomizers),
        }
    }

    pub fn pad(&mut self) {
        while !self.table.matrix.len().is_power_of_two() {
            let last_row = self.table.matrix.last().unwrap();
            let mut new_row = vec![E::zero(); self.table.base_width];
            new_row[Self::CYCLE] = last_row[Self::CYCLE] + E::one();
            new_row[Self::IP] = last_row[Self::IP];
            new_row[Self::CURR_INSTR] = E::zero();
            new_row[Self::NEXT_INSTR] = E::zero();
            new_row[Self::MP] = last_row[Self::MP];
            new_row[Self::MEM_VAL] = last_row[Self::MEM_VAL];
            new_row[Self::MEM_VAL_INV] = last_row[Self::MEM_VAL_INV];
            self.table.matrix.push(new_row);
        }
    }

    fn derive_matrix(processor_matrix: Vec<[E; 7]>) {}

    fn base_boundary_constraints(&self) -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(self.table.base_width);
        // All registers except CURR_INSTR and NEXT_INSTR should be zero
        vec![
            variables[Self::CYCLE].clone(),
            variables[Self::IP].clone(),
            variables[Self::MP].clone(),
            variables[Self::MEM_VAL].clone(),
            variables[Self::MEM_VAL_INV].clone(),
        ]
    }

    fn extension_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(11);
        vec![
            variables[Self::CYCLE].clone(),
            variables[Self::IP].clone(),
            variables[Self::MP].clone(),
            variables[Self::MEM_VAL].clone(),
            variables[Self::MEM_VAL_INV].clone(),
            variables[Self::INPUT_EVALUATION].clone(),
            variables[Self::OUTPUT_EVALUATION].clone(),
        ]
    }

    fn if_instr(instr: &OpCode, indeterminate: &Multivariate<E>) -> Multivariate<E> {
        Multivariate::constant(Into::<usize>::into(instr.clone()).into()) - indeterminate.clone()
    }

    /// returns a polynomial in X that evaluates to 0 in all instructions except
    /// for one provided
    fn if_not_instr(instr: &OpCode, indeterminate: &Multivariate<E>) -> Multivariate<E> {
        let mut accumulator = Multivariate::one();
        for opcode in OpCode::iterator() {
            if opcode != instr {
                let factor = indeterminate.clone() - E::from(Into::<usize>::into(opcode.clone()));
                accumulator = accumulator * factor;
            }
        }
        accumulator
    }

    fn transition_constraints(
        cycle: &Multivariate<E>,
        ip: &Multivariate<E>,
        curr_instr: &Multivariate<E>,
        next_instr: &Multivariate<E>,
        mp: &Multivariate<E>,
        mem_val: &Multivariate<E>,
        mem_val_inv: &Multivariate<E>,
        cycle_next: &Multivariate<E>,
        ip_next: &Multivariate<E>,
        curr_instr_next: &Multivariate<E>,
        next_instr_next: &Multivariate<E>,
        mp_next: &Multivariate<E>,
        mem_val_next: &Multivariate<E>,
        mem_val_inv_next: &Multivariate<E>,
    ) -> Vec<Multivariate<E>> {
        let zero = E::zero();
        let one = E::one();
        let two = one + one;
        let mem_val_is_zero = mem_val.clone() * mem_val_inv.clone() - one;
        let mut polynomials = vec![
            Multivariate::zero(),
            Multivariate::zero(),
            Multivariate::zero(),
        ];

        use OpCode::*;
        for instr in OpCode::iterator() {
            // max degree: 4
            let mut instr_polynomials = [
                Multivariate::zero(),
                Multivariate::zero(),
                Multivariate::zero(),
            ];

            match instr {
                IncrementPointer => {
                    instr_polynomials[0] = ip_next.clone() - ip.clone() - one;
                    instr_polynomials[1] = mp_next.clone() - mp.clone() + one;
                }
                DecrementPointer => {
                    instr_polynomials[0] = ip_next.clone() - ip.clone() - one;
                    instr_polynomials[1] = mp_next.clone() - mp.clone() - one;
                }
                Increment => {
                    instr_polynomials[0] = ip_next.clone() - ip.clone() - one;
                    instr_polynomials[1] = mp_next.clone() - mp.clone();
                    instr_polynomials[2] = mem_val_next.clone() - mem_val.clone() - one;
                }
                Decrement => {
                    instr_polynomials[0] = ip_next.clone() - ip.clone() - one;
                    instr_polynomials[1] = mp_next.clone() - mp.clone();
                    instr_polynomials[2] = mem_val_next.clone() - mem_val.clone() + one;
                }
                Write => {
                    instr_polynomials[0] = ip_next.clone() - ip.clone() - one;
                    instr_polynomials[1] = mp_next.clone() - mp.clone();
                }
                Read => {
                    instr_polynomials[0] = ip_next.clone() - ip.clone() - one;
                    instr_polynomials[1] = mp_next.clone() - mp.clone();
                    instr_polynomials[2] = mem_val_next.clone() - mem_val.clone();
                }
                LoopBegin => {
                    instr_polynomials[0] = mem_val.clone() * (ip_next.clone() - ip.clone() - two)
                        + mem_val_is_zero.clone() * (ip_next.clone() - next_instr.clone());
                    instr_polynomials[1] = mp_next.clone() - mp.clone();
                    instr_polynomials[2] = mem_val_next.clone() - mem_val.clone();
                }
                LoopEnd => {
                    instr_polynomials[0] = mem_val_is_zero.clone()
                        * (ip_next.clone() - ip.clone() - two)
                        + mem_val.clone() * (ip_next.clone() - next_instr.clone());
                    instr_polynomials[1] = mp_next.clone() - mp.clone();
                    instr_polynomials[2] = mem_val_next.clone() - mem_val.clone();
                }
            }

            // max degree: 7
            let deselector = Self::if_not_instr(instr, &curr_instr);

            for (polynomial, instr_polynomial) in polynomials.iter_mut().zip(instr_polynomials) {
                // max degree: 7 + 4 = 11
                *polynomial = polynomial.clone() + deselector.clone() * instr_polynomial;
            }
        }

        // cycle independent polynomials
        polynomials.push(cycle_next.clone() - cycle.clone() - one);
        polynomials.push(mem_val.clone() * mem_val_is_zero.clone());
        polynomials.push(mem_val_inv.clone() * mem_val_is_zero.clone());

        polynomials // max degree: 11
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(14);
        // current cycle
        let cycle = variables[Self::CYCLE].clone();
        let ip = variables[Self::IP].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let next_instr = variables[Self::NEXT_INSTR].clone();
        let mp = variables[Self::MP].clone();
        let mem_val = variables[Self::MEM_VAL].clone();
        let mem_val_inv = variables[Self::MEM_VAL_INV].clone();
        // next cycle
        let cycle_next = variables[7 + Self::CYCLE].clone();
        let ip_next = variables[7 + Self::IP].clone();
        let curr_instr_next = variables[7 + Self::CURR_INSTR].clone();
        let next_instr_next = variables[7 + Self::NEXT_INSTR].clone();
        let mp_next = variables[7 + Self::MP].clone();
        let mem_val_next = variables[7 + Self::MEM_VAL].clone();
        let mem_val_inv_next = variables[7 + Self::MEM_VAL_INV].clone();

        Self::transition_constraints(
            &cycle,
            &ip,
            &curr_instr,
            &next_instr,
            &mp,
            &mem_val,
            &mem_val_inv,
            &cycle_next,
            &ip_next,
            &curr_instr_next,
            &next_instr_next,
            &mp_next,
            &mem_val_next,
            &mem_val_inv_next,
        )
    }

    fn extension_transition_constraints(&self, challenges: &[E]) -> Vec<Multivariate<E>> {
        let mut challenges_iter = challenges.iter().copied();
        let a = challenges_iter.next().unwrap();
        let b = challenges_iter.next().unwrap();
        let c = challenges_iter.next().unwrap();
        let d = challenges_iter.next().unwrap();
        let e = challenges_iter.next().unwrap();
        let f = challenges_iter.next().unwrap();
        let alpha = challenges_iter.next().unwrap();
        let beta = challenges_iter.next().unwrap();
        let gamma = challenges_iter.next().unwrap();
        let delta = challenges_iter.next().unwrap();
        let eta = challenges_iter.next().unwrap();

        let variables = Multivariate::<E>::variables(22);
        // current cycle
        let cycle = variables[Self::CYCLE].clone();
        let ip = variables[Self::IP].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let next_instr = variables[Self::NEXT_INSTR].clone();
        let mp = variables[Self::MP].clone();
        let mem_val = variables[Self::MEM_VAL].clone();
        let mem_val_inv = variables[Self::MEM_VAL_INV].clone();
        let instruction_permutation = variables[Self::INSTRUCTION_PERMUTATION].clone();
        let memory_permutation = variables[Self::MEMORY_PERMUTATION].clone();
        let input_evaluation = variables[Self::INPUT_EVALUATION].clone();
        let output_evaluation = variables[Self::OUTPUT_EVALUATION].clone();
        // next cycle
        let cycle_next = variables[11 + Self::CYCLE].clone();
        let ip_next = variables[11 + Self::IP].clone();
        let curr_instr_next = variables[11 + Self::CURR_INSTR].clone();
        let next_instr_next = variables[11 + Self::NEXT_INSTR].clone();
        let mp_next = variables[11 + Self::MP].clone();
        let mem_val_next = variables[11 + Self::MEM_VAL].clone();
        let mem_val_inv_next = variables[11 + Self::MEM_VAL_INV].clone();
        let instruction_permutation_next = variables[11 + Self::INSTRUCTION_PERMUTATION].clone();
        let memory_permutation_next = variables[11 + Self::MEMORY_PERMUTATION].clone();
        let input_evaluation_next = variables[11 + Self::INPUT_EVALUATION].clone();
        let output_evaluation_next = variables[11 + Self::OUTPUT_EVALUATION].clone();

        // Base AIR polynomials
        let mut polynomials = Self::transition_constraints(
            &cycle,
            &ip,
            &curr_instr,
            &next_instr,
            &mp,
            &mem_val,
            &mem_val_inv,
            &cycle_next,
            &ip_next,
            &curr_instr_next,
            &next_instr_next,
            &mp_next,
            &mem_val_next,
            &mem_val_inv_next,
        );
        assert_eq!(polynomials.len(), 6, "expected 6 transition constraints");

        // running product for instruction table permutation
        let instruction_permutation_constraint = curr_instr.clone()
            * (instruction_permutation.clone()
                * (Multivariate::constant(alpha) - ip.clone() * a
                    + curr_instr.clone() * b
                    + next_instr.clone() * c)
                - instruction_permutation_next.clone())
            + instr_zerofier(&curr_instr)
                * (instruction_permutation_next.clone() - instruction_permutation.clone());
        polynomials.push(instruction_permutation_constraint);

        // running product for memory table permutation
        let memory_permutation_constraint = memory_permutation.clone()
            * (Multivariate::constant(beta) - cycle.clone() * d
                + mp.clone() * e
                + mem_val.clone() * f)
            - memory_permutation_next.clone();
        polynomials.push(memory_permutation_constraint);

        // running evaluation for input tape
        // TODO: think can remove curr_instr.clone() from the start here.
        let input_evaluation_constraint = curr_instr.clone()
            * Self::if_not_instr(&OpCode::Read, &curr_instr)
            * (input_evaluation_next.clone()
                - input_evaluation.clone() * gamma
                - mem_val_next.clone())
            + (input_evaluation_next.clone() - input_evaluation.clone())
                * Self::if_instr(&OpCode::Read, &curr_instr);
        polynomials.push(input_evaluation_constraint);

        // running evaluation for output tape
        let output_evaluation_constraint = curr_instr.clone()
            * Self::if_not_instr(&OpCode::Write, &curr_instr)
            * (output_evaluation_next.clone()
                - output_evaluation.clone() * delta
                - mem_val.clone())
            + (output_evaluation_next.clone() - output_evaluation.clone())
                * Self::if_instr(&OpCode::Write, &curr_instr);
        polynomials.push(output_evaluation_constraint);

        assert_eq!(polynomials.len(), 10);
        polynomials
    }

    fn extension_terminal_constraints(challenges: &[E], terminals: &[E]) -> Vec<Multivariate<E>> {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let d = challenges_iter.next().unwrap();
        let e = challenges_iter.next().unwrap();
        let f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let beta = challenges_iter.next().unwrap();
        let _gamma = challenges_iter.next().unwrap();
        let _delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();

        let mut terminal_iter = terminals.iter().copied();
        let processor_instruction_permutation_terminal = terminal_iter.next().unwrap();
        let processor_memory_permutation_terminal = terminal_iter.next().unwrap();
        let processor_input_evaluation_terminal = terminal_iter.next().unwrap();
        let processor_output_evaluation_terminal = terminal_iter.next().unwrap();
        let instruction_evaluation_terminal = terminal_iter.next().unwrap();

        let variables = Multivariate::<E>::variables(22);
        let cycle = variables[Self::CYCLE].clone();
        let mp = variables[Self::MP].clone();
        let mem_val = variables[Self::MEM_VAL].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let instruction_permutation = variables[Self::INSTRUCTION_PERMUTATION].clone();
        let memory_permutation = variables[Self::MEMORY_PERMUTATION].clone();
        let input_evaluation = variables[Self::INPUT_EVALUATION].clone();
        let output_evaluation = variables[Self::OUTPUT_EVALUATION].clone();
        vec![
            // running product for instruction table permutation
            instruction_permutation.clone() - processor_instruction_permutation_terminal,
            // running product for memory table permutation
            // TODO: this is so strange. Can't explain this terminal constraint
            // ...think it's to do with the padding and that the terminal constraints are checked
            // at the end
            (memory_permutation.clone()
                * (Multivariate::constant(beta)
                    - cycle.clone() * d
                    - mp.clone() * e
                    - mem_val.clone() * f)
                - processor_memory_permutation_terminal)
                * curr_instr.clone()
                + (memory_permutation.clone() - processor_memory_permutation_terminal)
                    * instr_zerofier(&curr_instr),
            // running evaluation for input table
            // TODO: why is this one so simple
            input_evaluation.clone() - processor_input_evaluation_terminal,
            // running evaluation for output table
            output_evaluation.clone() - processor_output_evaluation_terminal,
        ]
    }
}
