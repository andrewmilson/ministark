use super::table::Table;
use crate::table::instr_zerofier;
use crate::OpCode;
use algebra::Felt;
use algebra::Multivariate;
use std::mem;
use std::ops::Mul;

pub struct InstructionTable<E> {
    table: Table<E>,
}

impl<E: Felt> InstructionTable<E> {
    // Column indices
    const IP: usize = 0;
    const CURR_INSTR: usize = 1;
    const NEXT_INSTR: usize = 2;
    // Extension columns
    const PROCESSOR_PERMUTATION: usize = 3;
    const PROGRAM_EVALUATION: usize = 4;

    fn new(length: usize, num_randomizers: usize) -> Self {
        InstructionTable {
            table: Table::new(3, 5, length, num_randomizers),
        }
    }

    fn pad(&mut self) {
        while !self.table.matrix.len().is_power_of_two() {
            let mut new_row = vec![E::zero(); self.table.base_width];
            new_row[Self::IP] = self.table.matrix.last().unwrap()[Self::IP];
            new_row[Self::CURR_INSTR] = E::zero();
            new_row[Self::NEXT_INSTR] = E::zero();
            self.table.matrix.push(new_row);
        }
    }

    fn transition_constraints(
        ip: &Multivariate<E>,
        curr_instr: &Multivariate<E>,
        next_instr: &Multivariate<E>,
        ip_next: &Multivariate<E>,
        curr_instr_next: &Multivariate<E>,
        next_instr_next: &Multivariate<E>,
    ) -> Vec<Multivariate<E>> {
        let one = E::one();
        vec![
            // instruction pointer increases by 0 or 1
            (ip_next.clone() - ip.clone() - one) * (ip_next.clone() - ip.clone()),
            // if address increases the next instruction in the current row must equal the current
            // instruction in the next row
            (ip_next.clone() - ip.clone()) * (next_instr.clone() - curr_instr_next.clone()),
            // if address is the same, then current instruction is also
            (ip_next.clone() - ip.clone() - one) * (curr_instr_next.clone() - curr_instr.clone()),
            // if address is the same, then next instruction is also
            (ip_next.clone() - ip.clone() - one) * (next_instr_next.clone() - next_instr.clone()),
        ]
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(14);
        let ip = variables[Self::IP].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let next_instr = variables[Self::NEXT_INSTR].clone();
        let ip_next = variables[3 + Self::IP].clone();
        let curr_instr_next = variables[3 + Self::CURR_INSTR].clone();
        let next_instr_next = variables[3 + Self::NEXT_INSTR].clone();
        Self::transition_constraints(
            &ip,
            &curr_instr,
            &next_instr,
            &ip_next,
            &curr_instr_next,
            &next_instr_next,
        )
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::variables(3);
        // address starts at zero
        vec![variables[Self::IP].clone()]
    }

    fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
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

        let variables = Multivariate::<E>::variables(10);
        let ip = variables[Self::IP].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let next_instr = variables[Self::NEXT_INSTR].clone();
        let processor_permutation = variables[Self::PROCESSOR_PERMUTATION].clone();
        let program_evaluation = variables[Self::PROGRAM_EVALUATION].clone();
        let ip_next = variables[5 + Self::IP].clone();
        let curr_instr_next = variables[5 + Self::CURR_INSTR].clone();
        let next_instr_next = variables[5 + Self::NEXT_INSTR].clone();
        let processor_permutation_next = variables[5 + Self::PROCESSOR_PERMUTATION].clone();
        let program_evaluation_next = variables[5 + Self::PROGRAM_EVALUATION].clone();

        let mut polynomials = Self::transition_constraints(
            &ip,
            &curr_instr,
            &next_instr,
            &ip_next,
            &curr_instr_next,
            &next_instr_next,
        );

        let one = E::one();

        // - processer permutation changes correctly if ip changes
        // - processor permutation doesn't change if `curr_instr=0` i.e. padding
        // - processor permutation doesn't change if `ip` stays the same
        let processor_permutation_constraint = curr_instr.clone()
            * (ip.clone() - ip_next.clone() - one)
            * (processor_permutation_next.clone()
                - processor_permutation.clone()
                    * (ip_next.clone() * a
                        + curr_instr_next.clone() * b
                        + next_instr_next.clone() * c
                        - alpha))
            + instr_zerofier(&curr_instr)
                * (processor_permutation_next.clone() - processor_permutation.clone())
            + (ip.clone() - ip_next.clone())
                * (processor_permutation.clone() - processor_permutation_next.clone());
        polynomials.push(processor_permutation_constraint);

        // - no evaluation change if `ip` remains the same
        // - evaluation change if `ip` changes
        let program_evaluation_constraint = (ip_next.clone() - ip.clone() - one)
            * (program_evaluation_next.clone() - program_evaluation.clone())
            + (ip_next.clone() - ip.clone())
                * (program_evaluation_next.clone()
                    - program_evaluation.clone() * eta
                    - ip_next.clone() * a
                    - curr_instr_next.clone() * b
                    - next_instr_next.clone() * c);
        polynomials.push(program_evaluation_constraint);

        polynomials
    }

    fn extension_terminal_constraints(challenges: &[E], terminals: &[E]) -> Vec<Multivariate<E>> {
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

        let mut terminal_iter = terminals.iter().copied();
        let processor_instruction_permutation_terminal = terminal_iter.next().unwrap();
        let processor_memory_permutation_terminal = terminal_iter.next().unwrap();
        let processor_input_evaluation_terminal = terminal_iter.next().unwrap();
        let processor_output_evaluation_terminal = terminal_iter.next().unwrap();
        let instruction_evaluation_terminal = terminal_iter.next().unwrap();

        let variables = Multivariate::<E>::variables(5);

        vec![
            variables[Self::PROCESSOR_PERMUTATION].clone()
                - processor_instruction_permutation_terminal,
            variables[Self::PROGRAM_EVALUATION].clone() - instruction_evaluation_terminal,
        ]
    }
}
