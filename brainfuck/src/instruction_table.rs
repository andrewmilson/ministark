use super::table::Table;
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

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(14);
        let ip = variables[Self::IP].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let next_instr = variables[Self::NEXT_INSTR].clone();
        let ip_next = variables[3 + Self::IP].clone();
        let curr_instr_next = variables[3 + Self::CURR_INSTR].clone();
        let next_instr_next = variables[3 + Self::NEXT_INSTR].clone();

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

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::variables(3);
        // address starts at zero
        vec![variables[Self::IP].clone()]
    }
}
