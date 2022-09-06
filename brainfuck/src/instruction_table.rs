use super::table::Table;
use crate::util::instr_zerofier;
use algebra::Felt;
use algebra::Multivariate;

const BASE_WIDTH: usize = 3;
const EXTENSION_WIDTH: usize = 5;

pub struct InstructionTable<E> {
    num_padded_rows: usize,
    num_randomizers: usize,
    matrix: Vec<[E; BASE_WIDTH]>,
}

impl<E: Felt> InstructionTable<E> {
    // base columns
    const IP: usize = 0;
    const CURR_INSTR: usize = 1;
    const NEXT_INSTR: usize = 2;
    // extension columns
    const PROCESSOR_PERMUTATION: usize = 3;
    const PROGRAM_EVALUATION: usize = 4;

    pub fn new(num_randomizers: usize) -> Self {
        InstructionTable {
            num_padded_rows: 0,
            num_randomizers,
            matrix: Vec::new(),
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
}

impl<E: Felt> Table<E> for InstructionTable<E> {
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        todo!()
    }

    fn pad(&mut self, n: usize) {
        while self.matrix.len() < n {
            let mut new_row = [E::zero(); BASE_WIDTH];
            new_row[Self::IP] = self.matrix.last().unwrap()[Self::IP];
            new_row[Self::CURR_INSTR] = E::zero();
            new_row[Self::NEXT_INSTR] = E::zero();
            self.matrix.push(new_row);
            self.num_padded_rows += 1;
        }
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::variables(3);
        // address starts at zero
        vec![variables[Self::IP].clone()]
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

    fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        todo!()
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
                    * (Multivariate::constant(alpha)
                        - ip_next.clone() * a
                        - curr_instr_next.clone() * b
                        - next_instr_next.clone() * c))
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

    fn max_degree(&self) -> usize {
        todo!()
    }

    fn set_matrix(&mut self, matrix: Vec<[E; Self::BASE_WIDTH]>) {
        self.num_padded_rows = 0;
        self.matrix = matrix;
    }
}
