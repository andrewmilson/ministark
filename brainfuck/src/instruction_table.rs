use super::table::Table;
use crate::util::instr_zerofier;
use crate::util::interpolate_columns;
use crate::util::lift;
use algebra::ExtensionOf;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use mini_stark::number_theory_transform::number_theory_transform;

const BASE_WIDTH: usize = 3;
const EXTENSION_WIDTH: usize = 5;

pub struct InstructionTable<F, E> {
    num_padded_rows: usize,
    num_randomizers: usize,
    matrix: Vec<[F; BASE_WIDTH]>,
    extended_matrix: Option<Vec<[E; EXTENSION_WIDTH]>>,
    pub permutation_terminal: Option<E>,
    pub evaluation_terminal: Option<E>,
}

impl<F: StarkFelt + PrimeFelt, E: Felt<BaseFelt = F> + ExtensionOf<F>> InstructionTable<F, E> {
    // base columns
    pub const IP: usize = 0;
    pub const CURR_INSTR: usize = 1;
    pub const NEXT_INSTR: usize = 2;
    // extension columns
    pub const PROCESSOR_PERMUTATION: usize = 3;
    pub const PROGRAM_EVALUATION: usize = 4;

    pub fn new(num_randomizers: usize) -> Self {
        InstructionTable {
            num_padded_rows: 0,
            num_randomizers,
            matrix: Vec::new(),
            extended_matrix: None,
            permutation_terminal: None,
            evaluation_terminal: None,
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

impl<F: StarkFelt + PrimeFelt, E: Felt<BaseFelt = F> + ExtensionOf<F>> Table<F, E>
    for InstructionTable<F, E>
{
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        self.matrix.len() - self.num_padded_rows
    }

    fn height(&self) -> usize {
        self.matrix.len()
    }

    fn pad(&mut self, n: usize) {
        while self.matrix.len() < n {
            let mut new_row = [F::zero(); BASE_WIDTH];
            new_row[Self::IP] = self.matrix.last().unwrap()[Self::IP];
            new_row[Self::CURR_INSTR] = F::zero();
            new_row[Self::NEXT_INSTR] = F::zero();
            self.matrix.push(new_row);
            self.num_padded_rows += 1;
        }
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::variables(BASE_WIDTH);
        // address starts at zero
        vec![variables[Self::IP].clone()]
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(BASE_WIDTH * 2);
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
        let mut challenges_iter = challenges.iter().copied();
        let a = challenges_iter.next().unwrap();
        let b = challenges_iter.next().unwrap();
        let c = challenges_iter.next().unwrap();
        let d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let _gamma = challenges_iter.next().unwrap();
        let _delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();

        let variables = Multivariate::<E>::variables(EXTENSION_WIDTH);
        let ip = variables[Self::IP].clone();
        let curr_instr = variables[Self::CURR_INSTR].clone();
        let next_instr = variables[Self::NEXT_INSTR].clone();
        let evaluation = variables[Self::PROGRAM_EVALUATION].clone();

        vec![
            ip.clone(),
            evaluation.clone() - ip.clone() * a - curr_instr.clone() * b - next_instr.clone() * c,
        ]
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

        let variables = Multivariate::<E>::variables(EXTENSION_WIDTH * 2);
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

    fn extension_terminal_constraints(
        &self,
        challenges: &[E],
        terminals: &[E],
    ) -> Vec<Multivariate<E>> {
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

        let variables = Multivariate::<E>::variables(EXTENSION_WIDTH);

        vec![
            variables[Self::PROCESSOR_PERMUTATION].clone()
                - processor_instruction_permutation_terminal,
            variables[Self::PROGRAM_EVALUATION].clone() - instruction_evaluation_terminal,
        ]
    }

    fn interpolant_degree(&self) -> usize {
        self.matrix.len() - self.num_padded_rows
    }

    fn set_matrix(&mut self, matrix: Vec<[F; Self::BASE_WIDTH]>) {
        self.num_padded_rows = 0;
        self.matrix = matrix;
    }

    fn extend(&mut self, challenges: &[E], initials: &[E]) {
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

        let instr_permutation_initial = initials[0];
        let mem_permutation_initial = initials[1];

        // prepare
        let mut permutation_running_product = instr_permutation_initial;
        let mut evaluation_running_sum = E::zero();
        let mut previous_address = -F::one();

        let mut extended_matrix = Vec::new();
        let mut num_padded_rows = 0usize;
        for i in 0..self.matrix.len() {
            let base_row = self.matrix[i];
            let mut extension_row = [E::zero(); EXTENSION_WIDTH];
            extension_row[..BASE_WIDTH].copy_from_slice(&base_row.map(|v| v.into()));
            if extension_row[Self::CURR_INSTR].is_zero() {
                num_padded_rows += 1;
            } else if i > 0 && base_row[Self::IP] == self.matrix[i - 1][Self::IP] {
                // permutation argument
                // update running product
                // make sure new row is not padding
                // and that the instruction address didn't just change
                permutation_running_product *= alpha
                    - extension_row[Self::IP] * a
                    - extension_row[Self::CURR_INSTR] * b
                    - extension_row[Self::NEXT_INSTR] * c;
            }
            extension_row[Self::PROCESSOR_PERMUTATION] = permutation_running_product;
            // evaluation argument
            if base_row[Self::IP] != previous_address {
                evaluation_running_sum = eta * evaluation_running_sum
                    + extension_row[Self::IP] * a
                    + extension_row[Self::CURR_INSTR] * b
                    + extension_row[Self::NEXT_INSTR] * c;
            }
            extension_row[Self::PROGRAM_EVALUATION] = evaluation_running_sum;
            previous_address = base_row[Self::IP];
            extended_matrix.push(extension_row);
        }

        self.extended_matrix = Some(extended_matrix);
        self.permutation_terminal = Some(permutation_running_product);
        self.evaluation_terminal = Some(evaluation_running_sum);
    }

    fn base_lde(&mut self, offset: F, codeword_len: usize) -> Vec<Vec<E>> {
        println!("instr_lde");
        let polynomials = interpolate_columns(&self.matrix, self.num_randomizers);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = poly.scale(offset).coefficients;
                coefficients.resize(codeword_len, F::zero());
                lift(number_theory_transform(&coefficients))
            })
            .collect()
    }

    fn extension_lde(&mut self, offset: F, codeword_len: usize) -> Vec<Vec<E>> {
        println!("instr_lde_ext");
        let extension_rows = self
            .extended_matrix
            .as_ref()
            .unwrap()
            .iter()
            .map(|row| {
                [
                    row[Self::PROCESSOR_PERMUTATION],
                    row[Self::PROGRAM_EVALUATION],
                ]
            })
            .collect::<Vec<[E; 2]>>();
        let polynomials = interpolate_columns(&extension_rows, self.num_randomizers);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = poly.scale(offset.into()).coefficients;
                coefficients.resize(codeword_len, E::zero());
                let coef = number_theory_transform(&coefficients);
                lift(coef)
            })
            .collect()
    }
}
