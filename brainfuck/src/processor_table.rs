use super::table::Table;
use crate::util::if_instr;
use crate::util::if_not_instr;
use crate::util::instr_zerofier;
use crate::util::interpolate_columns;
use crate::util::lift;
use crate::OpCode;
use algebra::fp_u128::BaseFelt;
use algebra::ExtensionOf;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use mini_stark::number_theory_transform::number_theory_transform;
use num_bigint::BigUint;
use std::convert::From;

const BASE_WIDTH: usize = 7;
const EXTENSION_WIDTH: usize = 11;

pub struct ProcessorTable<F, E> {
    num_padded_rows: usize,
    num_randomizers: usize,
    matrix: Vec<[F; BASE_WIDTH]>,
    extended_matrix: Option<Vec<[E; EXTENSION_WIDTH]>>,
    pub instr_permutation_terminal: Option<E>,
    pub memory_permutation_terminal: Option<E>,
    pub input_evaluation_terminal: Option<E>,
    pub output_evaluation_terminal: Option<E>,
}

impl<F: StarkFelt + PrimeFelt, E: Felt<BaseFelt = F> + ExtensionOf<F>> ProcessorTable<F, E> {
    // base columns
    pub const CYCLE: usize = 0;
    pub const IP: usize = 1;
    pub const CURR_INSTR: usize = 2;
    pub const NEXT_INSTR: usize = 3;
    pub const MP: usize = 4;
    pub const MEM_VAL: usize = 5;
    pub const MEM_VAL_INV: usize = 6;
    // extension columns
    pub const INSTRUCTION_PERMUTATION: usize = 7;
    pub const MEMORY_PERMUTATION: usize = 8;
    pub const INPUT_EVALUATION: usize = 9;
    pub const OUTPUT_EVALUATION: usize = 10;

    pub fn new(num_randomizers: usize) -> Self {
        ProcessorTable {
            num_padded_rows: 0,
            num_randomizers,
            matrix: Vec::new(),
            extended_matrix: None,
            instr_permutation_terminal: None,
            memory_permutation_terminal: None,
            input_evaluation_terminal: None,
            output_evaluation_terminal: None,
        }
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
            let deselector = if_not_instr(instr, &curr_instr);

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
}

impl<F: StarkFelt + PrimeFelt, E: Felt<BaseFelt = F> + ExtensionOf<F>> Table<F, E>
    for ProcessorTable<F, E>
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
            let last_row = self.matrix.last().unwrap();
            let mut new_row = [F::zero(); BASE_WIDTH];
            new_row[Self::CYCLE] = last_row[Self::CYCLE] + F::one();
            new_row[Self::IP] = last_row[Self::IP];
            // TODO: nit. may be too verbose. remove
            new_row[Self::CURR_INSTR] = F::zero();
            // TODO: nit. may be too verbose. remove
            new_row[Self::NEXT_INSTR] = F::zero();
            new_row[Self::MP] = last_row[Self::MP];
            new_row[Self::MEM_VAL] = last_row[Self::MEM_VAL];
            new_row[Self::MEM_VAL_INV] = last_row[Self::MEM_VAL_INV];
            self.matrix.push(new_row);
            self.num_padded_rows += 1;
        }
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(BASE_WIDTH);
        // All registers except CURR_INSTR and NEXT_INSTR should be zero
        vec![
            variables[Self::CYCLE].clone(),
            variables[Self::IP].clone(),
            variables[Self::MP].clone(),
            variables[Self::MEM_VAL].clone(),
            variables[Self::MEM_VAL_INV].clone(),
        ]
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

    fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
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
        assert_eq!(polynomials.len(), 6, "unexpected transition constraints");

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
            * if_not_instr(&OpCode::Read, &curr_instr)
            * (input_evaluation_next.clone()
                - input_evaluation.clone() * gamma
                - mem_val_next.clone())
            + (input_evaluation_next.clone() - input_evaluation.clone())
                * if_instr(&OpCode::Read, &curr_instr);
        polynomials.push(input_evaluation_constraint);

        // running evaluation for output tape
        let output_evaluation_constraint = curr_instr.clone()
            * if_not_instr(&OpCode::Write, &curr_instr)
            * (output_evaluation_next.clone()
                - output_evaluation.clone() * delta
                - mem_val.clone())
            + (output_evaluation_next.clone() - output_evaluation.clone())
                * if_instr(&OpCode::Write, &curr_instr);
        polynomials.push(output_evaluation_constraint);

        assert_eq!(polynomials.len(), 10);
        polynomials
    }

    fn extension_terminal_constraints(
        &self,
        challenges: &[E],
        terminals: &[E],
    ) -> Vec<Multivariate<E>> {
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
        let processor_instr_permutation_terminal = terminal_iter.next().unwrap();
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
            instruction_permutation.clone() - processor_instr_permutation_terminal,
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

    fn interpolant_degree(&self) -> usize {
        self.height() + self.num_randomizers
    }

    fn set_matrix(&mut self, matrix: Vec<[F; BASE_WIDTH]>) {
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
        let _eta = challenges_iter.next().unwrap();

        let instr_permutation_initial = initials[0];
        let mem_permutation_initial = initials[1];

        // prepare
        let mut instr_permutation_running_product = instr_permutation_initial;
        let mut mem_permutation_running_product = mem_permutation_initial;
        let mut input_running_evaluation = E::zero();
        let mut output_running_evaluation = E::zero();

        // loop over all rows
        let mut extended_matrix = Vec::new();
        for i in 0..self.matrix.len() {
            let base_row = self.matrix[i];
            let mut extension_row = [E::zero(); EXTENSION_WIDTH];
            extension_row[..BASE_WIDTH].copy_from_slice(&base_row.map(|v| v.into()));

            // Permutations columns
            extension_row[Self::INSTRUCTION_PERMUTATION] = instr_permutation_running_product;
            extension_row[Self::MEMORY_PERMUTATION] = mem_permutation_running_product;
            // if not padding
            if !extension_row[Self::CURR_INSTR].is_zero() {
                instr_permutation_running_product *= alpha
                    - a * extension_row[Self::IP]
                    - b * extension_row[Self::CURR_INSTR]
                    - c * extension_row[Self::NEXT_INSTR];
                mem_permutation_running_product *= beta
                    - d * extension_row[Self::CYCLE]
                    - e * extension_row[Self::MP]
                    - f * extension_row[Self::MEM_VAL];
            }

            // Evaluation columns
            extension_row[Self::INPUT_EVALUATION] = input_running_evaluation;
            extension_row[Self::OUTPUT_EVALUATION] = output_running_evaluation;
            let curr_instr = Into::<BigUint>::into(base_row[Self::CURR_INSTR]);
            let read_instr = BigUint::from(Into::<usize>::into(OpCode::Read));
            let write_instr = BigUint::from(Into::<usize>::into(OpCode::Write));
            if curr_instr == read_instr {
                let next_row = self.matrix[i + 1];
                let input_val: E = next_row[Self::MEM_VAL].into();
                input_running_evaluation = input_running_evaluation * gamma + input_val;
            } else if curr_instr == write_instr {
                let next_row = self.matrix[i + 1];
                let output_val: E = next_row[Self::MEM_VAL].into();
                output_running_evaluation = output_running_evaluation * delta + output_val;
            }

            extended_matrix.push(extension_row);
        }

        self.extended_matrix = Some(extended_matrix);
        self.instr_permutation_terminal = Some(instr_permutation_running_product);
        self.memory_permutation_terminal = Some(mem_permutation_running_product);
        self.input_evaluation_terminal = Some(input_running_evaluation);
        self.output_evaluation_terminal = Some(output_running_evaluation);
    }

    fn base_lde(&mut self, offset: F, codeword_len: usize) -> Vec<Vec<E>> {
        println!("proc_lde");
        let polynomials = interpolate_columns(&self.matrix, self.num_randomizers);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = poly.scale(offset).coefficients;
                coefficients.resize(codeword_len, F::zero());
                let coef = number_theory_transform(&coefficients);
                lift(coef)
            })
            .collect()
    }

    fn extension_lde(&mut self, offset: F, codeword_len: usize) -> Vec<Vec<E>> {
        println!("proc_lde_ext");
        let extension_rows = self
            .extended_matrix
            .as_ref()
            .unwrap()
            .iter()
            .map(|row| {
                [
                    row[Self::INSTRUCTION_PERMUTATION],
                    row[Self::MEMORY_PERMUTATION],
                    row[Self::INPUT_EVALUATION],
                    row[Self::OUTPUT_EVALUATION],
                ]
            })
            .collect::<Vec<[E; 4]>>();
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
