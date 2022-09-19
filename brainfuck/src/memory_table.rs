use super::table::Table;
use crate::processor_table::ProcessorTable;
use crate::util::interpolate_columns;
use crate::util::lift;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use legacy_algebra::number_theory_transform::number_theory_transform;
use legacy_algebra::Multivariate;

const BASE_WIDTH: usize = 4;
const EXTENSION_WIDTH: usize = 5;

pub struct MemoryTable<F: Field> {
    num_padded_rows: usize,
    num_randomizers: usize,
    matrix: Vec<[F::BasePrimeField; BASE_WIDTH]>,
    extended_matrix: Option<Vec<[F; EXTENSION_WIDTH]>>,
}

impl<F> MemoryTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    // base columns
    pub const CYCLE: usize = 0;
    pub const MP: usize = 1;
    pub const MEM_VAL: usize = 2;
    pub const DUMMY: usize = 3;
    // extension columns
    pub const PERMUTATION: usize = 4;

    pub fn new(num_randomizers: usize) -> Self {
        MemoryTable {
            num_padded_rows: 0,
            num_randomizers,
            matrix: Vec::new(),
            extended_matrix: None,
        }
    }

    fn transition_constraints(
        cycle: &Multivariate<F>,
        mp: &Multivariate<F>,
        mem_val: &Multivariate<F>,
        dummy: &Multivariate<F>,
        cycle_next: &Multivariate<F>,
        mp_next: &Multivariate<F>,
        mem_val_next: &Multivariate<F>,
        dummy_next: &Multivariate<F>,
    ) -> Vec<Multivariate<F>> {
        let one = F::one();
        vec![
            // 1. memory pointer increases by one or zero
            // note: remember table is sorted by memory address
            (mp_next.clone() - mp.clone() - one) * (mp_next.clone() - mp.clone()),
            //
            // 2. the memory value changes only if (a.) the memory pointer does not increase or
            // (b.) the cycle count increases by one.These constraints are implied by 3.
            //
            // 3. if the memory pointer increases by one, then the memory value must be set to zero
            (mp_next.clone() - mp.clone()) * mem_val_next.clone(),
            // 4. dummy has to be zero or one
            (dummy_next.clone() - one) * dummy_next.clone(),
            // 5. if dummy is set the memory pointer can not change
            (mp_next.clone() - mp.clone()) * dummy.clone(),
            // 6. if dummy is set the memory value can not change
            (mem_val_next.clone() - mem_val.clone()) * dummy.clone(),
            // 7. if the memory pointer remains the same, then the cycle has to increase by one
            (mp_next.clone() - mp.clone() - one) * (cycle_next.clone() - cycle.clone() - one),
        ]
    }

    /// Outputs an unpadded but interweaved matrix
    pub fn derive_matrix(
        processor_matrix: &[[F::BasePrimeField; 7]],
    ) -> Vec<[F::BasePrimeField; BASE_WIDTH]> {
        // copy unpadded rows and sort
        // TODO: sorted by IP and then CYCLE. Check to see if processor table sorts by
        // cycle.
        let mut matrix = processor_matrix
            .iter()
            .filter_map(|row| {
                if row[ProcessorTable::<F>::CURR_INSTR].is_zero() {
                    None
                } else {
                    Some([
                        row[ProcessorTable::<F>::CYCLE],
                        row[ProcessorTable::<F>::MP],
                        row[ProcessorTable::<F>::MEM_VAL],
                        F::BasePrimeField::zero(), // dummy=no
                    ])
                }
            })
            .collect::<Vec<[F::BasePrimeField; 4]>>();
        matrix.sort_by_key(|row| row[Self::MP]);

        // insert dummy rows for smooth clk jumps
        for i in 0..matrix.len() - 1 {
            let curr_row = &matrix[i];
            let next_row = &matrix[i + 1];
            if curr_row[Self::MP] == next_row[Self::MP]
                && curr_row[Self::CYCLE] + F::BasePrimeField::one() != next_row[Self::CYCLE]
            {
                matrix.insert(
                    i + 1,
                    [
                        curr_row[Self::CYCLE] + F::BasePrimeField::one(),
                        curr_row[Self::MP],
                        curr_row[Self::MEM_VAL],
                        F::BasePrimeField::one(), // dummy=yes
                    ],
                )
            }
        }

        matrix
    }
}

impl<F> Table<F> for MemoryTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
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
            self.matrix.push([
                last_row[Self::CYCLE] + F::BasePrimeField::one(),
                last_row[Self::MP],
                last_row[Self::MEM_VAL],
                F::BasePrimeField::one(), // dummy=yes
            ]);
            self.num_padded_rows += 1;
        }
    }

    fn base_boundary_constraints() -> Vec<Multivariate<F>> {
        let variables = Multivariate::<F>::variables(BASE_WIDTH);
        vec![
            variables[Self::CYCLE].clone(),
            variables[Self::MP].clone(),
            variables[Self::MEM_VAL].clone(),
        ]
    }

    fn base_transition_constraints() -> Vec<Multivariate<F>> {
        let variables = Multivariate::<F>::variables(BASE_WIDTH * 2);
        let cycle = variables[Self::CYCLE].clone();
        let mp = variables[Self::MP].clone();
        let mem_val = variables[Self::MEM_VAL].clone();
        let dummy = variables[Self::DUMMY].clone();
        let cycle_next = variables[BASE_WIDTH + Self::CYCLE].clone();
        let mp_next = variables[BASE_WIDTH + Self::MP].clone();
        let mem_val_next = variables[BASE_WIDTH + Self::MEM_VAL].clone();
        let dummy_next = variables[BASE_WIDTH + Self::DUMMY].clone();
        Self::transition_constraints(
            &cycle,
            &mp,
            &mem_val,
            &dummy,
            &cycle_next,
            &mp_next,
            &mem_val_next,
            &dummy_next,
        )
    }

    fn extension_boundary_constraints(challenges: &[F]) -> Vec<Multivariate<F>> {
        let variables = Multivariate::<F>::variables(EXTENSION_WIDTH);
        vec![
            variables[Self::CYCLE].clone(),
            variables[Self::MP].clone(),
            variables[Self::MEM_VAL].clone(),
            // TODO: why is this not included?
            // variables[Self::PERMUTATION].clone() - F::one(),
        ]
    }

    fn extension_transition_constraints(challenges: &[F]) -> Vec<Multivariate<F>> {
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

        let variables = Multivariate::<F>::variables(EXTENSION_WIDTH * 2);
        let cycle = variables[Self::CYCLE].clone();
        let mp = variables[Self::MP].clone();
        let mem_val = variables[Self::MEM_VAL].clone();
        let dummy = variables[Self::DUMMY].clone();
        let permutation = variables[Self::PERMUTATION].clone();
        let cycle_next = variables[EXTENSION_WIDTH + Self::CYCLE].clone();
        let mp_next = variables[EXTENSION_WIDTH + Self::MP].clone();
        let mem_val_next = variables[EXTENSION_WIDTH + Self::MEM_VAL].clone();
        let dummy_next = variables[EXTENSION_WIDTH + Self::DUMMY].clone();
        let permutation_next = variables[EXTENSION_WIDTH + Self::PERMUTATION].clone();

        let mut polynomials = Self::transition_constraints(
            &cycle,
            &mp,
            &mem_val,
            &dummy,
            &cycle_next,
            &mp_next,
            &mem_val_next,
            &dummy_next,
        );

        let permutation_constraint = (permutation_next.clone()
            - permutation.clone()
                * (Multivariate::constant(beta)
                    - cycle.clone() * d
                    - mp.clone() * e
                    - mem_val.clone() * f))
            * (dummy.clone() - F::one())
            + (permutation_next.clone() - permutation.clone()) * dummy.clone();
        polynomials.push(permutation_constraint);

        polynomials
    }

    fn extension_terminal_constraints(
        &self,
        challenges: &[F],
        terminals: &[F],
    ) -> Vec<Multivariate<F>> {
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

        let variables = Multivariate::<F>::variables(EXTENSION_WIDTH);
        let cycle = variables[Self::CYCLE].clone();
        let mp = variables[Self::MP].clone();
        let mem_val = variables[Self::MEM_VAL].clone();
        let dummy = variables[Self::DUMMY].clone();
        let permutation = variables[Self::PERMUTATION].clone();

        vec![
            (permutation.clone()
                * (Multivariate::constant(beta)
                    - cycle.clone() * d
                    - mp.clone() * e
                    - mem_val.clone() * f)
                - processor_memory_permutation_terminal)
                * (dummy.clone() - F::one())
                + (permutation.clone() - processor_memory_permutation_terminal) * dummy.clone(),
        ]
    }

    fn interpolant_degree(&self) -> usize {
        self.matrix.len() + self.num_randomizers
    }

    fn set_matrix(&mut self, matrix: Vec<[F::BasePrimeField; BASE_WIDTH]>) {
        self.num_padded_rows = 0;
        self.matrix = matrix;
    }

    fn extend(&mut self, challenges: &[F], initials: &[F]) {
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
        let mut mem_permutation_running_product = mem_permutation_initial;

        // loop over all rows
        let mut extended_matrix = Vec::new();
        for base_row in &self.matrix {
            let mut extension_row = [F::zero(); EXTENSION_WIDTH];
            extension_row[..BASE_WIDTH].copy_from_slice(&base_row.map(F::from_base_prime_field));
            extension_row[Self::PERMUTATION] = mem_permutation_running_product;
            if extension_row[Self::DUMMY].is_zero() {
                mem_permutation_running_product *= beta
                    - extension_row[Self::CYCLE] * d
                    - extension_row[Self::MP] * e
                    - extension_row[Self::MEM_VAL] * f;
            }
            extended_matrix.push(extension_row);
        }

        self.extended_matrix = Some(extended_matrix);
        // TODO: terminal
        // self.permutation_terminal = mem_permutation_running_product;
    }

    fn base_lde(&mut self, offset: F::BasePrimeField, codeword_len: usize) -> Vec<Vec<F>> {
        println!("mem_lde");
        let polynomials = interpolate_columns(&self.matrix, self.num_randomizers);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = poly.scale(offset).coefficients;
                coefficients.resize(codeword_len, F::BasePrimeField::zero());
                lift(number_theory_transform(&coefficients))
            })
            .collect()
    }

    fn extension_lde(&mut self, offset: F::BasePrimeField, codeword_len: usize) -> Vec<Vec<F>> {
        println!("mem_lde_ext");
        let extension_rows = self
            .extended_matrix
            .as_ref()
            .unwrap()
            .iter()
            .map(|row| [row[Self::PERMUTATION]])
            .collect::<Vec<[F; 1]>>();
        let polynomials = interpolate_columns(&extension_rows, self.num_randomizers);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = poly.scale(F::from_base_prime_field(offset)).coefficients;
                coefficients.resize(codeword_len, F::zero());
                number_theory_transform(&coefficients)
            })
            .collect()
    }
}
