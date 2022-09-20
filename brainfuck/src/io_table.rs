use super::table::Table;
use crate::util::interpolate_columns;
use crate::util::lift;
use ark_ff::FftField;
use ark_ff::Field;
use legacy_algebra::number_theory_transform::number_theory_transform;
use legacy_algebra::scale_poly;
use legacy_algebra::Multivariate;
use num_traits::Zero;

const BASE_WIDTH: usize = 1;
const EXTENSION_WIDTH: usize = 2;

struct IoTable<F: Field> {
    num_padded_rows: usize,
    height: Option<usize>,
    matrix: Vec<[F::BasePrimeField; BASE_WIDTH]>,
    extended_matrix: Option<Vec<[F; EXTENSION_WIDTH]>>,
}

impl<F> IoTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    // base column
    pub const VALUE: usize = 0;
    // extension column
    pub const EVALUATION: usize = 1;

    pub fn new() -> Self {
        IoTable {
            num_padded_rows: 0,
            matrix: Vec::new(),
            height: None,
            extended_matrix: None,
        }
    }

    pub fn len(&self) -> usize {
        self.matrix.len() - self.num_padded_rows
    }

    fn height(&self) -> usize {
        self.height.unwrap_or(self.matrix.len())
    }

    fn set_height(&mut self, height: usize) {
        self.height = Some(height);
    }

    pub fn pad(&mut self, n: usize) {
        // TODO: seting length here seems kind of strange
        while self.matrix.len() < n {
            self.matrix.push([F::BasePrimeField::zero()]);
            self.num_padded_rows += 1;
        }
    }

    fn base_boundary_constraints() -> Vec<Multivariate<F>> {
        Vec::new()
    }

    fn extension_boundary_constraints() -> Vec<Multivariate<F>> {
        let variables = Multivariate::variables(EXTENSION_WIDTH);
        vec![variables[Self::EVALUATION].clone() - variables[Self::VALUE].clone()]
    }

    fn base_transition_constraints() -> Vec<Multivariate<F>> {
        Vec::new()
    }

    fn extension_transition_constraints(challenge: F) -> Vec<Multivariate<F>> {
        let variables = Multivariate::<F>::variables(EXTENSION_WIDTH * 2);
        let value = variables[Self::VALUE].clone();
        let evaluation = variables[Self::EVALUATION].clone();
        let value_next = variables[Self::VALUE].clone();
        let evaluation_next = variables[Self::EVALUATION].clone();
        vec![evaluation.clone() * challenge + value_next.clone() - evaluation_next.clone()]
    }

    fn extension_terminal_constraints(&self, challenge: F, terminal: F) -> Vec<Multivariate<F>> {
        let variables = Multivariate::<F>::variables(EXTENSION_WIDTH);
        let offset = challenge.pow([self.num_padded_rows as u64]);
        // In every padded row the running evaluation variable is multiplied by another
        // factor `challenge`. We need to multiply `challenge ^ padding_length` to get
        // the value of the evaluation terminal after all `2^k` rows.
        let actual_terminal = terminal * offset;
        vec![variables[Self::EVALUATION].clone() - actual_terminal]
    }

    fn set_matrix(&mut self, matrix: Vec<[F::BasePrimeField; BASE_WIDTH]>) {
        self.num_padded_rows = 0;
        self.matrix = matrix;
    }

    fn interpolant_degree(&self) -> usize {
        self.matrix.len()
    }

    fn extend(&mut self, challenge: F) {
        // prepare
        let mut extended_matrix = Vec::new();
        let mut io_running_evaluation = F::zero();
        let mut evaluation_terminal = F::zero();

        // loop over all rows
        for i in 0..self.matrix.len() {
            let base_row = self.matrix[i];
            let mut extension_row = [F::zero(); EXTENSION_WIDTH];
            extension_row[..BASE_WIDTH]
                .copy_from_slice(&base_row.map(|v| F::from_base_prime_field(v)));
            io_running_evaluation = io_running_evaluation * challenge + extension_row[Self::VALUE];
            extension_row[Self::EVALUATION] = io_running_evaluation;
            if !self.len().is_zero() && i == self.len() - 1 {
                evaluation_terminal = io_running_evaluation;
            }
            extended_matrix.push(extension_row);
        }

        self.extended_matrix = Some(extended_matrix);
        // TODO: terminal
        // evaluation_terminal
    }

    fn base_lde(&mut self, offset: F::BasePrimeField, codeword_len: usize) -> Vec<Vec<F>> {
        let polynomials = interpolate_columns(&self.matrix, 0);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = scale_poly(&poly, offset).coeffs;
                coefficients.resize(codeword_len, F::BasePrimeField::zero());
                lift(number_theory_transform(&coefficients))
            })
            .collect()
    }

    fn extension_lde(&mut self, offset: F::BasePrimeField, codeword_len: usize) -> Vec<Vec<F>> {
        let extension_rows = self
            .extended_matrix
            .as_ref()
            .unwrap()
            .iter()
            .map(|row| [row[Self::EVALUATION]])
            .collect::<Vec<[F; 1]>>();
        let polynomials = interpolate_columns(&extension_rows, 0);
        // return the codewords
        polynomials
            .into_iter()
            .map(|poly| {
                let mut coefficients = scale_poly(&poly, F::from_base_prime_field(offset)).coeffs;
                coefficients.resize(codeword_len, F::zero());
                number_theory_transform(&coefficients)
            })
            .collect()
    }
}

pub struct OutputTable<F: Field>(IoTable<F>);

impl<F> OutputTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    pub fn new() -> Self {
        OutputTable(IoTable::new())
    }
}

impl<F> Table<F> for OutputTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        self.0.len()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn set_height(&mut self, height: usize) {
        self.0.set_height(height)
    }

    fn pad(&mut self, n: usize) {
        self.0.pad(n)
    }

    fn base_boundary_constraints() -> Vec<Multivariate<F>> {
        IoTable::base_boundary_constraints()
    }

    fn base_transition_constraints() -> Vec<Multivariate<F>> {
        IoTable::base_transition_constraints()
    }

    fn extension_boundary_constraints(_challenges: &[F]) -> Vec<Multivariate<F>> {
        IoTable::extension_boundary_constraints()
    }

    fn extension_transition_constraints(challenges: &[F]) -> Vec<Multivariate<F>> {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let _d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let _gamma = challenges_iter.next().unwrap();
        let delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();
        IoTable::extension_transition_constraints(delta)
    }

    fn extension_terminal_constraints(
        &self,
        challenges: &[F],
        terminals: &[F],
    ) -> Vec<Multivariate<F>> {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let _d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let _gamma = challenges_iter.next().unwrap();
        let delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();

        let mut terminal_iter = terminals.iter().copied();
        let _processor_instruction_permutation_terminal = terminal_iter.next().unwrap();
        let _processor_memory_permutation_terminal = terminal_iter.next().unwrap();
        let _processor_input_evaluation_terminal = terminal_iter.next().unwrap();
        let processor_output_evaluation_terminal = terminal_iter.next().unwrap();
        let _instruction_evaluation_terminal = terminal_iter.next().unwrap();

        self.0
            .extension_terminal_constraints(delta, processor_output_evaluation_terminal)
    }

    fn interpolant_degree(&self) -> usize {
        self.0.interpolant_degree()
    }

    fn set_matrix(&mut self, matrix: Vec<[F::BasePrimeField; BASE_WIDTH]>) {
        self.0.set_matrix(matrix)
    }

    fn extend(&mut self, challenges: &[F], initials: &[F]) {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let _d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let _gamma = challenges_iter.next().unwrap();
        let delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();

        self.0.extend(delta);
    }

    fn base_lde(&mut self, offset: F::BasePrimeField, codeword_len: usize) -> Vec<Vec<F>> {
        println!("output_lde");
        self.0.base_lde(offset, codeword_len)
    }

    fn extension_lde(&mut self, offset: F::BasePrimeField, expansion_factor: usize) -> Vec<Vec<F>> {
        println!("output_lde_ext");
        self.0.extension_lde(offset, expansion_factor)
    }
}

pub struct InputTable<F: Field>(IoTable<F>);

impl<F> InputTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    pub fn new() -> Self {
        InputTable(IoTable::new())
    }
}

impl<F> Table<F> for InputTable<F>
where
    F: Field,
    F::BasePrimeField: FftField,
{
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        self.0.len()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn set_height(&mut self, height: usize) {
        self.0.set_height(height)
    }

    fn pad(&mut self, n: usize) {
        self.0.pad(n)
    }

    fn base_boundary_constraints() -> Vec<Multivariate<F>> {
        IoTable::base_boundary_constraints()
    }

    fn base_transition_constraints() -> Vec<Multivariate<F>> {
        IoTable::base_transition_constraints()
    }

    fn extension_boundary_constraints(_challenges: &[F]) -> Vec<Multivariate<F>> {
        IoTable::extension_boundary_constraints()
    }

    fn extension_transition_constraints(challenges: &[F]) -> Vec<Multivariate<F>> {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let _d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let gamma = challenges_iter.next().unwrap();
        let _delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();
        IoTable::extension_transition_constraints(gamma)
    }

    fn extension_terminal_constraints(
        &self,
        challenges: &[F],
        terminals: &[F],
    ) -> Vec<Multivariate<F>> {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let _d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let gamma = challenges_iter.next().unwrap();
        let _delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();

        let mut terminal_iter = terminals.iter().copied();
        let _processor_instruction_permutation_terminal = terminal_iter.next().unwrap();
        let _processor_memory_permutation_terminal = terminal_iter.next().unwrap();
        let processor_input_evaluation_terminal = terminal_iter.next().unwrap();
        let _processor_output_evaluation_terminal = terminal_iter.next().unwrap();
        let _instruction_evaluation_terminal = terminal_iter.next().unwrap();

        self.0
            .extension_terminal_constraints(gamma, processor_input_evaluation_terminal)
    }

    fn interpolant_degree(&self) -> usize {
        self.0.interpolant_degree()
    }

    fn set_matrix(&mut self, matrix: Vec<[F::BasePrimeField; BASE_WIDTH]>) {
        self.0.set_matrix(matrix)
    }

    fn extend(&mut self, challenges: &[F], _initials: &[F]) {
        let mut challenges_iter = challenges.iter().copied();
        let _a = challenges_iter.next().unwrap();
        let _b = challenges_iter.next().unwrap();
        let _c = challenges_iter.next().unwrap();
        let _d = challenges_iter.next().unwrap();
        let _e = challenges_iter.next().unwrap();
        let _f = challenges_iter.next().unwrap();
        let _alpha = challenges_iter.next().unwrap();
        let _beta = challenges_iter.next().unwrap();
        let gamma = challenges_iter.next().unwrap();
        let _delta = challenges_iter.next().unwrap();
        let _eta = challenges_iter.next().unwrap();

        self.0.extend(gamma)
    }

    fn base_lde(&mut self, offset: F::BasePrimeField, codeword_len: usize) -> Vec<Vec<F>> {
        println!("input_lde");
        self.0.base_lde(offset, codeword_len)
    }

    fn extension_lde(&mut self, offset: F::BasePrimeField, expansion_factor: usize) -> Vec<Vec<F>> {
        println!("input_lde_ext");
        self.0.extension_lde(offset, expansion_factor)
    }
}
