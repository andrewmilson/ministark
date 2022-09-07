use super::table::Table;
use algebra::Multivariate;
use algebra::PrimeFelt;

const BASE_WIDTH: usize = 1;
const EXTENSION_WIDTH: usize = 2;

struct IoTable<E> {
    num_padded_rows: usize,
    matrix: Vec<[E; BASE_WIDTH]>,
}

impl<E: PrimeFelt> IoTable<E> {
    // base column
    const VALUE: usize = 0;
    // extension column
    const EVALUATION: usize = 1;

    pub fn new() -> Self {
        IoTable {
            num_padded_rows: 0,
            matrix: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.matrix.len() - self.num_padded_rows
    }

    fn height(&self) -> usize {
        self.matrix.len()
    }

    pub fn pad(&mut self, n: usize) {
        // TODO: seting length here seems kind of strange
        while self.matrix.len() < n {
            self.matrix.push([E::zero()]);
            self.num_padded_rows += 1;
        }
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        Vec::new()
    }

    fn extension_boundary_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::variables(2);
        vec![variables[Self::EVALUATION].clone() - variables[Self::VALUE].clone()]
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        Vec::new()
    }

    fn extension_transition_constraints(challenge: E) -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(4);
        let value = variables[Self::VALUE].clone();
        let evaluation = variables[Self::EVALUATION].clone();
        let value_next = variables[Self::VALUE].clone();
        let evaluation_next = variables[Self::EVALUATION].clone();
        vec![evaluation.clone() * challenge + value_next.clone() - evaluation_next.clone()]
    }

    fn extension_terminal_constraints(&self, challenge: E, terminal: E) -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(2);
        let offset = challenge.pow(&[self.num_padded_rows as u64]);
        // In every padded row the running evaluation variable is multiplied by another
        // factor `challenge`. We need to multiply `challenge ^ padding_length` to get
        // the value of the evaluation terminal after all `2^k` rows.
        let actual_terminal = terminal * offset;
        vec![variables[Self::EVALUATION].clone() - actual_terminal]
    }

    fn set_matrix(&mut self, matrix: Vec<[E; BASE_WIDTH]>) {
        self.num_padded_rows = 0;
        self.matrix = matrix;
    }

    fn interpolant_degree(&self) -> usize {
        self.matrix.len()
    }
}

pub struct OutputTable<E>(IoTable<E>);

impl<E: PrimeFelt> OutputTable<E> {
    pub fn new() -> Self {
        OutputTable(IoTable::new())
    }
}

impl<E: PrimeFelt> Table<E> for OutputTable<E> {
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        self.0.len()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn pad(&mut self, n: usize) {
        self.0.pad(n)
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_boundary_constraints()
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_transition_constraints()
    }

    fn extension_boundary_constraints(_challenges: &[E]) -> Vec<Multivariate<E>> {
        IoTable::<E>::extension_boundary_constraints()
    }

    fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
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
        IoTable::<E>::extension_transition_constraints(delta)
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

    fn set_matrix(&mut self, matrix: Vec<[E; BASE_WIDTH]>) {
        self.0.set_matrix(matrix)
    }
}

pub struct InputTable<E>(IoTable<E>);

impl<E: PrimeFelt> InputTable<E> {
    pub fn new() -> Self {
        InputTable(IoTable::new())
    }
}

impl<E: PrimeFelt> Table<E> for InputTable<E> {
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        self.0.len()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn pad(&mut self, n: usize) {
        self.0.pad(n)
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_boundary_constraints()
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_transition_constraints()
    }

    fn extension_boundary_constraints(_challenges: &[E]) -> Vec<Multivariate<E>> {
        IoTable::<E>::extension_boundary_constraints()
    }

    fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
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
        IoTable::<E>::extension_transition_constraints(gamma)
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

    fn set_matrix(&mut self, matrix: Vec<[E; BASE_WIDTH]>) {
        self.0.set_matrix(matrix)
    }
}
