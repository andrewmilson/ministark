use super::table::Table;
use crate::processor_table::ProcessorTable;
use algebra::Multivariate;
use algebra::PrimeFelt;

pub struct IoTable<E> {
    table: Table<E>,
}

impl<E: PrimeFelt> IoTable<E> {
    const VALUE: usize = 0;
    const EVALUATION: usize = 1;

    pub fn new(length: usize) -> Self {
        IoTable {
            table: Table::new(1, 2, length, 0),
        }
    }

    pub fn pad(&mut self) {
        // TODO: seting length here seems kind of strange
        self.table.length = self.table.matrix.len();
        while !self.table.matrix.len().is_power_of_two() {
            self.table.matrix.push(vec![E::zero()]);
        }
        self.table.height = self.table.matrix.len();
    }

    pub fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        Vec::new()
    }

    pub fn base_transition_constraints() -> Vec<Multivariate<E>> {
        let variables = Multivariate::variables(2);
        vec![variables[Self::EVALUATION].clone() - variables[Self::VALUE].clone()]
    }

    pub fn extension_boundary_constraints() -> Vec<Multivariate<E>> {
        vec![]
    }

    pub fn extension_transition_constraints(challenge: E) -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(4);
        let value = variables[Self::VALUE].clone();
        let evaluation = variables[Self::EVALUATION].clone();
        let value_next = variables[Self::VALUE].clone();
        let evaluation_next = variables[Self::EVALUATION].clone();
        vec![evaluation.clone() * challenge + value_next.clone() - evaluation_next.clone()]
    }

    pub fn extension_terminal_constraints(
        &self,
        challenge: E,
        terminal: E,
    ) -> Vec<Multivariate<E>> {
        let variables = Multivariate::<E>::variables(2);
        let offset = challenge.pow(&[(self.table.height - self.table.length) as u64]);
        // In every padded row the running evaluation variable is multiplied by another
        // factor `challenge`. We need to multiply `challenge ^ padding_length` to get
        // the value of the evaluation terminal after all `2^k` rows.
        let actual_terminal = terminal * offset;
        vec![variables[Self::EVALUATION].clone() - actual_terminal]
    }
}

struct OutputTable<E> {
    io_table: IoTable<E>,
}

impl<E: PrimeFelt> OutputTable<E> {
    pub fn new(length: usize) -> Self {
        OutputTable {
            io_table: IoTable::new(length),
        }
    }

    pub fn pad(&mut self) {
        self.io_table.pad();
    }

    pub fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_boundary_constraints()
    }

    pub fn base_transition_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_transition_constraints()
    }

    pub fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        IoTable::<E>::extension_boundary_constraints()
    }

    pub fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
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

    pub fn extension_terminal_constraints(
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
        let processor_instruction_permutation_terminal = terminal_iter.next().unwrap();
        let processor_memory_permutation_terminal = terminal_iter.next().unwrap();
        let processor_input_evaluation_terminal = terminal_iter.next().unwrap();
        let processor_output_evaluation_terminal = terminal_iter.next().unwrap();
        let instruction_evaluation_terminal = terminal_iter.next().unwrap();

        self.io_table
            .extension_terminal_constraints(delta, processor_output_evaluation_terminal)
    }
}

struct InputTable<E> {
    io_table: IoTable<E>,
}

impl<E: PrimeFelt> InputTable<E> {
    pub fn new(length: usize) -> Self {
        InputTable {
            io_table: IoTable::new(length),
        }
    }

    pub fn pad(&mut self) {
        self.io_table.pad();
    }

    pub fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_boundary_constraints()
    }

    pub fn base_transition_constraints() -> Vec<Multivariate<E>> {
        IoTable::<E>::base_transition_constraints()
    }

    pub fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        IoTable::<E>::extension_boundary_constraints()
    }

    pub fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
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

    pub fn extension_terminal_constraints(
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
        let processor_instruction_permutation_terminal = terminal_iter.next().unwrap();
        let processor_memory_permutation_terminal = terminal_iter.next().unwrap();
        let processor_input_evaluation_terminal = terminal_iter.next().unwrap();
        let processor_output_evaluation_terminal = terminal_iter.next().unwrap();
        let instruction_evaluation_terminal = terminal_iter.next().unwrap();

        self.io_table
            .extension_terminal_constraints(gamma, processor_input_evaluation_terminal)
    }
}
