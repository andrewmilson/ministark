use super::table::Table;
use crate::processor_table::ProcessorTable;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;

const BASE_WIDTH: usize = 1;
const EXTENSION_WIDTH: usize = 2;

pub struct OutputTable<E> {
    num_padded_rows: usize,
    num_randomizers: usize,
    matrix: Vec<[E; BASE_WIDTH]>,
}

impl<E: Felt> Table<E> for OutputTable<E> {
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        todo!()
    }

    fn pad(&mut self, n: usize) {
        todo!()
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        todo!()
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        todo!()
    }

    fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        todo!()
    }

    fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        todo!()
    }

    fn extension_terminal_constraints(challenges: &[E], terminals: &[E]) -> Vec<Multivariate<E>> {
        todo!()
    }

    fn max_degree(&self) -> usize {
        todo!()
    }

    fn set_matrix(&mut self, matrix: Vec<[E; Self::BASE_WIDTH]>) {
        todo!()
    }
}

pub struct InputTable<E> {
    num_padded_rows: usize,
    num_randomizers: usize,
    matrix: Vec<[E; BASE_WIDTH]>,
}

impl<E: Felt> Table<E> for InputTable<E> {
    const BASE_WIDTH: usize = BASE_WIDTH;
    const EXTENSION_WIDTH: usize = EXTENSION_WIDTH;

    fn len(&self) -> usize {
        todo!()
    }

    fn pad(&mut self, n: usize) {
        todo!()
    }

    fn base_boundary_constraints() -> Vec<Multivariate<E>> {
        todo!()
    }

    fn base_transition_constraints() -> Vec<Multivariate<E>> {
        todo!()
    }

    fn extension_boundary_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        todo!()
    }

    fn extension_transition_constraints(challenges: &[E]) -> Vec<Multivariate<E>> {
        todo!()
    }

    fn extension_terminal_constraints(challenges: &[E], terminals: &[E]) -> Vec<Multivariate<E>> {
        todo!()
    }

    fn max_degree(&self) -> usize {
        todo!()
    }

    fn set_matrix(&mut self, matrix: Vec<[E; Self::BASE_WIDTH]>) {
        todo!()
    }
}

fn base_boundary_constraints() {}

// pub struct IoTable<E> {
//     pub matrix: Vec<[E; BASE_WIDTH]>,
// }

// impl<E: PrimeFelt> IoTable<E> {
//     // Base columns
//     const VALUE: usize = 0;
//     // Extension columns
//     const EVALUATION: usize = 1;

//     pub fn new() -> Self {
//         IoTable { matrix: Vec::new() }
//     }

//     /// Sets the non padded matrix
//     pub fn set_matrix(&mut self, matrix: Vec<[E; 1]>) {
//         self.table.length = matrix.len();
//         self.table.matrix = matrix.into_iter().map(|row|
// row.to_vec()).collect();     }

//     pub fn pad(&mut self, n: usize) {
//         // TODO: seting length here seems kind of strange
//         self.table.length = self.table.matrix.len();
//         while self.table.matrix.len() < n {
//             self.table.matrix.push(vec![E::zero()]);
//         }
//         self.table.height = self.table.matrix.len();
//     }

//     pub fn base_boundary_constraints() -> Vec<Multivariate<E>> {
//         Vec::new()
//     }

//     pub fn base_transition_constraints() -> Vec<Multivariate<E>> {
//         let variables = Multivariate::variables(2);
//         vec![variables[Self::EVALUATION].clone() -
// variables[Self::VALUE].clone()]     }

//     pub fn extension_boundary_constraints() -> Vec<Multivariate<E>> {
//         vec![]
//     }

//     pub fn extension_transition_constraints(challenge: E) ->
// Vec<Multivariate<E>> {         let variables =
// Multivariate::<E>::variables(4);         let value =
// variables[Self::VALUE].clone();         let evaluation =
// variables[Self::EVALUATION].clone();         let value_next =
// variables[Self::VALUE].clone();         let evaluation_next =
// variables[Self::EVALUATION].clone();         vec![evaluation.clone() *
// challenge + value_next.clone() - evaluation_next.clone()]     }

//     pub fn extension_terminal_constraints(
//         &self,
//         challenge: E,
//         terminal: E,
//     ) -> Vec<Multivariate<E>> {
//         let variables = Multivariate::<E>::variables(2);
//         let offset = challenge.pow(&[(self.table.height - self.table.length)
// as u64]);         // In every padded row the running evaluation variable is
// multiplied by another         // factor `challenge`. We need to multiply
// `challenge ^ padding_length` to get         // the value of the evaluation
// terminal after all `2^k` rows.         let actual_terminal = terminal *
// offset;         vec![variables[Self::EVALUATION].clone() - actual_terminal]
//     }
// }

// impl<E: PrimeFelt> OutputTable<E> {
//     pub fn new(length: usize) -> Self {
//         OutputTable {
//             io_table: IoTable::new(length),
//         }
//     }

//     pub fn pad(&mut self, n: usize) {
//         self.io_table.pad(n);
//     }

//     pub fn base_boundary_constraints() -> Vec<Multivariate<E>> {
//         IoTable::<E>::base_boundary_constraints()
//     }

//     pub fn base_transition_constraints() -> Vec<Multivariate<E>> {
//         IoTable::<E>::base_transition_constraints()
//     }

//     pub fn extension_boundary_constraints(challenges: &[E]) ->
// Vec<Multivariate<E>> {         IoTable::<E>::extension_boundary_constraints()
//     }

//     pub fn extension_transition_constraints(challenges: &[E]) ->
// Vec<Multivariate<E>> {         let mut challenges_iter =
// challenges.iter().copied();         let _a = challenges_iter.next().unwrap();
//         let _b = challenges_iter.next().unwrap();
//         let _c = challenges_iter.next().unwrap();
//         let _d = challenges_iter.next().unwrap();
//         let _e = challenges_iter.next().unwrap();
//         let _f = challenges_iter.next().unwrap();
//         let _alpha = challenges_iter.next().unwrap();
//         let _beta = challenges_iter.next().unwrap();
//         let _gamma = challenges_iter.next().unwrap();
//         let delta = challenges_iter.next().unwrap();
//         let _eta = challenges_iter.next().unwrap();
//         IoTable::<E>::extension_transition_constraints(delta)
//     }

//     /// Sets the non padded matrix
//     pub fn set_matrix(&mut self, matrix: Vec<[E; 1]>) {
//         self.io_table.set_matrix(matrix)
//     }

//     pub fn extension_terminal_constraints(
//         &self,
//         challenges: &[E],
//         terminals: &[E],
//     ) -> Vec<Multivariate<E>> {
//         let mut challenges_iter = challenges.iter().copied();
//         let _a = challenges_iter.next().unwrap();
//         let _b = challenges_iter.next().unwrap();
//         let _c = challenges_iter.next().unwrap();
//         let _d = challenges_iter.next().unwrap();
//         let _e = challenges_iter.next().unwrap();
//         let _f = challenges_iter.next().unwrap();
//         let _alpha = challenges_iter.next().unwrap();
//         let _beta = challenges_iter.next().unwrap();
//         let _gamma = challenges_iter.next().unwrap();
//         let delta = challenges_iter.next().unwrap();
//         let _eta = challenges_iter.next().unwrap();

//         let mut terminal_iter = terminals.iter().copied();
//         let processor_instruction_permutation_terminal =
// terminal_iter.next().unwrap();         let
// processor_memory_permutation_terminal = terminal_iter.next().unwrap();
//         let processor_input_evaluation_terminal =
// terminal_iter.next().unwrap();         let
// processor_output_evaluation_terminal = terminal_iter.next().unwrap();
//         let instruction_evaluation_terminal = terminal_iter.next().unwrap();

//         self.io_table
//             .extension_terminal_constraints(delta,
// processor_output_evaluation_terminal)     }
// }

// impl<E: PrimeFelt> InputTable<E> {
//     pub fn new(length: usize) -> Self {
//         InputTable {
//             io_table: IoTable::new(length),
//         }
//     }

//     pub fn pad(&mut self, n: usize) {
//         self.io_table.pad(n);
//     }

//     pub fn base_boundary_constraints() -> Vec<Multivariate<E>> {
//         IoTable::<E>::base_boundary_constraints()
//     }

//     pub fn base_transition_constraints() -> Vec<Multivariate<E>> {
//         IoTable::<E>::base_transition_constraints()
//     }

//     pub fn extension_boundary_constraints(challenges: &[E]) ->
// Vec<Multivariate<E>> {         IoTable::<E>::extension_boundary_constraints()
//     }

//     pub fn extension_transition_constraints(challenges: &[E]) ->
// Vec<Multivariate<E>> {         let mut challenges_iter =
// challenges.iter().copied();         let _a = challenges_iter.next().unwrap();
//         let _b = challenges_iter.next().unwrap();
//         let _c = challenges_iter.next().unwrap();
//         let _d = challenges_iter.next().unwrap();
//         let _e = challenges_iter.next().unwrap();
//         let _f = challenges_iter.next().unwrap();
//         let _alpha = challenges_iter.next().unwrap();
//         let _beta = challenges_iter.next().unwrap();
//         let gamma = challenges_iter.next().unwrap();
//         let _delta = challenges_iter.next().unwrap();
//         let _eta = challenges_iter.next().unwrap();
//         IoTable::<E>::extension_transition_constraints(gamma)
//     }

//     /// Sets the non padded matrix
//     pub fn set_matrix(&mut self, matrix: Vec<[E; 1]>) {
//         self.io_table.set_matrix(matrix)
//     }

//     pub fn extension_terminal_constraints(
//         &self,
//         challenges: &[E],
//         terminals: &[E],
//     ) -> Vec<Multivariate<E>> {
//         let mut challenges_iter = challenges.iter().copied();
//         let _a = challenges_iter.next().unwrap();
//         let _b = challenges_iter.next().unwrap();
//         let _c = challenges_iter.next().unwrap();
//         let _d = challenges_iter.next().unwrap();
//         let _e = challenges_iter.next().unwrap();
//         let _f = challenges_iter.next().unwrap();
//         let _alpha = challenges_iter.next().unwrap();
//         let _beta = challenges_iter.next().unwrap();
//         let gamma = challenges_iter.next().unwrap();
//         let _delta = challenges_iter.next().unwrap();
//         let _eta = challenges_iter.next().unwrap();

//         let mut terminal_iter = terminals.iter().copied();
//         let processor_instruction_permutation_terminal =
// terminal_iter.next().unwrap();         let
// processor_memory_permutation_terminal = terminal_iter.next().unwrap();
//         let processor_input_evaluation_terminal =
// terminal_iter.next().unwrap();         let
// processor_output_evaluation_terminal = terminal_iter.next().unwrap();
//         let instruction_evaluation_terminal = terminal_iter.next().unwrap();

//         self.io_table
//             .extension_terminal_constraints(gamma,
// processor_input_evaluation_terminal)     }
// }
