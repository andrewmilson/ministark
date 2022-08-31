use super::table::Table;
use algebra::Felt;

struct MemoryTable<E> {
    cycle: usize,
    mp: usize,
    memory_value: usize,
    table: Table<E>,
}

impl<E: Felt> MemoryTable<E> {
    fn derive_matrix(processor_matrix: Vec<[E; 7]>) {}
}
