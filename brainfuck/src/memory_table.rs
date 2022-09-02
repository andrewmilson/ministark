use super::table::Table;
use crate::processor_table::ProcessorTable;
use algebra::PrimeFelt;

pub struct MemoryTable<E> {
    table: Table<E>,
}

impl<E: PrimeFelt> MemoryTable<E> {
    const CYCLE: usize = 0;
    const MP: usize = 1;
    const MEM_VAL: usize = 2;

    /// Outputs an unpadded but interweaved matrix
    pub fn derive_matrix(processor_matrix: &[[E; 7]]) -> Vec<[E; 4]> {
        // copy unpadded rows and sort
        let mut matrix = processor_matrix
            .iter()
            .map(|row| {
                [
                    row[ProcessorTable::<E>::CYCLE],
                    row[ProcessorTable::<E>::MP],
                    row[ProcessorTable::<E>::MEM_VAL],
                    E::zero(),
                ]
            })
            .collect::<Vec<[E; 4]>>();
        matrix.sort_by_key(|row| row[Self::MP].into_bigint());

        // insert dummy rows for smooth clk jumps
        // for i in 0..

        todo!()
    }
}
