use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use brainfuck::InputTable;
use brainfuck::InstructionTable;
use brainfuck::MemoryTable;
use brainfuck::OutputTable;
use brainfuck::ProcessorTable;
use protocol::ProofStream;
use std::cmp::max;
use std::marker::PhantomData;
use std::vec;

mod protocol;

pub struct StarkParams {
    /// power of 2 expansion factor
    expansion_factor: usize,
    /// security level of generated proofs
    security_level: usize,
    // TODO: fri params. folding factor, queries, etc.
}

impl StarkParams {
    pub fn new(expansion_factor: usize, security_level: usize) -> StarkParams {
        assert!(expansion_factor >= 4, "must be 4 or greater");
        assert!(expansion_factor.is_power_of_two(), "not a power of two");
        StarkParams {
            expansion_factor,
            security_level,
        }
    }

    pub fn num_randomizers(&self) -> usize {
        self.security_level
    }

    pub fn security_level(&self) -> usize {
        self.security_level
    }

    pub fn expansion_factor(&self) -> usize {
        self.expansion_factor
    }
}

pub struct BrainFuckStark<E> {
    params: StarkParams,
    processor_table: ProcessorTable<E>,
    memory_table: MemoryTable<E>,
    instruction_table: InstructionTable<E>,
    input_table: InputTable<E>,
    output_table: OutputTable<E>,
}

impl<E: PrimeFelt + StarkFelt> BrainFuckStark<E> {
    pub fn new(params: StarkParams) -> BrainFuckStark<E> {
        let num_randomizers = params.num_randomizers();
        BrainFuckStark {
            params,
            processor_table: ProcessorTable::new(num_randomizers),
            memory_table: MemoryTable::new(num_randomizers),
            instruction_table: InstructionTable::new(num_randomizers),
            input_table: InputTable::new(num_randomizers),
            output_table: OutputTable::new(num_randomizers),
        }
    }

    fn max_degree(&self) {
        let max_degree = self
            .processor_table
            .max_degree()
            .max(self.memory_table.max_degree())
            .max(self.instruction_table.max_degree())
            .max(self.input_table.max_degree())
            .max(self.output_table.max_degree());
        if max_length.is_power_of_two() {
            max_length
        } else {
            max_length.next_power_of_two()
        }
    }

    pub fn prove<T: ProofStream<E>>(
        &mut self,
        processor_matrix: Vec<[E; 7]>,
        memory_matrix: Vec<[E; 4]>,
        instruction_matrix: Vec<[E; 3]>,
        input_matrix: Vec<[E; 1]>,
        output_matrix: Vec<[E; 1]>,
        proof_stream: &mut T,
    ) -> Vec<u8> {
        let padding_length = {
            let max_length = processor_matrix
                .len()
                .max(memory_matrix.len())
                .max(instruction_matrix.len())
                .max(input_matrix.len())
                .max(output_matrix.len());
            if max_length.is_power_of_two() {
                max_length
            } else {
                max_length.next_power_of_two()
            }
        };

        let Self {
            processor_table,
            memory_table,
            instruction_table,
            input_table,
            output_table,
            ..
        } = self;

        processor_table.set_matrix(processor_matrix);
        memory_table.set_matrix(memory_matrix);
        instruction_table.set_matrix(instruction_matrix);
        input_table.set_matrix(input_matrix);
        output_table.set_matrix(output_matrix);

        // pad tables to height 2^k
        processor_table.pad(padding_length);
        memory_table.pad(padding_length);
        instruction_table.pad(padding_length);
        input_table.pad(padding_length);
        output_table.pad(padding_length);

        // let tables = vec![
        //     &processor_table as &dyn Table<E>,
        //     &memory_table as &dyn Table<E>,
        //     &instruction_table as &dyn Table<E>,
        //     &input_table as &dyn Table<E>,
        //     &output_table as &dyn Table<E>,
        // ];

        // let max_degree = tables.iter().map(|table|
        // table.max_degree()).max().unwrap(); let fri_domain_length =

        Vec::new()
    }
}
