#![feature(generic_const_exprs, int_log)]

use algebra::ExtensionOf;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use algebra::Univariate;
use brainfuck::InputTable;
use brainfuck::InstructionTable;
use brainfuck::MemoryTable;
use brainfuck::OutputTable;
use brainfuck::ProcessorTable;
use brainfuck::Table;
use mini_stark::number_theory_transform::number_theory_transform;
use mini_stark::polynomial::Polynomial;
use protocol::ProofStream;
use rand::Rng;
use salted_merkle::SaltedMerkle;
use std::cmp::max;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::vec;

pub mod protocol;
mod salted_merkle;

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

pub struct BrainFuckStark<F, E> {
    params: StarkParams,
    processor_table: ProcessorTable<F, E>,
    memory_table: MemoryTable<F, E>,
    instruction_table: InstructionTable<F, E>,
    input_table: InputTable<F, E>,
    output_table: OutputTable<F, E>,
}

impl<F: PrimeFelt + StarkFelt, E: Felt + ExtensionOf<F>> BrainFuckStark<F, E>
where
    F: PrimeFelt + StarkFelt,
    E: Felt<BaseFelt = F> + ExtensionOf<F>,
{
    pub fn new(params: StarkParams) -> Self {
        let num_randomizers = params.num_randomizers();
        BrainFuckStark {
            params,
            processor_table: ProcessorTable::new(num_randomizers),
            memory_table: MemoryTable::new(num_randomizers),
            instruction_table: InstructionTable::new(num_randomizers),
            input_table: InputTable::new(),
            output_table: OutputTable::new(),
        }
    }

    fn max_degree(&self) -> usize {
        assert!(!self.processor_table.is_empty(), "tables not populated");
        // TODO: could be a bug here... Instead of rounding up to the power of two it
        // should be the next power of two.
        let max_degree = self
            .processor_table
            .max_degree()
            .max(self.memory_table.max_degree())
            .max(self.instruction_table.max_degree())
            .max(self.input_table.max_degree())
            .max(self.output_table.max_degree());
        ceil_power_of_two(max_degree) - 1
    }

    fn fri_codeword_length(&self) -> usize {
        (self.max_degree() + 1) * self.params.expansion_factor
    }

    fn sample_weights(&self, n: usize, randomness: u64) -> Vec<E> {
        (0..n)
            .map(|i| {
                let mut hash = DefaultHasher::new();
                (randomness + i as u64).hash(&mut hash);
                E::from(hash.finish())
            })
            .collect()
    }

    pub fn prove<T: ProofStream<F>>(
        &mut self,
        processor_matrix: Vec<[F; ProcessorTable::<F, E>::BASE_WIDTH]>,
        memory_matrix: Vec<[F; MemoryTable::<F, E>::BASE_WIDTH]>,
        instruction_matrix: Vec<[F; InstructionTable::<F, E>::BASE_WIDTH]>,
        input_matrix: Vec<[F; InputTable::<F, E>::BASE_WIDTH]>,
        output_matrix: Vec<[F; OutputTable::<F, E>::BASE_WIDTH]>,
        proof_stream: &mut T,
    ) -> Vec<u8> {
        let padding_length = {
            let max_length = processor_matrix
                .len()
                .max(memory_matrix.len())
                .max(instruction_matrix.len())
                .max(input_matrix.len())
                .max(output_matrix.len());
            ceil_power_of_two(max_length)
        };

        self.processor_table.set_matrix(processor_matrix);
        self.memory_table.set_matrix(memory_matrix);
        self.instruction_table.set_matrix(instruction_matrix);
        self.input_table.set_matrix(input_matrix);
        self.output_table.set_matrix(output_matrix);

        // pad tables to height 2^k
        self.processor_table.pad(padding_length);
        self.memory_table.pad(padding_length);
        self.instruction_table.pad(padding_length);
        self.input_table.pad(padding_length);
        self.output_table.pad(padding_length);

        let codeword_len = self.fri_codeword_length();

        let randomizer_codewords = {
            let mut rng = rand::thread_rng();
            let n = ceil_power_of_two(self.max_degree());
            let coefficients = (0..n).map(|_| E::rand(&mut rng)).collect();
            let polynomial = Univariate::new(coefficients);
            polynomial.scale(F::GENERATOR.into());
            let mut coefficients = polynomial.coefficients;
            coefficients.resize(codeword_len, E::zero());
            vec![number_theory_transform(&coefficients)]
        };

        let base_codewords = vec![
            self.processor_table.base_lde(F::GENERATOR, codeword_len),
            self.memory_table.base_lde(F::GENERATOR, codeword_len),
            self.instruction_table.base_lde(F::GENERATOR, codeword_len),
            self.input_table.base_lde(F::GENERATOR, codeword_len),
            self.output_table.base_lde(F::GENERATOR, codeword_len),
            randomizer_codewords,
        ]
        .concat();
        let zipped_codeword = (0..codeword_len)
            .map(|i| base_codewords.iter().map(|codeword| codeword[i]).collect())
            .collect::<Vec<Vec<E>>>();
        let base_tree = SaltedMerkle::new(&zipped_codeword);
        proof_stream.push(protocol::ProofObject::MerkleRoot(base_tree.root()));

        // get coefficients for table extensions
        let challenges = self.sample_weights(11, proof_stream.prover_fiat_shamir());

        // println!("{}", base_tree.root());

        Vec::new()
    }
}

/// Rounds the input value up the the nearest power of two
fn ceil_power_of_two(value: usize) -> usize {
    if value.is_power_of_two() {
        value
    } else {
        value.next_power_of_two()
    }
}
