#![feature(generic_const_exprs, int_log)]

use algebra::ExtensionOf;
use algebra::Felt;
use algebra::Multivariate;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use algebra::UniformRand;
use algebra::Univariate;
use brainfuck::permutation_argument;
use brainfuck::InputTable;
use brainfuck::InstructionTable;
use brainfuck::MemoryTable;
use brainfuck::OutputTable;
use brainfuck::ProcessorTable;
use brainfuck::Table;
use fri::Fri;
use mini_stark::number_theory_transform::number_theory_transform;
use mini_stark::polynomial::Polynomial;
use num_traits::Zero;
use protocol::ProofStream;
use rand::Rng;
use salted_merkle::SaltedMerkle;
use std::cmp::max;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;
use std::iter::empty;
use std::marker::PhantomData;
use std::vec;

mod fri;
mod merkle;
pub mod protocol;
mod salted_merkle;

pub trait Config {
    type BaseFelt: PrimeFelt + StarkFelt;
    type ExtensionFelt: Felt<BaseFelt = Self::BaseFelt> + ExtensionOf<Self::BaseFelt>;
    // type MerkleParams: ark_crypto_primitives::merkle_tree::Config;

    const EXPANSION_FACTOR: usize;
    const SECURITY_LEVEL: usize;
    const NUM_RANDOMIZERS: usize = Self::SECURITY_LEVEL;
}

// Generate fri config from stark config
impl<P: Config> fri::Config for P {
    type BaseFelt = P::BaseFelt;
    type ExtensionFelt = P::ExtensionFelt;
    // type MerkleParams = P::MerkleParams;

    const EXPANSION_FACTOR: usize = P::EXPANSION_FACTOR;
    const SECURITY_LEVEL: usize = P::SECURITY_LEVEL;

    // assert!(EXPANSION_FACTOR >= 4, "must be 4 or greater");
    // assert!(EXPANSION_FACTOR.is_power_of_two(), "not a power of two");
}

// pub struct StarkParams {
//     /// power of 2 expansion factor
//     EXPANSION_FACTOR: usize,
//     /// security level of generated proofs
//     SECURITY_LEVEL: usize,
//     // TODO: fri params. folding factor, queries, etc.
// }

// impl StarkParams {
//     pub fn new(EXPANSION_FACTOR: usize, SECURITY_LEVEL: usize) -> StarkParams
// {
//         StarkParams {
//             EXPANSION_FACTOR,
//             SECURITY_LEVEL,
//         }
//     }

//     pub fn NUM_RANDOMIZERS(&self) -> usize {
//         0 //self.SECURITY_LEVEL
//     }

//     pub fn SECURITY_LEVEL(&self) -> usize {
//         self.SECURITY_LEVEL
//     }

//     pub fn EXPANSION_FACTOR(&self) -> usize {
//         self.EXPANSION_FACTOR
//     }

//     pub fn fri_params(&self) -> FriParams {
//         let expension_factor = self.EXPANSION_FACTOR();
//         let num_colinearity_checks = self.SECURITY_LEVEL() /
// expension_factor.ilog2() as usize;         FriParams::new(expension_factor,
// num_colinearity_checks)     }
// }

pub struct BrainFuckStark<P: Config> {
    fri: Fri<P>,
    processor_table: ProcessorTable<P::BaseFelt, P::ExtensionFelt>,
    memory_table: MemoryTable<P::BaseFelt, P::ExtensionFelt>,
    instruction_table: InstructionTable<P::BaseFelt, P::ExtensionFelt>,
    input_table: InputTable<P::BaseFelt, P::ExtensionFelt>,
    output_table: OutputTable<P::BaseFelt, P::ExtensionFelt>,
}

impl<P: Config> BrainFuckStark<P> {
    pub fn new(params: P) -> Self {
        BrainFuckStark {
            fri: Fri::new(params),
            processor_table: ProcessorTable::new(P::NUM_RANDOMIZERS),
            memory_table: MemoryTable::new(P::NUM_RANDOMIZERS),
            instruction_table: InstructionTable::new(P::NUM_RANDOMIZERS),
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
        (self.max_degree() + 1) * P::EXPANSION_FACTOR
    }

    fn sample_weights(&self, n: usize, randomness: u64) -> Vec<P::ExtensionFelt> {
        (0..n)
            .map(|i| {
                let mut hash = DefaultHasher::new();
                (randomness as u128 + i as u128).hash(&mut hash);
                P::ExtensionFelt::from(hash.finish())
            })
            .collect()
    }

    fn sample_indices(n: usize, randomness: u64, bound: usize) -> Vec<usize> {
        let mut indices = vec![];
        let mut hasher = DefaultHasher::new();
        for _ in 0..n {
            randomness.hash(&mut hasher);
            indices.push((hasher.finish() as usize) % bound)
        }
        indices
    }

    pub fn prove<T: ProofStream<P::ExtensionFelt>>(
        &mut self,
        processor_matrix: Vec<
            [P::BaseFelt; ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH],
        >,
        memory_matrix: Vec<[P::BaseFelt; MemoryTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH]>,
        instruction_matrix: Vec<
            [P::BaseFelt; InstructionTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH],
        >,
        input_matrix: Vec<[P::BaseFelt; InputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH]>,
        output_matrix: Vec<[P::BaseFelt; OutputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH]>,
        proof_stream: &mut T,
    ) -> Vec<u8> {
        let mut rng = rand::thread_rng();

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
        let offset = P::BaseFelt::GENERATOR;
        let omega = P::BaseFelt::get_root_of_unity(codeword_len.ilog2());

        let randomizer_codewords = {
            let n = ceil_power_of_two(self.max_degree());
            let coefficients = (0..n).map(|_| P::ExtensionFelt::rand(&mut rng)).collect();
            let polynomial = Univariate::new(coefficients);
            polynomial.scale(offset.into());
            let mut coefficients = polynomial.coefficients;
            coefficients.resize(codeword_len, P::ExtensionFelt::zero());
            vec![number_theory_transform(&coefficients)]
        };

        let base_processor_lde = self.processor_table.base_lde(offset, codeword_len);
        let base_memory_lde = self.memory_table.base_lde(offset, codeword_len);
        let base_instruction_lde = self.instruction_table.base_lde(offset, codeword_len);
        let base_input_lde = self.input_table.base_lde(offset, codeword_len);
        let base_output_lde = self.output_table.base_lde(offset, codeword_len);

        let base_codewords = vec![
            base_processor_lde.clone(),
            base_memory_lde.clone(),
            base_instruction_lde.clone(),
            base_input_lde.clone(),
            base_output_lde.clone(),
        ]
        .concat();

        let base_degree_bounds = vec![
            vec![
                self.processor_table.interpolant_degree();
                ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.memory_table.interpolant_degree();
                MemoryTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.instruction_table.interpolant_degree();
                InstructionTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.input_table.interpolant_degree();
                InputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.output_table.interpolant_degree();
                OutputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
        ]
        .concat();

        let base_zipped_codeword = (0..codeword_len)
            .map(|i| base_codewords.iter().map(|codeword| codeword[i]).collect())
            .collect::<Vec<Vec<P::ExtensionFelt>>>();
        let base_tree = SaltedMerkle::new(&base_zipped_codeword);
        proof_stream.push(protocol::ProofObject::MerkleRoot(base_tree.root()));
        // TODO: base degree bounds

        // get coefficients for table extensions
        let challenges = self.sample_weights(11, proof_stream.prover_fiat_shamir());
        let initials = vec![
            P::ExtensionFelt::rand(&mut rng),
            P::ExtensionFelt::rand(&mut rng),
        ];

        self.processor_table.extend(&challenges, &initials);
        self.memory_table.extend(&challenges, &initials);
        self.instruction_table.extend(&challenges, &initials);
        self.input_table.extend(&challenges, &initials);
        self.output_table.extend(&challenges, &initials);

        let terminals = vec![
            self.processor_table.instr_permutation_terminal.unwrap(),
            self.processor_table.memory_permutation_terminal.unwrap(),
            self.processor_table.input_evaluation_terminal.unwrap(),
            self.processor_table.output_evaluation_terminal.unwrap(),
            self.instruction_table.evaluation_terminal.unwrap(),
        ];

        let ext_processor_lde = self.processor_table.extension_lde(offset, codeword_len);
        let ext_memory_lde = self.memory_table.extension_lde(offset, codeword_len);
        let ext_instruction_lde = self.instruction_table.extension_lde(offset, codeword_len);
        let ext_input_lde = self.input_table.extension_lde(offset, codeword_len);
        let ext_output_lde = self.output_table.extension_lde(offset, codeword_len);

        let extension_codewords = vec![
            ext_processor_lde.clone(),
            ext_memory_lde.clone(),
            ext_instruction_lde.clone(),
            ext_input_lde.clone(),
            ext_output_lde.clone(),
        ]
        .concat();

        let extension_degree_bounds = vec![
            vec![
                self.processor_table.interpolant_degree();
                ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                    - ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.memory_table.interpolant_degree();
                MemoryTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                    - MemoryTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.instruction_table.interpolant_degree();
                InstructionTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                    - InstructionTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.input_table.interpolant_degree();
                InputTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                    - InputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
            vec![
                self.output_table.interpolant_degree();
                OutputTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                    - OutputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            ],
        ]
        .concat();

        let ext_zipped_codeword = (0..codeword_len)
            .map(|i| {
                extension_codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect()
            })
            .collect::<Vec<Vec<P::ExtensionFelt>>>();
        let extension_tree = SaltedMerkle::new(&ext_zipped_codeword);
        proof_stream.push(protocol::ProofObject::MerkleRoot(extension_tree.root()));
        // TODO: extension degree bounds

        let processor_codewords = vec![base_processor_lde, ext_processor_lde].concat();
        let memory_codewords = vec![base_memory_lde, ext_memory_lde].concat();
        let instruction_codewords = vec![base_instruction_lde, ext_instruction_lde].concat();
        let input_codewords = vec![base_input_lde, ext_input_lde].concat();
        let output_codewords = vec![base_output_lde, ext_output_lde].concat();

        let mut quotient_codewords = vec![
            self.processor_table.all_quotients(
                codeword_len,
                &processor_codewords,
                &challenges,
                &terminals,
            ),
            self.memory_table.all_quotients(
                codeword_len,
                &memory_codewords,
                &challenges,
                &terminals,
            ),
            self.instruction_table.all_quotients(
                codeword_len,
                &instruction_codewords,
                &challenges,
                &terminals,
            ),
            self.input_table
                .all_quotients(codeword_len, &input_codewords, &challenges, &terminals),
            self.output_table.all_quotients(
                codeword_len,
                &output_codewords,
                &challenges,
                &terminals,
            ),
        ]
        .concat();

        let mut quotient_degree_bounds = vec![
            self.processor_table
                .all_quotient_degree_bounds(&challenges, &terminals),
            self.memory_table
                .all_quotient_degree_bounds(&challenges, &terminals),
            self.instruction_table
                .all_quotient_degree_bounds(&challenges, &terminals),
            self.input_table
                .all_quotient_degree_bounds(&challenges, &terminals),
            self.output_table
                .all_quotient_degree_bounds(&challenges, &terminals),
        ]
        .concat();

        // Instruction permutation
        quotient_codewords.push(permutation_argument::quotient(
            &processor_codewords
                [ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::INSTRUCTION_PERMUTATION],
            &instruction_codewords
                [InstructionTable::<P::BaseFelt, P::ExtensionFelt>::PROCESSOR_PERMUTATION],
        ));
        quotient_degree_bounds.push(
            usize::max(
                self.processor_table.interpolant_degree(),
                self.instruction_table.interpolant_degree(),
            ) - 1,
        );

        // Memory permutation
        quotient_codewords.push(permutation_argument::quotient(
            &processor_codewords
                [ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::MEMORY_PERMUTATION],
            &memory_codewords[MemoryTable::<P::BaseFelt, P::ExtensionFelt>::PERMUTATION],
        ));
        quotient_degree_bounds.push(
            usize::max(
                self.processor_table.interpolant_degree(),
                self.memory_table.interpolant_degree(),
            ) - 1,
        );

        for &terminal in &terminals {
            proof_stream.push(protocol::ProofObject::Terminal(terminal));
        }

        // get weights for non-linear combination
        // - 1 for randomizer polynomial
        // - 2 for every other polynomial (base, extension, quotients)
        let num_base_polynomials = ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            + MemoryTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            + InstructionTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            + InputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH
            + OutputTable::<P::BaseFelt, P::ExtensionFelt>::BASE_WIDTH;
        let num_extension_polynomials =
            ProcessorTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                + MemoryTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                + InstructionTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                + InputTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                + OutputTable::<P::BaseFelt, P::ExtensionFelt>::EXTENSION_WIDTH
                - num_base_polynomials;
        let num_randomizer_polynomials = randomizer_codewords.len();
        let num_quotient_polynomials = quotient_degree_bounds.len();
        let weights_seed = proof_stream.prover_fiat_shamir();
        let weights = self.sample_weights(
            num_randomizer_polynomials
                + 2 * (num_base_polynomials + num_extension_polynomials + num_quotient_polynomials),
            weights_seed,
        );

        // compute terms of non-linear combination polynomial
        let mut terms = randomizer_codewords.clone();
        assert_eq!(base_codewords.len(), num_base_polynomials);
        assert_eq!(extension_codewords.len(), num_extension_polynomials);
        assert_eq!(quotient_codewords.len(), num_quotient_polynomials);

        for (codeword, degree_bound) in empty()
            .chain(base_codewords.into_iter().zip(base_degree_bounds))
            .chain(extension_codewords.into_iter().zip(extension_degree_bounds))
            .chain(quotient_codewords.into_iter().zip(quotient_degree_bounds))
        {
            let shift = (self.max_degree() - degree_bound) as u64;
            // dot product of codeword and evaluation of the polynomial `x^shift`
            let shifted_codeword = (0..codeword_len)
                .map(|i| {
                    P::ExtensionFelt::from(offset * omega.pow(&[i as u64])).pow(&[shift])
                        * codeword[i]
                })
                .collect::<Vec<P::ExtensionFelt>>();
            terms.push(codeword);
            terms.push(shifted_codeword);
        }

        assert_eq!(terms.len(), weights.len());
        let combination_codeword = terms
            .into_iter()
            .zip(weights)
            .map(|(term, weight)| term.into_iter().map(|v| v * weight).collect())
            .fold(vec![P::ExtensionFelt::zero(); codeword_len], sum);
        let combination_tree = SaltedMerkle::new(&combination_codeword);
        proof_stream.push(protocol::ProofObject::MerkleRoot(combination_tree.root()));

        // get indices of leafs to prove non-linear combination
        let indices_seed = proof_stream.prover_fiat_shamir();
        let indices = Self::sample_indices(P::SECURITY_LEVEL, indices_seed, codeword_len);

        let row_step = P::EXPANSION_FACTOR;
        assert!(row_step.is_power_of_two());

        // open leafs of zipped codewords at indicated positions
        for &index in &indices {
            let idx = (index + row_step) % codeword_len;
            let (element, salt, path) = base_tree.open(idx);
            assert!(SaltedMerkle::verify(
                base_tree.root(),
                idx,
                salt,
                &path,
                &element
            ));
            proof_stream.push(protocol::ProofObject::LeafItems(element));
            proof_stream.push(protocol::ProofObject::MerklePathWithSalt((salt, path)));
            let (element, salt, path) = extension_tree.open(idx);
            assert!(SaltedMerkle::verify(
                extension_tree.root(),
                idx,
                salt,
                &path,
                &element
            ));
            proof_stream.push(protocol::ProofObject::LeafItems(element));
            proof_stream.push(protocol::ProofObject::MerklePathWithSalt((salt, path)));
        }

        // TODO: merge these arrays

        // open combination codewords at same positions
        for index in indices {
            let (element, salt, path) = combination_tree.open(index);
            assert!(SaltedMerkle::verify(
                combination_tree.root(),
                index,
                salt,
                &path,
                &element
            ));
            proof_stream.push(protocol::ProofObject::LeafItem(element));
            proof_stream.push(protocol::ProofObject::MerklePathWithSalt((salt, path)));
        }

        // prove low degree of combination polynomial and collect indices
        self.fri.prove(proof_stream, &combination_codeword);

        // the final proof is just the serialized stream
        proof_stream.serialize()
    }
}

/// Rounds the input value up the the nearest power of two
pub(crate) fn ceil_power_of_two(value: usize) -> usize {
    if value.is_power_of_two() {
        value
    } else {
        value.next_power_of_two()
    }
}

fn sum<E: Felt>(a: Vec<E>, b: Vec<E>) -> Vec<E> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    (0..n).map(|i| a[i] + b[i]).collect()
}
