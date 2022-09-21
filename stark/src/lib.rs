#![feature(generic_const_exprs, int_log)]

use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::PrimeField;
use ark_ff::UniformRand;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use brainfuck::permutation_argument;
use brainfuck::InputTable;
use brainfuck::InstructionTable;
use brainfuck::MemoryTable;
use brainfuck::OutputTable;
use brainfuck::ProcessorTable;
use brainfuck::Table;
use fri::Fri;
use legacy_algebra::number_theory_transform::number_theory_transform;
use legacy_algebra::scale_poly;
use num_traits::Zero;
use protocol::ProofObject;
use protocol::ProofStream;
use salted_merkle::SaltedMerkle;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::Hash;
use std::hash::Hasher;
use std::ops::Add;
use std::vec;

mod fri;
mod merkle;
pub mod protocol;
mod salted_merkle;

/// The domain separator, used when proving statements on gemini.
pub(crate) const PROTOCOL_NAME: &[u8] = b"ZK-WASM-v0";

pub trait Config {
    /// Base prime field
    type Fp: PrimeField + FftField;
    /// Extension field element
    type Fx: Field<BasePrimeField = Self::Fp>;

    // type MTParams: ark_crypto_primitives::merkle_tree::Config;
    // type SaltedMTParams: ark_crypto_primitives::merkle_tree::Config;

    const EXPANSION_FACTOR: usize;
    const SECURITY_LEVEL: usize;
    const NUM_RANDOMIZERS: usize = Self::SECURITY_LEVEL;
}

// Generate fri config from stark config
impl<P: Config> fri::Config for P {
    type Fp = P::Fp;
    type Fx = P::Fx;

    const EXPANSION_FACTOR: usize = P::EXPANSION_FACTOR;
    const SECURITY_LEVEL: usize = P::SECURITY_LEVEL;

    // assert!(EXPANSION_FACTOR >= 4, "must be 4 or greater");
    // assert!(EXPANSION_FACTOR.is_power_of_two(), "not a power of two");
}

pub struct BrainFuckStark<P: Config> {
    fri: Fri<P>,
    processor_table: ProcessorTable<P::Fx>,
    memory_table: MemoryTable<P::Fx>,
    instruction_table: InstructionTable<P::Fx>,
    input_table: InputTable<P::Fx>,
    output_table: OutputTable<P::Fx>,
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
        // assert!(!self.processor_table.is_empty(), "tables not populated");
        // TODO: could be a bug here... Instead of rounding up to the power of two it
        // should be the next power of two.
        let max_degree = [
            self.processor_table.max_degree(),
            self.memory_table.max_degree(),
            self.instruction_table.max_degree(),
            self.input_table.max_degree(),
            self.output_table.max_degree(),
        ]
        .into_iter()
        .max()
        .unwrap();

        println!(
            "{} {} {} {} {}",
            self.processor_table.max_degree(),
            self.memory_table.max_degree(),
            self.instruction_table.max_degree(),
            self.input_table.max_degree(),
            self.output_table.max_degree(),
        );
        ceil_power_of_two(max_degree) - 1
    }

    fn fri_codeword_length(&self) -> usize {
        (self.max_degree() + 1) * P::EXPANSION_FACTOR
    }

    fn sample_weights(&self, n: usize, randomness: u64) -> Vec<P::Fx> {
        (0..n)
            .map(|i| {
                let mut hash = DefaultHasher::new();
                (randomness as u128 + i as u128).hash(&mut hash);
                P::Fx::from(hash.finish())
            })
            .collect()
    }

    fn sample_indices(n: usize, randomness: u64, bound: usize) -> Vec<usize> {
        println!("N:{n}, randomness:{randomness}, bound:{bound}");
        let mut indices = vec![];
        let mut hasher = DefaultHasher::new();
        for _ in 0..n {
            randomness.hash(&mut hasher);
            indices.push((hasher.finish() as usize) % bound)
        }
        indices
    }

    pub fn prove(
        &mut self,
        processor_matrix: Vec<[P::Fp; ProcessorTable::<P::Fx>::BASE_WIDTH]>,
        memory_matrix: Vec<[P::Fp; MemoryTable::<P::Fx>::BASE_WIDTH]>,
        instruction_matrix: Vec<[P::Fp; InstructionTable::<P::Fx>::BASE_WIDTH]>,
        input_matrix: Vec<[P::Fp; InputTable::<P::Fx>::BASE_WIDTH]>,
        output_matrix: Vec<[P::Fp; OutputTable::<P::Fx>::BASE_WIDTH]>,
        proof_stream: &mut impl ProofStream<P::Fx>,
    ) -> Vec<u8> {
        let mut rng = rand::thread_rng();

        let padding_length = {
            let max_length = [
                processor_matrix.len(),
                memory_matrix.len(),
                instruction_matrix.len(),
                input_matrix.len(),
                output_matrix.len(),
            ]
            .into_iter()
            .max()
            .unwrap();
            ceil_power_of_two(max_length)
        };

        println!("padding len: {padding_length}");

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

        proof_stream.push(ProofObject::TraceLen(padding_length as u64));

        let codeword_len = self.fri_codeword_length();
        println!("Codeword len {}", codeword_len);
        let offset = P::Fp::GENERATOR;
        let omega = P::Fp::get_root_of_unity(codeword_len as u64).unwrap();

        let randomizer_codewords = {
            let n = ceil_power_of_two(self.max_degree());
            let coefficients = (0..n).map(|_| P::Fx::rand(&mut rng)).collect::<Vec<_>>();
            let polynomial = DensePolynomial::from_coefficients_vec(coefficients);
            let polynomial = scale_poly(&polynomial, P::Fx::from_base_prime_field(offset));
            let mut coefficients = polynomial.coeffs;
            coefficients.resize(codeword_len, P::Fx::zero());
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
            vec![self.processor_table.interpolant_degree(); ProcessorTable::<P::Fx>::BASE_WIDTH],
            vec![self.memory_table.interpolant_degree(); MemoryTable::<P::Fx>::BASE_WIDTH],
            vec![
                self.instruction_table.interpolant_degree();
                InstructionTable::<P::Fx>::BASE_WIDTH
            ],
            vec![self.input_table.interpolant_degree(); InputTable::<P::Fx>::BASE_WIDTH],
            vec![self.output_table.interpolant_degree(); OutputTable::<P::Fx>::BASE_WIDTH],
        ]
        .concat();

        let base_zipped_codeword = (0..codeword_len)
            .map(|i| base_codewords.iter().map(|codeword| codeword[i]).collect())
            .collect::<Vec<Vec<P::Fx>>>();
        let base_tree = SaltedMerkle::new(&base_zipped_codeword);
        proof_stream.push(ProofObject::MerkleRoot(base_tree.root()));
        // TODO: base degree bounds

        // get coefficients for table extensions
        let challenges = self.sample_weights(11, proof_stream.prover_fiat_shamir());
        let initials = vec![P::Fx::rand(&mut rng), P::Fx::rand(&mut rng)];

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
                ProcessorTable::<P::Fx>::EXTENSION_WIDTH - ProcessorTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.memory_table.interpolant_degree();
                MemoryTable::<P::Fx>::EXTENSION_WIDTH - MemoryTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.instruction_table.interpolant_degree();
                InstructionTable::<P::Fx>::EXTENSION_WIDTH - InstructionTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.input_table.interpolant_degree();
                InputTable::<P::Fx>::EXTENSION_WIDTH - InputTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.output_table.interpolant_degree();
                OutputTable::<P::Fx>::EXTENSION_WIDTH - OutputTable::<P::Fx>::BASE_WIDTH
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
            .collect::<Vec<Vec<P::Fx>>>();
        let extension_tree = SaltedMerkle::new(&ext_zipped_codeword);
        proof_stream.push(ProofObject::MerkleRoot(extension_tree.root()));
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
            &processor_codewords[ProcessorTable::<P::Fx>::INSTRUCTION_PERMUTATION],
            &instruction_codewords[InstructionTable::<P::Fx>::PROCESSOR_PERMUTATION],
        ));
        quotient_degree_bounds.push(
            std::cmp::max(
                self.processor_table.interpolant_degree(),
                self.instruction_table.interpolant_degree(),
            ) - 1,
        );

        // Memory permutation
        quotient_codewords.push(permutation_argument::quotient(
            &processor_codewords[ProcessorTable::<P::Fx>::MEMORY_PERMUTATION],
            &memory_codewords[MemoryTable::<P::Fx>::PERMUTATION],
        ));
        quotient_degree_bounds.push(
            std::cmp::max(
                self.processor_table.interpolant_degree(),
                self.memory_table.interpolant_degree(),
            ) - 1,
        );

        for &terminal in &terminals {
            proof_stream.push(ProofObject::Terminal(terminal));
        }

        // get weights for non-linear combination
        // - 1 for randomizer polynomial
        // - 2 for every other polynomial (base, extension, quotients)
        let num_base_polynomials = ProcessorTable::<P::Fx>::BASE_WIDTH
            + MemoryTable::<P::Fx>::BASE_WIDTH
            + InstructionTable::<P::Fx>::BASE_WIDTH
            + InputTable::<P::Fx>::BASE_WIDTH
            + OutputTable::<P::Fx>::BASE_WIDTH;
        let num_extension_polynomials = ProcessorTable::<P::Fx>::EXTENSION_WIDTH
            + MemoryTable::<P::Fx>::EXTENSION_WIDTH
            + InstructionTable::<P::Fx>::EXTENSION_WIDTH
            + InputTable::<P::Fx>::EXTENSION_WIDTH
            + OutputTable::<P::Fx>::EXTENSION_WIDTH
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
        let mut terms = randomizer_codewords;
        assert_eq!(base_codewords.len(), num_base_polynomials);
        assert_eq!(extension_codewords.len(), num_extension_polynomials);
        assert_eq!(quotient_codewords.len(), num_quotient_polynomials);

        for (codeword, degree_bound) in std::iter::empty()
            .chain(base_codewords.into_iter().zip(base_degree_bounds))
            .chain(extension_codewords.into_iter().zip(extension_degree_bounds))
            .chain(quotient_codewords.into_iter().zip(quotient_degree_bounds))
        {
            let shift = (self.max_degree() - degree_bound) as u64;
            // dot product of codeword and evaluation of the polynomial `x^shift`
            let shifted_codeword = (0..codeword_len)
                .map(|i| {
                    P::Fx::from_base_prime_field(offset * omega.pow([i as u64])).pow([shift])
                        * codeword[i]
                })
                .collect::<Vec<P::Fx>>();
            terms.push(codeword);
            terms.push(shifted_codeword);
        }

        assert_eq!(terms.len(), weights.len());
        let combination_codeword = terms
            .into_iter()
            .zip(weights)
            .map(|(term, weight)| term.into_iter().map(|v| v * weight).collect())
            .fold(vec![P::Fx::zero(); codeword_len], sum);
        let combination_tree = SaltedMerkle::new(&combination_codeword);
        proof_stream.push(ProofObject::MerkleRoot(combination_tree.root()));

        // get indices of leafs to prove non-linear combination
        let indices_seed = proof_stream.prover_fiat_shamir();
        let indices = Self::sample_indices(P::SECURITY_LEVEL, indices_seed, codeword_len);

        let row_step = P::EXPANSION_FACTOR;
        assert!(row_step.is_power_of_two());

        // open leafs of zipped codewords at indicated positions
        for &index in &indices {
            for distance in [0, row_step] {
                let idx = (index + distance) % codeword_len;
                let (element, salt, path) = base_tree.open(idx);
                assert!(SaltedMerkle::verify(
                    base_tree.root(),
                    idx,
                    salt,
                    &path,
                    &element
                ));
                proof_stream.push(ProofObject::LeafItems(element));
                proof_stream.push(ProofObject::MerklePathWithSalt((salt, path)));
                let (element, salt, path) = extension_tree.open(idx);
                assert!(SaltedMerkle::verify(
                    extension_tree.root(),
                    idx,
                    salt,
                    &path,
                    &element
                ));
                proof_stream.push(ProofObject::LeafItems(element));
                proof_stream.push(ProofObject::MerklePathWithSalt((salt, path)));
            }
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
            proof_stream.push(ProofObject::LeafItem(element));
            proof_stream.push(ProofObject::MerklePathWithSalt((salt, path)));
        }

        // prove low degree of combination polynomial and collect indices
        self.fri.prove(proof_stream, &combination_codeword);

        // the final proof is just the serialized stream
        proof_stream.serialize()
    }

    pub fn verify(
        &mut self,
        proof: &[u8],
        proof_stream: &mut impl ProofStream<P::Fx>,
    ) -> Result<(), &str> {
        let mut proof_stream = proof_stream.deserialize(proof);

        // Get the trace length
        let trace_len = match proof_stream.pull() {
            ProofObject::TraceLen(len) => len as usize,
            _ => return Err("Expected the trace length"),
        };

        self.processor_table.set_height(trace_len);
        self.memory_table.set_height(trace_len);
        self.instruction_table.set_height(trace_len);
        self.input_table.set_height(trace_len);
        self.output_table.set_height(trace_len);

        let codeword_len = self.fri_codeword_length();

        // get the merkle root of the base tables
        let base_root = match proof_stream.pull() {
            ProofObject::MerkleRoot(root) => root,
            _ => return Err("Expected to recieve base columns merkle root"),
        };

        // Get coefficients for table extension
        let challenges = self.sample_weights(11, proof_stream.verifier_fiat_shamir());

        // get root of table extensions
        let extension_root = match proof_stream.pull() {
            ProofObject::MerkleRoot(root) => root,
            _ => return Err("Expected to recieve extension columns root"),
        };

        // get terminals
        let mut terminals = Vec::new();
        for i in 0..5 {
            match proof_stream.pull() {
                ProofObject::Terminal(terminal) => terminals.push(terminal),
                _ => return Err("Expected to receive terminal"),
            }
        }

        let base_degree_bounds = vec![
            vec![self.processor_table.interpolant_degree(); ProcessorTable::<P::Fx>::BASE_WIDTH],
            vec![self.memory_table.interpolant_degree(); MemoryTable::<P::Fx>::BASE_WIDTH],
            vec![
                self.instruction_table.interpolant_degree();
                InstructionTable::<P::Fx>::BASE_WIDTH
            ],
            vec![self.input_table.interpolant_degree(); InputTable::<P::Fx>::BASE_WIDTH],
            vec![self.output_table.interpolant_degree(); OutputTable::<P::Fx>::BASE_WIDTH],
        ]
        .concat();

        let extension_degree_bounds = vec![
            vec![
                self.processor_table.interpolant_degree();
                ProcessorTable::<P::Fx>::EXTENSION_WIDTH - ProcessorTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.memory_table.interpolant_degree();
                MemoryTable::<P::Fx>::EXTENSION_WIDTH - MemoryTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.instruction_table.interpolant_degree();
                InstructionTable::<P::Fx>::EXTENSION_WIDTH - InstructionTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.input_table.interpolant_degree();
                InputTable::<P::Fx>::EXTENSION_WIDTH - InputTable::<P::Fx>::BASE_WIDTH
            ],
            vec![
                self.output_table.interpolant_degree();
                OutputTable::<P::Fx>::EXTENSION_WIDTH - OutputTable::<P::Fx>::BASE_WIDTH
            ],
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

        // get weights for nonlinear combination
        // - 1 randomizer
        // - 2 for every other polynomial
        let num_base_polynomials = ProcessorTable::<P::Fx>::BASE_WIDTH
            + MemoryTable::<P::Fx>::BASE_WIDTH
            + InstructionTable::<P::Fx>::BASE_WIDTH
            + InputTable::<P::Fx>::BASE_WIDTH
            + OutputTable::<P::Fx>::BASE_WIDTH;
        let num_extension_polynomials = ProcessorTable::<P::Fx>::EXTENSION_WIDTH
            + MemoryTable::<P::Fx>::EXTENSION_WIDTH
            + InstructionTable::<P::Fx>::EXTENSION_WIDTH
            + InputTable::<P::Fx>::EXTENSION_WIDTH
            + OutputTable::<P::Fx>::EXTENSION_WIDTH
            - num_base_polynomials;
        let num_randomizer_polynomials = 1;
        // +2 for instr and memory permutations
        let num_quotient_polynomials = quotient_degree_bounds.iter().sum::<usize>() + 2;

        let weights_seed = proof_stream.verifier_fiat_shamir();
        let weights = self.sample_weights(
            num_randomizer_polynomials
                + 2 * (num_base_polynomials + num_extension_polynomials + num_quotient_polynomials),
            weights_seed,
        );

        // pull Merkle root of combination codeword
        let combination_root = match proof_stream.pull() {
            ProofObject::MerkleRoot(root) => root,
            _ => return Err("Expected to receive combination codeword root"),
        };

        // get indices of leafs to verify non-linear combination
        let indices_seed = proof_stream.verifier_fiat_shamir();
        let indices = Self::sample_indices(P::SECURITY_LEVEL, indices_seed, codeword_len);

        let row_step = P::EXPANSION_FACTOR;
        let mut values: HashMap<usize, Vec<P::Fx>> = HashMap::new();
        // open leafs of zipped codewords at indicated positions
        for &index in &indices {
            for distance in [0, row_step] {
                let idx = (index + distance) % codeword_len;
                let entry = values.entry(idx).or_default();
                assert!(entry.is_empty());

                // base codewords
                let elements = match proof_stream.pull() {
                    ProofObject::LeafItems(elements) => elements,
                    _ => return Err("Expected to receive leaf items"),
                };
                let (salt, path) = match proof_stream.pull() {
                    ProofObject::MerklePathWithSalt(v) => v,
                    _ => return Err("Expected to receive a merkle path with a salt"),
                };
                if !SaltedMerkle::verify(base_root, idx, salt, &path, &elements) {
                    return Err("Invalid base codeword path");
                }
                entry.extend_from_slice(&elements);

                // extension codewords
                let elements = match proof_stream.pull() {
                    ProofObject::LeafItems(elements) => elements,
                    _ => return Err("Expected to receive leaf items"),
                };
                let (salt, path) = match proof_stream.pull() {
                    ProofObject::MerklePathWithSalt(v) => v,
                    _ => return Err("Expected to receive a merkle path with a salt"),
                };
                if !SaltedMerkle::verify(extension_root, idx, salt, &path, &elements) {
                    return Err("Invalid extension codeword path");
                }
                entry.extend_from_slice(&elements);
            }
        }

        Ok(())
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

fn sum<F: Field>(a: Vec<F>, b: Vec<F>) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    (0..n).map(|i| a[i] + b[i]).collect()
}

// N:128, randomness:18379656959775521745, bound:65536

// [41632, 4154, 32355, 58925, 35723, 11523, 50548, 10183, 55426, 30236, 9332,
// 57351, 55176, 7459, 6162, 45442, 17456, 49953, 57238, 5900, 60643, 64444,
// 8760, 25804, 16733, 23567, 47786, 36823, 15199, 45294, 32288, 10850, 54885,
// 54791, 32419, 29928, 52087, 10145, 11106, 10331, 26702, 15512, 54250, 13846,
// 18270, 11058, 54121, 14267, 22768, 45353, 20153, 26066, 42700, 53711, 47453,
// 43590, 17467, 32822, 49982, 27828, 42153, 6915, 39507, 2672, 12254, 37317,
// 25340, 22264, 36643, 60181, 42489, 47171, 37530, 30558, 65430, 3833, 6941,
// 55384, 24549, 21018, 1506, 24639, 6026, 25051, 50825, 44697, 27220, 23299,
// 13208, 7293, 56223, 60358, 20707, 53506, 54058, 63002, 35065, 34235, 18013,
// 36328, 51715, 12971, 43149, 44677, 62391, 49205, 18746, 15327, 52197, 19130,
// 23144, 16484, 65084, 27609, 3764, 27627, 18376, 25400, 17676, 18728, 13350,
// 15083, 10647, 11866, 31200, 43080, 8611, 15518]

// [672, 58, 3683, 1581, 2955, 3331, 1396, 1991, 2178, 1564, 1140, 7, 1928,
// 3363, 2066, 386, 1072, 801, 3990, 1804, 3299, 3004, 568, 1228, 349, 3087,
// 2730, 4055, 2911, 238, 3616, 2658, 1637, 1543, 3747, 1256, 2935, 1953, 2914,
// 2139, 2126, 3224, 1002, 1558, 1886, 2866, 873, 1979, 2288, 297, 3769, 1490,
// 1740, 463, 2397, 2630, 1083, 54, 830, 3252, 1193, 2819, 2643, 2672, 4062,
// 453, 764, 1784, 3875, 2837, 1529, 2115, 666, 1886, 3990, 3833, 2845, 2136,
// 4069, 538, 1506, 63, 1930, 475, 1673, 3737, 2644, 2819, 920, 3197, 2975,
// 3014, 227, 258, 810, 1562, 2297, 1467, 1629, 3560, 2563, 683, 2189, 3717,
// 951, 53, 2362, 3039, 3045, 2746, 2664, 100, 3644, 3033, 3764, 3051, 1992,
// 824, 1292, 2344, 1062, 2795, 2455, 3674, 2528, 2120, 419, 3230]
