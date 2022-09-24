#![feature(generic_const_exprs, int_log)]

use crate::merkle::Merkle;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::UniformRand;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use brainfuck::permutation_argument;
use brainfuck::InputTable;
use brainfuck::InstructionTable;
use brainfuck::MemoryTable;
use brainfuck::OutputTable;
use brainfuck::ProcessorTable;
use brainfuck::Table;
use fri::Fri;
use legacy_algebra::number_theory_transform::inverse_number_theory_transform;
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
use std::iter::zip;
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
    program: Vec<usize>,
    fri: Fri<P>,
    processor_table: ProcessorTable<P::Fx>,
    memory_table: MemoryTable<P::Fx>,
    instruction_table: InstructionTable<P::Fx>,
    input_table: InputTable<P::Fx>,
    output_table: OutputTable<P::Fx>,
}

impl<P: Config> BrainFuckStark<P> {
    pub fn new(params: P, program: Vec<usize>) -> Self {
        BrainFuckStark {
            program,
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

    pub fn sample_indices(n: usize, randomness: u64, bound: usize) -> Vec<usize> {
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

        proof_stream.push(ProofObject::PublicInputs(input_matrix.clone()));
        proof_stream.push(ProofObject::PublicOutputs(output_matrix.clone()));

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

        let all_base_codewords =
            vec![randomizer_codewords.clone(), base_codewords.clone()].concat();

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
            .map(|i| {
                all_base_codewords
                    .iter()
                    .map(|codeword| codeword[i])
                    .collect()
            })
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
        // let combination_codeword = terms
        //     .into_iter()
        //     .zip(weights)
        //     .map(|(term, weight)| term.into_iter().map(|v| v * weight).collect())
        //     .fold(vec![P::Fx::zero(); codeword_len], sum);
        // let combination_codeword = terms
        //     .iter()
        //     .map(|term| P::Fx::sum_of_products(term, &weights))
        //     .fold(vec![P::Fx::zero(); codeword_len], P::Fx::sum);
        //     .collect::<Vec<P::Fx>>();
        let combination_codeword = terms
            .iter()
            .zip(weights)
            .map(|(term, weight)| term.iter().copied().map(|v| v * weight).collect())
            .fold(vec![P::Fx::zero(); codeword_len], sum);

        {
            println!("TRYING TO INTERPOLATE");
            let coefficients = inverse_number_theory_transform(&combination_codeword);
            let my_poly = DensePolynomial::from_coefficients_vec(coefficients);
            println!("MY POLY DEGREE {}", my_poly.degree());
            println!("CODEWORD LEN {}", codeword_len);
        }

        let combination_tree = Merkle::new(&combination_codeword);
        proof_stream.push(ProofObject::MerkleRoot(combination_tree.root()));

        // get indices of leafs to prove non-linear combination
        let indices_seed = proof_stream.prover_fiat_shamir();
        let indices = Self::sample_indices(P::SECURITY_LEVEL, indices_seed, codeword_len);

        let row_step = codeword_len / padding_length;
        println!("Row step is: {}", row_step);
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
        let mut is_first = true;
        for index in indices {
            if is_first {
                let individual_terms = terms.iter().map(|term| term[index]).collect::<Vec<P::Fx>>();
                println!("{:?}", individual_terms);
                is_first = false;
            }
            let (element, path) = combination_tree.open(index);
            assert!(Merkle::verify(
                combination_tree.root(),
                index,
                &path,
                &element
            ));
            // Not needed. Verifier can re-compute the leaf
            // proof_stream.push(ProofObject::LeafItem(element));
            proof_stream.push(ProofObject::MerklePath(path));
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
        input: &[usize],
        output: &[usize],
    ) -> Result<(), &str>
    where
        [(); InputTable::<P::Fx>::BASE_WIDTH]: Sized,
    {
        let mut proof_stream = proof_stream.deserialize(proof);

        // TODO: refactor the how these get used by the verifier
        //      makes V have to allocate so much memory defeats the point
        // Get the public inputs
        let public_inputs = match proof_stream.pull() {
            ProofObject::PublicInputs(inputs) => inputs,
            _ => return Err("Expected the trace length"),
        };

        // Get the public outputs
        let public_outputs = match proof_stream.pull() {
            ProofObject::PublicOutputs(outputs) => outputs,
            _ => return Err("Expected the trace length"),
        };

        // Get the trace length
        let trace_len = match proof_stream.pull() {
            ProofObject::TraceLen(len) => len as usize,
            _ => return Err("Expected the trace length"),
        };

        self.input_table.set_matrix(public_inputs);
        self.input_table.pad(trace_len);
        self.output_table.set_matrix(public_outputs);
        self.output_table.pad(trace_len);

        self.processor_table.set_height(trace_len);
        self.memory_table.set_height(trace_len);
        self.instruction_table.set_height(trace_len);

        let codeword_len = self.fri_codeword_length();

        println!("Codeword shlength: {}", codeword_len);

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
        for _ in 0..5 {
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

        let quotient_degree_bounds = vec![
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
        // *difference quotients TODO learn more about these
        let num_quotient_polynomials = quotient_degree_bounds.len() + 2;

        let weights_seed = proof_stream.verifier_fiat_shamir();
        println!("{num_randomizer_polynomials} + 2 * ({num_base_polynomials} + {num_extension_polynomials} + {num_quotient_polynomials})");
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

        let row_step = codeword_len / trace_len;
        let mut values: HashMap<usize, Vec<P::Fx>> = HashMap::new();
        // open leafs of zipped codewords at indicated positions
        for &index in &indices {
            for distance in [0, row_step] {
                let idx = (index + distance) % codeword_len;
                let entry = values.entry(idx).or_default();
                // assert!(entry.is_empty());
                // BUG: sample indices should not sample row_step
                let is_empty = entry.is_empty();

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
                if !is_empty {
                    assert_eq!(elements[0], entry[0]);
                } else {
                    entry.extend_from_slice(&elements);
                }

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
                if is_empty {
                    entry.extend_from_slice(&elements);
                }
            }
        }

        assert_eq!(num_base_polynomials, base_degree_bounds.len());

        let offset = P::Fp::GENERATOR;
        let omega = P::Fp::get_root_of_unity(codeword_len as u64).unwrap();
        let omicron = P::Fp::get_root_of_unity(trace_len as u64).unwrap();

        let base_widths = vec![
            ProcessorTable::<P::Fx>::BASE_WIDTH,
            MemoryTable::<P::Fx>::BASE_WIDTH,
            InstructionTable::<P::Fx>::BASE_WIDTH,
            InputTable::<P::Fx>::BASE_WIDTH,
            OutputTable::<P::Fx>::BASE_WIDTH,
        ];

        // verify non-linear combination
        for index in &indices {
            let entry = values.remove(index).unwrap();

            assert_eq!(
                entry.len(),
                num_base_polynomials + num_extension_polynomials + 1
            );

            // collect terms: randomizer
            let mut terms = entry[0..num_randomizer_polynomials].to_vec();

            // collect terms: base
            for i in num_randomizer_polynomials..num_randomizer_polynomials + num_base_polynomials {
                terms.push(entry[i]);
                let shift = self.max_degree() - base_degree_bounds[i - num_randomizer_polynomials];
                terms.push(
                    entry[i]
                        * P::Fx::from_base_prime_field(
                            (offset * omega.pow([*index as u64])).pow([shift as u64]),
                        ),
                );
            }

            // collect terms: extension
            let extension_offset = num_randomizer_polynomials + num_base_polynomials;

            assert_eq!(
                terms.len(),
                2 * extension_offset - num_randomizer_polynomials
            );

            assert_eq!(
                extension_offset,
                base_widths.iter().sum::<usize>() + num_randomizer_polynomials
            );

            for i in 0..num_extension_polynomials {
                terms.push(entry[extension_offset + i]);
                let shift = self.max_degree() - extension_degree_bounds[i];
                terms.push(
                    entry[extension_offset + i]
                        * P::Fx::from_base_prime_field(
                            (offset * omega.pow([*index as u64])).pow([shift as u64]),
                        ),
                );
            }

            // collect terms: quotients
            // quotients need to be computed
            // skip randomizers
            let mut entry_iter = entry.into_iter().skip(num_randomizer_polynomials);
            let mut points: Vec<Vec<P::Fx>> = vec![
                (0..ProcessorTable::<P::Fx>::BASE_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect(),
                (0..MemoryTable::<P::Fx>::BASE_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect(),
                (0..InstructionTable::<P::Fx>::BASE_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect(),
                (0..InputTable::<P::Fx>::BASE_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect(),
                (0..OutputTable::<P::Fx>::BASE_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect(),
            ];

            points[0].extend_from_slice(
                &(ProcessorTable::<P::Fx>::BASE_WIDTH..ProcessorTable::<P::Fx>::EXTENSION_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect::<Vec<P::Fx>>(),
            );

            points[1].extend_from_slice(
                &(MemoryTable::<P::Fx>::BASE_WIDTH..MemoryTable::<P::Fx>::EXTENSION_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect::<Vec<P::Fx>>(),
            );

            points[2].extend_from_slice(
                &(InstructionTable::<P::Fx>::BASE_WIDTH
                    ..InstructionTable::<P::Fx>::EXTENSION_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect::<Vec<P::Fx>>(),
            );

            points[3].extend_from_slice(
                &(InputTable::<P::Fx>::BASE_WIDTH..InputTable::<P::Fx>::EXTENSION_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect::<Vec<P::Fx>>(),
            );

            points[4].extend_from_slice(
                &(OutputTable::<P::Fx>::BASE_WIDTH..OutputTable::<P::Fx>::EXTENSION_WIDTH)
                    .map(|_| entry_iter.next().unwrap())
                    .collect::<Vec<P::Fx>>(),
            );

            let boundary_constraints = vec![
                ProcessorTable::<P::Fx>::extension_boundary_constraints(&challenges),
                MemoryTable::<P::Fx>::extension_boundary_constraints(&challenges),
                InstructionTable::<P::Fx>::extension_boundary_constraints(&challenges),
                InputTable::<P::Fx>::extension_boundary_constraints(&challenges),
                OutputTable::<P::Fx>::extension_boundary_constraints(&challenges),
            ];

            let boundary_quotient_degree_bounds = vec![
                self.processor_table
                    .boundary_quotient_degree_bounds(&challenges),
                self.memory_table
                    .boundary_quotient_degree_bounds(&challenges),
                self.instruction_table
                    .boundary_quotient_degree_bounds(&challenges),
                self.input_table
                    .boundary_quotient_degree_bounds(&challenges),
                self.output_table
                    .boundary_quotient_degree_bounds(&challenges),
            ];

            let transition_constraints = vec![
                ProcessorTable::<P::Fx>::extension_transition_constraints(&challenges),
                MemoryTable::<P::Fx>::extension_transition_constraints(&challenges),
                InstructionTable::<P::Fx>::extension_transition_constraints(&challenges),
                InputTable::<P::Fx>::extension_transition_constraints(&challenges),
                OutputTable::<P::Fx>::extension_transition_constraints(&challenges),
            ];

            let transition_quotient_degree_bounds = vec![
                self.processor_table
                    .transition_quotient_degree_bounds(&challenges),
                self.memory_table
                    .transition_quotient_degree_bounds(&challenges),
                self.instruction_table
                    .transition_quotient_degree_bounds(&challenges),
                self.input_table
                    .transition_quotient_degree_bounds(&challenges),
                self.output_table
                    .transition_quotient_degree_bounds(&challenges),
            ];

            let terminal_constraints_ext = vec![
                self.processor_table
                    .extension_terminal_constraints(&challenges, &terminals),
                self.memory_table
                    .extension_terminal_constraints(&challenges, &terminals),
                self.instruction_table
                    .extension_terminal_constraints(&challenges, &terminals),
                self.input_table
                    .extension_terminal_constraints(&challenges, &terminals),
                self.output_table
                    .extension_terminal_constraints(&challenges, &terminals),
            ];

            let terminal_quotient_degree_bounds = vec![
                self.processor_table
                    .terminal_quotient_degree_bounds(&challenges, &terminals),
                self.memory_table
                    .terminal_quotient_degree_bounds(&challenges, &terminals),
                self.instruction_table
                    .terminal_quotient_degree_bounds(&challenges, &terminals),
                self.input_table
                    .terminal_quotient_degree_bounds(&challenges, &terminals),
                self.output_table
                    .terminal_quotient_degree_bounds(&challenges, &terminals),
            ];

            let extension_widths = vec![
                ProcessorTable::<P::Fx>::EXTENSION_WIDTH,
                MemoryTable::<P::Fx>::EXTENSION_WIDTH,
                InstructionTable::<P::Fx>::EXTENSION_WIDTH,
                InputTable::<P::Fx>::EXTENSION_WIDTH,
                OutputTable::<P::Fx>::EXTENSION_WIDTH,
            ];

            let mut base_acc_index = num_randomizer_polynomials;
            let mut ext_acc_index = extension_offset;
            for (
                (
                    (
                        (
                            (
                                (
                                    (
                                        (point, boundary_constraints),
                                        boundary_quotient_degree_bounds,
                                    ),
                                    transition_constraints,
                                ),
                                transition_quotient_degree_bounds,
                            ),
                            terminal_constraints_ext,
                        ),
                        terminal_quotient_degree_bounds,
                    ),
                    base_width,
                ),
                extension_width,
            ) in points
                .iter()
                .zip(&boundary_constraints)
                .zip(&boundary_quotient_degree_bounds)
                .zip(&transition_constraints)
                .zip(&transition_quotient_degree_bounds)
                .zip(&terminal_constraints_ext)
                .zip(&terminal_quotient_degree_bounds)
                .zip(&base_widths)
                .zip(&extension_widths)
            {
                for (constraint, bound) in
                    zip(boundary_constraints, boundary_quotient_degree_bounds)
                {
                    let eval = constraint.evaluate(point);
                    let quotient = eval
                        / (P::Fx::from_base_prime_field(offset * omega.pow([*index as u64]))
                            - P::Fx::one());
                    terms.push(quotient);
                    let shift = self.max_degree() - bound;
                    terms.push(
                        quotient
                            * P::Fx::from_base_prime_field(
                                (offset * omega.pow([*index as u64])).pow([shift as u64]),
                            ),
                    );
                }

                // transition
                let next_index = (index + row_step) % codeword_len;
                let mut next_point = values.get(&next_index).unwrap()
                    [base_acc_index..base_acc_index + base_width]
                    .to_vec();
                next_point.extend_from_slice(
                    &values.get(&next_index).unwrap()
                        [ext_acc_index..ext_acc_index + extension_width - base_width],
                );
                base_acc_index += base_width;
                ext_acc_index += extension_width - base_width;
                for (constraint, bound) in
                    zip(transition_constraints, transition_quotient_degree_bounds)
                {
                    let eval_point = vec![point.clone(), next_point.clone()].concat();
                    let eval = constraint.evaluate(&eval_point);
                    // If the trace length is 0 then there is no subgroup where the transition
                    // polynomials should be zero. The fast zerofier (based on
                    // group theory) needs a non-empty group. Forcing it on an
                    // empty group generates a division by zero error.
                    let quotient = if trace_len == 0 {
                        P::Fx::zero()
                    } else {
                        // transition quotients apply to all
                        let omicron_inv = omicron.inverse().unwrap();
                        let evaluation_domain = offset * omega.pow([*index as u64]);
                        // TODO: get a grasp on codomain intricacies
                        // To get the quotient we evaluate the transition constraints at the point
                        // and divide out the zerofier i.e. = eval_result /
                        // ((ω^i - 1)(ω^i - o)(ω^i - o^2)...(ω^i - o^(n - 2)))
                        // Observe the group theory based zerofier gives us:
                        // ((ω^i - o^0)...(ω^i - o^(n - 1))) = (ω^i)^n - 1
                        // using this the formula is now much easier to compute
                        // = eval_result * (ω^i - o^(n - 1)) / ((ω^i)^n - 1)
                        // = eval_result * (ω^i - 1/o)) / ((ω^i)^n - 1)
                        eval * P::Fx::from_base_prime_field(evaluation_domain - omicron_inv)
                            / (P::Fx::from_base_prime_field(
                                evaluation_domain.pow([trace_len as u64]),
                            ) - P::Fx::one())
                    };
                    terms.push(quotient);
                    let shift = self.max_degree() - bound;
                    terms.push(
                        quotient
                            * P::Fx::from_base_prime_field(
                                (offset * omega.pow([*index as u64])).pow([shift as u64]),
                            ),
                    )
                }

                // terminal
                for (constraint, bound) in
                    zip(terminal_constraints_ext, terminal_quotient_degree_bounds)
                {
                    let eval = constraint.evaluate(point);
                    let quotient = eval
                        / (P::Fx::from_base_prime_field(
                            offset * omega.pow([*index as u64]) - omicron.inverse().unwrap(),
                        ));
                    terms.push(quotient);
                    let shift = self.max_degree() - bound;
                    terms.push(
                        quotient
                            * P::Fx::from_base_prime_field(
                                (offset * omega.pow([*index as u64])).pow([shift as u64]),
                            ),
                    );
                }
            }

            // Instruction permutation
            let quotient = (points[0][ProcessorTable::<P::Fx>::INSTRUCTION_PERMUTATION]
                - points[2][InstructionTable::<P::Fx>::PROCESSOR_PERMUTATION])
                / (P::Fx::from_base_prime_field(offset * omega.pow([*index as u64]))
                    - P::Fx::one());
            terms.push(quotient);
            let degree_bound = std::cmp::max(
                self.processor_table.interpolant_degree(),
                self.instruction_table.interpolant_degree(),
            ) - 1;
            let shift = self.max_degree() - degree_bound;
            terms.push(
                quotient
                    * P::Fx::from_base_prime_field(
                        (offset * omega.pow([*index as u64])).pow([shift as u64]),
                    ),
            );

            // Memory permutation
            let quotient = (points[0][ProcessorTable::<P::Fx>::MEMORY_PERMUTATION]
                - points[1][MemoryTable::<P::Fx>::PERMUTATION])
                / (P::Fx::from_base_prime_field(offset * omega.pow([*index as u64]))
                    - P::Fx::one());
            terms.push(quotient);
            let degree_bound = std::cmp::max(
                self.processor_table.interpolant_degree(),
                self.memory_table.interpolant_degree(),
            ) - 1;
            let shift = self.max_degree() - degree_bound;
            terms.push(
                quotient
                    * P::Fx::from_base_prime_field(
                        (offset * omega.pow([*index as u64])).pow([shift as u64]),
                    ),
            );

            assert_eq!(terms.len(), weights.len());

            // compute inner product of weights and terms
            let inner_product = P::Fx::sum_of_products(&weights, &terms);

            // compare inner product against the combination codeword value
            // Not needed. Verifier can re-compute the leaf
            // let combination_leaf = match proof_stream.pull() {
            //     ProofObject::LeafItem(leaf) => leaf,
            //     _ => return Err("Expected to receive combination codeword leaf"),
            // };
            let combination_path = match proof_stream.pull() {
                ProofObject::MerklePath(path) => path,
                _ => return Err("Expected to receive combination codeword path"),
            };

            if !Merkle::verify(combination_root, *index, &combination_path, &inner_product) {
                return Err("Failed combination codeword verification");
            }
        }

        // verify low degree of combination polynomial
        self.fri
            .verify(&mut proof_stream, codeword_len, combination_root)?;

        // verify external terminals:
        // for ea in self.evaluation_arguments:
        // ea.select_terminal(terminals) == ea.compute_terminal(challenges)

        // verifier the external terminals
        let mut terminals = terminals.into_iter();
        let _processor_instruction_permutation_terminal = terminals.next().unwrap();
        let _processor_memory_permutation_terminal = terminals.next().unwrap();
        let processor_input_evaluation_terminal = terminals.next().unwrap();
        let processor_output_evaluation_terminal = terminals.next().unwrap();
        let instruction_evaluation_terminal = terminals.next().unwrap();

        if processor_input_evaluation_terminal
            != brainfuck::evaluation_argument::compute_io_terminal(input, challenges[8])
        {
            return Err("Failed to verify input evaluation argument");
        }

        if processor_output_evaluation_terminal
            != brainfuck::evaluation_argument::compute_io_terminal(output, challenges[9])
        {
            return Err("Failed to verify output evaluation argument");
        }

        if instruction_evaluation_terminal
            != brainfuck::evaluation_argument::compute_program_terminal(&self.program, &challenges)
        {
            return Err("Failed to verify program evaluation argument");
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
