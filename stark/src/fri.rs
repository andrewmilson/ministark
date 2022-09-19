use crate::ceil_power_of_two;
use crate::merkle::Merkle;
use crate::protocol::ProofStream;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::PrimeField;
use num_traits::One;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;
use std::iter::zip;

pub trait Config {
    /// Base prime field
    type Fp: PrimeField + FftField;
    /// Extension field element
    type Fx: Field<BasePrimeField = Self::Fp>;

    const EXPANSION_FACTOR: usize;
    const SECURITY_LEVEL: usize;
    const NUM_COLINEARITY_CHECKS: usize =
        Self::SECURITY_LEVEL / Self::EXPANSION_FACTOR.ilog2() as usize;
}

pub struct Fri<P: Config> {
    params: P,
}

impl<P: Config> Fri<P> {
    pub fn new(params: P) -> Self {
        Fri { params }
    }

    fn sample_indices(
        &self,
        n: usize,
        reduced_size: usize,
        max: usize,
        randomness: u64,
    ) -> Vec<usize> {
        assert!(n <= reduced_size);
        let mut indices = Vec::new();
        let mut reduced_indices = vec![false; reduced_size];
        let mut counter = 0;
        while indices.len() < n {
            let mut hasher = DefaultHasher::new();
            randomness.hash(&mut hasher);
            counter.hash(&mut hasher);
            let hash = hasher.finish();
            let index = hash as usize % max;
            let reduced_index = index % reduced_size;
            if !reduced_indices[reduced_index] {
                indices.push(index);
                reduced_indices[reduced_index] = true;
            }
            counter += 1;
        }
        indices
    }

    pub fn commit(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        codeword: &[P::Fx],
    ) -> (Vec<Vec<P::Fx>>, Vec<Merkle<P::Fx>>) {
        let one = P::Fx::one();
        let two = one + one;

        let mut codeword = codeword.to_vec();
        let mut omega = P::Fp::get_root_of_unity(codeword.len() as u64).unwrap();
        let mut offset = P::Fp::GENERATOR;

        let mut codewords = Vec::new();
        let mut trees = Vec::new();

        while codeword.len() >= ceil_power_of_two(P::NUM_COLINEARITY_CHECKS)
            && codeword.len() >= P::EXPANSION_FACTOR
        {
            let tree = Merkle::new(&codeword);
            let root = tree.root();

            // Skip the first round
            if !trees.is_empty() {
                // TODO: HELP: is this needed on the last round?
                proof_stream.push(crate::protocol::ProofObject::MerkleRoot(root));
            }

            // only prepare next round if necessary
            if codeword.len() == P::EXPANSION_FACTOR {
                break;
            }

            trees.push(tree);
            codewords.push(codeword.clone());

            // get challenge for split and fold
            let alpha = P::Fx::from(proof_stream.prover_fiat_shamir());
            let (lhs, rhs) = codeword.split_at(codeword.len() / 2);
            let n = codeword.len();
            codeword = zip(lhs, rhs)
                .enumerate()
                .map(|(i, (&l, &r))| {
                    (one + alpha / P::Fx::from_base_prime_field(offset * omega.pow([i as u64])) * l
                        + (one
                            - alpha
                                / P::Fx::from_base_prime_field(
                                    offset * omega.pow([(n / 2 + i) as u64]),
                                ))
                            * r)
                        / two
                })
                .collect();

            omega.square_in_place();
            offset.square_in_place();
        }

        // send last codeword
        proof_stream.push(crate::protocol::ProofObject::Codeword(codeword.clone()));
        codewords.push(codeword);

        (codewords, trees)
    }

    pub fn query(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        curr_tree: &Merkle<P::Fx>,
        next_tree: &Merkle<P::Fx>,
        indices: &[usize],
    ) {
        let lhs_indices = indices.to_vec();
        let rhs_indices = indices
            .iter()
            .map(|i| i + curr_tree.leafs.len() / 2)
            .collect::<Vec<usize>>();

        // reveal leafs
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            proof_stream.push(crate::protocol::ProofObject::FriLeafs((
                curr_tree.leafs[lhs_indices[i]],
                curr_tree.leafs[rhs_indices[i]],
                next_tree.leafs[lhs_indices[i]],
            )));
        }

        // reveal authentication paths
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            let mp = crate::protocol::ProofObject::MerklePath(curr_tree.open(lhs_indices[i]).1);
            proof_stream.push(mp);
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(rhs_indices[i]).1,
            ));
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                next_tree.open(lhs_indices[i]).1,
            ));
        }
    }

    pub fn query_last(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        curr_tree: &Merkle<P::Fx>,
        last_codeword: &[P::Fx],
        indices: &[usize],
    ) {
        let lhs_indices = indices.to_vec();
        let rhs_indices = indices
            .iter()
            .map(|i| i + curr_tree.leafs.len() / 2)
            .collect::<Vec<usize>>();

        // reveal leafs
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            proof_stream.push(crate::protocol::ProofObject::FriLeafs((
                curr_tree.leafs[lhs_indices[i]],
                curr_tree.leafs[rhs_indices[i]],
                last_codeword[lhs_indices[i]],
            )));
        }

        // reveal authentication paths
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(lhs_indices[i]).1,
            ));
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(rhs_indices[i]).1,
            ));
        }
    }

    pub fn prove(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        codeword: &[P::Fx],
    ) -> Vec<usize> {
        // commit phase
        let (codewords, trees) = self.commit(proof_stream, codeword);

        // query phase
        let last_codeword = codewords.last().unwrap();
        println!("Codewords: {}", codewords.len());
        println!("Last codeword len: {}", last_codeword.len());
        let top_level_indices = self.sample_indices(
            P::NUM_COLINEARITY_CHECKS,
            ceil_power_of_two(P::NUM_COLINEARITY_CHECKS),
            codeword.len() / 2,
            proof_stream.prover_fiat_shamir(),
        );
        println!("Hello!");
        for i in 0..trees.len() - 1 {
            let indices = top_level_indices
                .iter()
                .map(|index| index % (codewords[i].len() / 2))
                .collect::<Vec<usize>>();
            self.query(proof_stream, &trees[i], &trees[i + 1], &indices);
        }
        let indices = top_level_indices
            .iter()
            .map(|index| index % last_codeword.len())
            .collect::<Vec<usize>>();
        self.query_last(proof_stream, trees.last().unwrap(), last_codeword, &indices);

        top_level_indices
    }
}
