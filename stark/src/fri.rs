use crate::ceil_power_of_two;
use crate::merkle::Merkle;
use crate::protocol::ProofStream;
use legacy_algebra::fp_u128::BaseFelt;
use legacy_algebra::ExtensionOf;
use legacy_algebra::Felt;
use legacy_algebra::PrimeFelt;
use legacy_algebra::StarkFelt;
use num_traits::One;
use std::collections::hash_map::DefaultHasher;
use std::hash;
use std::hash::Hash;
use std::hash::Hasher;
use std::iter::zip;
use std::marker::PhantomData;

pub trait Config {
    type BaseFelt: PrimeFelt + StarkFelt;
    type ExtensionFelt: Felt<BaseFelt = Self::BaseFelt> + ExtensionOf<Self::BaseFelt>;
    // type MerkleParams: ark_crypto_primitives::merkle_tree::Config;

    const EXPANSION_FACTOR: usize;
    const SECURITY_LEVEL: usize;
    const NUM_COLINEARITY_CHECKS: usize =
        Self::SECURITY_LEVEL / Self::EXPANSION_FACTOR.ilog2() as usize;
}

// pub struct FriParams {
//     EXPANSION_FACTOR: usize,
//     num_colinearity_tests: usize,
// }

// impl FriParams {
//     pub fn new(EXPANSION_FACTOR: usize, num_colinearity_tests: usize) ->
// FriParams {         Self {
//             EXPANSION_FACTOR,
//             num_colinearity_tests,
//         }
//     }

//     pub fn EXPANSION_FACTOR(&self) -> usize {
//         self.EXPANSION_FACTOR
//     }

//     pub fn num_colinearity_tests(&self) -> usize {
//         self.num_colinearity_tests
//     }
// }

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
        proof_stream: &mut impl ProofStream<P::ExtensionFelt>,
        codeword: &[P::ExtensionFelt],
    ) -> (Vec<Vec<P::ExtensionFelt>>, Vec<Merkle<P::ExtensionFelt>>) {
        let one = P::ExtensionFelt::one();
        let two = one + one;

        let mut codeword = codeword.to_vec();
        let mut omega = P::BaseFelt::get_root_of_unity(codeword.len().ilog2());
        let mut offset = P::BaseFelt::GENERATOR;

        let mut codewords = Vec::new();
        let mut trees = Vec::new();

        while codeword.len() >= ceil_power_of_two(P::NUM_COLINEARITY_CHECKS)
            && codeword.len() >= P::EXPANSION_FACTOR
        {
            let tree = Merkle::new(&codeword);
            let root = tree.root();

            // Skip the first round
            if trees.len() != 0 {
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
            let alpha = P::ExtensionFelt::from(proof_stream.prover_fiat_shamir());
            let (lhs, rhs) = codeword.split_at(codeword.len() / 2);
            codeword = zip(lhs, rhs)
                .enumerate()
                .map(|(i, (&l, &r))| {
                    (one + alpha / P::ExtensionFelt::from(offset * omega.pow(&[i as u64])) * l
                        + (one - alpha / P::ExtensionFelt::from(offset * omega.pow(&[i as u64])))
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
        proof_stream: &mut impl ProofStream<P::ExtensionFelt>,
        curr_tree: &Merkle<P::ExtensionFelt>,
        next_tree: &Merkle<P::ExtensionFelt>,
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
        proof_stream: &mut impl ProofStream<P::ExtensionFelt>,
        curr_tree: &Merkle<P::ExtensionFelt>,
        last_codeword: &[P::ExtensionFelt],
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
        proof_stream: &mut impl ProofStream<P::ExtensionFelt>,
        codeword: &[P::ExtensionFelt],
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
