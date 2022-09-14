use crate::ceil_power_of_two;
use crate::merkle::Merkle;
use crate::protocol::ProofStream;
use algebra::ExtensionOf;
use algebra::Felt;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use std::collections::hash_map::DefaultHasher;
use std::hash;
use std::hash::Hash;
use std::hash::Hasher;
use std::iter::zip;
use std::marker::PhantomData;

pub struct FriParams {
    expansion_factor: usize,
    num_colinearity_tests: usize,
}

impl FriParams {
    pub fn new(expansion_factor: usize, num_colinearity_tests: usize) -> FriParams {
        Self {
            expansion_factor,
            num_colinearity_tests,
        }
    }

    pub fn expansion_factor(&self) -> usize {
        self.expansion_factor
    }

    pub fn num_colinearity_tests(&self) -> usize {
        self.num_colinearity_tests
    }
}

pub struct Fri<F> {
    params: FriParams,
    _phantom: PhantomData<F>,
}

impl<E> Fri<E>
where
    E: Felt,
    E::BaseFelt: StarkFelt,
{
    pub fn new(params: FriParams) -> Self {
        Fri {
            params,
            _phantom: PhantomData,
        }
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
        proof_stream: &mut impl ProofStream<E>,
        codeword: &[E],
    ) -> (Vec<Vec<E>>, Vec<Merkle<E>>) {
        let one = E::one();
        let two = one + one;

        let mut codeword = codeword.to_vec();
        let mut omega = E::BaseFelt::get_root_of_unity(codeword.len().ilog2());
        let mut offset = E::BaseFelt::GENERATOR;

        let mut codewords = Vec::new();
        let mut trees = Vec::new();

        while codeword.len() >= ceil_power_of_two(self.params.num_colinearity_tests())
            && codeword.len() >= self.params.expansion_factor()
        {
            let tree = Merkle::new(&codeword);
            let root = tree.root();

            // Skip the first round
            if trees.len() != 0 {
                // TODO: HELP: is this needed on the last round?
                proof_stream.push(crate::protocol::ProofObject::MerkleRoot(root));
            }

            // only prepare next round if necessary
            if codeword.len() == self.params.expansion_factor() {
                break;
            }

            trees.push(tree);
            codewords.push(codeword.clone());

            // get challenge for split and fold
            let alpha = E::from(proof_stream.prover_fiat_shamir());
            let (lhs, rhs) = codeword.split_at(codeword.len() / 2);
            codeword = zip(lhs, rhs)
                .enumerate()
                .map(|(i, (&l, &r))| {
                    (one + alpha / E::from(offset * omega.pow(&[i as u64])) * l
                        + (one - alpha / E::from(offset * omega.pow(&[i as u64]))) * r)
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
        proof_stream: &mut impl ProofStream<E>,
        curr_tree: &Merkle<E>,
        next_tree: &Merkle<E>,
        indices: &[usize],
    ) {
        let lhs_indices = indices.to_vec();
        let rhs_indices = indices
            .iter()
            .map(|i| i + curr_tree.leafs.len() / 2)
            .collect::<Vec<usize>>();

        // reveal leafs
        for i in 0..self.params.num_colinearity_tests() {
            proof_stream.push(crate::protocol::ProofObject::FriLeafs((
                curr_tree.leafs[lhs_indices[i]],
                curr_tree.leafs[rhs_indices[i]],
                next_tree.leafs[lhs_indices[i]],
            )));
        }

        // reveal authentication paths
        for i in 0..self.params.num_colinearity_tests() {
            println!("WE ARE HERE");
            let mp = crate::protocol::ProofObject::MerklePath(curr_tree.open(lhs_indices[i]).1);
            println!("MERKLE PATH: {}", serde_json::to_string(&mp).unwrap());
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
        proof_stream: &mut impl ProofStream<E>,
        curr_tree: &Merkle<E>,
        last_codeword: &[E],
        indices: &[usize],
    ) {
        let lhs_indices = indices.to_vec();
        let rhs_indices = indices
            .iter()
            .map(|i| i + curr_tree.leafs.len() / 2)
            .collect::<Vec<usize>>();

        // reveal leafs
        for i in 0..self.params.num_colinearity_tests() {
            proof_stream.push(crate::protocol::ProofObject::FriLeafs((
                curr_tree.leafs[lhs_indices[i]],
                curr_tree.leafs[rhs_indices[i]],
                last_codeword[lhs_indices[i]],
            )));
        }

        // reveal authentication paths
        for i in 0..self.params.num_colinearity_tests() {
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(lhs_indices[i]).1,
            ));
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(rhs_indices[i]).1,
            ));
        }
    }

    pub fn prove(&self, proof_stream: &mut impl ProofStream<E>, codeword: &[E]) -> Vec<usize> {
        // commit phase
        let (codewords, trees) = self.commit(proof_stream, codeword);

        // query phase
        let last_codeword = codewords.last().unwrap();
        println!("Last codeword len: {}", last_codeword.len());
        let top_level_indices = self.sample_indices(
            self.params.num_colinearity_tests(),
            ceil_power_of_two(self.params.num_colinearity_tests()),
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
