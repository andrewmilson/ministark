use crate::hash::Digest;
use crate::hash::ElementHashFn;
use crate::hash::HashFn;
use crate::Matrix;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use snafu::Snafu;
use std::collections::VecDeque;
use std::fmt::Debug;
use std::iter::zip;
use std::marker::PhantomData;

/// Merkle tree error
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(display("tree must contain `{min}` leaves, but `{actual}` were provided"))]
    TooFewLeaves { min: usize, actual: usize },
    #[snafu(display("number of leaves must be a power of two, but `{n}` were provided"))]
    NumberOfLeavesNotPowerOfTwo { n: usize },
    #[snafu(display("leaf index `{i}` cannot exceed the number of leaves (`{n}`)"))]
    LeafIndexOutOfBounds { i: usize, n: usize },
    #[snafu(display("proof is invalid"))]
    InvalidProof,
}

pub trait MerkleTree: Sized + Send + Sync + Clone {
    type Proof: CanonicalSerialize + CanonicalDeserialize + Clone + Send + Sync;
    type Root: Digest;

    /// Returns the root of the merkle tree
    fn root(&self) -> Self::Root;

    /// Generates a merkle proof
    ///
    /// # Errors
    ///
    /// Returns an error if the leaf index is out of bounds.
    fn prove(&self, indices: &[usize]) -> Result<Self::Proof, Error>;

    /// Verifies a merkle proof
    ///
    /// # Errors
    ///
    /// This function returns an error if the proof fails verification.
    fn verify(root: &Self::Root, proof: Self::Proof, indices: &[usize]) -> Result<(), Error>;

    /// Returns the number of security bits
    fn security_level_bits() -> u32;
}

// TODO: all these merkle tree abstractions are way out of control. need to
// refactor
pub trait MerkleTreeConfig: Send + Sync + Sized + 'static {
    type Digest: Digest;
    type Leaf: CanonicalDeserialize + CanonicalSerialize + Clone + Send + Sync + Sized + 'static;

    fn hash_leaves(depth: u32, l0: &Self::Leaf, l1: &Self::Leaf) -> Self::Digest;

    fn hash_nodes(depth: u32, n0: &Self::Digest, n1: &Self::Digest) -> Self::Digest;

    fn security_level_bits() -> u32;
}

/// Merkle View contains information needed to verify multiple Merkle paths.
///
/// Inspired by Starkware's Solidity verifier
/// <https://etherscan.io/address/0xe9664D230490d5A515ef7Ef30033d8075a8D0E24#code#F24#L1>
#[derive(Debug, Clone, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub struct MerkleView<
    N: CanonicalDeserialize + CanonicalSerialize + Clone,
    L: CanonicalDeserialize + CanonicalSerialize + Clone,
> {
    pub nodes: Vec<N>,
    pub initial_leaves: Vec<L>,
    pub sibling_leaves: Vec<L>,
    pub height: u32,
}

/// Merkle tree implemented as a full power-of-two arity tree.
///
/// ```text
///       #        <- root node
///     /   \
///   #       #    <- nodes
///  / \     / \
/// +   +   +   +  <- leaves
/// ```
pub struct MerkleTreeImpl<C: MerkleTreeConfig> {
    pub nodes: Vec<C::Digest>,
    pub leaves: Vec<C::Leaf>,
}

impl<C: MerkleTreeConfig> Clone for MerkleTreeImpl<C> {
    fn clone(&self) -> Self {
        Self {
            nodes: self.nodes.clone(),
            leaves: self.leaves.clone(),
        }
    }
}

impl<C: MerkleTreeConfig> MerkleTreeImpl<C> {
    /// # Errors
    ///
    /// This function will return an error if:
    /// * there are less than two leaves
    /// * the number of leaves is not a power of two
    pub fn new(leaves: Vec<C::Leaf>) -> Result<Self, Error> {
        const MIN_LEAVES: usize = 2;

        let n = leaves.len();
        if n < MIN_LEAVES {
            return Err(Error::TooFewLeaves {
                min: MIN_LEAVES,
                actual: n,
            });
        } else if !n.is_power_of_two() {
            return Err(Error::NumberOfLeavesNotPowerOfTwo { n });
        }

        let nodes = build_merkle_nodes::<C>(&leaves);
        Ok(Self { nodes, leaves })
    }

    /// Returns the height of the merkle tree
    /// i.e. for the merkle tree below `height=1`
    /// ```text
    ///   +
    ///  / \
    /// +   +
    /// ```
    fn height(&self) -> u32 {
        self.leaves.len().ilog2()
    }
}

impl<C: MerkleTreeConfig> MerkleTree for MerkleTreeImpl<C> {
    type Proof = MerkleView<C::Digest, C::Leaf>;
    type Root = C::Digest;

    fn root(&self) -> C::Digest {
        self.nodes[1].clone()
    }

    fn prove(&self, indices: &[usize]) -> Result<MerkleView<C::Digest, C::Leaf>, Error> {
        let num_leaves = self.leaves.len();
        for &i in indices {
            if i >= num_leaves {
                return Err(Error::LeafIndexOutOfBounds { i, n: num_leaves });
            }
        }

        let mut indices = indices.to_vec();
        indices.sort_unstable();
        indices.dedup();

        let mut initial_leaves = Vec::new();
        let mut sibling_leaves = Vec::new();

        // handle leaves and specify the internal node indices
        let mut node_queue = VecDeque::new();
        let mut leaf_queue = VecDeque::from_iter(indices);
        while let Some(index) = leaf_queue.pop_front() {
            initial_leaves.push(self.leaves[index].clone());
            node_queue.push_back((num_leaves + index) >> 1);

            if let Some(&next_index) = leaf_queue.front() {
                let are_siblings = index ^ 1 == next_index;
                if are_siblings {
                    initial_leaves.push(self.leaves[next_index].clone());
                    leaf_queue.pop_front();
                    continue;
                }
            }

            sibling_leaves.push(self.leaves[index ^ 1].clone());
        }

        // handle internal nodes
        let mut nodes = Vec::new();
        while let Some(index) = node_queue.pop_front() {
            if index > 2 {
                node_queue.push_back(index >> 1);
            }

            if let Some(next_index) = node_queue.front() {
                let are_siblings = index ^ 1 == *next_index;
                if are_siblings {
                    node_queue.pop_front();
                    continue;
                }
            }

            nodes.push(self.nodes[index ^ 1].clone());
        }

        Ok(MerkleView {
            nodes,
            initial_leaves,
            sibling_leaves,
            height: self.height(),
        })
    }

    fn verify(
        root: &C::Digest,
        proof: MerkleView<C::Digest, C::Leaf>,
        indices: &[usize],
    ) -> Result<(), Error> {
        let height = proof.height;
        let num_leaves = 1 << height;
        for &i in indices {
            if i >= num_leaves {
                return Err(Error::LeafIndexOutOfBounds { i, n: num_leaves });
            }
        }

        let mut indices = indices.to_vec();
        indices.sort_unstable();
        indices.dedup();

        // handle leaves and specify the internal node indices
        let mut node_queue = VecDeque::new();
        let mut siblings = VecDeque::from_iter(proof.sibling_leaves);
        let mut leaf_queue = zip(indices, proof.initial_leaves).collect::<VecDeque<_>>();
        while let Some((index, leaf)) = leaf_queue.pop_front() {
            let node_index = (num_leaves + index) >> 1;

            if let Some((next_index, next_leaf)) = leaf_queue.front() {
                let are_siblings = index ^ 1 == *next_index;
                if are_siblings {
                    let running_hash = C::hash_leaves(height - 1, &leaf, next_leaf);
                    node_queue.push_back((node_index, running_hash));
                    leaf_queue.pop_front();
                    continue;
                }
            }

            let sibling = siblings.pop_front().unwrap();
            let running_hash = if index % 2 == 0 {
                C::hash_leaves(height - 1, &leaf, &sibling)
            } else {
                C::hash_leaves(height - 1, &sibling, &leaf)
            };
            node_queue.push_back((node_index, running_hash));
        }
        assert!(siblings.is_empty());

        // handle internal nodes
        let mut nodes = VecDeque::from_iter(proof.nodes);
        while let Some((index, hash)) = node_queue.pop_front() {
            let depth = index.ilog2();

            if depth == 0 {
                assert!(node_queue.is_empty());
                // compare against the root
                return if *root == hash {
                    Ok(())
                } else {
                    Err(Error::InvalidProof)
                };
            }

            if let Some((next_index, next_hash)) = node_queue.front() {
                let are_siblings = index ^ 1 == *next_index;
                if are_siblings {
                    let running_hash = C::hash_nodes(depth - 1, &hash, next_hash);
                    node_queue.push_back((index >> 1, running_hash));
                    node_queue.pop_front();
                    continue;
                }
            }

            let sibling = nodes.pop_front().unwrap();
            let running_hash = if index % 2 == 0 {
                C::hash_nodes(depth - 1, &hash, &sibling)
            } else {
                C::hash_nodes(depth - 1, &sibling, &hash)
            };
            node_queue.push_back((index >> 1, running_hash));
        }

        Ok(())
    }

    fn security_level_bits() -> u32 {
        C::security_level_bits()
    }
}

/// Merkle tree that supports proving/verifying rows of a matrix
///
/// Inspired by plonky3's MMCS
/// <https://github.com/Plonky3/Plonky3/blob/main/commit/src/mmcs.rs>
pub trait MatrixMerkleTree<T>: MerkleTree + Sized {
    fn from_matrix(m: &Matrix<T>) -> Self;

    fn prove_rows(&self, row_ids: &[usize]) -> Result<Self::Proof, Error> {
        self.prove(row_ids)
    }

    fn verify_rows(
        root: &Self::Root,
        row_ids: &[usize],
        rows: &[impl AsRef<[T]>],
        proof: Self::Proof,
    ) -> Result<(), Error>;
}

pub struct MatrixMerkleTreeImpl<H: HashFn> {
    merkle_tree: MerkleTreeImpl<HashedLeafConfig<H>>,
}

impl<H: HashFn> Clone for MatrixMerkleTreeImpl<H> {
    fn clone(&self) -> Self {
        Self {
            merkle_tree: self.merkle_tree.clone(),
        }
    }
}

impl<H: HashFn> MatrixMerkleTreeImpl<H> {
    fn new(leaves: Vec<H::Digest>) -> Result<Self, Error> {
        assert!(leaves.len().is_power_of_two());
        Ok(Self {
            merkle_tree: MerkleTreeImpl::new(leaves)?,
        })
    }
}

impl<H: HashFn> MerkleTree for MatrixMerkleTreeImpl<H> {
    type Proof = MerkleView<H::Digest, H::Digest>;
    type Root = H::Digest;

    fn root(&self) -> Self::Root {
        self.merkle_tree.root()
    }

    fn prove(&self, indices: &[usize]) -> Result<Self::Proof, Error> {
        self.merkle_tree.prove(indices)
    }

    fn verify(root: &Self::Root, proof: Self::Proof, indices: &[usize]) -> Result<(), Error> {
        MerkleTreeImpl::<HashedLeafConfig<H>>::verify(root, proof, indices)
    }

    fn security_level_bits() -> u32 {
        H::COLLISION_RESISTANCE
    }
}

impl<F: Field, H: ElementHashFn<F> + Send + Sync + 'static> MatrixMerkleTree<F>
    for MatrixMerkleTreeImpl<H>
{
    fn from_matrix(m: &Matrix<F>) -> Self {
        Self::new(hash_rows::<F, H>(m)).unwrap()
    }

    fn verify_rows(
        root: &Self::Root,
        row_ids: &[usize],
        rows: &[impl AsRef<[F]>],
        proof: Self::Proof,
    ) -> Result<(), Error> {
        // remove duplicates and sort
        let mut instances = zip(row_ids, rows).collect::<Vec<_>>();
        instances.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
        instances.dedup_by(|(a, _), (b, _)| a == b);

        let (indices, rows): (Vec<_>, Vec<_>) = instances.into_iter().unzip();
        let initial_leaves = rows
            .iter()
            .map(|r| H::hash_elements(r.as_ref().iter().copied()))
            .collect::<Vec<_>>();
        if proof.initial_leaves == initial_leaves {
            Self::verify(root, proof, &indices)
        } else {
            Err(Error::InvalidProof)
        }
    }
}

pub struct HashedLeafConfig<H: HashFn>(PhantomData<H>);

impl<H: HashFn> Clone for HashedLeafConfig<H> {
    fn clone(&self) -> Self {
        Self(PhantomData)
    }
}

impl<H: HashFn> MerkleTreeConfig for HashedLeafConfig<H> {
    type Digest = H::Digest;
    type Leaf = H::Digest;

    fn hash_leaves(_: u32, l0: &H::Digest, l1: &H::Digest) -> H::Digest {
        H::merge(l0, l1)
    }

    fn hash_nodes(_: u32, n0: &Self::Digest, n1: &Self::Digest) -> Self::Digest {
        H::merge(n0, n1)
    }

    fn security_level_bits() -> u32 {
        H::COLLISION_RESISTANCE
    }
}

pub fn hash_rows<F: Field, H: ElementHashFn<F>>(matrix: &Matrix<F>) -> Vec<H::Digest> {
    let num_rows = matrix.num_rows();
    let mut row_hashes = vec![H::Digest::default(); num_rows];

    #[cfg(not(feature = "parallel"))]
    let chunk_size = row_hashes.len();
    #[cfg(feature = "parallel")]
    let chunk_size = core::cmp::max(
        row_hashes.len() / rayon::current_num_threads().next_power_of_two(),
        128,
    );

    ark_std::cfg_chunks_mut!(row_hashes, chunk_size)
        .enumerate()
        .for_each(|(chunk_offset, chunk)| {
            let offset = chunk_size * chunk_offset;
            let mut row_buffer = vec![F::zero(); matrix.num_cols()];
            for (i, row_hash) in chunk.iter_mut().enumerate() {
                matrix.read_row(offset + i, &mut row_buffer);
                *row_hash = H::hash_elements(row_buffer.iter().copied());
            }
        });

    row_hashes
}

#[cfg(feature = "parallel")]
pub fn build_merkle_nodes<C: MerkleTreeConfig>(leaves: &[C::Leaf]) -> Vec<C::Digest> {
    let n = leaves.len();
    let num_subtrees = core::cmp::min(rayon::current_num_threads().next_power_of_two(), n / 2);
    let mut nodes = vec![C::Digest::default(); n];

    // code adapted from winterfell
    // https://github.com/facebook/winterfell
    rayon::scope(|s| {
        for i in 0..num_subtrees {
            let nodes = unsafe { &mut *core::ptr::addr_of_mut!(nodes[..]) };
            s.spawn(move |_| {
                // generate first layer of nodes from leaf nodes
                let batch_size = n / num_subtrees;
                let leaf_offset = batch_size * i;
                let depth = (n / 2).ilog2();
                for j in (0..batch_size).step_by(2) {
                    let lhs = &leaves[leaf_offset + j];
                    let rhs = &leaves[leaf_offset + j + 1];
                    nodes[(n + leaf_offset + j) / 2] = C::hash_leaves(depth, lhs, rhs);
                }

                // generate remaining nodes
                let mut batch_size = n / num_subtrees / 4;
                let mut start_idx = n / 4 + batch_size * i;
                let mut depth = (n / 4).ilog2();
                while start_idx >= num_subtrees {
                    for k in (start_idx..(start_idx + batch_size)).rev() {
                        nodes[k] = C::hash_nodes(depth, &nodes[k * 2], &nodes[k * 2 + 1]);
                    }
                    start_idx /= 2;
                    batch_size /= 2;
                    depth -= 1;
                }
            });
        }
    });

    // finish the tip of the tree
    for i in (1..num_subtrees).rev() {
        let layer = i.ilog2();
        nodes[i] = C::hash_nodes(layer, &nodes[i * 2], &nodes[i * 2 + 1]);
    }

    nodes
}

#[cfg(not(feature = "parallel"))]
pub fn build_merkle_nodes<C: MerkleTreeConfig>(leaves: &[C::Leaf]) -> Vec<C::Digest> {
    let n = leaves.len();
    assert!(n.is_power_of_two());
    let mut nodes = vec![C::Digest::default(); n];

    // generate first layer of nodes from leaf nodes
    let depth = (n / 2).ilog2();
    for i in 0..n / 2 {
        nodes[n / 2 + i] = C::hash_leaves(depth, &leaves[i * 2], &leaves[i * 2 + 1]);
    }

    // generate remaining nodes
    let num_remaining_layers = (n / 2).ilog2();
    for depth in (0..num_remaining_layers).rev() {
        let size = 1 << depth;
        let offset = size;
        for i in offset..offset + size {
            nodes[i] = C::hash_nodes(depth, &nodes[i * 2], &nodes[i * 2 + 1]);
        }
    }

    nodes
}

#[cfg(test)]
mod tests {
    use super::Error;
    use super::MatrixMerkleTree;
    use super::MatrixMerkleTreeImpl;
    use super::MerkleTree;
    use super::MerkleTreeConfig;
    use super::MerkleTreeImpl;
    use crate::hash::HashFn;
    use crate::hash::Sha256HashFn;
    use crate::utils::GpuAllocator;
    use crate::utils::SerdeOutput;
    use crate::Matrix;
    use ark_ff::MontFp as Fp;
    use digest::Digest;
    use ministark_gpu::fields::p3618502788666131213697322783095070105623107215331596699973092056135872020481::ark::Fp;
    use sha2::Sha256;

    #[test]
    fn verify() -> Result<(), Error> {
        let leaves = vec![1u32, 2, 3, 4, 5, 6, 7, 8];
        let tree = MerkleTreeImpl::<UnhashedLeafConfig>::new(leaves)?;
        let commitment = tree.root();
        let i = 3;

        let proof = tree.prove(&[i])?;

        MerkleTreeImpl::<UnhashedLeafConfig>::verify(&commitment, proof, &[i])
    }

    #[test]
    fn prove_all_leaves() -> Result<(), Error> {
        let column: &[Fp] = &[Fp!("1"), Fp!("2"), Fp!("3"), Fp!("4")];
        let matrix = Matrix::new(vec![column.to_vec_in(GpuAllocator)]);
        let tree = MatrixMerkleTreeImpl::<Sha256HashFn>::from_matrix(&matrix);
        let commitment = tree.root();
        let row_ids = [0, 1, 2, 3];
        let rows = row_ids.map(|i| [column[i]]);

        let proof = MatrixMerkleTree::<Fp>::prove_rows(&tree, &row_ids)?;

        MatrixMerkleTreeImpl::<Sha256HashFn>::verify_rows(&commitment, &row_ids, &rows, proof)
    }

    #[test]
    fn verify_hashed_leaves() -> Result<(), Error> {
        let leaves = [1u32, 2, 3, 4, 5, 6, 7, 8];
        let hashed_leaves = leaves
            .iter()
            .map(|&v| Sha256::digest(v.to_be_bytes()))
            .map(SerdeOutput::new)
            .collect();
        let tree = MerkleTreeImpl::<HashedLeafConfig>::new(hashed_leaves)?;
        let commitment = tree.root();
        let i = 3;

        let proof = tree.prove(&[i])?;

        MerkleTreeImpl::<HashedLeafConfig>::verify(&commitment, proof, &[i])
    }

    #[test]
    fn verify_large_tree() -> Result<(), Error> {
        let leaves = (0..1 << 10).collect::<Vec<u32>>();
        let tree = MerkleTreeImpl::<UnhashedLeafConfig>::new(leaves)?;
        let commitment = tree.root();
        let i = 378;

        let proof = tree.prove(&[i])?;

        MerkleTreeImpl::<UnhashedLeafConfig>::verify(&commitment, proof, &[i])
    }

    struct HashedLeafConfig;

    impl MerkleTreeConfig for HashedLeafConfig {
        type Digest = SerdeOutput<Sha256>;
        type Leaf = SerdeOutput<Sha256>;

        fn hash_leaves(
            _: u32,
            l0: &SerdeOutput<Sha256>,
            l1: &SerdeOutput<Sha256>,
        ) -> SerdeOutput<Sha256> {
            Sha256HashFn::merge(l0, l1)
        }

        fn hash_nodes(_: u32, n0: &Self::Digest, n1: &Self::Digest) -> Self::Digest {
            Sha256HashFn::merge(n0, n1)
        }

        fn security_level_bits() -> u32 {
            Sha256HashFn::COLLISION_RESISTANCE
        }
    }

    struct UnhashedLeafConfig;

    impl MerkleTreeConfig for UnhashedLeafConfig {
        type Digest = SerdeOutput<Sha256>;
        type Leaf = u32;

        fn hash_leaves(_: u32, l0: &u32, l1: &u32) -> SerdeOutput<Sha256> {
            let l0_bytes = l0.to_be_bytes();
            let l1_bytes = l1.to_be_bytes();
            Sha256HashFn::hash_chunks([&l0_bytes[..], &l1_bytes[..]])
        }

        fn hash_nodes(_: u32, n0: &Self::Digest, n1: &Self::Digest) -> Self::Digest {
            Sha256HashFn::merge(n0, n1)
        }

        fn security_level_bits() -> u32 {
            Sha256HashFn::COLLISION_RESISTANCE
        }
    }
}
