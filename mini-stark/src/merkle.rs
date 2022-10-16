//! Use arkwork_rs or re make this. Just used for personal education.
use anyhow::Result;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::marker::PhantomData;
use thiserror::Error;

/// MerkleTree tree error
#[derive(Error, Debug)]
pub enum MerkleTreeError {
    #[error("tree must contain at least `{0}` leaves, but `{1}` were provided")]
    TooFewLeaves(usize, usize),
    #[error("number of leaves must be a power of two, but `{0}` were provided")]
    NumberOfLeavesNotPowerOfTwo(usize),
    #[error("leaf index `{0}` cannot exceed the number of leaves (`{1}`)")]
    LeafIndexOutOfBounds(usize, usize),
    #[error("proof is invalid")]
    InvalidProof,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct MerkleProof(Vec<u8>);

impl MerkleProof {
    pub fn new<D: Digest>(merkle_path: Vec<Output<D>>) -> Self {
        MerkleProof(merkle_path.into_iter().flatten().collect())
    }

    pub fn parse<D: Digest>(&self) -> Vec<Output<D>> {
        // TODO: would be great if this whole thing could be better.
        let chunk_size = <D as digest::OutputSizeUser>::output_size();
        let chunks = self.0.chunks(chunk_size);
        chunks
            .map(|chunk| Output::<D>::from_slice(chunk).clone())
            .collect()
    }
}

/// Merkle tree implemented as a full power-of-two arity tree.
///
/// ```text
///       #        <- root node
///     /   \
///   #       #    <- nodes
///  / \     / \
/// #   #   #   #  <- leaf nodes
/// |   |   |   |
/// +   +   +   +  <- values
/// ```
pub struct MerkleTree<D: Digest> {
    nodes: Vec<Output<D>>,
    leaf_nodes: Vec<Output<D>>,
}

impl<D: Digest> MerkleTree<D> {
    pub fn new(leaf_nodes: Vec<Output<D>>) -> Result<Self, MerkleTreeError> {
        let n = leaf_nodes.len();
        if n < 2 {
            return Err(MerkleTreeError::TooFewLeaves(2, n));
        } else if !n.is_power_of_two() {
            return Err(MerkleTreeError::NumberOfLeavesNotPowerOfTwo(n));
        }

        let nodes = build_merkle_nodes::<D>(&leaf_nodes);
        Ok(MerkleTree { nodes, leaf_nodes })
    }

    pub fn root(&self) -> &Output<D> {
        &self.nodes[1]
    }

    pub fn prove(&self, index: usize) -> Result<MerkleProof, MerkleTreeError> {
        if index >= self.leaf_nodes.len() {
            return Err(MerkleTreeError::LeafIndexOutOfBounds(
                self.leaf_nodes.len(),
                index,
            ));
        }

        // TODO: batch proofs
        // TODO: could omit leaf_nodes[index]
        let mut path = vec![
            self.leaf_nodes[index].clone(),
            self.leaf_nodes[index ^ 1].clone(),
        ];

        let mut index = (index + self.nodes.len()) >> 1;
        while index > 1 {
            path.push(self.nodes[index ^ 1].clone());
            index >>= 1;
        }

        Ok(MerkleProof::new::<D>(path))
    }
}

fn build_merkle_nodes<D: Digest>(leaf_nodes: &[Output<D>]) -> Vec<Output<D>> {
    let n = leaf_nodes.len();
    let mut nodes = vec![Output::<D>::default(); n];
    // generate first row of nodes (parents of leaves)
    ark_std::cfg_iter_mut!(nodes[n / 2..])
        .enumerate()
        .for_each(|(i, node)| {
            let mut hasher = D::new();
            hasher.update(&leaf_nodes[2 * i]);
            hasher.update(&leaf_nodes[2 * i + 1]);
            *node = hasher.finalize();
        });
    // generate remainding nodes
    for i in (1..n / 2).rev() {
        let mut hasher = D::new();
        hasher.update(&nodes[2 * i]);
        hasher.update(&leaf_nodes[2 * i + 1]);
        nodes[i] = hasher.finalize();
    }
    nodes
}
