//! Use arkwork_rs or re make this. Just used for personal education.
use anyhow::Result;
use digest::Digest;
use digest::Output;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;
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

/// Merkle tree implemented as a full (power of 2) arity tree.
///
/// ```text
///       #        <- root node
///     /   \
///   #       #    <- nodes
///  / \     / \
/// #   #   #   #  <- leaves
/// |   |   |   |
/// +   +   +   +  <- values
/// ```
pub struct MerkleTree<D: Digest> {
    nodes: Vec<Output<D>>,
    leafs: Vec<Output<D>>,
}

impl<D: Digest> MerkleTree<D> {
    pub fn new(leafs: Vec<Output<D>>) -> Result<Self, MerkleTreeError> {
        let n = leafs.len();
        if n < 2 {
            return Err(MerkleTreeError::TooFewLeaves(2, n));
        } else if !n.is_power_of_two() {
            return Err(MerkleTreeError::NumberOfLeavesNotPowerOfTwo(n));
        }

        let nodes = build_merkle_nodes::<D>(&leafs);
        Ok(MerkleTree { nodes, leafs })
    }

    pub fn root(&self) -> &Output<D> {
        &self.nodes[1]
    }
}

fn build_merkle_nodes<D: Digest>(leafs: &[Output<D>]) -> Vec<Output<D>> {
    let n = leafs.len();
    let mut nodes = vec![Output::<D>::default(); n];
    // generate first row of nodes (parents of leaves)
    for i in 0..n / 2 {
        let mut hasher = D::new();
        hasher.update(&leafs[2 * i]);
        hasher.update(&leafs[2 * i + 1]);
        nodes[n / 2 + i] = hasher.finalize();
    }
    // generate remainding nodes
    for i in (1..n / 2).rev() {
        let mut hasher = D::new();
        hasher.update(&nodes[2 * i]);
        hasher.update(&leafs[2 * i + 1]);
        nodes[i] = hasher.finalize();
    }
    nodes
}

// #[cfg(test)]
// mod tests {
//     use super::MerkleTree;

//     #[test]
//     fn test_verify() {
//         let data = vec![584395, 344, 2, 543, 5435, 343, 76, 88];
//         let tree = MerkleTree::new(&data);
//         let open_index = 4;
//         let (element, path) = tree.open(open_index);

//         assert!(MerkleTree::verify(tree.root(), open_index, &path, &element))
//     }
// }
