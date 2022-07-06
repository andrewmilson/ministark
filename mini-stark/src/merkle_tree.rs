use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

pub struct MerkleTree;

impl MerkleTree {
    pub fn commit<T: Hash + Clone>(leafs: &[T]) -> u64 {
        assert!(leafs.len().is_power_of_two(), "length must be power of two");

        let mut hash = DefaultHasher::new();

        if leafs.len() == 1 {
            leafs[0].hash(&mut hash);
        } else {
            Self::commit(&leafs[..(leafs.len() / 2)].to_vec()).hash(&mut hash);
            Self::commit(&leafs[(leafs.len() / 2)..].to_vec()).hash(&mut hash);
        }

        hash.finish()
    }

    pub fn open<T: Hash + Clone>(index: usize, leafs: &[T]) -> Vec<u64> {
        assert!(leafs.len().is_power_of_two(), "length must be power of two");
        assert!(index < leafs.len(), "cannot open invalid index");

        if leafs.len() == 2 {
            vec![Self::commit(&[leafs[1 - index].clone()])]
        } else if index < leafs.len() / 2 {
            let mut path = Self::open(index, &leafs[..(leafs.len() / 2)].to_vec());
            path.push(Self::commit(&leafs[(leafs.len() / 2)..].to_vec()));
            path
        } else {
            let mut path = Self::open(
                index - leafs.len() / 2,
                &leafs[(leafs.len() / 2)..].to_vec(),
            );
            path.push(Self::commit(&leafs[..(leafs.len() / 2)].to_vec()));
            path
        }
    }

    fn verify_path(root: u64, index: usize, path: &[u64], leaf_hash: u64) -> bool {
        assert!(index < (1 << path.len()), "cannot verify invalid index");

        if path.len() == 1 {
            let mut hasher = DefaultHasher::new();

            if index == 0 {
                leaf_hash.hash(&mut hasher);
                path[0].hash(&mut hasher);
            } else {
                path[0].hash(&mut hasher);
                leaf_hash.hash(&mut hasher);
            }

            hasher.finish() == root
        } else {
            let mut hasher = DefaultHasher::new();

            if index % 2 == 0 {
                leaf_hash.hash(&mut hasher);
                path[0].hash(&mut hasher);
            } else {
                path[0].hash(&mut hasher);
                leaf_hash.hash(&mut hasher);
            }

            Self::verify_path(root, index / 2, &path[1..], hasher.finish())
        }
    }

    pub fn verify<T: Hash + Clone>(root: u64, index: usize, path: &[u64], leaf: T) -> bool {
        let mut hasher = DefaultHasher::new();
        leaf.hash(&mut hasher);
        let leaf_hash = hasher.finish();

        Self::verify_path(root, index, path, leaf_hash)
    }
}
