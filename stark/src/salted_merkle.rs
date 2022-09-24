//! Use arkwork_rs or re make this. Just used for personal education.
use rand::Rng;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;

// TODO: Using default hasher to make prototyping easier. Fix later
pub struct SaltedMerkle<T> {
    pub leafs: Vec<(T, u64)>,
    nodes: Vec<u64>,
}

impl<T: Hash + Clone> SaltedMerkle<T> {
    pub fn new(values: &[T]) -> Self {
        let n = values.len();
        assert!(n.is_power_of_two());
        let mut rng = rand::thread_rng();
        // append salt to leafs
        let leafs = values
            .iter()
            .map(|v| (v.clone(), rng.gen()))
            .collect::<Vec<(T, u64)>>();
        // make room for nodes
        let mut nodes = vec![0; n * 2];
        for (i, leaf) in leafs.iter().enumerate() {
            let mut hasher = DefaultHasher::new();
            leaf.0.hash(&mut hasher);
            leaf.1.hash(&mut hasher);
            nodes[n + i] = hasher.finish();
        }
        for i in (0..n).rev() {
            let mut hasher = DefaultHasher::new();
            // let state = ((nodes[i * 2] as u128) << 64) + nodes[i + 1] as u128;
            // state
            nodes[i * 2].hash(&mut hasher);
            nodes[i * 2 + 1].hash(&mut hasher);
            nodes[i] = hasher.finish();
        }
        SaltedMerkle { leafs, nodes }
    }

    pub fn root(&self) -> u64 {
        self.nodes[1]
    }

    // Returns element, salt and authentication path
    pub fn open(&self, i: usize) -> (T, u64, Vec<u64>) {
        let (element, salt) = self.leafs[i].clone();
        let mut authentication_path = vec![];
        let mut i = i | self.leafs.len();
        assert!(self.leafs.len().is_power_of_two());
        while i > 1 {
            authentication_path.push(self.nodes[i ^ 1]);
            i >>= 1;
        }
        (element, salt, authentication_path)
    }

    pub fn verify(root: u64, i: usize, salt: u64, path: &[u64], element: &T) -> bool {
        let mut i = i;
        let mut running_hath = DefaultHasher::new();
        element.hash(&mut running_hath);
        salt.hash(&mut running_hath);
        for &node in path {
            let hash = running_hath.finish();
            running_hath = DefaultHasher::new();
            if i % 2 == 0 {
                hash.hash(&mut running_hath);
                node.hash(&mut running_hath);
            } else {
                node.hash(&mut running_hath);
                hash.hash(&mut running_hath);
            }
            i >>= 1;
        }
        running_hath.finish() == root
    }
}

#[cfg(test)]
mod tests {
    use super::SaltedMerkle;

    #[test]
    fn test_verify() {
        let data = vec![584395, 344, 2, 543, 5435, 343, 76, 88];
        let tree = SaltedMerkle::new(&data);
        let open_index = 4;
        let (element, salt, path) = tree.open(open_index);

        assert!(SaltedMerkle::verify(
            tree.root(),
            open_index,
            salt,
            &path,
            &element
        ))
    }
}
