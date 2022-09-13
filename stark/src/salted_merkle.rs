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
        let depth = n.ilog2();
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
            leaf.hash(&mut hasher);
            nodes[n + i] = hasher.finish();
        }
        for i in (0..n).rev() {
            let mut hasher = DefaultHasher::new();
            let state = ((nodes[i * 2] as u128) << 64) + nodes[i + 1] as u128;
            state.hash(&mut hasher);
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
        let mut i = i + self.leafs.len();
        assert!(self.leafs.len().is_power_of_two());
        while i > 0 {
            authentication_path.push(self.nodes[i]);
            i >>= 1;
        }
        (element, salt, authentication_path)
    }
}
