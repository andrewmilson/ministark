use crate::polynomial::*;
use algebra::Felt;
use serde::Deserialize;
use serde::Serialize;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;

#[derive(Serialize, Deserialize, Clone)]
pub enum ProofObject<E> {
    Polynomial(Polynomial<E>),
    MerkleRoot(u64),
    Codeword(Vec<E>),
    MerklePath(Vec<u64>),
    YValue(E),
    YValues((E, E, E)),
}

pub trait ProofStream<E> {
    fn push(&mut self, object: ProofObject<E>);
    fn pull(&mut self) -> ProofObject<E>;
    fn serialize(&self) -> Vec<u8>;
    fn deserialize(&self, bytes: &[u8]) -> Self;
    fn prover_fiat_shamir(&self) -> u64;
    fn verifier_fiat_shamir(&self) -> u64;
}
pub struct StandardProofStream<E> {
    pub objects: Vec<ProofObject<E>>,
    pub read_index: usize,
}

impl<E: Felt> StandardProofStream<E> {
    pub fn new() -> Self {
        Default::default()
    }
}

impl<E: Felt> Default for StandardProofStream<E> {
    fn default() -> Self {
        StandardProofStream {
            objects: vec![],
            read_index: 0,
        }
    }
}

impl<E: Felt> ProofStream<E> for StandardProofStream<E> {
    fn push(&mut self, object: ProofObject<E>) {
        self.objects.push(object);
    }

    fn pull(&mut self) -> ProofObject<E> {
        assert!(
            self.read_index < self.objects.len(),
            "ProofStream: cannot pull object; queue empty."
        );
        let object = self.objects[self.read_index].clone();
        self.read_index += 1;
        object
    }

    fn serialize(&self) -> Vec<u8> {
        serde_json::to_vec(&self.objects).unwrap()
    }

    fn deserialize(&self, bytes: &[u8]) -> StandardProofStream<E> {
        let objects = serde_json::from_slice(bytes).unwrap();
        StandardProofStream {
            objects,
            read_index: 0,
        }
    }

    fn prover_fiat_shamir(&self) -> u64 {
        let mut hash = DefaultHasher::new();
        self.serialize().hash(&mut hash);
        hash.finish()
    }

    fn verifier_fiat_shamir(&self) -> u64 {
        let serialized_objects =
            serde_json::to_vec(&self.objects[..self.read_index].to_vec()).unwrap();
        let mut hash = DefaultHasher::new();
        serialized_objects.hash(&mut hash);
        hash.finish()
    }
}
