use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_serialize::CanonicalSerialize;
use rand::Rng;
use serde::Deserialize;
use serde::Serialize;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;

// impl<T: CanonicalSerialize> CanonicalSerialize for UncompressedUnchecked<T> {
// 	// always invoke `T::serialize_uncompressed()`
// }
// impl<T: CanonicalDeserialize> CanonicalDeserialize for
// UncompressedUnchecked<T> { 	// always invoke
// `T::deserialize_uncompressed_unchecked()` }
// impl<T: CanonicalSerialize> serde::Serialize for UncompressedUnchecked<T> {
// 	// always invoke `CanonicalSerialize::serialize_uncompressed()`
// }
// impl<T: CanonicalDeserialize> serde::DeserializeOwned for
// UncompressedUnchecked<T> { 	// always invoke
// `CanonicalDeserialize::deserialize_uncompressed_unchecked()` }

// #[derive(Serialize, Deserialize, Clone)]
pub enum ProofObject<F: Field> {
    MerkleRoot(u64),
    Codeword(Vec<F>),
    MerklePath(Vec<u64>),
    YValue(F),
    YValues((F, F, F)),

    // new vals
    Terminal(F),
    LeafItems(Vec<F>),
    LeafItem(F),
    MerklePathWithSalt((u64, Vec<u64>)),
    FriLeafs((F, F, F)),
}

pub trait ProofStream<F: Field> {
    fn push(&mut self, object: ProofObject<F>);
    fn pull(&mut self) -> ProofObject<F>;
    fn serialize(&self) -> Vec<u8>;
    fn deserialize(&self, bytes: &[u8]) -> Self;
    fn prover_fiat_shamir(&self) -> u64;
    fn verifier_fiat_shamir(&self) -> u64;
}

// #[derive(Serialize, Deserialize)]
pub struct StandardProofStream<F: Field> {
    pub objects: Vec<ProofObject<F>>,
    pub read_index: usize,
}

impl<F: Field> StandardProofStream<F> {
    pub fn new() -> Self {
        Default::default()
    }
}

impl<F: Field> Default for StandardProofStream<F> {
    fn default() -> Self {
        StandardProofStream {
            objects: vec![],
            read_index: 0,
        }
    }
}

impl<F: Field> ProofStream<F> for StandardProofStream<F> {
    fn push(&mut self, object: ProofObject<F>) {
        // self.objects.push(object);
        // todo!()
    }

    fn pull(&mut self) -> ProofObject<F> {
        // assert!(
        //     self.read_index < self.objects.len(),
        //     "ProofStream: cannot pull object; queue empty."
        // );
        // let object = self.objects[self.read_index].clone();
        // self.read_index += 1;
        // object
        todo!()
    }

    fn serialize(&self) -> Vec<u8> {
        // serde_json::to_vec(&self.objects.iter().map(|object| match object {

        // })).unwrap()
        todo!()
    }

    fn deserialize(&self, bytes: &[u8]) -> StandardProofStream<F> {
        // let objects = serde_json::from_slice(bytes).unwrap();
        // StandardProofStream {
        //     objects,
        //     read_index: 0,
        // }
        todo!()
    }

    fn prover_fiat_shamir(&self) -> u64 {
        // let serialized_objects = serde_json::to_vec(self).unwrap();
        // let mut hash = DefaultHasher::new();
        // serialized_objects.hash(&mut hash);
        // hash.finish()
        // todo!()
        let mut rng = rand::thread_rng();
        rng.gen()
    }

    fn verifier_fiat_shamir(&self) -> u64 {
        // let serialized_objects =
        //     serde_json::to_vec(&self.objects[..self.read_index].to_vec()).unwrap();
        // let mut hash = DefaultHasher::new();
        // serialized_objects.hash(&mut hash);
        // hash.finish()
        todo!()
    }
}
