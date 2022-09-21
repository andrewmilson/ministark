use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ark_serialize::Valid;
use rand::Rng;
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

#[derive(Debug, Clone)]
pub enum ProofObject<F: Field> {
    MerkleRoot(u64),
    Codeword(Vec<F>),
    MerklePath(Vec<u64>),
    YValue(F),
    YValues((F, F, F)),

    // new vals
    TraceLen(u64),
    Terminal(F),
    LeafItems(Vec<F>),
    LeafItem(F),
    MerklePathWithSalt((u64, Vec<u64>)),
    FriLeafs((F, F, F)),
}

impl<F: Field> CanonicalSerialize for ProofObject<F> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        match self {
            Self::MerkleRoot(v) => {
                0u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::Codeword(v) => {
                1u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::MerklePath(v) => {
                2u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::YValue(v) => {
                3u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::YValues(v) => {
                4u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::TraceLen(v) => {
                5u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::Terminal(v) => {
                6u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::LeafItems(v) => {
                7u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::LeafItem(v) => {
                8u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::MerklePathWithSalt(v) => {
                9u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
            Self::FriLeafs(v) => {
                10u8.serialize_compressed(&mut writer)?;
                v.serialize_compressed(&mut writer)?;
            }
        }
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        1 + match self {
            Self::MerkleRoot(v) => v.serialized_size(compress),
            Self::Codeword(v) => v.serialized_size(compress),
            Self::MerklePath(v) => v.serialized_size(compress),
            Self::YValue(v) => v.serialized_size(compress),
            Self::YValues(v) => v.serialized_size(compress),
            Self::TraceLen(v) => v.serialized_size(compress),
            Self::Terminal(v) => v.serialized_size(compress),
            Self::LeafItems(v) => v.serialized_size(compress),
            Self::LeafItem(v) => v.serialized_size(compress),
            Self::MerklePathWithSalt(v) => v.serialized_size(compress),
            Self::FriLeafs(v) => v.serialized_size(compress),
        }
    }
}

impl<F: Field> Valid for ProofObject<F> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

impl<F: Field> CanonicalDeserialize for ProofObject<F> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let mut buf = [0u8];
        reader.read_exact(&mut buf)?;
        let object_type = buf[0];
        Ok(match object_type {
            0 => Self::MerkleRoot(u64::deserialize_with_mode(reader, compress, validate)?),
            1 => Self::Codeword(Vec::<F>::deserialize_with_mode(reader, compress, validate)?),
            2 => Self::MerklePath(Vec::<u64>::deserialize_with_mode(
                reader, compress, validate,
            )?),
            3 => Self::YValue(F::deserialize_with_mode(reader, compress, validate)?),
            4 => Self::YValues(<(F, F, F)>::deserialize_with_mode(
                reader, compress, validate,
            )?),
            5 => Self::TraceLen(u64::deserialize_with_mode(reader, compress, validate)?),
            6 => Self::Terminal(F::deserialize_with_mode(reader, compress, validate)?),
            7 => Self::LeafItems(Vec::<F>::deserialize_with_mode(reader, compress, validate)?),
            8 => Self::LeafItem(F::deserialize_with_mode(reader, compress, validate)?),
            9 => Self::MerklePathWithSalt(<(u64, Vec<u64>)>::deserialize_with_mode(
                reader, compress, validate,
            )?),
            10 => Self::FriLeafs(<(F, F, F)>::deserialize_with_mode(
                reader, compress, validate,
            )?),
            _ => return Err(ark_serialize::SerializationError::UnexpectedFlags),
        })
    }
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
        self.objects.push(object);
    }

    fn pull(&mut self) -> ProofObject<F> {
        assert!(
            self.read_index < self.objects.len(),
            "ProofStream: cannot pull object; queue empty."
        );
        let object = self.objects[self.read_index].clone();
        self.read_index += 1;
        object
    }

    fn serialize(&self) -> Vec<u8> {
        // serde_json::to_vec(&self.objects.iter().map(|object| match object {

        // })).unwrap()
        // todo!()
        // Compress::Yes
        let mut buf = Vec::<u8>::new();
        self.objects.serialize_compressed(&mut buf).unwrap();
        buf
    }

    fn deserialize(&self, bytes: &[u8]) -> StandardProofStream<F> {
        // let mut buf = Vec::<u8>::new();
        // self.objects(&mut buf).unwrap();
        // buf
        let objects = Vec::<ProofObject<F>>::deserialize_compressed(bytes).unwrap();
        StandardProofStream {
            objects,
            read_index: 0,
        }
    }

    fn prover_fiat_shamir(&self) -> u64 {
        let serialized_objects = self.serialize();
        let mut hash = DefaultHasher::new();
        serialized_objects.hash(&mut hash);
        hash.finish()
    }

    fn verifier_fiat_shamir(&self) -> u64 {
        let mut serialized_objects = Vec::<u8>::new();
        self.objects[..self.read_index]
            .serialize_compressed(&mut serialized_objects)
            .unwrap();
        let mut hash = DefaultHasher::new();
        serialized_objects.hash(&mut hash);
        hash.finish()
    }
}
