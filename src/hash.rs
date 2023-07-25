use crate::utils::SerdeOutput;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use core::fmt::Debug;
use digest::Digest as _;
use num_traits::ops::bytes;
use sha2::Sha256;

/// NOTE: ever so slightly modified Hasher trait from Winterfell
pub trait HashFn: Send + Sync + 'static {
    /// Specifies a digest type returned by this hasher.
    type Digest: Digest;

    /// Collision resistance of the hash function measured in bits.
    const COLLISION_RESISTANCE: u32;

    /// Returns a hash of the provided sequence of bytes.
    fn hash(bytes: impl IntoIterator<Item = u8>) -> Self::Digest;

    /// Returns a hash of two digests. This method is intended for use in
    /// construction of Merkle trees.
    fn merge(v0: &Self::Digest, v1: &Self::Digest) -> Self::Digest;

    /// Returns hash(`seed` || `value`). This method is intended for use in PRNG
    /// and PoW contexts.
    fn merge_with_int(seed: &Self::Digest, value: u64) -> Self::Digest;
}

/// NOTE: ever so slightly modified ElementHasher trait from Winterfell
/// Defines a cryptographic hash function for hashing field elements.
///
/// This trait defines a hash procedure for a sequence of field elements. The
/// elements can be either in the base field specified for this hasher, or in an
/// extension of the base field.
pub trait ElementHashFn<F: Field>: HashFn {
    /// Returns a hash of the provided field elements.
    fn hash_elements(elements: impl IntoIterator<Item = F>) -> Self::Digest;
}

/// NOTE: Digest trait from winterfell
/// Defines output type for a cryptographic hash function.
pub trait Digest:
    Debug + Default + Clone + Eq + PartialEq + Send + Sync + CanonicalSerialize + CanonicalDeserialize
{
    /// Returns this digest serialized into an array of bytes.
    ///
    /// Ideally, the length of the returned array should be defined by an
    /// associated constant, but using associated constants in const
    /// generics is not supported by Rust yet. Thus, we put an upper limit
    /// on the possible digest size. For digests which are smaller than 32
    /// bytes, the unused bytes should be set to 0.
    fn as_bytes(&self) -> [u8; 32];
}

pub struct Sha256HashFn;

impl HashFn for Sha256HashFn {
    type Digest = SerdeOutput<Sha256>;

    const COLLISION_RESISTANCE: u32 = 128;

    fn hash(bytes: impl IntoIterator<Item = u8>) -> SerdeOutput<Sha256> {
        let mut hasher = Sha256::new();
        bytes.into_iter().for_each(|b| hasher.update([b]));
        SerdeOutput::new(hasher.finalize())
    }

    fn merge(v0: &SerdeOutput<Sha256>, v1: &SerdeOutput<Sha256>) -> SerdeOutput<Sha256> {
        let mut hasher = Sha256::new();
        hasher.update(**v0);
        hasher.update(**v1);
        SerdeOutput::new(hasher.finalize())
    }

    fn merge_with_int(seed: &SerdeOutput<Sha256>, value: u64) -> SerdeOutput<Sha256> {
        let mut hasher = Sha256::new();
        hasher.update(**seed);
        hasher.update(value.to_be_bytes());
        SerdeOutput::new(hasher.finalize())
    }
}

impl<F: Field> ElementHashFn<F> for Sha256HashFn {
    fn hash_elements(elements: impl IntoIterator<Item = F>) -> Self::Digest {
        let mut byte_buffer = Vec::new();
        for element in elements {
            element.serialize_uncompressed(&mut byte_buffer).unwrap();
        }
        Self::hash(byte_buffer)
    }
}
