use crate::challenges::Challenges;
use crate::merkle::MatrixMerkleTree;
use crate::merkle::MerkleTree;
use crate::stark::Stark;
use crate::Matrix;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ark_serialize::Valid;

/// STARK execution trace
#[allow(clippy::len_without_is_empty)]
pub trait Trace: Send + Sync {
    type Fp: FftField;
    type Fq: Field<BasePrimeField = Self::Fp>;

    /// Returns the number of rows in this execution trace.
    fn len(&self) -> usize {
        self.base_columns().num_rows()
    }

    /// Returns a reference to the base trace columns.
    fn base_columns(&self) -> &Matrix<Self::Fp>;

    /// Builds and returns the extension trace columns
    /// These columns require auxiliary random elements to be constructed.
    /// Returns None if there are no columns that require this.
    fn build_extension_columns(
        &self,
        _challenges: &Challenges<Self::Fq>,
    ) -> Option<Matrix<Self::Fq>> {
        None
    }
}

pub struct Queries<C: Stark> {
    pub base_trace_values: Vec<C::Fp>,
    pub extension_trace_values: Vec<C::Fq>,
    pub composition_trace_values: Vec<C::Fq>,
    pub base_trace_proof: <C::MerkleTree as MerkleTree>::Proof,
    pub extension_trace_proof: Option<<C::MerkleTree as MerkleTree>::Proof>,
    pub composition_trace_proof: <C::MerkleTree as MerkleTree>::Proof,
}

impl<C: Stark> CanonicalSerialize for Queries<C> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.base_trace_values
            .serialize_with_mode(&mut writer, compress)?;
        self.extension_trace_values
            .serialize_with_mode(&mut writer, compress)?;
        self.composition_trace_values
            .serialize_with_mode(&mut writer, compress)?;
        self.base_trace_proof
            .serialize_with_mode(&mut writer, compress)?;
        self.extension_trace_proof
            .serialize_with_mode(&mut writer, compress)?;
        self.composition_trace_proof
            .serialize_with_mode(&mut writer, compress)?;
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.base_trace_values.serialized_size(compress)
            + self.extension_trace_values.serialized_size(compress)
            + self.composition_trace_values.serialized_size(compress)
            + self.base_trace_proof.serialized_size(compress)
            + self.extension_trace_proof.serialized_size(compress)
            + self.composition_trace_proof.serialized_size(compress)
    }
}

impl<C: Stark> Valid for Queries<C> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

impl<C: Stark> CanonicalDeserialize for Queries<C> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        Ok(Self {
            base_trace_values: <_>::deserialize_with_mode(&mut reader, compress, validate)?,
            extension_trace_values: <_>::deserialize_with_mode(&mut reader, compress, validate)?,
            composition_trace_values: <_>::deserialize_with_mode(&mut reader, compress, validate)?,
            base_trace_proof: <_>::deserialize_with_mode(&mut reader, compress, validate)?,
            extension_trace_proof: <_>::deserialize_with_mode(&mut reader, compress, validate)?,
            composition_trace_proof: <_>::deserialize_with_mode(&mut reader, compress, validate)?,
        })
    }
}

impl<C: Stark> Clone for Queries<C> {
    fn clone(&self) -> Self {
        Self {
            base_trace_values: self.base_trace_values.clone(),
            extension_trace_values: self.extension_trace_values.clone(),
            composition_trace_values: self.composition_trace_values.clone(),
            base_trace_proof: self.base_trace_proof.clone(),
            extension_trace_proof: self.extension_trace_proof.clone(),
            composition_trace_proof: self.composition_trace_proof.clone(),
        }
    }
}

impl<C: Stark> Queries<C> {
    pub fn new(
        base_trace_lde: &Matrix<C::Fp>,
        extension_trace_lde: Option<&Matrix<C::Fq>>,
        composition_trace_lde: &Matrix<C::Fq>,
        base_tree: &C::MerkleTree,
        extension_tree: Option<&C::MerkleTree>,
        composition_tree: &C::MerkleTree,
        positions: &[usize],
    ) -> Self {
        let base_trace_proof = MatrixMerkleTree::<C::Fp>::prove_rows(base_tree, positions).unwrap();
        let extension_trace_proof = extension_tree.map(|extension_tree| {
            MatrixMerkleTree::<C::Fq>::prove_rows(extension_tree, positions).unwrap()
        });
        let composition_trace_proof =
            MatrixMerkleTree::<C::Fq>::prove_rows(composition_tree, positions).unwrap();

        let mut base_trace_values = Vec::new();
        let mut extension_trace_values = Vec::new();
        let mut composition_trace_values = Vec::new();
        for &position in positions {
            // execution trace
            let base_trace_row = base_trace_lde.get_row(position).unwrap();
            base_trace_values.extend(base_trace_row);

            if let Some(extension_trace_lde) = extension_trace_lde {
                // TODO: suport ark DomainCoeff on evaluate_at
                let extension_trace_row = extension_trace_lde.get_row(position).unwrap();
                extension_trace_values.extend(extension_trace_row);
            }

            // composition trace
            let composition_trace_row = composition_trace_lde.get_row(position).unwrap();
            composition_trace_values.extend(composition_trace_row);
        }
        Self {
            base_trace_values,
            extension_trace_values,
            composition_trace_values,
            base_trace_proof,
            extension_trace_proof,
            composition_trace_proof,
        }
    }
}
