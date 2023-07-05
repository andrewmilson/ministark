use crate::challenges::Challenges;
use crate::merkle;
use crate::merkle::MerkleTree;
use crate::Matrix;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use digest::Digest;

/// STARK witness
pub trait Witness {
    type Fp: FftField;
    type Fq: Field<BasePrimeField = Self::Fp>;

    /// Returns the number of rows in this execution trace.
    fn trace_len(&self) -> usize {
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

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Queries<Fp: Field, Fq: Field> {
    pub base_trace_values: Vec<Fp>,
    pub extension_trace_values: Vec<Fq>,
    pub composition_trace_values: Vec<Fq>,
    pub base_trace_proofs: Vec<merkle::Proof>,
    pub extension_trace_proofs: Vec<merkle::Proof>,
    pub composition_trace_proofs: Vec<merkle::Proof>,
}

impl<Fp: Field, Fq: Field> Queries<Fp, Fq> {
    pub fn new<D: Digest>(
        base_trace_lde: &Matrix<Fp>,
        extension_trace_lde: Option<&Matrix<Fq>>,
        composition_trace_lde: &Matrix<Fq>,
        base_commitment: &MerkleTree<D>,
        extension_commitment: Option<&MerkleTree<D>>,
        composition_commitment: &MerkleTree<D>,
        positions: &[usize],
    ) -> Self {
        let mut base_trace_values = Vec::new();
        let mut extension_trace_values = Vec::new();
        let mut composition_trace_values = Vec::new();
        let mut base_trace_proofs = Vec::new();
        let mut extension_trace_proofs = Vec::new();
        let mut composition_trace_proofs = Vec::new();
        for &position in positions {
            // execution trace
            let base_trace_row = base_trace_lde.get_row(position).unwrap();
            base_trace_values.extend(base_trace_row);
            let base_proof = base_commitment.prove(position).unwrap();
            base_trace_proofs.push(base_proof);

            if let Some(extension_trace_lde) = extension_trace_lde {
                // TODO: suport ark DomainCoeff on evaluate_at
                let extension_trace_row = extension_trace_lde.get_row(position).unwrap();
                extension_trace_values.extend(extension_trace_row);
                let extension_proof = extension_commitment
                    .as_ref()
                    .unwrap()
                    .prove(position)
                    .unwrap();
                extension_trace_proofs.push(extension_proof);
            }

            // composition trace
            let composition_trace_row = composition_trace_lde.get_row(position).unwrap();
            composition_trace_values.extend(composition_trace_row);
            let composition_proof = composition_commitment.prove(position).unwrap();
            composition_trace_proofs.push(composition_proof);
        }
        Self {
            base_trace_values,
            extension_trace_values,
            composition_trace_values,
            base_trace_proofs,
            extension_trace_proofs,
            composition_trace_proofs,
        }
    }
}
