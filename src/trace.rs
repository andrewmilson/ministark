use crate::challenges::Challenges;
use crate::merkle::MatrixMerkleTree;
use crate::merkle::MerkleTree;
use crate::Matrix;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;

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

// #[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
// pub struct Queries<S: Stark> {
//     pub base_trace_values: Vec<S::Fp>,
//     pub extension_trace_values: Vec<S::Fq>,
//     pub composition_trace_values: Vec<S::Fq>,
//     pub base_trace_proofs: Vec<<S::MerkleTree as MerkleTree>::Proof>,
//     pub extension_trace_proofs: Vec<<S::MerkleTree as MerkleTree>::Proof>,
//     pub composition_trace_proofs: Vec<<S::MerkleTree as MerkleTree>::Proof>,
// }

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Queries<
    Fp: Field,
    Fq: Field,
    MerkleProof: CanonicalDeserialize + CanonicalSerialize + Clone + Send + Sync,
> {
    pub base_trace_values: Vec<Fp>,
    pub extension_trace_values: Vec<Fq>,
    pub composition_trace_values: Vec<Fq>,
    pub base_trace_proofs: Vec<MerkleProof>,
    pub extension_trace_proofs: Vec<MerkleProof>,
    pub composition_trace_proofs: Vec<MerkleProof>,
}

impl<
        Fp: Field,
        Fq: Field,
        MerkleProof: CanonicalDeserialize + CanonicalSerialize + Clone + Send + Sync,
    > Queries<Fp, Fq, MerkleProof>
{
    pub fn new<M: MerkleTree<Proof = MerkleProof> + MatrixMerkleTree<Fp> + MatrixMerkleTree<Fq>>(
        base_trace_lde: &Matrix<Fp>,
        extension_trace_lde: Option<&Matrix<Fq>>,
        composition_trace_lde: &Matrix<Fq>,
        base_tree: &M,
        extension_tree: Option<&M>,
        composition_tree: &M,
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
            let base_proof = MatrixMerkleTree::<Fp>::prove_row(base_tree, position).unwrap();
            base_trace_proofs.push(base_proof);

            if let Some(extension_trace_lde) = extension_trace_lde {
                // TODO: suport ark DomainCoeff on evaluate_at
                let extension_trace_row = extension_trace_lde.get_row(position).unwrap();
                extension_trace_values.extend(extension_trace_row);
                let extension_tree = extension_tree.unwrap();
                let extension_proof =
                    MatrixMerkleTree::<Fp>::prove_row(extension_tree, position).unwrap();
                extension_trace_proofs.push(extension_proof);
            }

            // composition trace
            let composition_trace_row = composition_trace_lde.get_row(position).unwrap();
            composition_trace_values.extend(composition_trace_row);
            let composition_proof =
                MatrixMerkleTree::<Fp>::prove_row(composition_tree, position).unwrap();
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
