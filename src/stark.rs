use crate::air::AirConfig;
use crate::challenges::Challenges;
use crate::composer::DeepCompositionCoeffs;
use crate::debug::default_validate_constraints;
use crate::hints::Hints;
use crate::merkle::MatrixMerkleTree;
use crate::merkle::MerkleTree;
use crate::prover::default_prove;
use crate::prover::ProvingError;
use crate::random::draw_multiple;
use crate::random::PublicCoin;
use crate::verifier::default_verify;
use crate::verifier::VerificationError;
use crate::Air;
use crate::Matrix;
use crate::Proof;
use crate::ProofOptions;
use crate::StarkExtensionOf;
use crate::Trace;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use ministark_gpu::GpuFftField;

pub trait Stark: Sized + Send + Sync {
    type Fp: GpuFftField + FftField;
    type Fq: StarkExtensionOf<Self::Fp>;
    type AirConfig: AirConfig<Fp = Self::Fp, Fq = Self::Fq>;
    type Digest: Digest;
    type PublicCoin: PublicCoin<Digest = Self::Digest, Field = Self::Fq>;
    type MerkleTree: MerkleTree<Root = Output<Self::Digest>>
        + MatrixMerkleTree<Self::Fp>
        + MatrixMerkleTree<Self::Fq>;
    type Trace: Trace<Fp = Self::Fp, Fq = Self::Fq>;
    type Witness;

    fn get_public_inputs(&self) -> <Self::AirConfig as AirConfig>::PublicInputs;

    fn gen_public_coin(&self, air: &Air<Self::AirConfig>) -> Self::PublicCoin {
        let mut seed = Vec::new();
        air.public_inputs().serialize_compressed(&mut seed).unwrap();
        air.trace_len().serialize_compressed(&mut seed).unwrap();
        air.options().serialize_compressed(&mut seed).unwrap();
        PublicCoin::new(Self::Digest::digest(&seed))
    }

    fn gen_deep_coeffs(
        &self,
        public_coin: &mut Self::PublicCoin,
        air: &Air<Self::AirConfig>,
    ) -> DeepCompositionCoeffs<Self::Fq> {
        let num_execution_trace = air.trace_arguments().len();
        let num_composition_trace = air.ce_blowup_factor();
        DeepCompositionCoeffs {
            execution_trace: draw_multiple(public_coin, num_execution_trace),
            composition_trace: draw_multiple(public_coin, num_composition_trace),
            degree: (public_coin.draw(), public_coin.draw()),
        }
    }

    fn generate_trace(&self, witness: Self::Witness) -> Self::Trace;

    async fn prove(
        &self,
        options: ProofOptions,
        witness: Self::Witness,
    ) -> Result<Proof<Self::Fp, Self::Fq, Self::Digest, Self::MerkleTree>, ProvingError> {
        default_prove(self, options, witness)
    }

    /// Check the AIR constraints are valid
    fn validate_constraints(
        &self,
        challenges: &Challenges<Self::Fq>,
        hints: &Hints<Self::Fq>,
        base_trace: &Matrix<Self::Fp>,
        extension_trace: Option<&Matrix<Self::Fq>>,
    ) {
        // TODO: this is unfinished
        default_validate_constraints(self, challenges, hints, base_trace, extension_trace);
    }

    #[allow(clippy::too_many_lines)]
    fn verify(
        &self,
        proof: Proof<Self::Fp, Self::Fq, Self::Digest, Self::MerkleTree>,
    ) -> Result<(), VerificationError> {
        default_verify(self, proof)
    }
}
