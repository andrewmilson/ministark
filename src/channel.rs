use crate::challenges::Challenges;
use crate::fri;
use crate::fri::FriProof;
use crate::hints::Hints;
use crate::random::PublicCoin;
use crate::stark::Stark;
use crate::trace::Queries;
use crate::Air;
use crate::Proof;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use std::collections::BTreeSet;

pub struct ProverChannel<'a, S: Stark> {
    air: &'a Air<S::AirConfig>,
    pub public_coin: S::PublicCoin,
    base_trace_commitment: S::Digest,
    extension_trace_commitment: Option<S::Digest>,
    composition_trace_commitment: S::Digest,
    fri_layer_commitments: Vec<S::Digest>,
    fri_remainder_coeffs: Vec<S::Fq>,
    execution_trace_ood_evals: Vec<S::Fq>,
    composition_trace_ood_evals: Vec<S::Fq>,
    pow_nonce: u64,
}

// impl<'a, A: Air, D: Digest> ProverChannel<'a, A, D> {
impl<'a, S: Stark> ProverChannel<'a, S> {
    pub fn new(air: &'a Air<S::AirConfig>, public_coin: S::PublicCoin) -> Self {
        ProverChannel {
            air,
            public_coin,
            extension_trace_commitment: None,
            base_trace_commitment: S::Digest::default(),
            composition_trace_commitment: S::Digest::default(),
            execution_trace_ood_evals: Vec::new(),
            composition_trace_ood_evals: Vec::new(),
            fri_layer_commitments: Vec::new(),
            fri_remainder_coeffs: Vec::new(),
            pow_nonce: 0,
        }
    }

    pub fn commit_base_trace(&mut self, commitment: S::Digest) {
        self.public_coin.reseed_with_digest(&commitment);
        self.base_trace_commitment = commitment;
    }

    pub fn commit_extension_trace(&mut self, commitment: S::Digest) {
        self.public_coin.reseed_with_digest(&commitment);
        self.extension_trace_commitment = Some(commitment);
    }

    pub fn commit_composition_trace(&mut self, commitment: S::Digest) {
        self.public_coin.reseed_with_digest(&commitment);
        self.composition_trace_commitment = commitment;
    }

    pub fn get_ood_point(&mut self) -> S::Fq {
        self.public_coin.draw()
    }

    pub fn send_ood_evals(
        &mut self,
        execution_trace_oods: Vec<S::Fq>,
        composition_trace_oods: Vec<S::Fq>,
    ) {
        let ood_evals = [execution_trace_oods.clone(), composition_trace_oods.clone()].concat();
        self.public_coin.reseed_with_field_elements(&ood_evals);
        self.execution_trace_ood_evals = execution_trace_oods;
        self.composition_trace_ood_evals = composition_trace_oods;
    }

    pub fn grind_fri_commitments(&mut self) {
        let grinding_factor = self.air.options().grinding_factor;
        if grinding_factor == 0 {
            // skip if there is no grinding required
            return;
        }

        let nonce = self
            .public_coin
            .grind_proof_of_work(grinding_factor)
            .expect("nonce not found");
        assert!(self
            .public_coin
            .verify_proof_of_work(grinding_factor, nonce));

        self.pow_nonce = nonce;
        self.public_coin.reseed_with_int(self.pow_nonce);
    }

    pub fn get_fri_query_positions(&mut self) -> BTreeSet<usize> {
        // TODO: voulnerability if multiple positions are the same
        let lde_domain_size = self.air.trace_len() * self.air.lde_blowup_factor();
        let num_queries = self.air.options().num_queries as usize;
        self.public_coin.draw_queries(num_queries, lde_domain_size)
    }

    pub fn build_proof(
        self,
        trace_queries: Queries<S>,
        fri_proof: FriProof<S::Fq, S::Digest, S::MerkleTree>,
    ) -> Proof<S> {
        Proof {
            options: self.air.options(),
            trace_len: self.air.trace_len(),
            base_trace_commitment: self.base_trace_commitment,
            extension_trace_commitment: self.extension_trace_commitment,
            composition_trace_commitment: self.composition_trace_commitment,
            execution_trace_ood_evals: self.execution_trace_ood_evals,
            composition_trace_ood_evals: self.composition_trace_ood_evals,
            pow_nonce: self.pow_nonce,
            fri_proof,
            trace_queries,
        }
    }
}

// FRI prover channel implementation
// Inspired by Winterfell: https://github.com/facebook/winterfell/blob/main/fri/src/prover/channel.rs
impl<'a, S: Stark> fri::ProverChannel for ProverChannel<'a, S> {
    type Digest = S::Digest;
    type Field = S::Fq;

    fn commit_fri_layer(&mut self, commitment: S::Digest) {
        self.public_coin.reseed_with_digest(&commitment);
        self.fri_layer_commitments.push(commitment);
    }

    fn commit_remainder(&mut self, remainder_coeffs: &[Self::Field]) {
        self.public_coin
            .reseed_with_field_element_vector(remainder_coeffs);
        self.fri_remainder_coeffs = remainder_coeffs.to_vec();
    }

    fn draw_fri_alpha(&mut self) -> S::Fq {
        self.public_coin.draw()
    }
}

// TODO: maybe just have a VerifierChannel
#[derive(Debug, Clone, CanonicalDeserialize, CanonicalSerialize)]
pub struct VerifierChannelArtifacts<F: Field> {
    pub air_challenges: Challenges<F>,
    pub air_hints: Hints<F>,
    pub fri_alphas: Vec<F>,
    pub query_positions: Vec<usize>,
}
