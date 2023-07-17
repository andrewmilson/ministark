use crate::fri;
use crate::fri::FriProof;
use crate::merkle::MerkleTree;
use crate::random::PublicCoin;
use crate::stark::Stark;
use crate::trace::Queries;
use crate::Air;
use crate::Proof;
use alloc::vec::Vec;
use digest::generic_array::GenericArray;
use digest::Output;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::BTreeSet;

pub struct ProverChannel<'a, S: Stark> {
    air: &'a Air<S::AirConfig>,
    pub public_coin: S::PublicCoin,
    base_trace_commitment: Output<S::Digest>,
    extension_trace_commitment: Option<Output<S::Digest>>,
    composition_trace_commitment: Output<S::Digest>,
    fri_layer_commitments: Vec<Output<S::Digest>>,
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
            base_trace_commitment: GenericArray::default(),
            composition_trace_commitment: GenericArray::default(),
            execution_trace_ood_evals: Vec::new(),
            composition_trace_ood_evals: Vec::new(),
            fri_layer_commitments: Vec::new(),
            fri_remainder_coeffs: Vec::new(),
            pow_nonce: 0,
        }
    }

    pub fn commit_base_trace(&mut self, commitment: &Output<S::Digest>) {
        self.public_coin.reseed_with_hash(commitment);
        self.base_trace_commitment = commitment.clone();
    }

    pub fn commit_extension_trace(&mut self, commitment: &Output<S::Digest>) {
        self.public_coin.reseed_with_hash(commitment);
        self.extension_trace_commitment = Some(commitment.clone());
    }

    pub fn commit_composition_trace(&mut self, commitment: &Output<S::Digest>) {
        self.public_coin.reseed_with_hash(commitment);
        self.composition_trace_commitment = commitment.clone();
    }

    pub fn get_ood_point(&mut self) -> S::Fq {
        self.public_coin.draw()
    }

    pub fn send_execution_trace_ood_evals(&mut self, evals: Vec<S::Fq>) {
        for eval in &evals {
            self.public_coin.reseed_with_field_element(eval);
        }
        self.execution_trace_ood_evals = evals;
    }

    pub fn send_composition_trace_ood_evals(&mut self, evals: Vec<S::Fq>) {
        for eval in &evals {
            self.public_coin.reseed_with_field_element(eval);
        }
        self.composition_trace_ood_evals = evals;
    }

    pub fn grind_fri_commitments(&mut self) {
        let grinding_factor = self.air.options().grinding_factor;
        if grinding_factor == 0 {
            // skip if there is no grinding required
            return;
        }

        #[cfg(not(feature = "parallel"))]
        let nonce = (1..u64::MAX).find(|&nonce| {
            self.public_coin
                .verify_proof_of_work(grinding_factor, nonce)
        });

        #[cfg(feature = "parallel")]
        let nonce = (1..u64::MAX).into_par_iter().find_any(|&nonce| {
            self.public_coin
                .verify_proof_of_work(grinding_factor, nonce)
        });

        self.pow_nonce = nonce.expect("nonce not found");
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
        trace_queries: Queries<S::Fp, S::Fq, <S::MerkleTree as MerkleTree>::Proof>,
        fri_proof: FriProof<S::Fq, S::Digest, S::MerkleTree>,
    ) -> Proof<S::Fp, S::Fq, S::Digest, S::MerkleTree> {
        Proof {
            options: self.air.options(),
            trace_len: self.air.trace_len(),
            base_trace_commitment: self.base_trace_commitment.into(),
            extension_trace_commitment: self.extension_trace_commitment.map(|o| o.into()),
            composition_trace_commitment: self.composition_trace_commitment.into(),
            execution_trace_ood_evals: self.execution_trace_ood_evals,
            composition_trace_ood_evals: self.composition_trace_ood_evals,
            pow_nonce: self.pow_nonce,
            fri_proof,
            trace_queries,
        }
    }
}

// FRI prover channel implementation
impl<'a, S: Stark> fri::ProverChannel for ProverChannel<'a, S> {
    type Digest = S::Digest;
    type Field = S::Fq;

    fn commit_fri_layer(&mut self, commitment: &Output<S::Digest>) {
        self.public_coin.reseed_with_hash(commitment);
        self.fri_layer_commitments.push(commitment.clone());
    }

    fn commit_remainder(&mut self, remainder_coeffs: &[Self::Field]) {
        self.public_coin
            .reseed_with_field_elements(remainder_coeffs);
        self.fri_remainder_coeffs = remainder_coeffs.to_vec();
    }

    fn draw_fri_alpha(&mut self) -> S::Fq {
        self.public_coin.draw()
    }
}
