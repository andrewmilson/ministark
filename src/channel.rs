use crate::air::AirConfig;
use crate::fri;
use crate::fri::FriProof;
use crate::random::PublicCoin;
use crate::trace::Queries;
use crate::Air;
use crate::Proof;
use alloc::vec::Vec;
use digest::generic_array::GenericArray;
use digest::Digest;
use digest::Output;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::BTreeSet;

pub struct ProverChannel<'a, A: AirConfig, D: Digest, C: PublicCoin> {
    air: &'a Air<A>,
    pub public_coin: C,
    base_trace_commitment: Output<D>,
    extension_trace_commitment: Option<Output<D>>,
    composition_trace_commitment: Output<D>,
    fri_layer_commitments: Vec<Output<D>>,
    execution_trace_ood_evals: Vec<A::Fq>,
    composition_trace_ood_evals: Vec<A::Fq>,
    pow_nonce: u64,
}

// impl<'a, A: Air, D: Digest> ProverChannel<'a, A, D> {
impl<'a, A: AirConfig, D: Digest, C: PublicCoin<Field = A::Fq, Digest = D>>
    ProverChannel<'a, A, D, C>
{
    pub fn new(air: &'a Air<A>, public_coin: C) -> Self {
        ProverChannel {
            air,
            public_coin,
            extension_trace_commitment: None,
            base_trace_commitment: GenericArray::default(),
            composition_trace_commitment: GenericArray::default(),
            execution_trace_ood_evals: Vec::default(),
            composition_trace_ood_evals: Vec::default(),
            fri_layer_commitments: Vec::default(),
            pow_nonce: 0,
        }
    }

    pub fn commit_base_trace(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed_with_hash(commitment);
        self.base_trace_commitment = commitment.clone();
    }

    pub fn commit_extension_trace(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed_with_hash(commitment);
        self.extension_trace_commitment = Some(commitment.clone());
    }

    pub fn commit_composition_trace(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed_with_hash(commitment);
        self.composition_trace_commitment = commitment.clone();
    }

    pub fn get_ood_point(&mut self) -> C::Field {
        self.public_coin.draw()
    }

    pub fn send_execution_trace_ood_evals(&mut self, evals: Vec<A::Fq>) {
        for eval in &evals {
            self.public_coin.reseed_with_field_element(eval);
        }
        self.execution_trace_ood_evals = evals;
    }

    pub fn send_composition_trace_ood_evals(&mut self, evals: Vec<A::Fq>) {
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
        trace_queries: Queries<A::Fp, A::Fq>,
        fri_proof: FriProof<A::Fq>,
    ) -> Proof<A::Fp, A::Fq> {
        Proof {
            options: self.air.options(),
            trace_len: self.air.trace_len(),
            base_trace_commitment: self.base_trace_commitment.to_vec(),
            extension_trace_commitment: self.extension_trace_commitment.map(|o| o.to_vec()),
            composition_trace_commitment: self.composition_trace_commitment.to_vec(),
            execution_trace_ood_evals: self.execution_trace_ood_evals,
            composition_trace_ood_evals: self.composition_trace_ood_evals,
            pow_nonce: self.pow_nonce,
            fri_proof,
            trace_queries,
        }
    }
}

// FRI prover channel implementation
impl<'a, A: AirConfig, D: Digest, C: PublicCoin<Field = A::Fq, Digest = D>> fri::ProverChannel
    for ProverChannel<'a, A, D, C>
{
    type Digest = D;
    type Field = A::Fq;

    fn commit_fri_layer(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed_with_hash(commitment);
        self.fri_layer_commitments.push(commitment.clone());
    }

    fn draw_fri_alpha(&mut self) -> A::Fq {
        self.public_coin.draw()
    }
}
