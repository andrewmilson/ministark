use crate::fri;
use crate::fri::FriProof;
use crate::random::PublicCoin;
use crate::trace::Queries;
use crate::Air;
use crate::Proof;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use digest::Digest;
use digest::Output;
use fast_poly::GpuField;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Deref;

pub struct ProverChannel<'a, A: Air, D: Digest> {
    air: &'a A,
    pub public_coin: PublicCoin<D>,
    base_trace_commitment: Output<D>,
    extension_trace_commitment: Option<Output<D>>,
    composition_trace_commitment: Output<D>,
    fri_layer_commitments: Vec<Output<D>>,
    ood_trace_states: (Vec<A::Fp>, Vec<A::Fp>),
    ood_constraint_evaluations: Vec<A::Fp>,
    pow_nonce: u64,
}

// impl<'a, A: Air, D: Digest> ProverChannel<'a, A, D> {
impl<'a, A: Air, D: Digest> ProverChannel<'a, A, D> {
    pub fn new(air: &'a A) -> Self {
        let mut seed = Vec::new();
        // Seed the public coin with:
        // 1. serialized public imputs
        air.pub_inputs().serialize_compressed(&mut seed).unwrap();
        // 2. various metadata about the air and proof
        // TODO: field bytes?
        air.trace_info().serialize_compressed(&mut seed).unwrap();
        air.options().serialize_compressed(&mut seed).unwrap();
        let public_coin = PublicCoin::<D>::new(&seed);
        ProverChannel {
            air,
            public_coin,
            extension_trace_commitment: None,
            base_trace_commitment: Default::default(),
            composition_trace_commitment: Default::default(),
            ood_trace_states: Default::default(),
            ood_constraint_evaluations: Default::default(),
            fri_layer_commitments: Default::default(),
            pow_nonce: 0,
        }
    }

    pub fn commit_base_trace(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed(&commitment.deref());
        self.base_trace_commitment = commitment.clone();
    }

    pub fn commit_extension_trace(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed(&commitment.deref());
        self.extension_trace_commitment = Some(commitment.clone());
    }

    pub fn commit_composition_trace(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed(&commitment.deref());
        self.composition_trace_commitment = commitment.clone();
    }

    pub fn get_ood_point<F: GpuField>(&mut self) -> F {
        self.public_coin.draw()
    }

    pub fn send_ood_trace_states(&mut self, evals: &[A::Fp], next_evals: &[A::Fp]) {
        assert_eq!(evals.len(), next_evals.len());
        self.public_coin.reseed(&evals);
        self.public_coin.reseed(&next_evals);
        self.ood_trace_states = (evals.to_vec(), next_evals.to_vec());
    }

    pub fn send_ood_constraint_evaluations(&mut self, evals: &[A::Fp]) {
        self.public_coin.reseed(&evals);
        self.ood_constraint_evaluations = evals.to_vec();
    }

    pub fn grind_fri_commitments(&mut self) {
        let grinding_factor = self.air.options().grinding_factor as u32;
        if grinding_factor == 0 {
            // skip if there is no grinding required
            return;
        }

        #[cfg(not(feature = "parallel"))]
        let nonce = (1..u64::MAX)
            .find(|&nonce| self.public_coin.check_leading_zeros(nonce) >= grinding_factor);

        #[cfg(feature = "parallel")]
        let nonce = (1..u64::MAX)
            .into_par_iter()
            .find_any(|&nonce| self.public_coin.check_leading_zeros(nonce) >= grinding_factor);

        self.pow_nonce = nonce.expect("nonce not found");
        self.public_coin.reseed(&self.pow_nonce);
    }

    pub fn get_fri_query_positions(&mut self) -> Vec<usize> {
        // TODO: voulnerability if multiple positions are the same
        let num_queries = self.air.options().num_queries;
        let lde_domain_size = self.air.trace_len() * self.air.lde_blowup_factor();
        let mut rng = self.public_coin.draw_rng();
        (0..num_queries)
            .map(|_| rng.gen_range(0..lde_domain_size))
            .collect()
    }

    pub fn build_proof(
        self,
        trace_queries: Queries<A::Fp>,
        fri_proof: FriProof<A::Fp>,
    ) -> Proof<A> {
        Proof {
            options: *self.air.options(),
            trace_info: self.air.trace_info().clone(),
            base_trace_commitment: self.base_trace_commitment.to_vec(),
            extension_trace_commitment: self.extension_trace_commitment.map(|o| o.to_vec()),
            composition_trace_commitment: self.composition_trace_commitment.to_vec(),
            public_inputs: self.air.pub_inputs().clone(),
            ood_trace_states: self.ood_trace_states,
            ood_constraint_evaluations: self.ood_constraint_evaluations,
            pow_nonce: self.pow_nonce,
            fri_proof,
            trace_queries,
        }
    }
}

// FRI prover channel implementation
impl<'a, A: Air, D: Digest> fri::ProverChannel<A::Fp> for ProverChannel<'a, A, D> {
    type Digest = D;

    fn commit_fri_layer(&mut self, commitment: &Output<D>) {
        self.public_coin.reseed(&commitment.deref());
        self.fri_layer_commitments.push(commitment.clone());
    }

    fn draw_fri_alpha(&mut self) -> A::Fp {
        self.public_coin.draw()
    }
}

// pub struct VerifierChannel<'a, A: Air, D: Digest> {
//     public_coin: PublicCoin<D>,
//     base_trace_commitment: Output<D>,
//     extension_trace_commitment: Option<Output<D>>,
//     composition_trace_commitment: Output<D>,
//     ood_trace_states: (Vec<A::Fp>, Vec<A::Fp>),
//     ood_constraint_evaluations: Vec<A::Fp>,
//     fri_pow_nonce: u64,
// }

// impl<'a, A: Air, D: Digest> VerifierChannel<'a, A, D> {
//     pub fn new(proof: &Proof<A>) -> Result<Self, VerificationError> {
//         let Proof {
//             base_trace_commitment,
//             extension_trace_commitment,
//             composition_trace_commitment,
//             fri_pow_nonce,
//             ..
//         } = proof;

//         let seed = Vec::new();
//         proof.public_inputs.serialize_compressed(&mut seed).unwrap();
//         proof.trace_info.serialize_compressed(&mut seed).unwrap();
//         proof.options.serialize_compressed(&mut seed).unwrap();
//         let public_coin = PublicCoin::<Sha256>::new(&seed);

//         let base_trace_commitment =
// Output::<D>::from_iter(base_trace_commitment);         let
// extension_trace_commitment =
// extension_trace_commitment.map(Output::<D>::from_iter);         let
// composition_trace_commitment =
// Output::<D>::from_iter(composition_trace_commitment);

//         Ok(VerifierChannel {
//             public_coin,
//             fri_pow_nonce,
//             base_trace_commitment,
//             extension_trace_commitment,
//             composition_trace_commitment,
//         })
//     }
// }
