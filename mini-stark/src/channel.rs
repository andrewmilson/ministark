use crate::challenges::Challenges;
use crate::random::PublicCoin;
use crate::Air;
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use digest::Digest;
use digest::Output;
use fast_poly::GpuField;
use rand_chacha::ChaCha20Rng;
use std::ops::Deref;

pub struct ProverChannel<'a, A: Air, D: Digest> {
    air: &'a A,
    pub public_coin: PublicCoin<D>,
    base_trace_commitment: Output<D>,
    extension_trace_commitment: Option<Output<D>>,
    composition_trace_commitment: Output<D>,
    ood_trace_states: (Vec<A::Fp>, Vec<A::Fp>),
    ood_constraint_evaluations: Vec<A::Fp>,
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

    // TODO: make this generic
    pub fn get_constraint_composition_coeffs(&mut self) -> Vec<(A::Fp, A::Fp)> {
        let mut rng = self.public_coin.draw_rng();
        (0..self.air.num_constraints())
            .map(|_| (A::Fp::rand(&mut rng), A::Fp::rand(&mut rng)))
            .collect()
    }

    pub fn send_ood_trace_states(&mut self, evals: Vec<A::Fp>, next_evals: Vec<A::Fp>) {
        self.public_coin.reseed(&evals);
        self.public_coin.reseed(&next_evals);
        self.ood_trace_states = (evals, next_evals);
    }

    pub fn send_ood_constraint_evaluations(&mut self, evals: Vec<A::Fp>) {
        self.public_coin.reseed(&evals);
        self.ood_constraint_evaluations = evals;
    }

    pub fn get_challenges<F: GpuField>(&mut self, num_challenges: usize) -> Challenges<F> {
        let mut rng = self.public_coin.draw_rng();
        Challenges::new(&mut rng, num_challenges)
    }
}
