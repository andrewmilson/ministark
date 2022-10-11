use crate::challenges::Challenges;
use crate::random::PublicCoin;
use crate::Air;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use fast_poly::GpuField;

pub struct ProverChannel<'a, A: Air, D: Digest> {
    air: &'a A,
    public_coin: PublicCoin<D>,
    base_trace_commitment: Output<D>,
    extension_trace_commitment: Option<Output<D>>,
    constraint_commitment: Output<D>,
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
            constraint_commitment: Default::default(),
        }
    }

    pub fn commit_base_trace(&mut self, commitment: &Output<D>) {
        self.base_trace_commitment = commitment.clone();
        self.public_coin.reseed(commitment.to_vec());
    }

    pub fn commit_extension_trace(&mut self, commitment: &Output<D>) {
        self.extension_trace_commitment = Some(commitment.clone());
        self.public_coin.reseed(commitment.to_vec());
    }

    pub fn commit_constraints(&mut self, commitment: &Output<D>) {
        self.constraint_commitment = commitment.clone();
        self.public_coin.reseed(commitment.to_vec());
    }

    pub fn get_ood_point<F: GpuField>(&mut self) -> F {
        self.public_coin.draw()
    }

    pub fn send_ood_trace_states<F: GpuField>(&mut self, states: [Vec<F>; 2]) -> F {}

    pub fn send_ood_constraint_states<F: GpuField>(&mut self, states: [Vec<F>; 2]) -> F {}

    pub fn get_challenges<F: GpuField>(&mut self, num_challenges: usize) -> Challenges<F> {
        let mut rng = self.public_coin.draw_rng();
        Challenges::new(&mut rng, num_challenges)
    }
}
