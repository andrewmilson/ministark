use crate::random::PublicCoin;
use crate::Air;
use anyhow::Result;
// use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use ark_serialize::SerializationError;
use digest::Digest;
use digest::Output;

pub struct ProverChannel<'a, A: Air, D: Digest> {
    air: &'a A,
    public_coin: PublicCoin<D>,
    trace_commitments: Vec<Output<D>>,
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
            trace_commitments: Vec::new(),
        }
    }

    pub fn commit_trace(&mut self, commitment: &Output<D>) {
        self.trace_commitments.push(commitment.clone());
        self.public_coin.reseed(commitment.to_vec());
    }
}
