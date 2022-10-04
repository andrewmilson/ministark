use crate::random::PublicCoin;
use crate::Air;
use anyhow::Result;
// use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use ark_serialize::SerializationError;
use digest::Digest;
use sha2::Sha256;

pub struct ProverChannel<'a, A: Air, D: Digest> {
    air: &'a A,
    public_coin: PublicCoin<D>,
}

// impl<'a, A: Air, D: Digest> ProverChannel<'a, A, D> {
impl<'a, A: Air, D: Digest> ProverChannel<'a, A, D> {
    pub fn new(air: &'a A) -> Result<Self, SerializationError> {
        let mut seed = Vec::new();
        // Seed the public coin with:
        // 1. serialized public imputs
        air.pub_inputs().serialize_compressed(&mut seed)?;
        // 2. various metadata about the air and proof
        // TODO: field bytes?
        air.trace_info().serialize_compressed(&mut seed)?;
        air.options().serialize_compressed(&mut seed)?;
        let public_coin = PublicCoin::<D>::new(&seed);
        Ok(ProverChannel { public_coin, air })
    }
}
