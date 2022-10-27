// use crate::channel::VerifierChannel;
use crate::prover::Proof;
use crate::random::PublicCoin;
use crate::Air;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use sha2::Sha256;
use std::ops::Deref;

/// Errors that are returned during verificatino of a STARK proof
#[derive(Debug)]
pub enum VerificationError {
    Fail,
    // TODO
}

impl<A: Air> Proof<A> {
    pub fn verify(self) -> Result<(), VerificationError> {
        let Proof {
            base_trace_commitment,
            extension_trace_commitment,
            composition_trace_commitment,
            trace_info,
            public_inputs,
            options,
            ..
        } = self;

        let mut seed = Vec::new();
        public_inputs.serialize_compressed(&mut seed).unwrap();
        trace_info.serialize_compressed(&mut seed).unwrap();
        options.serialize_compressed(&mut seed).unwrap();
        let mut public_coin = PublicCoin::<Sha256>::new(&seed);

        let air = A::new(trace_info, public_inputs, options);

        let base_trace_root = Output::<Sha256>::from_iter(base_trace_commitment);
        public_coin.reseed(&base_trace_root.deref());
        let challenges = air.get_challenges(&mut public_coin);

        if let Some(extension_trace_commitment) = extension_trace_commitment {
            let extension_trace_root = Output::<Sha256>::from_iter(extension_trace_commitment);
            public_coin.reseed(&extension_trace_root.deref());
        }

        let composition_trace_root = Output::<Sha256>::from_iter(composition_trace_commitment);

        Ok(())
    }
}
