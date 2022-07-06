use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
};

use rand::Rng;

use crate::rescue_prime::RescuePrime;

use mini_stark::{
    polynomial::Polynomial, prime_field_u128, ProofObject, ProofStream, StandardProofStream, Stark,
    StarkElement,
};

pub struct SignatureProofStream<E: StarkElement> {
    standard_proof_stream: StandardProofStream<E>,
    document: String,
    prefix: u64,
}

impl<E: StarkElement> SignatureProofStream<E> {
    pub fn new(document: String) -> SignatureProofStream<E> {
        let mut hash = DefaultHasher::new();
        document.hash(&mut hash);
        let prefix = hash.finish();

        SignatureProofStream {
            standard_proof_stream: StandardProofStream::new(),
            document,
            prefix,
        }
    }
}

impl<E: StarkElement> ProofStream<E> for SignatureProofStream<E> {
    fn prover_fiat_shamir(&self) -> u64 {
        let mut hash = DefaultHasher::new();
        let prefix_bytes = self.prefix.to_be_bytes();
        prefix_bytes
            .into_iter()
            .chain(self.standard_proof_stream.serialize())
            .collect::<Vec<u8>>()
            .hash(&mut hash);

        hash.finish()
    }

    fn verifier_fiat_shamir(&self) -> u64 {
        let prefix_bytes = self.prefix.to_be_bytes();
        let serialized_objects = serde_json::to_vec(
            &self.standard_proof_stream.objects[..self.standard_proof_stream.read_index].to_vec(),
        )
        .unwrap();
        let mut hash = DefaultHasher::new();
        prefix_bytes
            .into_iter()
            .chain(serialized_objects)
            .collect::<Vec<u8>>()
            .hash(&mut hash);

        hash.finish()
    }

    fn deserialize(&self, bytes: &[u8]) -> SignatureProofStream<E> {
        let mut signature_proof_stream = SignatureProofStream::new(self.document.clone());
        signature_proof_stream.standard_proof_stream.objects =
            serde_json::from_slice(bytes).unwrap();
        signature_proof_stream
    }

    fn push(&mut self, object: ProofObject<E>) {
        self.standard_proof_stream.push(object);
    }

    fn pull(&mut self) -> ProofObject<E> {
        self.standard_proof_stream.pull()
    }

    fn serialize(&self) -> Vec<u8> {
        self.standard_proof_stream.serialize()
    }
}

pub struct StarkSignatureScheme {
    rescue_prime: RescuePrime,
    stark: Stark<prime_field_u128::BaseElement>,
    transition_zerofier: Polynomial<prime_field_u128::BaseElement>,
    transition_zerofier_codeword: Vec<prime_field_u128::BaseElement>,
    transition_zerofier_root: u64,
}

impl StarkSignatureScheme {
    pub fn new() -> StarkSignatureScheme {
        let expansion_factor = 4;
        // let num_colinearity_checks = 64;
        let num_colinearity_checks = 64;
        let security_level = num_colinearity_checks * 2;

        let rescue_prime = RescuePrime::new();
        let num_cycles = rescue_prime.N + 1;
        let state_width = rescue_prime.m;

        let stark = Stark::new(
            expansion_factor,
            num_colinearity_checks,
            security_level,
            state_width,
            num_cycles,
            3,
        );

        let (transition_zerofier, transition_zerofier_codeword, transition_zerofier_root) =
            stark.preprocess();

        StarkSignatureScheme {
            rescue_prime,
            stark,
            transition_zerofier,
            transition_zerofier_codeword,
            transition_zerofier_root,
        }
    }

    pub fn stark_prove<T: ProofStream<prime_field_u128::BaseElement>>(
        &self,
        input_element: prime_field_u128::BaseElement,
        proof_stream: &mut T,
    ) -> Vec<u8> {
        let output_element = self.rescue_prime.hash(input_element);
        let trace = self.rescue_prime.trace(input_element);
        let transition_constraints = self.rescue_prime.transition_constraints(self.stark.omicron);
        let boundary_constraints = self.rescue_prime.boundary_constraints(output_element);
        let proof = self.stark.prove(
            trace,
            transition_constraints,
            &boundary_constraints,
            &self.transition_zerofier,
            &self.transition_zerofier_codeword,
            proof_stream,
        );

        proof
    }

    pub fn stark_verify(
        &self,
        output_element: prime_field_u128::BaseElement,
        stark_proof: &[u8],
        proof_stream: &mut SignatureProofStream<prime_field_u128::BaseElement>,
    ) -> Result<(), &str> {
        let boundary_constraints = self.rescue_prime.boundary_constraints(output_element);
        let transition_constraints = self.rescue_prime.transition_constraints(self.stark.omicron);
        self.stark.verify(
            stark_proof,
            transition_constraints,
            &boundary_constraints,
            self.transition_zerofier_root,
            proof_stream,
        )
    }

    pub fn keygen(&self) -> (prime_field_u128::BaseElement, prime_field_u128::BaseElement) {
        let mut rng = rand::thread_rng();
        let secret_key = prime_field_u128::BaseElement::from(rng.gen::<u64>());
        let public_key = self.rescue_prime.hash(secret_key);
        (secret_key, public_key)
    }

    pub fn sign(&self, secret_key: prime_field_u128::BaseElement, document: String) -> Vec<u8> {
        let mut signature_proof_stream = SignatureProofStream::new(document);

        self.stark_prove(secret_key, &mut signature_proof_stream)
    }

    pub fn verify(
        &self,
        public_key: prime_field_u128::BaseElement,
        document: String,
        signature: &[u8],
    ) -> Result<(), &str> {
        let mut signature_proof_stream = SignatureProofStream::new(document);
        self.stark_verify(public_key, signature, &mut signature_proof_stream)
    }
}
