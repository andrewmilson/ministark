use crate::challenges::Challenges;
use crate::composer::DeepCompositionCoeffs;
use crate::fri;
use crate::fri::FriProof;
use crate::prover::Proof;
use crate::random::PublicCoin;
use crate::trace::Queries;
use crate::Air;
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use digest::Digest;
use digest::Output;
use fast_poly::GpuField;
use std::ops::Deref;

pub struct ProverChannel<'a, A: Air, D: Digest> {
    air: &'a A,
    pub public_coin: PublicCoin<D>,
    base_trace_commitment: Output<D>,
    extension_trace_commitment: Option<Output<D>>,
    composition_trace_commitment: Output<D>,
    fri_layer_commitments: Vec<Output<D>>,
    fri_pow_nonce: u64,
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
            fri_layer_commitments: Default::default(),
            fri_pow_nonce: 0,
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

    // TODO: make this generic
    /// Output is of the form `(trace_coeffs, composition_coeffs,
    /// degree_adjustment_coeffs)`
    pub fn get_deep_composition_coeffs(&mut self) -> DeepCompositionCoeffs<A::Fp> {
        let mut rng = self.public_coin.draw_rng();

        // execution trace coeffs
        let trace_info = self.air.trace_info();
        let num_execution_trace_cols =
            trace_info.num_base_columns + trace_info.num_extension_columns;
        let mut execution_trace_coeffs = Vec::new();
        for _ in 0..num_execution_trace_cols {
            execution_trace_coeffs.push((
                A::Fp::rand(&mut rng),
                A::Fp::rand(&mut rng),
                A::Fp::rand(&mut rng),
            ));
        }

        // composition trace coeffs
        let num_composition_trace_cols = self.air.ce_blowup_factor();
        let mut composition_trace_coeffs = Vec::new();
        for _ in 0..num_composition_trace_cols {
            composition_trace_coeffs.push(A::Fp::rand(&mut rng));
        }

        DeepCompositionCoeffs {
            trace: execution_trace_coeffs,
            constraints: composition_trace_coeffs,
            degree: (A::Fp::rand(&mut rng), A::Fp::rand(&mut rng)),
        }
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

    pub fn get_challenges<F: GpuField>(&mut self, num_challenges: usize) -> Challenges<F> {
        let mut rng = self.public_coin.draw_rng();
        Challenges::new(&mut rng, num_challenges)
    }

    pub fn grind_fri_commitments(&mut self) {
        let grinding_factor = self.air.options().fri_grinding_factor as u32;
        if grinding_factor == 0 {
            // skip if there is no grinding required
            return;
        }

        // ark_std::cfg_into_iter!(1..u64::MAX).find(predicate)
        // TODO: make parallel
        let nonce = (1..u64::MAX)
            .find(|nonce| {
                let mut hasher = D::new();
                hasher.update(&self.public_coin.seed);
                hasher.update(nonce.to_le_bytes());
                let bytes = hasher.finalize();
                let seed_head = u32::from_le_bytes(bytes[..4].try_into().unwrap());
                let num_leading_zeros = seed_head.trailing_zeros();
                num_leading_zeros >= grinding_factor
            })
            .expect("nonce not found");

        self.fri_pow_nonce = nonce;
        self.public_coin.reseed(&nonce);
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
    ) -> Proof<A::Fp> {
        Proof {
            options: *self.air.options(),
            trace_info: self.air.trace_info().clone(),
            base_trace_commitment: self.base_trace_commitment.to_vec(),
            extension_trace_commitment: self.extension_trace_commitment.map(|o| o.to_vec()),
            composition_trace_commitment: self.composition_trace_commitment.to_vec(),
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
