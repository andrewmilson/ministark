use crate::challenges::Challenges;
// use crate::channel::VerifierChannel;
use crate::prover::Proof;
use crate::random::PublicCoin;
use crate::utils::Timer;
use crate::Air;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
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
        let _timer = Timer::new("Verification");

        let Proof {
            base_trace_commitment,
            extension_trace_commitment,
            composition_trace_commitment,
            ood_constraint_evaluations,
            ood_trace_states,
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

        let composition_coeffs = air.get_constraint_composition_coeffs(&mut public_coin);
        let composition_trace_root = Output::<Sha256>::from_iter(composition_trace_commitment);
        public_coin.reseed(&composition_trace_root.deref());

        let z = public_coin.draw::<A::Fp>();
        println!("{z}");
        public_coin.reseed(&ood_trace_states.0);
        public_coin.reseed(&ood_trace_states.1);
        let calculated_ood_constraint_evaluation = ood_constraint_evaluation(
            composition_coeffs,
            &challenges,
            &ood_trace_states.0,
            &ood_trace_states.1,
            &air,
            z,
        );

        // TODO: why not proof only include single value for constraint eval?
        let mut acc = A::Fp::one();
        let provided_ood_constraint_evaluation =
            ood_constraint_evaluations
                .iter()
                .fold(A::Fp::zero(), |mut res, value| {
                    res += acc * value;
                    acc *= z;
                    res
                });

        if calculated_ood_constraint_evaluation != provided_ood_constraint_evaluation {
            println!("NOOOOOOOO");
        }

        println!("BOOM!");

        Ok(())
    }
}

fn ood_constraint_evaluation<A: Air>(
    mut composition_coefficients: Vec<(A::Fp, A::Fp)>,
    challenges: &Challenges<A::Fp>,
    curr_trace_evals: &[A::Fp],
    next_trace_evals: &[A::Fp],
    air: &A,
    x: A::Fp,
) -> A::Fp {
    // TODO: refactor constraint and their divisors so they are grouped together
    let boundary_constraints = air.boundary_constraints();
    let transition_constraints = air.transition_constraints();
    let terminal_constraints = air.terminal_constraints();

    let boundary_divisor_degree = 1;
    let transition_divisor_degree = air.trace_len() - 1;
    let terminal_divisor_degree = 1;

    let trace_domain = air.trace_domain();
    let first_trace_x = A::Fp::one();
    let last_trace_x = trace_domain.group_gen_inv;
    // TODO docs
    let boundary_divisor = (x - first_trace_x).inverse().unwrap();
    let terminal_divisor = (x - last_trace_x).inverse().unwrap();
    let transition_divisor = (x - last_trace_x)
        * trace_domain
            .evaluate_vanishing_polynomial(x)
            .inverse()
            .unwrap();

    let boundary_iter = boundary_constraints
        .iter()
        .map(|constraint| (constraint, boundary_divisor, boundary_divisor_degree));
    let transition_iter = transition_constraints
        .iter()
        .map(|constraint| (constraint, transition_divisor, transition_divisor_degree));
    let terminal_iter = terminal_constraints
        .iter()
        .map(|constraint| (constraint, terminal_divisor, terminal_divisor_degree));

    let mut result = A::Fp::zero();
    let trace_degree = air.trace_len() - 1;
    let composition_degree = air.composition_degree();
    for (constraint, divisor, divisor_degree) in
        boundary_iter.chain(transition_iter).chain(terminal_iter)
    {
        // TODO: proper errors
        let (alpha, beta) = composition_coefficients.pop().unwrap();
        let evaluation = constraint.evaluate(challenges, curr_trace_evals, next_trace_evals);
        let quotient = evaluation * divisor;

        // TODO: handle case when degree is 0?
        let evaluation_degree = constraint.degree() * trace_degree - divisor_degree;
        assert!(evaluation_degree <= composition_degree);
        let degree_adjustment = (composition_degree - evaluation_degree) as u64;

        result += quotient * (alpha * x.pow([degree_adjustment]) + beta)
    }

    result
}
