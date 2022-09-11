use crate::number_theory_transform::*;
use crate::polynomial::MultivariatePolynomial;
use crate::polynomial::Polynomial;
use crate::Fri;
use crate::MerkleTree;
use crate::ProofObject;
use crate::ProofStream;
use algebra::StarkFelt;
use fast_poly::allocator::PageAlignedAllocator;
use rand::Rng;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::Hash;
use std::hash::Hasher;
use std::iter::once;
use std::iter::repeat;
use std::time::Duration;
use std::time::Instant;

pub struct Stark<E> {
    expansion_factor: usize,
    num_randomizers: usize,
    randomized_trace_length: usize,
    num_registers: usize,
    original_trace_length: usize,
    pub omega: E,
    pub omicron: E,
    omicron_domain: Vec<E>,
    pub fri: Fri<E>,
}

impl<E: StarkFelt<BaseFelt = E>> Stark<E> {
    pub fn new(
        expansion_factor: usize,
        num_colinearity_checks: usize,
        security_level: usize,
        num_registers: usize,
        num_cycles: usize,
        transition_constraints_degree: usize,
    ) -> Self {
        // TODO: fix
        // assert!(
        //     E::FIELD_ORDER_BITS >= security_level as u32,
        //     "the bits needed to represent the order of the field must at least be the
        // security level" );
        assert!(
            expansion_factor.is_power_of_two(),
            "expansion factor must be a power of 2"
        );
        assert!(
            num_colinearity_checks * 2 >= security_level,
            "number of colinearity checks must be at least half of security level"
        );

        let original_trace_length = num_cycles;
        let num_randomizers = 4 * num_colinearity_checks;
        let randomized_trace_length = original_trace_length + num_randomizers;
        let omicron_domain_length =
            (randomized_trace_length * transition_constraints_degree + 1).next_power_of_two();

        let fri_domain_length = omicron_domain_length * expansion_factor;

        let omega = E::get_root_of_unity(fri_domain_length.ilog2());
        let omicron = E::get_root_of_unity(omicron_domain_length.ilog2());

        let mut accumulator = E::one();
        let omicron_domain = (0..omicron_domain_length)
            .map(|_| {
                let res = accumulator;
                accumulator *= omicron;
                res
            })
            .collect::<Vec<E>>();

        let fri = Fri::new(
            E::GENERATOR,
            omega,
            fri_domain_length,
            expansion_factor,
            num_colinearity_checks,
        );

        Stark {
            expansion_factor,
            num_randomizers,
            randomized_trace_length,
            num_registers,
            original_trace_length,
            omega,
            omicron,
            omicron_domain,
            fri,
        }
    }

    // Returns: zerofier
    pub fn preprocess(&self) -> (Polynomial<E>, Vec<E>, u64) {
        let transition_zerofier = fast_zerofier(
            &self.omicron_domain[..(self.original_trace_length - 1)],
            // self.omicron,
            // self.omicron_domain.len(),
        );
        let transition_zerofier_codeword =
            fast_coset_evaluate(&transition_zerofier, E::GENERATOR, self.fri.domain_length);
        let transition_zerofier_root = MerkleTree::commit(&transition_zerofier_codeword);
        (
            transition_zerofier,
            transition_zerofier_codeword,
            transition_zerofier_root,
        )
    }

    fn transition_degree_bounds(
        &self,
        transition_constraints: &[MultivariatePolynomial<E>],
    ) -> Vec<u128> {
        let point_degrees = once(1).chain(
            repeat((self.original_trace_length + self.num_randomizers - 1) as u128)
                .take(2 * self.num_registers),
        );

        transition_constraints
            .iter()
            .map(|transition_constraint| {
                transition_constraint
                    .powers
                    .iter()
                    .map(|pad| {
                        pad.iter()
                            .zip(point_degrees.clone())
                            .map(|(power, point_degree)| power * point_degree)
                            .sum()
                    })
                    .max()
                    .unwrap_or(0)
            })
            .collect()
    }

    // Degrees of each transition quotient
    // Transition quotients are the result of dividing out the zerofier.
    fn transition_quotient_degree_bounds(
        &self,
        transition_constraints: &[MultivariatePolynomial<E>],
    ) -> Vec<u128> {
        self.transition_degree_bounds(transition_constraints)
            .into_iter()
            .map(|degree| degree - (self.original_trace_length - 1) as u128)
            .collect()
    }

    fn max_degree(&self, transition_constraints: &[MultivariatePolynomial<E>]) -> u128 {
        let max_degree = self
            .transition_quotient_degree_bounds(transition_constraints)
            .into_iter()
            .max()
            .unwrap();
        (1 << (max_degree.ilog2() + 1)) - 1
    }

    // fn transition_zerofier(&self) -> Polynomial<E> {
    //     let domain = &self.omicron_domain[0..(self.original_trace_length - 1)];
    //     Polynomial::zerofier_domain(domain)
    // }

    fn boundary_zerofiers(&self, boundary: &[(usize, usize, E)]) -> Vec<Polynomial<E>> {
        (0..self.num_registers)
            .map(|register| {
                Polynomial::zerofier_domain(
                    &boundary
                        .iter()
                        .copied()
                        .filter_map(|(cycle, boundary_register, _)| {
                            if boundary_register == register {
                                Some(self.omicron.pow(&[cycle as u64]))
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<E>>(),
                )
            })
            .collect()
    }

    fn boundary_interpolants(&self, boundary: &[(usize, usize, E)]) -> Vec<Polynomial<E>> {
        (0..self.num_registers)
            .map(|register| {
                let (domain, values): (Vec<E>, Vec<E>) = boundary
                    .iter()
                    .copied()
                    .filter_map(|(cycle, boundary_register, value)| {
                        if boundary_register == register {
                            Some((self.omicron.pow(&[cycle as u64]), value))
                        } else {
                            None
                        }
                    })
                    .unzip();

                Polynomial::interpolate(&domain, &values)
            })
            .collect()
    }

    fn boundary_quotient_degree_bounds(
        &self,
        randomized_trace_length: usize,
        boundary: &[(usize, usize, E)],
    ) -> Vec<usize> {
        let randomized_trace_degree = randomized_trace_length - 1;
        self.boundary_zerofiers(boundary)
            .iter()
            .map(|boundary_zerofier| randomized_trace_degree - boundary_zerofier.degree() as usize)
            .collect()
    }

    fn sample_weights(&self, num_weights: usize, randomness: u64) -> Vec<E> {
        (0..num_weights)
            .map(|i| {
                let mut hash = DefaultHasher::new();
                (randomness + i as u64).hash(&mut hash);
                E::from(hash.finish())
            })
            .collect()
    }

    pub fn prove<T: ProofStream<E>>(
        &self,
        trace: Vec<Vec<E, PageAlignedAllocator>>,
        transition_constraints: Vec<MultivariatePolynomial<E>>,
        boundary: &[(usize, usize, E)],
        transition_zerofier: &Polynomial<E>,
        transition_zerofier_codeword: &[E],
        proof_stream: &mut T,
    ) -> Vec<u8> {
        let proof_start = Instant::now();
        // Vec<(Section name, duration)>
        let mut profile: Vec<(&str, Duration)> = vec![];

        let mut rng = rand::thread_rng();

        // concatenate randomizers
        // let trace = trace
        //     .into_iter()
        //     .chain((0..self.num_randomizers).map(|_| {
        //         (0..self.num_registers)
        //             .map(|_| E::from(rng.gen::<u64>()))
        //             .collect::<Vec<E>>()
        //     }))
        //     .collect::<Vec<Vec<E>>>();

        // concatenate randomizers
        let trace = trace
            .into_iter()
            .chain((0..self.num_randomizers).map(|_| {
                let mut res = Vec::new_in(PageAlignedAllocator);
                res.resize(self.num_registers, E::zero());
                res.iter_mut().for_each(|v| *v = E::from(rng.gen::<u64>()));
                res
            }))
            .collect::<Vec<Vec<E, PageAlignedAllocator>>>();

        // interpolate
        let now = Instant::now();
        let trace_domain = (0..trace.len())
            .map(|i| self.omicron.pow(&[i as u64]))
            .collect::<Vec<E>>();
        let trace_polynomials = (0..self.num_registers)
            .map(|register| {
                let single_trace = trace
                    .iter()
                    .map(|state| state[register])
                    .collect::<Vec<E>>();
                fast_interpolate(
                    &trace_domain,
                    &single_trace,
                    // self.omicron,
                    // self.omicron_domain.len(),
                )
            })
            .collect::<Vec<Polynomial<E>>>();
        profile.push(("Interpolation of trace", now.elapsed()));

        // subtract boundary interpolants and divide out boundary zerofiers
        let boundary_interpolants = self.boundary_interpolants(boundary);
        let boundary_zerofiers = self.boundary_zerofiers(boundary);
        let boundary_quotients = trace_polynomials
            .iter()
            .cloned()
            .zip(boundary_interpolants.into_iter())
            .zip(boundary_zerofiers.into_iter())
            .map(|((trace_polynomial, interpolant), zerofier)| {
                (trace_polynomial - interpolant) / zerofier
            })
            .collect::<Vec<Polynomial<E>>>();

        // commit to the boundary quotients
        let boundary_quotient_codewords = (0..self.num_registers)
            .map(|register| {
                fast_coset_evaluate(
                    &boundary_quotients[register],
                    E::GENERATOR,
                    self.fri.domain_length,
                )
            })
            .collect::<Vec<Vec<E>>>();

        for codeword in &boundary_quotient_codewords {
            let merkle_root = MerkleTree::commit(codeword);
            proof_stream.push(ProofObject::MerkleRoot(merkle_root));
        }

        // symbolically evaluate transition constraints
        let now = Instant::now();
        let point = once(Polynomial::x())
            .chain(trace_polynomials.clone())
            .chain(
                trace_polynomials
                    .into_iter()
                    .map(|trace_polynomial| trace_polynomial.scale(self.omicron)),
            )
            .collect::<Vec<Polynomial<E>>>();

        // Slow
        // ====
        let transition_polynomials = transition_constraints
            .iter()
            .map(|transition_constraint| transition_constraint.evaluate_symbolic(&point))
            .collect::<Vec<Polynomial<E>>>();
        profile.push(("Symbolic transition polynomial evaluation", now.elapsed()));

        // This performs poorly
        // ====================
        // let transition_polynomials = transition_constraints
        //     .iter()
        //     .map(|transition_constraint| {
        //         fast_evaluate_symbolic(
        //             transition_constraint,
        //             &point,
        //             self.omicron,
        //             self.omicron_domain.len(),
        //         )
        //     })
        //     .collect::<Vec<Polynomial<E>>>();

        // divide out the zerofier
        let now = Instant::now();
        let transition_quotients = transition_polynomials
            .iter()
            .map(|transition_polynomial| {
                fast_coset_divide(
                    transition_polynomial,
                    transition_zerofier,
                    E::GENERATOR,
                    self.omicron,
                    self.omicron_domain.len(),
                )
            })
            .collect::<Vec<Polynomial<E>>>();
        profile.push(("Transition quotient generation", now.elapsed()));

        // commit to the randomizer polynomial
        let randomizer_polynomial = Polynomial::new(
            (0..(self.max_degree(&transition_constraints) + 1))
                .map(|_| E::from(rng.gen::<u64>()))
                .collect(),
        );
        let randomizer_codeword =
            fast_coset_evaluate(&randomizer_polynomial, E::GENERATOR, self.fri.domain_length);
        let randomizer_root = MerkleTree::commit(&randomizer_codeword);
        proof_stream.push(ProofObject::MerkleRoot(randomizer_root));

        // get weights for nonlinear combination
        // - 1 randomizer
        // - 2 for every transition quotient
        // - 3 for every boundary quotient
        let weights = self.sample_weights(
            1 + 2 * (transition_quotients.len() + boundary_quotients.len()),
            proof_stream.prover_fiat_shamir(),
        );

        assert_eq!(
            transition_quotients
                .clone()
                .into_iter()
                .map(|transition_quotient| transition_quotient.degree() as u128)
                .collect::<Vec<u128>>(),
            self.transition_quotient_degree_bounds(&transition_constraints),
            "transition quotient degrees do not match with expectation"
        );

        let now = Instant::now();
        let transition_constraints_max_degree = self.max_degree(&transition_constraints);
        let transition_quotient_degree_bounds =
            self.transition_quotient_degree_bounds(&transition_constraints);
        let boundary_quotiend_degree_bounds =
            self.boundary_quotient_degree_bounds(trace.len(), boundary);
        // println!("Randomizer degree: {}", randomizer_polynomial.degree());
        // println!(
        //     "Transition quotient degree: {:?}",
        //     transition_quotients
        //         .iter()
        //         .map(|tq| tq.degree())
        //         .collect::<Vec<_>>()
        // );
        let terms = once(randomizer_polynomial)
            .chain(transition_quotients.into_iter().enumerate().flat_map(
                |(i, transition_quotient)| {
                    let shift =
                        transition_constraints_max_degree - transition_quotient_degree_bounds[i];
                    let result = (Polynomial::x() ^ shift) * transition_quotient.clone();
                    vec![transition_quotient, result]
                },
            ))
            .chain((0..self.num_registers).flat_map(|register| {
                let shift = transition_constraints_max_degree
                    .max(boundary_quotiend_degree_bounds[register] as u128)
                    - transition_constraints_max_degree
                        .min(boundary_quotiend_degree_bounds[register] as u128);
                vec![
                    boundary_quotients[register].clone(),
                    (Polynomial::x() ^ shift) * boundary_quotients[register].clone(),
                ]
            }));
        profile.push(("Term generation", now.elapsed()));

        // println!("Terms length: {}", terms.clone().collect::<Vec<_>>().len());

        // take a weighted sum
        let combination = terms
            .zip(weights)
            .map(|(term, weight)| term * Polynomial::new(vec![weight]))
            .reduce(|a, b| a + b)
            .unwrap();

        // compute matching codeword
        let now = Instant::now();
        let combined_codeword =
            fast_coset_evaluate(&combination, E::GENERATOR, self.fri.domain_length);
        profile.push(("Combination codeword evaluation", now.elapsed()));

        // prove low degree of combination polynomial
        let now = Instant::now();
        let indices = self.fri.prove(combined_codeword, proof_stream);
        profile.push(("Fri proof generation", now.elapsed()));

        // process indices
        let now = Instant::now();
        let duplicated_indices = indices
            .iter()
            .cloned()
            .chain(
                indices
                    .iter()
                    .cloned()
                    .map(|index| (index + self.expansion_factor) % self.fri.domain_length),
            )
            .collect::<Vec<usize>>();
        let mut quadrupled_indices = duplicated_indices
            .iter()
            .cloned()
            .chain(
                duplicated_indices
                    .iter()
                    .cloned()
                    .map(|index| (index + self.fri.domain_length / 2) % self.fri.domain_length),
            )
            .collect::<Vec<usize>>();
        quadrupled_indices.sort_unstable();
        profile.push(("Index generation", now.elapsed()));

        // open indicated positions in the boundary quotient codewords
        let now = Instant::now();
        for codeword in boundary_quotient_codewords.iter() {
            for index in quadrupled_indices.iter().cloned() {
                proof_stream.push(ProofObject::YValue(codeword[index]));
                let path = MerkleTree::open(index, codeword);
                proof_stream.push(ProofObject::MerklePath(path));
            }
        }

        // as well as the randomizer
        for index in quadrupled_indices.iter().copied() {
            proof_stream.push(ProofObject::YValue(randomizer_codeword[index]));
            let path = MerkleTree::open(index, &randomizer_codeword);
            proof_stream.push(ProofObject::MerklePath(path));
        }

        // as well as the zerofier
        for index in quadrupled_indices {
            proof_stream.push(ProofObject::YValue(transition_zerofier_codeword[index]));
            let path = MerkleTree::open(index, transition_zerofier_codeword);
            proof_stream.push(ProofObject::MerklePath(path));
        }
        profile.push(("Commitment generation", now.elapsed()));

        let serialized_proof = proof_stream.serialize();

        profile.push((
            "Other",
            proof_start.elapsed() - profile.iter().map(|item| item.1).sum::<Duration>(),
        ));

        profile.sort_by(|a, b| b.1.cmp(&a.1));

        println!("============================");
        println!("STARK proof generation profile:");
        for (name, duration) in profile {
            println!("- {}: {:.1?}", name, duration);
        }
        println!("============================");

        serialized_proof
    }

    pub fn verify<T: ProofStream<E>>(
        &self,
        proof: &[u8],
        transition_constraints: Vec<MultivariatePolynomial<E>>,
        boundary: &[(usize, usize, E)],
        transition_zerofier_root: u64,
        proof_stream: &mut T,
    ) -> Result<(), &str> {
        let verification_start = Instant::now();
        // Vec<(section, duration)>
        let mut profile: Vec<(&str, Duration)> = vec![];

        // deserialize proof
        let now = Instant::now();
        let mut proof_stream = proof_stream.deserialize(proof);
        profile.push(("Deserialization of proof", now.elapsed()));

        // get merkle roots of boundary quotient codewords
        let mut boundary_quotient_roots = vec![];
        for _ in 0..self.num_registers {
            match proof_stream.pull() {
                ProofObject::MerkleRoot(root) => boundary_quotient_roots.push(root),
                _ => return Err("Expected to recieve boundary quotient merkle root"),
            }
        }

        // get merkle root of randomizer polynomial
        let randomizer_root = match proof_stream.pull() {
            ProofObject::MerkleRoot(root) => root,
            _ => return Err("Expected to receive randomizer merkle root"),
        };

        // get weights for nonlinear combination
        let weights = self.sample_weights(
            1 + 2 * (transition_constraints.len() + self.boundary_interpolants(boundary).len()),
            proof_stream.verifier_fiat_shamir(),
        );

        // verify low degree of combination polynomial
        let mut polynomial_values = vec![];

        let now = Instant::now();
        self.fri.verify(&mut proof_stream, &mut polynomial_values)?;
        profile.push(("Fri verification", now.elapsed()));

        polynomial_values.sort_by(|a, b| a.0.cmp(&b.0));

        let (indices, values): (Vec<usize>, Vec<E>) = polynomial_values.into_iter().unzip();

        // read and verify leafs, which are elements of boundary quotient codewords
        let mut duplicated_indices = indices
            .iter()
            .cloned()
            .chain(
                indices
                    .iter()
                    .cloned()
                    .map(|index| (index + self.expansion_factor) % self.fri.domain_length),
            )
            .collect::<Vec<usize>>();
        duplicated_indices.sort_unstable();

        let mut leafs = vec![];

        let now = Instant::now();
        for root in &boundary_quotient_roots {
            let mut leaf = HashMap::new();

            for index in duplicated_indices.iter().cloned() {
                leaf.insert(
                    index,
                    match proof_stream.pull() {
                        ProofObject::YValue(value) => value,
                        _ => return Err("Expected to receive a value in a codeword"),
                    },
                );

                let path = match proof_stream.pull() {
                    ProofObject::MerklePath(path) => path,
                    _ => return Err("Expected to receive merkle path"),
                };

                if !MerkleTree::verify(*root, index, &path, leaf.get(&index).unwrap()) {
                    return Err("Failed to verify merkle path");
                }
            }

            leafs.push(leaf);
        }
        profile.push(("Decommitment of boundary conditions", now.elapsed()));

        let now = Instant::now();
        let mut randomizer = HashMap::new();
        for index in duplicated_indices.iter().cloned() {
            let randomizer_value = match proof_stream.pull() {
                ProofObject::YValue(value) => value,
                _ => return Err("Expected randomizer value"),
            };

            randomizer.insert(index, randomizer_value);

            let path = match proof_stream.pull() {
                ProofObject::MerklePath(path) => path,
                _ => return Err("Expected to receive merkle path"),
            };

            if !MerkleTree::verify(randomizer_root, index, &path, randomizer_value) {
                return Err("Failed to verify merkle path for randomizer");
            }
        }
        profile.push(("Decommitment of randomizer values", now.elapsed()));

        let now = Instant::now();
        let mut transition_zerofier = HashMap::new();
        for index in duplicated_indices.iter().cloned() {
            let zerofier_value = match proof_stream.pull() {
                ProofObject::YValue(value) => value,
                _ => return Err("Expected zerofier value"),
            };

            transition_zerofier.insert(index, zerofier_value);

            let path = match proof_stream.pull() {
                ProofObject::MerklePath(path) => path,
                _ => return Err("Expected to receive merkle path"),
            };

            if !MerkleTree::verify(transition_zerofier_root, index, &path, zerofier_value) {
                return Err("Failed to verify merkle path for randomizer");
            }
        }
        profile.push(("Decommitment of zerofier values", now.elapsed()));

        let boundary_zerofiers = self.boundary_zerofiers(boundary);
        let boundary_interpolants = self.boundary_interpolants(boundary);
        let boundary_quotient_degree_bounds =
            self.boundary_quotient_degree_bounds(self.randomized_trace_length, boundary);
        let transition_quotient_degree_bounds =
            self.transition_quotient_degree_bounds(&transition_constraints);
        let transition_constraints_max_degree = self.max_degree(&transition_constraints);

        let mut transition_constraint_duration = Duration::ZERO;
        let mut combination_duration = Duration::ZERO;

        // verify leafs of combination polynomial
        for (i, current_index) in indices.into_iter().enumerate() {
            // get trace values by applying a correction to the boundary quotient values
            // (which are the leafs)
            let domain_current_index = E::GENERATOR * self.omega.pow(&[current_index as u64]);
            let next_index = (current_index + self.expansion_factor) % self.fri.domain_length;
            let domain_next_index = E::GENERATOR * self.omega.pow(&[next_index as u64]);

            let (current_trace, next_trace): (Vec<E>, Vec<E>) = (0..self.num_registers)
                .map(|register| {
                    let current = *leafs[register].get(&current_index).unwrap()
                        * boundary_zerofiers[register].evaluate(domain_current_index)
                        + boundary_interpolants[register].evaluate(domain_current_index);
                    let next = *leafs[register].get(&next_index).unwrap()
                        * boundary_zerofiers[register].evaluate(domain_next_index)
                        + boundary_interpolants[register].evaluate(domain_next_index);

                    (current, next)
                })
                .unzip();

            let point = once(domain_current_index)
                .chain(current_trace)
                .chain(next_trace)
                .collect::<Vec<E>>();

            let now = Instant::now();
            let transition_constraints_values = transition_constraints
                .iter()
                .map(|transition_constraint| transition_constraint.evaluate(&point))
                .collect::<Vec<E>>();
            transition_constraint_duration += now.elapsed();

            // compute non-linear combination
            let mut terms = vec![*randomizer.get(&current_index).unwrap()];

            for register in 0..transition_constraints_values.len() {
                let transition_constraint_value = transition_constraints_values[register];
                let quotient = transition_constraint_value / transition_zerofier[&current_index];
                terms.push(quotient);
                let shift =
                    transition_constraints_max_degree - transition_quotient_degree_bounds[register];
                // TODO: Been a while since I looked at this. Double check u64 cast is safe here
                terms.push(quotient * domain_current_index.pow(&[shift as u64]));
            }

            for register in 0..self.num_registers {
                let boundary_quotient_value = *leafs[register].get(&current_index).unwrap();
                terms.push(boundary_quotient_value);
                let shift = transition_constraints_max_degree
                    .max(boundary_quotient_degree_bounds[register] as u128)
                    - transition_constraints_max_degree
                        .min(boundary_quotient_degree_bounds[register] as u128);
                // TODO: Same as above ^^. Here too.
                terms.push(boundary_quotient_value * domain_current_index.pow(&[shift as u64]));
            }

            let now = Instant::now();
            let combination = terms
                .into_iter()
                .zip(weights.iter().cloned())
                .map(|(term, weight)| term * weight)
                .reduce(|a, b| a + b)
                .unwrap();
            combination_duration += now.elapsed();

            if combination != values[i] {
                println!(
                    "does it match on itteration {}? {}, {}",
                    i, combination, values[i]
                );
                return Err("Calculated value does not equal value from proover");
            }
        }

        profile.push((
            "Other",
            verification_start.elapsed() - profile.iter().map(|item| item.1).sum::<Duration>(),
        ));

        profile.sort_by(|a, b| b.1.cmp(&a.1));

        println!("============================");
        println!("STARK proof verification profile:");
        for (name, duration) in profile {
            println!("- {}: {:.1?}", name, duration);
        }
        println!("============================");

        Ok(())
    }
}

// #[cfg(test)]
// mod tests {
//     use crate::{
//         prime_field_u128::BaseElement, protocol::StandardProofStream,
// rescue_prime::RescuePrime,     };
//     use num_traits::One;
//     use rand::Rng;

//     // use tests::Bencher;
//     use super::*;

//     // #[bench]
//     // fn bench_interpolate_500_points(b: &mut Bencher) {
//     //     b.iter(|| {

//     //     })
//     // }

//     #[ignore]
//     #[test]
//     fn test_stark() {
//         let expansion_factor = 4;
//         let num_colinearity_checks = 2;
//         let security_level = 2;

//         let rp = RescuePrime::new();
//         let output_element = BaseElement::new(0xdeadbeef);

//         for trial in 0..20 {
//             let input_element = output_element;
//             println!("running trial with input: {}", input_element);
//             let output_element = rp.hash(input_element);
//             let num_cycles = rp.N + 1;
//             let state_width = rp.m;

//             let stark = Stark::new(
//                 expansion_factor,
//                 num_colinearity_checks,
//                 security_level,
//                 state_width,
//                 num_cycles,
//                 2,
//             );

//             let (transition_zerofier, transition_zerofier_codeword,
// transition_zerofier_root) =                 stark.preprocess();

//             // prove honestly
//             println!("honest proof generation ...");

//             // prove
//             let mut trace = rp.trace(input_element);
//             let air = rp.transition_constraints(stark.omicron);
//             let boundary = rp.boundary_constraints(output_element);
//             let mut proover_proof_stream = StandardProofStream::new();

//             // self.stark.prove(
//             //     trace,
//             //     transition_constraints,
//             //     &boundary_constraints,
//             //     &self.transition_zerofier,
//             //     &self.transition_zerofier_codeword,
//             //     proof_stream,
//             // )

//             let proof = stark.prove(
//                 trace.clone(),
//                 air.clone(),
//                 &boundary,
//                 &transition_zerofier,
//                 &transition_zerofier_codeword,
//                 &mut proover_proof_stream,
//             );

//             // verify
//             let mut verifier_proof_stream = StandardProofStream::new();
//             let verdict = stark.verify(
//                 &proof,
//                 air.clone(),
//                 &boundary,
//                 transition_zerofier_root,
//                 &mut verifier_proof_stream,
//             );

//             println!("{:?}", verdict.unwrap());
//             assert!(verdict.is_ok(), "valid stark proof fails to verify");
//             println!("success \\o/");

//             println!("verifying false claim ...");
//             // verify false claim
//             let output_element_ = output_element + BaseElement::one();
//             let boundary_ = rp.boundary_constraints(output_element_);
//             let mut verifier_proof_stream = StandardProofStream::new();
//             let verdict = stark.verify(
//                 &proof,
//                 air.clone(),
//                 &boundary_,
//                 transition_zerofier_root,
//                 &mut verifier_proof_stream,
//             );

//             assert!(verdict.is_err(), "invalid stark proof verifies");
//             println!("proof rejected! \\o/");

//             if trial == 19 {
//                 // verify with false witness
//                 print!("attempting to prove with false witness (should fail)
// ...");                 let mut rng = rand::thread_rng();
//                 let cycle = rng.gen::<usize>() % trace.len();
//                 let register = rng.gen::<usize>() % state_width;
//                 let error = BaseElement::from(rng.gen::<u64>());

//                 trace[cycle][register] += error;

//                 let mut proover_proof_stream = StandardProofStream::new();
//                 // should fail
//                 stark.prove(
//                     trace,
//                     air.clone(),
//                     &boundary,
//                     &transition_zerofier,
//                     &transition_zerofier_codeword,
//                     &mut proover_proof_stream,
//                 );
//             }
//         }
//     }
// }
