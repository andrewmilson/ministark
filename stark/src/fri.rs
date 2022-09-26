use crate::ceil_power_of_two;
use crate::merkle::Merkle;
use crate::protocol::ProofObject;
use crate::protocol::ProofStream;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Evaluations;
use ark_poly::GeneralEvaluationDomain;
use ark_poly::Polynomial;
use brainfuck::InputTable;
use brainfuck::Table;
use num_traits::One;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;
use std::iter::zip;

pub trait Config {
    /// Base prime field
    type Fp: PrimeField + FftField;
    /// Extension field element
    type Fx: FftField + Field<BasePrimeField = Self::Fp>;

    const EXPANSION_FACTOR: usize;
    const SECURITY_LEVEL: usize;
    const NUM_COLINEARITY_CHECKS: usize =
        Self::SECURITY_LEVEL / Self::EXPANSION_FACTOR.ilog2() as usize;
}

pub struct Fri<P: Config> {
    _params: P,
}

impl<P: Config> Fri<P> {
    pub fn new(_params: P) -> Self {
        Fri { _params }
    }

    fn sample_indices(
        &self,
        n: usize,
        reduced_size: usize,
        max: usize,
        randomness: u64,
    ) -> Vec<usize> {
        assert!(n <= reduced_size);
        let mut indices = Vec::new();
        let mut reduced_indices = vec![false; reduced_size];
        let mut counter = 0;
        while indices.len() < n {
            let mut hasher = DefaultHasher::new();
            randomness.hash(&mut hasher);
            counter.hash(&mut hasher);
            let hash = hasher.finish();
            let index = hash as usize % max;
            let reduced_index = index % reduced_size;
            if !reduced_indices[reduced_index] {
                indices.push(index);
                reduced_indices[reduced_index] = true;
            }
            counter += 1;
        }
        indices
    }

    pub fn commit(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        codeword: &[P::Fx],
    ) -> (Vec<Vec<P::Fx>>, Vec<Merkle<P::Fx>>)
    where
        [(); InputTable::<P::Fx>::BASE_WIDTH]: Sized,
    {
        let one = P::Fx::one();
        assert_ne!(P::Fx::characteristic(), [2]);
        let two = one + one;

        let mut codeword = codeword.to_vec();
        let mut omega = P::Fp::get_root_of_unity(codeword.len() as u64).unwrap();
        let mut offset = P::Fp::GENERATOR;

        let mut codewords = Vec::new();
        let mut trees = Vec::new();

        while codeword.len() >= ceil_power_of_two(P::NUM_COLINEARITY_CHECKS)
            && codeword.len() >= P::EXPANSION_FACTOR
        {
            let tree = Merkle::new(&codeword);
            let root = tree.root();

            // Skip the first round
            if !trees.is_empty() {
                // TODO: HELP: is this needed on the last round?
                proof_stream.push(crate::protocol::ProofObject::MerkleRoot(root));
            }

            // only prepare next round if necessary
            if codeword.len() == ceil_power_of_two(P::NUM_COLINEARITY_CHECKS)
                || codeword.len() == P::EXPANSION_FACTOR
            {
                break;
            }

            trees.push(tree);
            codewords.push(codeword.clone());

            // get challenge for split and fold
            let alpha = P::Fx::from(proof_stream.prover_fiat_shamir());
            println!("Alpha: {}", alpha);
            let (lhs, rhs) = codeword.split_at(codeword.len() / 2);
            codeword = zip(lhs, rhs)
                .enumerate()
                .map(|(i, (&l, &r))| {
                    ((one + alpha / P::Fx::from_base_prime_field(offset * omega.pow([i as u64])))
                        * l
                        + (one
                            - alpha / P::Fx::from_base_prime_field(offset * omega.pow([i as u64])))
                            * r)
                        / two
                })
                .collect();

            // println!("NEW CODEWORD DEGREE");
            // determine_codeword_degree(&codeword);

            omega.square_in_place();
            offset.square_in_place();
        }

        // send last codeword
        proof_stream.push(crate::protocol::ProofObject::Codeword(codeword.clone()));
        codewords.push(codeword);

        (codewords, trees)
    }

    pub fn query(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        curr_tree: &Merkle<P::Fx>,
        next_tree: &Merkle<P::Fx>,
        indices: &[usize],
    ) where
        [(); InputTable::<P::Fx>::BASE_WIDTH]: Sized,
    {
        let lhs_indices = indices.to_vec();
        let rhs_indices = indices
            .iter()
            .map(|i| i + curr_tree.leafs.len() / 2)
            .collect::<Vec<usize>>();

        // reveal leafs
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            proof_stream.push(crate::protocol::ProofObject::FriLeafs((
                curr_tree.leafs[lhs_indices[i]],
                curr_tree.leafs[rhs_indices[i]],
                next_tree.leafs[lhs_indices[i]],
            )));
        }

        // reveal authentication paths
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            let mp = crate::protocol::ProofObject::MerklePath(curr_tree.open(lhs_indices[i]).1);
            proof_stream.push(mp);
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(rhs_indices[i]).1,
            ));
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                next_tree.open(lhs_indices[i]).1,
            ));
        }
    }

    pub fn query_last(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        curr_tree: &Merkle<P::Fx>,
        last_codeword: &[P::Fx],
        indices: &[usize],
    ) where
        [(); InputTable::<P::Fx>::BASE_WIDTH]: Sized,
    {
        let lhs_indices = indices.to_vec();
        let rhs_indices = indices
            .iter()
            .map(|i| i + curr_tree.leafs.len() / 2)
            .collect::<Vec<usize>>();

        // reveal leafs
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            proof_stream.push(crate::protocol::ProofObject::FriLeafs((
                curr_tree.leafs[lhs_indices[i]],
                curr_tree.leafs[rhs_indices[i]],
                last_codeword[lhs_indices[i]],
            )));
        }

        // reveal authentication paths
        for i in 0..P::NUM_COLINEARITY_CHECKS {
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(lhs_indices[i]).1,
            ));
            proof_stream.push(crate::protocol::ProofObject::MerklePath(
                curr_tree.open(rhs_indices[i]).1,
            ));
        }
    }

    pub fn prove(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        codeword: &[P::Fx],
    ) -> Vec<usize>
    where
        [(); InputTable::<P::Fx>::BASE_WIDTH]: Sized,
    {
        // commit phase
        let (codewords, trees) = self.commit(proof_stream, codeword);

        // query phase
        let last_codeword = codewords.last().unwrap();
        println!("Codewords: {}", codewords.len());
        println!("Last codeword len: {}", last_codeword.len());
        let top_level_indices = self.sample_indices(
            P::NUM_COLINEARITY_CHECKS,
            ceil_power_of_two(P::NUM_COLINEARITY_CHECKS),
            codeword.len() / 2,
            proof_stream.prover_fiat_shamir(),
        );
        println!("Hello!");
        for i in 0..trees.len() - 1 {
            let indices = top_level_indices
                .iter()
                .map(|index| index % (codewords[i].len() / 2))
                .collect::<Vec<usize>>();
            self.query(proof_stream, &trees[i], &trees[i + 1], &indices);
        }
        let indices = top_level_indices
            .iter()
            .map(|index| index % last_codeword.len())
            .collect::<Vec<usize>>();
        self.query_last(proof_stream, trees.last().unwrap(), last_codeword, &indices);

        top_level_indices
    }

    pub fn verify(
        &self,
        proof_stream: &mut impl ProofStream<P::Fx>,
        codeword_len: usize,
        combination_root: u64,
    ) -> Result<(), &str>
    where
        [(); InputTable::<P::Fx>::BASE_WIDTH]: Sized,
    {
        let mut offset = P::Fp::GENERATOR;
        let mut omega = P::Fp::get_root_of_unity(codeword_len as u64).unwrap();

        // extract all roots and alphas
        let mut roots = vec![combination_root];
        let mut alphas = Vec::new();
        let mut round_len = codeword_len;
        while round_len >= ceil_power_of_two(P::NUM_COLINEARITY_CHECKS)
            && round_len >= P::EXPANSION_FACTOR
        {
            if round_len != codeword_len {
                roots.push(match proof_stream.pull() {
                    ProofObject::MerkleRoot(root) => root,
                    _ => return Err("Expected root"),
                })
            }
            let alpha = P::Fx::from(proof_stream.verifier_fiat_shamir());
            alphas.push(alpha);
            round_len /= 2;
        }

        // extract the last codeword
        let last_codeword = match proof_stream.pull() {
            ProofObject::Codeword(codeword) => codeword,
            _ => return Err("Expected last codeword"),
        };
        let last_root = roots.last().unwrap();

        // check if it matches the given root
        // TODO: why check? no point
        assert_eq!(
            *last_root,
            Merkle::new(&last_codeword).root(),
            "last codeword is not well formed"
        );

        // check if the last codeword is low degree
        let degree = last_codeword.len() / P::EXPANSION_FACTOR;
        let last_coset_domain =
            GeneralEvaluationDomain::new_coset(last_codeword.len(), P::Fx::GENERATOR).unwrap();
        let last_evals = Evaluations::from_vec_and_domain(last_codeword.clone(), last_coset_domain);
        let poly = last_evals.clone().interpolate();

        // let poly = legacy_algebra::scale_poly(
        //     &DensePolynomial::from_coefficients_vec(inverse_number_theory_transform(
        //         &last_codeword,
        //     )),
        //     P::Fx::from_base_prime_field(offset.inverse().unwrap()),
        // );

        assert!(poly.degree() <= degree);

        // verify by evaluating
        assert_eq!(
            last_evals.evals,
            poly.evaluate_over_domain(last_coset_domain).evals
        );

        // get indices
        let top_level_indices = self.sample_indices(
            P::NUM_COLINEARITY_CHECKS,
            ceil_power_of_two(P::NUM_COLINEARITY_CHECKS),
            codeword_len / 2,
            proof_stream.verifier_fiat_shamir(),
        );

        // let codeword_from = last_codeword;

        let mut n = codeword_len / 2;
        let mut round = 0;
        let mut omega = P::Fp::get_root_of_unity(codeword_len as u64).unwrap();
        let mut offset = P::Fp::GENERATOR;

        // for every pair of consecutive rounds check consistency of subsequent layers
        while n >= ceil_power_of_two(P::NUM_COLINEARITY_CHECKS) && n >= P::EXPANSION_FACTOR {
            let is_last_round =
                n == ceil_power_of_two(P::NUM_COLINEARITY_CHECKS) || n == P::EXPANSION_FACTOR;

            // TODO: rename from c
            // fold c indices
            let c_indices = top_level_indices
                .iter()
                .map(|index| index % n)
                .collect::<Vec<usize>>();

            // infer a and b indices
            let a_indices = c_indices.clone();
            let b_indices = c_indices
                .iter()
                .map(|index| index + n)
                .collect::<Vec<usize>>();

            // read values and check colinearity
            let mut aa = Vec::new();
            let mut bb = Vec::new();
            let mut cc = Vec::new();

            for s in 0..P::NUM_COLINEARITY_CHECKS {
                // println!("Check: {s}");
                let (ay, by, cy) = match proof_stream.pull() {
                    ProofObject::FriLeafs(v) => v,
                    a => return Err("expected FRI y values"),
                };

                aa.push(ay);
                bb.push(by);
                cc.push(cy);

                // colinearity check
                let ax = P::Fx::from_base_prime_field(offset * omega.pow([a_indices[s] as u64]));
                let bx = P::Fx::from_base_prime_field(offset * omega.pow([b_indices[s] as u64]));
                let cx = alphas[round];

                println!(
                    "Degree is: {}",
                    interpolate(&[ax, bx, cx], &[ay, by, cy]).degree()
                );
                // if interpolate(&[ax, bx, cx], &[ay, by, cy]).degree() != 1 {
                //     return Err("colinearity check failed");
                // }
            }

            // verify authentication paths
            for i in 0..P::NUM_COLINEARITY_CHECKS {
                let path = match proof_stream.pull() {
                    ProofObject::MerklePath(path) => path,
                    _ => return Err("Expected merkle path for aa"),
                };
                if !Merkle::verify(roots[round], a_indices[i], &path, &aa[i]) {
                    return Err("Failed to verify merkle path for aa");
                }
                let path = match proof_stream.pull() {
                    ProofObject::MerklePath(path) => path,
                    _ => return Err("Expected merkle path for bb"),
                };
                if !Merkle::verify(roots[round], b_indices[i], &path, &bb[i]) {
                    return Err("Failed to verify merkle path for bb");
                }
                if is_last_round {
                    if cc[i] != last_codeword[c_indices[i]] {
                        return Err("Failed to verify last codeword for cc");
                    }
                } else {
                    let path = match proof_stream.pull() {
                        ProofObject::MerklePath(path) => path,
                        _ => return Err("Expected merkle path for cc"),
                    };
                    if !Merkle::verify(roots[round + 1], b_indices[i], &path, &cc[i]) {
                        return Err("Failed to verify merkle path for cc");
                    }
                }
            }

            round += 1;
            n /= 2;
            omega.square_in_place();
            offset.square_in_place();
        }

        Ok(())
    }
}

pub fn interpolate<E: Field>(domain: &[E], values: &[E]) -> DensePolynomial<E> {
    assert_eq!(
        domain.len(),
        values.len(),
        "number of elements in domain does not match number of values -- cannot interpolate"
    );
    // Generate master numerator polynomial: (x - domain[0]) * (x - domain[1]) *
    // ....
    let root = zerofier_domain(domain);

    // Generate the numerator for each item in the domain
    let numerators = domain
        .iter()
        .copied()
        // root / (x - domain[i])
        .map(|d| &root / &DensePolynomial::from_coefficients_vec(vec![-d, E::one()]))
        .collect::<Vec<DensePolynomial<E>>>();

    // Generate denominators by evaluating numerator polys at each x
    let mut inverse_denominators = numerators
        .iter()
        .zip(domain)
        .map(|(numerator, d)| numerator.evaluate(d))
        .collect::<Vec<E>>();
    ark_ff::batch_inversion(&mut inverse_denominators);

    // Generate output polynomial
    let mut output_coefficients = vec![E::zero(); values.len()];

    for ((y, numerator), inverse_denominator) in values
        .iter()
        .copied()
        .zip(numerators)
        .zip(inverse_denominators)
    {
        let y_slice = y * inverse_denominator;
        for (j, coefficient) in numerator.coeffs.into_iter().enumerate() {
            output_coefficients[j] += coefficient * y_slice;
        }
    }

    DensePolynomial::from_coefficients_vec(output_coefficients)
}

fn zerofier_domain<E: Field>(domain: &[E]) -> DensePolynomial<E> {
    let x = DensePolynomial::from_coefficients_vec(vec![E::zero(), E::one()]);
    let mut accumulator = DensePolynomial::from_coefficients_vec(vec![E::one()]);
    for element in domain.iter() {
        let subtraction = &x - &DensePolynomial::from_coefficients_vec(vec![*element]);
        accumulator = accumulator.naive_mul(&subtraction);
    }
    accumulator
}

// fn determine_codeword_degree<F>(evaluations: &[F])
// where
//     F: Field,
//     F::BasePrimeField: FftField,
// {
//     let poly = DensePolynomial::from_coefficients_vec(
//         legacy_algebra::number_theory_transform::inverse_number_theory_transform(evaluations),
//     );
//     println!("Degree is: {}", poly.degree());
// }

// Alpha: 15280113618990651205
// Alpha: 11196160987341373762
// Alpha: 13946831793190648870
// Alpha: 2876196548255908070
// Alpha: 6134970219535787624
// Alpha: 13842569450404052519
// Alpha: 8853488135211430262
// Alpha: 11952022665098211853
// Alpha: 1191585754718550683
// Alpha: 11093522932086459321
// Alpha: 16409664712104079389
