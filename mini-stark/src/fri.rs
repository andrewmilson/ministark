use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use crate::{polynomial::*, FieldElement, MerkleTree, ProofObject, ProofStream};

pub struct Fri<E> {
    offset: E,
    omega: E,
    pub domain_length: usize,
    expansion_factor: usize,
    num_colinearity_tests: usize,
}

impl<E: FieldElement> Fri<E> {
    pub fn new(
        offset: E,
        omega: E,
        domain_length: usize,
        expansion_factor: usize,
        num_colinearity_tests: usize,
    ) -> Self {
        Fri {
            offset,
            omega,
            domain_length,
            expansion_factor,
            num_colinearity_tests,
        }
    }

    fn num_rounds(&self) -> usize {
        let mut codeword_length = self.domain_length;
        let mut num_rounds = 0;
        while codeword_length > self.expansion_factor
            && 4 * self.num_colinearity_tests < codeword_length
        {
            codeword_length /= 2;
            num_rounds += 1;
        }
        num_rounds
    }

    pub fn eval_domain(&self) -> Vec<E> {
        (0..self.domain_length)
            .map(|i| self.offset * self.omega.pow((i as u128).into()))
            .collect()
    }

    pub fn prove<T: ProofStream<E>>(&self, codeword: Vec<E>, proof_stream: &mut T) -> Vec<usize> {
        assert_eq!(
            self.domain_length,
            codeword.len(),
            "domain length does not match length of initial codeword"
        );

        // commit phase
        let codewords = self.commit(codeword, proof_stream);

        // get indices
        let top_level_indices = self.sample_indices(
            proof_stream.prover_fiat_shamir(),
            codewords[1].len(),
            codewords.last().unwrap().len(),
            self.num_colinearity_tests,
        );
        let mut indices = top_level_indices.clone();

        for i in 0..(codewords.len() - 1) {
            // fold
            indices = indices
                .into_iter()
                .map(|index| index % (codewords[i].len() / 2))
                .collect();
            self.query(&codewords[i], &codewords[i + 1], &indices, proof_stream);
        }

        top_level_indices
    }

    fn sample_indices(
        &self,
        seed: u64,
        size: usize,
        reduced_size: usize,
        number: usize,
    ) -> Vec<usize> {
        assert!(number <= reduced_size, "cannot sample more indices than available in last codeword; requested: {}, available: {}", number, reduced_size);
        assert!(
            number <= 2 * reduced_size,
            "not enough entropy in indices wrt last codeword"
        );

        let mut indices = vec![];
        let mut reduced_indices = vec![];
        let mut counter = 0;

        while indices.len() < number {
            let mut hash = DefaultHasher::new();
            (seed + counter).hash(&mut hash);
            let index = hash.finish() as usize % size;
            let reduced_index = index % reduced_size;
            counter += 1;
            if !reduced_indices.contains(&reduced_index) {
                indices.push(index);
                reduced_indices.push(reduced_index);
            }
        }

        indices
    }

    fn commit<T: ProofStream<E>>(&self, codeword: Vec<E>, proof_stream: &mut T) -> Vec<Vec<E>> {
        let one = E::one();
        // TODO: fix for galois fields
        let two = one.double();
        let mut omega = self.omega;
        let mut offset = self.offset;
        let mut codeword = codeword;
        let mut codewords = vec![];

        for r in 0..self.num_rounds() {
            // make sure omega has the right order
            assert!(
                omega.pow((codeword.len() as u128 - 1).into()) == omega.inverse().unwrap(),
                "error in commit: omega does not have the right order!",
            );

            let root = MerkleTree::commit(&codeword);
            proof_stream.push(ProofObject::MerkleRoot(root));

            if r == self.num_rounds() - 1 {
                break;
            }

            // get challenge
            let alpha = E::from(proof_stream.prover_fiat_shamir());

            // split and fold
            let mut new_codeword = vec![];
            for i in 0..(codeword.len() / 2) {
                new_codeword.push(
                    two.inverse().unwrap()
                        * ((one + alpha / (offset * (omega.pow((i as u128).into()))))
                            * codeword[i]
                            + (one - alpha / (offset * (omega.pow((i as u128).into()))))
                                * codeword[codeword.len() / 2 + i]),
                );
            }

            // collect codeword
            codewords.push(codeword);
            codeword = new_codeword;

            omega.square_in_place();
            offset.square_in_place();
        }

        // send last codeword
        proof_stream.push(ProofObject::Codeword(codeword.clone()));

        // collect last codeword
        codewords.push(codeword);

        codewords
    }

    fn query<T: ProofStream<E>>(
        &self,
        current_codeword: &[E],
        next_codeword: &[E],
        c_indices: &[usize],
        proof_stream: &mut T,
    ) {
        let a_indices = c_indices;
        let b_indices: Vec<usize> = c_indices
            .iter()
            .map(|index| index + current_codeword.len() / 2)
            .collect();

        for s in 0..self.num_colinearity_tests {
            proof_stream.push(ProofObject::YValues((
                current_codeword[a_indices[s]],
                current_codeword[b_indices[s]],
                next_codeword[c_indices[s]],
            )));
        }

        for s in 0..self.num_colinearity_tests {
            proof_stream.push(ProofObject::MerklePath(MerkleTree::open(
                a_indices[s],
                current_codeword,
            )));
            proof_stream.push(ProofObject::MerklePath(MerkleTree::open(
                b_indices[s],
                current_codeword,
            )));
            proof_stream.push(ProofObject::MerklePath(MerkleTree::open(
                c_indices[s],
                next_codeword,
            )));
        }
    }

    pub fn verify<T: ProofStream<E>>(
        &self,
        proof_stream: &mut T,
        polynomial_values: &mut Vec<(usize, E)>,
    ) -> Result<(), &str> {
        let mut omega = self.omega;
        let mut offset = self.offset;

        // extract all roots and alphas
        let mut roots = vec![];
        let mut alphas = vec![];

        for _ in 0..self.num_rounds() {
            match proof_stream.pull() {
                ProofObject::MerkleRoot(root) => {
                    roots.push(root);
                    alphas.push(E::from(proof_stream.verifier_fiat_shamir()));
                }
                _ => return Err("Didn't receive expected number of merkle roots"),
            }
        }

        // extract last codeword
        let last_codeword = match proof_stream.pull() {
            ProofObject::Codeword(codeword) => codeword,
            _ => return Err("Expected to receive the last codeword"),
        };

        // check if it matches the given root
        if *roots.last().unwrap() != MerkleTree::commit(&last_codeword) {
            return Err("Last codeword is not well formed");
        }

        // check if it's low degree
        let mut last_omega = omega;
        let mut last_offset = offset;
        for _ in 0..(self.num_rounds() - 1) {
            last_omega.square_in_place();
            last_offset.square_in_place();
        }

        // assert that last_omega has the right order
        assert!(
            last_omega.inverse().unwrap()
                == last_omega.pow((last_codeword.len() as u32 - 1).into()),
            "Omega does not have right order"
        );

        // Slow version
        // ============
        // compute interpolant
        let last_domain = (0..last_codeword.len())
            .map(|i| last_offset * (last_omega.pow((i as u32).into())))
            .collect::<Vec<E>>();
        let poly = Polynomial::interpolate(&last_domain, &last_codeword);

        // // Fast version
        // // ============
        // let last_domain = (0..last_codeword.len())
        //     .map(|i| last_offset * (last_omega.pow((i as u32).into())))
        //     .collect::<Vec<E>>();
        // let coefficients = inverse_number_theory_transfer(last_omega, &last_codeword);
        // let poly = Polynomial::new(coefficients).scale(last_offset.inverse());
        // println!("Fri 35");

        assert_eq!(
            poly.evaluate_domain(&last_domain),
            last_codeword,
            "re-evaluated codeword does not match original!"
        );

        let expected_degree = last_codeword.len() / self.expansion_factor - 1;
        if poly.degree() > expected_degree.try_into().unwrap() {
            return Err("last codeword does not correspond to polynomial of low enough degree");
        }

        // get indices
        let top_level_indices = self.sample_indices(
            proof_stream.verifier_fiat_shamir(),
            self.domain_length >> 1,
            self.domain_length >> (self.num_rounds() - 1),
            self.num_colinearity_tests,
        );

        // for every round, check consistency of subsequent layers
        for r in 0..(self.num_rounds() - 1) {
            // fold c indices
            let c_indices: Vec<usize> = top_level_indices
                .iter()
                .map(|index| index % (self.domain_length >> (r + 1)))
                .collect();

            // infer a and b indices
            let a_indices = c_indices.clone();
            let b_indices: Vec<usize> = c_indices
                .iter()
                .map(|index| index + (self.domain_length >> (r + 1)))
                .collect();

            // read values and check colinearity
            let mut aa = vec![];
            let mut bb = vec![];
            let mut cc = vec![];

            for s in 0..self.num_colinearity_tests {
                match proof_stream.pull() {
                    ProofObject::YValues((ay, by, cy)) => {
                        aa.push(ay);
                        bb.push(by);
                        cc.push(cy);

                        // record top-layer values for later verification
                        if r == 0 {
                            polynomial_values.push((a_indices[s], ay));
                            polynomial_values.push((b_indices[s], by));
                        }

                        // colinearity check
                        let ax = offset * omega.pow((a_indices[s] as u32).into());
                        let bx = offset * omega.pow((b_indices[s] as u32).into());
                        let cx = alphas[r];

                        if !Polynomial::test_colinearity(vec![(ax, ay), (bx, by), (cx, cy)]) {
                            return Err("colinearity check failure");
                        }
                    }
                    _ => return Err("Expected to y-values"),
                }
            }

            // verify authentication paths
            for i in 0..self.num_colinearity_tests {
                match proof_stream.pull() {
                    ProofObject::MerklePath(path) => {
                        if !MerkleTree::verify(roots[r], a_indices[i], &path, aa[i]) {
                            return Err("merkle authentication path verification fails for aa");
                        }
                    }
                    _ => return Err("Expected merkle path for aa"),
                }

                match proof_stream.pull() {
                    ProofObject::MerklePath(path) => {
                        if !MerkleTree::verify(roots[r], b_indices[i], &path, bb[i]) {
                            return Err("merkle authentication path verification fails for bb");
                        }
                    }
                    _ => return Err("Expected merkle path for bb"),
                }

                match proof_stream.pull() {
                    ProofObject::MerklePath(path) => {
                        if !MerkleTree::verify(roots[r + 1], c_indices[i], &path, cc[i]) {
                            return Err("merkle authentication path verification fails for cc");
                        }
                    }
                    _ => return Err("Expected merkle path for cc"),
                }
            }

            omega.square_in_place();
            offset.square_in_place();
        }

        Ok(())
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use rand::Rng;

//     #[test]
//     fn test_merkle() {
//         let n = 64;
//         let mut rng = rand::thread_rng();
//         let leafs = (0..n).map(|_| rng.gen()).collect::<Vec<u128>>();
//         let root = MerkleTree::commit(&leafs);

//         // opening any leaf should work
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             assert!(MerkleTree::verify(root, i, &path, leafs[i]));
//         }

//         // opening non-leafs should not work
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             assert!(!MerkleTree::verify(root, i, &path, u128::MAX));
//         }

//         // opening wrong leafs should not work
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             let j = (i + 1 + rng.gen_range(0..(n - 1))) % n;
//             assert!(!MerkleTree::verify(root, i, &path, leafs[j]));
//         }

//         // opening leafs with the wrong index should not work
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             let j = (i + 1 + rng.gen_range(0..(n - 1))) % n;
//             assert!(!MerkleTree::verify(root, j, &path, leafs[i]))
//         }

//         // opening leafs to a false root should not work
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             assert!(!MerkleTree::verify(rng.gen(), i, &path, leafs[i]));
//         }

//         // opening leafs with even one falsehood in the path should not work
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             for j in 0..path.len() {
//                 let mut fake_path = path.clone();
//                 fake_path[j] = rng.gen();
//                 assert!(!MerkleTree::verify(root, i, &fake_path, leafs[i]));
//             }
//         }

//         // opening leafs to a different root should not work
//         let fake_root = MerkleTree::commit(&(0..n).map(|_| rng.gen::<u128>()).collect());
//         for i in 0..n {
//             let path = MerkleTree::open(i, &leafs);
//             assert!(!MerkleTree::verify(fake_root, i, &path, leafs[i]));
//         }
//     }

//     #[test]
//     fn test_fri() {
//         let field = Field::main();
//         let degree = 63;
//         let expansion_factor: usize = 4;
//         let num_colinearity_tests = 17;

//         let initial_codeword_length = (degree + 1) * expansion_factor;
//         let mut log_codeword_length = 0;
//         let mut codeword_length = initial_codeword_length;
//         while codeword_length > 1 {
//             codeword_length /= 2;
//             log_codeword_length += 1;
//         }

//         assert_eq!(
//             1 << log_codeword_length,
//             initial_codeword_length,
//             "log not computed correctly"
//         );

//         let omega = field.primitive_nth_root(initial_codeword_length as u128);
//         let generator = field.generator();

//         assert!(
//             omega ^ (1 << log_codeword_length) == field.one(),
//             "omega not nth root of unity"
//         );
//         assert!(
//             omega ^ (1 << (log_codeword_length - 1)) != field.one(),
//             "omega not primitive"
//         );

//         let fri = Fri::new(
//             generator,
//             omega,
//             initial_codeword_length,
//             expansion_factor,
//             num_colinearity_tests,
//         );

//         let polynomial = Polynomial::new(
//             (0..degree + 1)
//                 .map(|i| FieldElement::new(i as u128, field))
//                 .collect(),
//         );
//         let domain = (0..initial_codeword_length)
//             .map(|i| omega ^ i as u128)
//             .collect();

//         let mut codeword = polynomial.evaluate_domain(&domain);

//         // test valid codeword
//         println!("testing valid codeword ...");
//         let mut proof_stream = StandardProofStream::new();

//         fri.prove(codeword.clone(), &mut proof_stream);
//         println!("");
//         let mut points = vec![];
//         let verdict = fri.verify(&mut proof_stream, &mut points);

//         assert!(
//             verdict.is_ok(),
//             "rejecting proof, but proof should be valid!"
//         );

//         for (x, y) in points {
//             assert_eq!(
//                 polynomial.evaluate(omega ^ x as u128),
//                 y,
//                 "polynomial evaluates to wrong value",
//             );
//         }
//         println!("success! \\o/");

//         // disturb then test for failure
//         println!("testing invalid codeword ...");
//         let mut proof_stream = StandardProofStream::new();

//         for i in 0..(degree / 3) {
//             codeword[i] = field.zero();
//         }

//         fri.prove(codeword, &mut proof_stream);
//         let mut points = vec![];
//         assert!(
//             fri.verify(&mut proof_stream, &mut points).is_err(),
//             "proof should fail, but is accepted ..."
//         );
//         println!("success! \\o/");
//     }
// }
