use crate::hash::Digest;
use crate::merkle;
use crate::merkle::MatrixMerkleTree;
use crate::merkle::MerkleTree;
use crate::random::PublicCoin;
use crate::utils::GpuAllocator;
use crate::utils::GpuVec;
use crate::Matrix;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use ark_poly::domain::DomainCoeff;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ministark_gpu::prelude::*;
use ministark_gpu::utils::bit_reverse;
use ministark_gpu::utils::bit_reverse_index;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use snafu::Snafu;
use std::iter::zip;
use std::marker::PhantomData;

#[derive(Clone, Copy)]
pub struct FriOptions {
    folding_factor: usize,
    max_remainder_coeffs: usize,
    blowup_factor: usize,
}

impl FriOptions {
    pub const fn new(
        blowup_factor: usize,
        folding_factor: usize,
        max_remainder_coeffs: usize,
    ) -> Self {
        Self {
            folding_factor,
            max_remainder_coeffs,
            blowup_factor,
        }
    }

    pub const fn num_layers(&self, mut domain_size: usize) -> usize {
        let mut num_layers = 0;
        while domain_size > self.max_remainder_coeffs * self.blowup_factor {
            domain_size /= self.folding_factor;
            num_layers += 1;
        }
        num_layers
    }

    pub const fn remainder_size(&self, mut domain_size: usize) -> usize {
        while domain_size > self.max_remainder_coeffs * self.blowup_factor {
            domain_size /= self.folding_factor;
        }
        domain_size
    }

    pub const fn domain_offset<F: GpuField>(&self) -> F::FftField
    where
        F::FftField: FftField,
    {
        // TODO: Can specify a domain offset in the air which might cause issues with
        // having this a constant
        F::FftField::GENERATOR
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct FriProof<F: Field, D: Digest, M: MatrixMerkleTree<F>> {
    pub layers: Vec<LayerProof<F, D, M>>,
    pub remainder_coeffs: Vec<F>,
}

impl<F: GpuField + Field, D: Digest, M: MatrixMerkleTree<F, Root = D>> FriProof<F, D, M>
where
    F::FftField: FftField,
{
    pub fn new(layers: Vec<LayerProof<F, D, M>>, remainder_coeffs: Vec<F>) -> Self {
        Self {
            layers,
            remainder_coeffs,
        }
    }
}

struct FriLayer<F: GpuField, M: MerkleTree> {
    merkle_tree: M,
    evaluations: Matrix<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct LayerProof<F: Field, D: Digest, M: MatrixMerkleTree<F>> {
    pub flattenend_rows: Vec<F>,
    pub merkle_proof: M::Proof,
    pub commitment: D,
}

impl<F: GpuField + Field, D: Digest, M: MatrixMerkleTree<F, Root = D>> LayerProof<F, D, M>
where
    F::FftField: FftField,
{
    pub fn new<const N: usize>(rows: Vec<[F; N]>, merkle_proof: M::Proof, commitment: D) -> Self {
        Self {
            flattenend_rows: rows.into_iter().flatten().collect(),
            merkle_proof,
            commitment,
        }
    }

    pub fn verify<const N: usize>(self, positions: &[usize]) -> Result<(), merkle::Error> {
        let commitment = &self.commitment;
        let (rows, _remainder) = &self.flattenend_rows.as_chunks::<N>();
        M::verify_rows(commitment, positions, rows, self.merkle_proof)
    }
}

pub struct FriProver<F: GpuField, D: Digest, M: MerkleTree> {
    options: FriOptions,
    layers: Vec<FriLayer<F, M>>,
    remainder_coeffs: Vec<F>,
    _phantom: PhantomData<D>,
}

impl<
        F: GpuField + Field + DomainCoeff<F::FftField>,
        D: Digest,
        M: MatrixMerkleTree<F, Root = D>,
    > FriProver<F, D, M>
where
    F::FftField: FftField,
{
    pub const fn new(options: FriOptions) -> Self {
        Self {
            options,
            layers: Vec::new(),
            remainder_coeffs: Vec::new(),
            _phantom: PhantomData,
        }
    }

    pub fn into_proof(self, positions: &[usize]) -> FriProof<F, D, M> {
        let folding_factor = self.options.folding_factor;
        // let (last_layer, initial_layers) = self.layers.split_last().unwrap();
        let mut proof_layers = Vec::new();
        let mut positions = positions.to_vec();
        for layer in &self.layers {
            positions = fold_positions(&positions, folding_factor);
            proof_layers.push(match folding_factor {
                2 => query_layer::<F, D, M, 2>(layer, &positions),
                4 => query_layer::<F, D, M, 4>(layer, &positions),
                6 => query_layer::<F, D, M, 6>(layer, &positions),
                8 => query_layer::<F, D, M, 8>(layer, &positions),
                16 => query_layer::<F, D, M, 16>(layer, &positions),
                _ => unimplemented!("folding factor {folding_factor} is not supported"),
            });
        }

        // // layers store interlaved evaluations so they need to be un-interleaved
        // let remainder_commitment = last_layer.tree.root().to_vec();
        // let last_evals = &last_layer.evaluations;
        // let mut remainder = vec![F::zero(); last_evals.len()];
        // let num_eval_chunks = last_evals.len() / folding_factor;
        // for i in 0..num_eval_chunks {
        //     for j in 0..folding_factor {
        //         remainder[i + num_eval_chunks * j] = last_evals[i * folding_factor +
        // j];     }
        // }

        FriProof::new(proof_layers, self.remainder_coeffs)
    }

    pub fn build_layers(
        &mut self,
        channel: &mut impl ProverChannel<Field = F, Digest = D>,
        mut evaluations: GpuVec<F>,
    ) {
        assert!(self.layers.is_empty());
        for _ in 0..self.options.num_layers(evaluations.len()) {
            evaluations = match self.options.folding_factor {
                2 => self.build_layer::<2>(channel, evaluations),
                4 => self.build_layer::<4>(channel, evaluations),
                8 => self.build_layer::<8>(channel, evaluations),
                16 => self.build_layer::<16>(channel, evaluations),
                folding_factor => unreachable!("folding factor {folding_factor} not supported"),
            }
        }
        self.set_remainder(channel, evaluations);
    }

    /// Builds a single layer of the FRI protocol
    /// Returns the evaluations for the next layer.
    fn build_layer<const N: usize>(
        &mut self,
        channel: &mut impl ProverChannel<Field = F, Digest = D>,
        evaluations: GpuVec<F>,
    ) -> GpuVec<F> {
        // Each layer requires decommitting to `folding_factor` many evaluations e.g.
        // `folding_factor = 2` decommits to an evaluation for LHS_i and RHS_i
        // (0 ≤ i < n/2) which requires two merkle paths if the evaluations are
        // committed to in their natural order. If we instead commit to bit-reversed
        // evaluations then LHS_i and RHS_i are next to each other and only require a
        // single merkle path for their decommitment.

        // TODO: update docs with bit reversed evals
        let (cosets, remainder) = evaluations.as_chunks::<N>();
        assert!(remainder.is_empty());

        let matrix = Matrix::from_arrays(cosets);
        let merkle_tree = M::from_matrix(&matrix);
        channel.commit_fri_layer(merkle_tree.root());

        self.layers.push(FriLayer {
            merkle_tree,
            evaluations: matrix,
        });

        // return the next evaluations
        apply_drp(
            evaluations,
            F::FftField::ONE,
            channel.draw_fri_alpha(),
            self.options.folding_factor,
        )
    }

    fn set_remainder(
        &mut self,
        channel: &mut impl ProverChannel<Field = F, Digest = D>,
        mut evaluations: GpuVec<F>,
    ) {
        let domain_size = evaluations.len();
        assert!(domain_size.is_power_of_two());
        assert!(domain_size <= self.options.max_remainder_coeffs * self.options.blowup_factor);
        let domain = Radix2EvaluationDomain::new(domain_size).unwrap();
        bit_reverse(&mut evaluations);
        let coeffs = domain.ifft(&evaluations);
        let max_degree = domain_size / self.options.blowup_factor - 1;
        let (remainder_coeffs, zero_coeffs) = coeffs.split_at(max_degree + 1);
        assert!(zero_coeffs.iter().all(F::is_zero));
        channel.commit_remainder(remainder_coeffs);
        self.remainder_coeffs = remainder_coeffs.to_vec();
    }
}

#[derive(Debug, Snafu)]
pub enum VerificationError {
    #[snafu(display("queries do not resolve to their commitment in layer {layer}"))]
    LayerCommitmentInvalid { layer: usize },
    #[snafu(display("degree respecting projection is invalid for layer {layer}"))]
    InvalidDegreeRespectingProjection { layer: usize },
    #[snafu(display("the number of query positions does not match the number of evaluations"))]
    NumPositionEvaluationMismatch,
    #[snafu(display("remainder is invalid"))]
    RemainderCommitmentInvalid,
    #[snafu(display("remainder is not a degree {degree} polynomial"))]
    RemainderDegreeMismatch { degree: usize },
    #[snafu(display("{size} can't be divided by {folding_factor} (layer {layer})"))]
    CodewordTruncation {
        size: usize,
        folding_factor: usize,
        layer: usize,
    },
}

/// Fri verifier adapted from Winterfell to match Starkware's verifier
/// <https://github.com/facebook/winterfell/blob/main/fri/src/verifier/mod.rs#L58>
pub struct FriVerifier<F: GpuField + Field, D: Digest, M: MatrixMerkleTree<F, Root = D>>
where
    F::FftField: FftField,
{
    options: FriOptions,
    layer_commitments: Vec<D>,
    pub layer_alphas: Vec<F>,
    proof: FriProof<F, D, M>,
    domain: Radix2EvaluationDomain<F::FftField>,
}

impl<
        F: GpuField + Field + DomainCoeff<F::FftField>,
        D: Digest,
        M: MatrixMerkleTree<F, Root = D>,
    > FriVerifier<F, D, M>
where
    F::FftField: FftField,
{
    pub fn new(
        public_coin: &mut impl PublicCoin<Field = F, Digest = D>,
        options: FriOptions,
        proof: FriProof<F, D, M>,
        max_poly_degree: usize,
    ) -> Result<Self, VerificationError> {
        let folding_factor = options.folding_factor;
        let domain_offset = options.domain_offset::<F>();
        let domain_size = max_poly_degree.next_power_of_two() * options.blowup_factor;
        let domain = Radix2EvaluationDomain::new_coset(domain_size, domain_offset).unwrap();

        let mut layer_alphas = Vec::new();
        let mut layer_commitments = Vec::new();
        let mut layer_codeword_len = domain_size;
        for (i, layer) in proof.layers.iter().enumerate() {
            // TODO: batch merkle tree proofs
            // get the merkle root from the first merkle path
            public_coin.reseed_with_digest(&layer.commitment);
            let alpha = public_coin.draw();
            layer_alphas.push(alpha);
            layer_commitments.push(layer.commitment.clone());

            if i != proof.layers.len() - 1 && layer_codeword_len % folding_factor != 0 {
                return Err(VerificationError::CodewordTruncation {
                    size: layer_codeword_len,
                    folding_factor,
                    layer: i,
                });
            }

            layer_codeword_len /= folding_factor;
        }

        public_coin.reseed_with_field_element_vector(&proof.remainder_coeffs);

        // TODO: add back in
        // let remainder_root =
        // Output::<D>::from_slice(&proof.remainder_commitment).clone();
        // let remainder_root =
        // Output::<D>::from_slice(&proof.remainder_commitment).clone();
        // public_coin.reseed_with_hash(&remainder_root);
        // let remainder_alpha = public_coin.draw();
        // layer_alphas.push(remainder_alpha);
        // layer_commitments.push(remainder_root);

        Ok(Self {
            options,
            layer_commitments,
            layer_alphas,
            proof,
            domain,
        })
    }

    pub fn verify_generic<const N: usize>(
        self,
        positions: &[usize],
        evaluations: &[F],
    ) -> Result<(), VerificationError> {
        let folding_domain = Radix2EvaluationDomain::new(N).unwrap();
        let mut layers = self.proof.layers.into_iter();
        let mut layer_alphas = self.layer_alphas.into_iter();
        let mut layer_commitments = self.layer_commitments.into_iter();
        let mut positions = positions.to_vec();
        let mut evaluations = evaluations.to_vec();
        let mut domain_size = self.domain.size();
        let mut domain_generator = self.domain.group_gen();

        // verify all layers except remainder
        for i in 0..self.options.num_layers(domain_size) {
            let folded_positions = fold_positions(&positions, N);
            let layer_alpha = layer_alphas.next().unwrap();
            let layer_commitment = layer_commitments.next().unwrap();

            // TODO: change assert to error. Check remainder
            let layer = layers.next().unwrap();
            let (rows, _) = &layer.flattenend_rows.as_chunks::<N>();
            assert_eq!(rows.len(), folded_positions.len());

            // verify the layer values against the layer's commitment
            M::verify_rows(
                &layer_commitment,
                &folded_positions,
                rows,
                layer.merkle_proof,
            )
            .map_err(|_| VerificationError::LayerCommitmentInvalid { layer: i })?;

            let query_values = get_query_values(rows, &positions, &folded_positions);
            // println!("evaluatinos: {:?}", evaluations);
            // println!("query vals: {:?}", query_values);
            assert_eq!(evaluations.len(), query_values.len());
            for (i, (eval, query_val)) in zip(&evaluations, &query_values).enumerate() {
                assert_eq!(eval, query_val, "mismatch at {i}");
            }

            if evaluations != query_values {
                return Err(VerificationError::InvalidDegreeRespectingProjection { layer: i });
            }

            let polys = rows
                .iter()
                .zip(&folded_positions)
                .map(move |(chunk, position)| {
                    let bit_rev_position = bit_reverse_index(domain_size / N, *position);
                    // let offset = domain_offset.pow([N.pow(i as u32) as u64])
                    //     * domain_generator.pow([bit_rev_position as u64]);
                    let offset = domain_generator.pow([bit_rev_position as u64]);
                    let domain = folding_domain.get_coset(offset).unwrap();
                    let mut chunk = *chunk;
                    bit_reverse(&mut chunk);
                    let mut coeffs = domain.ifft(&chunk);
                    for coeff in &mut coeffs {
                        *coeff *= F::from(N as u64);
                    }
                    DensePolynomial::from_coefficients_vec(coeffs)
                });

            // prepare for next layer
            evaluations = polys.map(|poly| poly.evaluate(&layer_alpha)).collect();
            positions = folded_positions;
            domain_generator = domain_generator.pow([N as u64]);
            domain_size /= N;
        }
        verify_remainder::<F>(
            self.proof.remainder_coeffs,
            &positions,
            &evaluations,
            domain_generator,
            domain_size,
            self.options.blowup_factor,
        )
    }

    pub fn verify(self, positions: &[usize], evaluations: &[F]) -> Result<(), VerificationError> {
        if positions.len() != evaluations.len() {
            return Err(VerificationError::NumPositionEvaluationMismatch);
        }

        match self.options.folding_factor {
            2 => self.verify_generic::<2>(positions, evaluations),
            4 => self.verify_generic::<4>(positions, evaluations),
            8 => self.verify_generic::<8>(positions, evaluations),
            16 => self.verify_generic::<16>(positions, evaluations),
            // TODO: move this to options
            folding_factor => unreachable!("folding factor {folding_factor} not supported"),
        }
    }
}

fn verify_remainder<F: GpuField + Field + DomainCoeff<F::FftField>>(
    remainder_coeffs: Vec<F>,
    positions: &[usize],
    expected_evaluations: &[F],
    domain_generator: F::FftField,
    domain_size: usize,
    blowup_factor: usize,
) -> Result<(), VerificationError>
where
    F::FftField: FftField,
{
    let remainder_poly = DensePolynomial::from_coefficients_vec(remainder_coeffs);
    let expected_degree = domain_size / blowup_factor - 1;
    if remainder_poly.degree() > expected_degree {
        return Err(VerificationError::RemainderDegreeMismatch {
            degree: expected_degree,
        });
    }
    let xs = positions
        .iter()
        .map(|p| domain_generator.pow([bit_reverse_index(domain_size, *p) as u64]));
    for (i, x) in xs.enumerate() {
        // TODO: fix types and just call horner_evaluate
        let y = remainder_poly
            .iter()
            .rfold(F::zero(), move |mut result, coeff| {
                result *= x;
                result + coeff
            });
        if expected_evaluations[i] != y {
            return Err(VerificationError::RemainderCommitmentInvalid);
        }
    }
    Ok(())
}

pub trait ProverChannel {
    type Digest: Digest;
    type Field: GpuField;

    fn commit_fri_layer(&mut self, layer_root: Self::Digest);

    fn commit_remainder(&mut self, remainder_coeffs: &[Self::Field]);

    fn draw_fri_alpha(&mut self) -> Self::Field;
}

/// Performs a degree respecting projection (drp) on polynomial evaluations.
// Example for `folding_factor = 2`:
// ```text
// 1. interpolate evals over the evaluation domain to obtain f(x):
//    ┌─────────┬───┬───┬───┬───┬───┬───┬───┬───┐
//    │ i       │ 0 │ 1 │ 2 │ 3 │ 4 │ 5 │ 6 │ 7 │
//    ├─────────┼───┼───┼───┼───┼───┼───┼───┼───┤
//    │ eval[i] │ 9 │ 2 │ 3 │ 5 │ 9 │ 2 │ 3 │ 5 │
//    └─────────┴───┴───┴───┴───┴───┴───┴───┴───┘
//    ┌──────┬───────┬───────┬───────┬───────┬───────┬───────┬───────┬───────┐
//    │ x    │ o*Ω^0 │ o*Ω^1 │ o*Ω^2 │ o*Ω^3 │ o*Ω^4 │ o*Ω^5 │ o*Ω^6 │ o*Ω^7 │
//    ├──────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
//    │ f(x) │ 9     │ 2     │ 3     │ 5     │ 9     │ 2     │ 3     │ 5     │
//    └──────┴───────┴───────┴───────┴───────┴───────┴───────┴───────┴───────┘
//      f(x) = c0 * x^0 + c1 * x^1 + c2 * x^2 + c3 * x^3 +
//             c4 * x^4 + c5 * x^5 + c6 * x^6 + c7 * x^7
//
// 2. perform a random linear combination of odd and even coefficients of f(x):
//    f_e(x) = c0 + c2 * x + c4 * x^2 + c6 * x^3
//    f_o(x) = c1 + c3 * x + c5 * x^2 + c7 * x^3
//    f(x)   = f_e(x) + x * f_o(x)
//    f'(x)  = f_e(x) + α * f_o(x)
//    α      = <random field element sent from verifier>
//
// 4. obtain the DRP by evaluating f'(x) over a new domain of half the size:
//    ┌───────┬───────────┬───────────┬───────────┬───────────┐
//    │ x     │ (o*Ω^0)^2 │ (o*Ω^1)^2 │ (o*Ω^2)^2 │ (o*Ω^3)^2 │
//    ├───────┼───────────┼───────────┼───────────┼───────────┤
//    │ f'(x) │ 82        │ 12        │ 57        │ 34        │
//    └───────┴───────────┴───────────┴───────────┴───────────┘
//    ┌────────┬────┬────┬────┬────┐
//    │ i      │ 0  │ 1  │ 2  │ 3  │
//    ├────────┼────┼────┼────┼────┤
//    │ drp[i] │ 82 │ 12 │ 57 │ 34 │
//    └────────┴────┴────┴────┴────┘
// ```
pub fn apply_drp<F: GpuField + Field + DomainCoeff<F::FftField>>(
    mut evals: GpuVec<F>,
    domain_offset: F::FftField,
    alpha: F,
    folding_factor: usize,
) -> GpuVec<F>
where
    F::FftField: FftField,
{
    let n = evals.len();
    let domain = Radix2EvaluationDomain::new_coset(n, domain_offset).unwrap();
    // TODO: integrate bit reverse into fft
    bit_reverse(&mut evals);
    let mut coeffs = ifft(evals, domain);
    let fold_fact = F::from(folding_factor as u64);
    for coeff in &mut coeffs {
        *coeff *= fold_fact;
    }

    let alpha_powers = (0..folding_factor)
        .map(|i| alpha.pow([i as u64]))
        .collect::<Vec<F>>();

    let drp_coeffs = ark_std::cfg_chunks!(coeffs, folding_factor)
        .map(|chunk| {
            chunk
                .iter()
                .zip(&alpha_powers)
                .map(|(v, alpha)| *v * alpha)
                .sum()
        })
        .collect::<Vec<F>>()
        .to_vec_in(GpuAllocator);

    let drp_offset = domain_offset.pow([folding_factor as u64]);
    let drp_domain = Radix2EvaluationDomain::new_coset(n / folding_factor, drp_offset).unwrap();

    // return the drp evals
    let mut evals = fft(drp_coeffs, drp_domain);
    bit_reverse(&mut evals);
    evals
}

// requires ownership when the gpu feature is enabled
#[allow(clippy::needless_pass_by_value)]
fn ifft<F: GpuField + Field + DomainCoeff<F::FftField>>(
    evals: GpuVec<F>,
    domain: Radix2EvaluationDomain<F::FftField>,
) -> GpuVec<F>
where
    F::FftField: FftField,
{
    #[cfg(feature = "gpu")]
    if domain.size() >= GpuFft::<F>::MIN_SIZE {
        let mut coeffs = evals;
        let mut ifft = GpuIfft::from(domain);
        ifft.encode(&mut coeffs);
        ifft.execute();
        return coeffs;
    }

    let coeffs = domain.ifft(&evals);
    coeffs.to_vec_in(GpuAllocator)
}

// requires ownership when the gpu feature is enabled
#[allow(clippy::needless_pass_by_value)]
fn fft<F: GpuField + Field + DomainCoeff<F::FftField>>(
    coeffs: GpuVec<F>,
    domain: Radix2EvaluationDomain<F::FftField>,
) -> GpuVec<F>
where
    F::FftField: FftField,
{
    #[cfg(feature = "gpu")]
    if domain.size() >= GpuFft::<F>::MIN_SIZE {
        let mut evals = coeffs;
        let mut fft = GpuFft::from(domain);
        fft.encode(&mut evals);
        fft.execute();
        return evals;
    }

    let evals = domain.fft(&coeffs);
    evals.to_vec_in(GpuAllocator)
}

/// # Panics
/// Panics is positions are not all unique and sorted
pub fn fold_positions(positions: &[usize], folding_factor: usize) -> Vec<usize> {
    assert!(positions.array_windows().all(|[a, b]| a < b));
    let mut res = positions
        .iter()
        .map(|p| p / folding_factor)
        .collect::<Vec<usize>>();
    res.dedup();
    res
    // let mut res = positions
    //     .iter()
    //     .map(|pos| pos % max)
    //     .collect::<Vec<usize>>();
    // res.sort_unstable();
    // res.dedup();
    // res
}

// from winterfell
pub fn get_query_values<F: Field, const N: usize>(
    chunks: &[[F; N]],
    positions: &[usize],
    folded_positions: &[usize],
) -> Vec<F> {
    positions
        .iter()
        .map(|position| {
            let i = folded_positions
                .iter()
                .position(|&v| v == position / N)
                .unwrap();
            chunks[i][position % N]
        })
        .collect()
}

fn query_layer<F: GpuField + Field, D: Digest, M: MatrixMerkleTree<F, Root = D>, const N: usize>(
    layer: &FriLayer<F, M>,
    positions: &[usize],
) -> LayerProof<F, D, M>
where
    F::FftField: FftField,
{
    let merkle_proof = layer.merkle_tree.prove_rows(positions).unwrap();
    let mut rows: Vec<[F; N]> = Vec::new();
    for &position in positions {
        let row = layer.evaluations.get_row(position).unwrap();
        rows.push(row.try_into().unwrap());
    }
    LayerProof::new(rows, merkle_proof, layer.merkle_tree.root())
}
