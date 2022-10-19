use crate::merkle::MerkleProof;
use crate::merkle::MerkleTree;
use crate::utils::interleave;
use crate::Matrix;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::GpuIfft;
use fast_poly::GpuField;
use fast_poly::GpuVec;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Clone, Copy)]
pub struct FriOptions {
    folding_factor: usize,
    max_remainder_size: usize,
    blowup_factor: usize,
}

impl FriOptions {
    pub fn new(blowup_factor: usize, folding_factor: usize, max_remainder_size: usize) -> Self {
        FriOptions {
            folding_factor,
            max_remainder_size,
            blowup_factor,
        }
    }

    pub fn num_layers(&self, mut domain_size: usize) -> usize {
        let mut num_layers = 0;
        while domain_size > self.max_remainder_size {
            domain_size /= self.folding_factor;
            num_layers += 1;
        }
        num_layers
    }

    pub fn remainder_size(&self, mut domain_size: usize) -> usize {
        while domain_size > self.max_remainder_size {
            domain_size /= self.folding_factor;
        }
        domain_size
    }

    pub fn domain_offset<F: GpuField>(&self) -> F {
        F::GENERATOR
    }
}

pub struct FriProver<F: GpuField, D: Digest> {
    options: FriOptions,
    layers: Vec<FriLayer<F, D>>,
}

struct FriLayer<F: GpuField, D: Digest> {
    tree: MerkleTree<D>,
    evaluations: Vec<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct FriProof<F: GpuField> {
    layers: Vec<FriProofLayer<F>>,
    remainder: Vec<F>,
}

impl<F: GpuField> FriProof<F> {
    pub fn new(layers: Vec<FriProofLayer<F>>, remainder: Vec<F>) -> Self {
        FriProof { layers, remainder }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct FriProofLayer<F: GpuField> {
    values: Vec<F>,
    proofs: Vec<MerkleProof>,
}

impl<F: GpuField> FriProofLayer<F> {
    pub fn new<const N: usize>(values: Vec<[F; N]>, proofs: Vec<MerkleProof>) -> Self {
        let values = values.into_iter().flatten().collect();
        FriProofLayer { values, proofs }
    }
}

impl<F: GpuField, D: Digest> FriProver<F, D> {
    pub fn new(options: FriOptions) -> Self {
        FriProver {
            options,
            layers: Vec::new(),
        }
    }

    pub fn into_proof(self, positions: &[usize]) -> FriProof<F> {
        let folding_factor = self.options.folding_factor;
        let (last_layer, initial_layers) = self.layers.split_last().unwrap();
        let mut domain_size = self.layers[0].evaluations.len();
        let mut proof_layers = Vec::new();
        let mut positions = positions.to_vec();
        for layer in initial_layers {
            let num_eval_chunks = domain_size / folding_factor;
            positions = fold_positions(&positions, num_eval_chunks);
            domain_size = num_eval_chunks;

            proof_layers.push(match folding_factor {
                2 => query_layer::<F, D, 2>(layer, &positions),
                4 => query_layer::<F, D, 4>(layer, &positions),
                6 => query_layer::<F, D, 6>(layer, &positions),
                8 => query_layer::<F, D, 8>(layer, &positions),
                16 => query_layer::<F, D, 16>(layer, &positions),
                _ => unimplemented!("folding factor {folding_factor} is not supported"),
            });
        }

        // layers store interlaved evaluations so they need to be un-interleaved
        let last_evals = &last_layer.evaluations;
        let mut remainder = vec![F::zero(); last_evals.len()];
        let num_eval_chunks = last_evals.len() / folding_factor;
        for i in 0..num_eval_chunks {
            for j in 0..folding_factor {
                remainder[i + num_eval_chunks * j] = last_evals[i * folding_factor + j];
            }
        }

        FriProof::new(proof_layers, remainder)
    }

    pub fn build_layers(
        &mut self,
        channel: &mut impl ProverChannel<F, Digest = D>,
        mut evaluations: Matrix<F>,
    ) {
        assert!(self.layers.is_empty());

        for _ in 0..self.options.num_layers(evaluations.num_rows()) + 1 {
            match self.options.folding_factor {
                2 => self.build_layer::<2>(channel, &mut evaluations.0[0]),
                4 => self.build_layer::<4>(channel, &mut evaluations.0[0]),
                8 => self.build_layer::<8>(channel, &mut evaluations.0[0]),
                16 => self.build_layer::<16>(channel, &mut evaluations.0[0]),
                folding_factor => unreachable!("folding factor {folding_factor} not supported"),
            }
        }
    }

    fn build_layer<const N: usize>(
        &mut self,
        channel: &mut impl ProverChannel<F, Digest = D>,
        evaluations: &mut GpuVec<F>,
    ) {
        // Each layer requires decommitting to `folding_factor` many evaluations e.g.
        // `folding_factor = 2` decommits to an evaluation for LHS_i and RHS_i
        // (0 ≤ i < n/2) which requires two merkle paths if the evaluations are
        // committed to in their natural order. If we instead commit to interleaved
        // evaluations i.e. [[LHS0, RHS0], [LHS1, RHS1], ...] LHS_i and RHS_i
        // only require a single merkle path for their decommitment.
        let interleaved_evals = interleave::<F, N>(evaluations);
        let hashed_evals = ark_std::cfg_iter!(interleaved_evals)
            .map(|chunk| {
                let mut buff = Vec::new();
                chunk.serialize_compressed(&mut buff).unwrap();
                let mut hasher = D::new();
                hasher.update(buff);
                hasher.finalize()
            })
            .collect();

        let evals_merkle_tree = MerkleTree::<D>::new(hashed_evals).unwrap();
        channel.commit_fri_layer(evals_merkle_tree.root());

        let alpha = channel.draw_fri_alpha();
        *evaluations = apply_drp(
            evaluations,
            self.options.domain_offset(),
            alpha,
            self.options.folding_factor,
        );

        self.layers.push(FriLayer {
            tree: evals_merkle_tree,
            evaluations: interleaved_evals.into_iter().flatten().collect(),
        })
    }
}

pub trait ProverChannel<F: GpuField> {
    type Digest: Digest;

    fn commit_fri_layer(&mut self, layer_root: &Output<Self::Digest>);

    fn draw_fri_alpha(&mut self) -> F;
}

/// Performs a degree respecting projection (drp) on polynomial evaluations.
/// Example for `folding_factor = 2`:
/// 1. interpolate evals over the evaluation domain to obtain f(x):
///    ┌─────────┬───┬───┬───┬───┬───┬───┬───┬───┐
///    │ i       │ 0 │ 1 │ 2 │ 3 │ 4 │ 5 │ 6 │ 7 │
///    ├─────────┼───┼───┼───┼───┼───┼───┼───┼───┤
///    │ eval[i] │ 9 │ 2 │ 3 │ 5 │ 9 │ 2 │ 3 │ 5 │
///    └─────────┴───┴───┴───┴───┴───┴───┴───┴───┘
///    ┌──────┬───────┬───────┬───────┬───────┬───────┬───────┬───────┬───────┐
///    │ x    │ o*Ω^0 │ o*Ω^1 │ o*Ω^2 │ o*Ω^3 │ o*Ω^4 │ o*Ω^5 │ o*Ω^6 │ o*Ω^7 │
///    ├──────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
///    │ f(x) │ 9     │ 2     │ 3     │ 5     │ 9     │ 2     │ 3     │ 5     │
///    └──────┴───────┴───────┴───────┴───────┴───────┴───────┴───────┴───────┘
///      f(x) = c0 * x^0 + c1 * x^1 + c2 * x^2 + c3 * x^3 +                   
///             c4 * x^4 + c5 * x^5 + c6 * x^6 + c7 * x^7                     
///    
/// 2. perform a random linear combination of odd and even coefficients of f(x):
///    f_e(x) = c0 + c2 * x + c4 * x^2 + c6 * x^3
///    f_o(x) = c1 + c3 * x + c5 * x^2 + c7 * x^3
///    f(x)   = f_e(x) + x * f_o(x)
///    f'(x)  = f_e(x) + α * f_o(x)
///   
/// 4. obtain the DRP by evaluating f'(x) over a new domain of half the size:
///    ┌───────┬───────────┬───────────┬───────────┬───────────┐
///    │ x     │ (o*Ω^0)^2 │ (o*Ω^1)^2 │ (o*Ω^2)^2 │ (o*Ω^3)^2 │
///    ├───────┼───────────┼───────────┼───────────┼───────────┤
///    │ f'(x) │ 82        │ 12        │ 57        │ 34        │
///    └───────┴───────────┴───────────┴───────────┴───────────┘
///    ┌────────┬────┬────┬────┬────┐
///    │ i      │ 0  │ 1  │ 2  │ 3  │
///    ├────────┼────┼────┼────┼────┤
///    │ drp[i] │ 82 │ 12 │ 57 │ 34 │
///    └────────┴────┴────┴────┴────┘
pub fn apply_drp<F: GpuField>(
    evals: &[F],
    domain_offset: F,
    alpha: F,
    folding_factor: usize,
) -> GpuVec<F> {
    let n = evals.len();
    let domain = Radix2EvaluationDomain::new_coset(n, domain_offset).unwrap();
    let mut coeffs = if n >= 2048 {
        let mut coeffs = evals.to_vec_in(PageAlignedAllocator);
        let mut ifft = GpuIfft::from(domain);
        ifft.encode(&mut coeffs);
        ifft.execute();
        coeffs
    } else {
        domain.ifft(&evals).to_vec_in(PageAlignedAllocator)
    };

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
        .collect::<Vec<F>>();

    let drp_offset = domain_offset.pow([folding_factor as u64]);
    let drp_domain = Radix2EvaluationDomain::new_coset(n / folding_factor, drp_offset).unwrap();
    drp_domain.fft(&drp_coeffs).to_vec_in(PageAlignedAllocator)
}

fn fold_positions(positions: &[usize], max: usize) -> Vec<usize> {
    let mut res = positions
        .iter()
        .map(|pos| pos % max)
        .collect::<Vec<usize>>();
    res.sort();
    res.dedup();
    res
}

fn query_layer<F: GpuField, D: Digest, const N: usize>(
    layer: &FriLayer<F, D>,
    positions: &[usize],
) -> FriProofLayer<F> {
    let proofs = positions
        .iter()
        .map(|pos| {
            layer
                .tree
                .prove(*pos)
                .expect("failed to generate Merkle proof")
        })
        .collect::<Vec<MerkleProof>>();
    // let chunked_evals = layer
    //     .evaluations
    //     .array_chunks::<N>()
    //     .collect::<Vec<&[F; N]>>();
    let mut values = Vec::<[F; N]>::new();
    for &position in positions {
        let i = position * N;
        let chunk = &layer.evaluations[i..i + N];
        values.push(chunk.try_into().unwrap());
    }
    FriProofLayer::new(values, proofs)
}
