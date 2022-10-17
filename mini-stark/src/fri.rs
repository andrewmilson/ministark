use crate::merkle::MerkleProof;
use crate::merkle::MerkleTree;
use crate::utils::interleave;
use crate::Matrix;
use ark_poly::evaluations;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use fast_poly::allocator::PageAlignedAllocator;
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
        let chunked_evals = interleave::<F, N>(&evaluations);
        let hashed_evals = ark_std::cfg_iter!(chunked_evals)
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
        *evaluations = apply_drp(&chunked_evals, self.options.domain_offset(), alpha);

        self.layers.push(FriLayer {
            tree: evals_merkle_tree,
            evaluations: chunked_evals.into_iter().flatten().collect(),
        })
    }
}

pub trait ProverChannel<F: GpuField> {
    type Digest: Digest;

    fn commit_fri_layer(&mut self, layer_root: &Output<Self::Digest>);

    fn draw_fri_alpha(&mut self) -> F;
}

// Apply degree respecting projection
pub fn apply_drp<F: GpuField, const N: usize>(
    values: &[[F; N]],
    domain_offset: F,
    alpha: F,
) -> GpuVec<F> {
    // TODO: whole thing can use optimization
    let n = values.len();
    let domain = Radix2EvaluationDomain::new_coset(N, domain_offset).unwrap();

    let mut drp = vec![F::zero(); n];

    ark_std::cfg_iter_mut!(drp)
        .zip(values)
        .for_each(|(drp, values)| {
            let mut coeffs = values.to_vec();
            domain.ifft_in_place(&mut coeffs);
            println!("<evals>");
            prety_print(values);
            println!("</evals>");
            println!("<poly>");
            prety_print(&coeffs);
            println!("</poly>");
            let poly = DensePolynomial::from_coefficients_vec(coeffs);
            *drp = poly.evaluate(&alpha);
        });

    drp.to_vec_in(PageAlignedAllocator)
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

fn prety_print(vals: &[impl GpuField]) {
    for val in vals {
        println!("{val},");
    }
}
