use crate::merkle::MerkleTree;
use crate::Matrix;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;
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

pub struct FriProof {
    layers: Vec<FriProofLayer>,
    remainder: Vec<u8>,
    // num_partitions: u8,
}

pub struct FriProofLayer {
    values: Vec<u8>,
    paths: Vec<u8>,
}

impl<F: GpuField, D: Digest> FriProver<F, D> {
    pub fn new(options: FriOptions) -> Self {
        FriProver {
            options,
            layers: Vec::new(),
        }
    }

    pub fn build_proof(&mut self, positions: &[usize]) -> FriProof {
        todo!()
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
        let chunked_evals = evaluations.array_chunks().copied().collect::<Vec<[F; N]>>();
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
fn apply_drp<F: GpuField, const N: usize>(
    values: &[[F; N]],
    domain_offset: F,
    alpha: F,
) -> GpuVec<F> {
    // TODO: whole thing can use optimization
    let n = values.len();
    let domain = Radix2EvaluationDomain::new_coset(n, domain_offset).unwrap();
    let mut drp = vec![F::zero(); n];

    ark_std::cfg_iter_mut!(drp)
        .zip(values)
        .for_each(|(drp, values)| {
            let mut coeffs = values.to_vec();
            domain.ifft_in_place(&mut coeffs);
            let poly = DensePolynomial::from_coefficients_vec(coeffs);
            *drp = poly.evaluate(&alpha);
        });

    drp.to_vec_in(PageAlignedAllocator)
}
