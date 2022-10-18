#![feature(allocator_api)]

use ark_ff::FftField;
use ark_ff::UniformRand;
use ark_ff_optimized::fp64::Fp;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use digest::Digest;
use digest::Output;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::GpuFft;
use fast_poly::plan::GpuIfft;
use fast_poly::GpuField;
use mini_stark::merkle::MerkleTree;
use sha2::Sha256;
use std::hash;

const BENCHMARK_TREE_DEPTH: [usize; 4] = [14, 15, 16, 17];

fn bench_commitment<F: GpuField, D: Digest>(c: &mut Criterion, name: &str) {
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group(name);
    group.sample_size(10);

    for d in BENCHMARK_TREE_DEPTH {
        let n = 1 << d;
        let leaves = (0..n).map(|_| Fp::rand(&mut rng)).collect::<Vec<Fp>>();
        let leaf_nodes = leaves
            .iter()
            .map(|leaf| {
                let mut hasher = D::new();
                let mut bytes = Vec::new();
                leaf.serialize_compressed(&mut bytes).unwrap();
                hasher.update(&bytes);
                hasher.finalize()
            })
            .collect::<Vec<Output<D>>>();

        group.bench_with_input(BenchmarkId::new("new", n), &n, |b, _| {
            b.iter(|| MerkleTree::<D>::new(leaf_nodes.clone()))
        });
    }
}

fn bench_commitments(c: &mut Criterion) {
    bench_commitment::<Fp, Sha256>(c, "Sha256 Merkle Tree");
}

criterion_group!(benches, bench_commitments);
criterion_main!(benches);
