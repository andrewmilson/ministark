#![feature(allocator_api)]

use ark_ff::Field;
use ark_ff_optimized::fp64::Fp;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use digest::Digest;
use ministark::merkle::MatrixMerkleTree;
use ministark::merkle::MatrixMerkleTreeImpl;
use ministark::utils::GpuAllocator;
use ministark::Matrix;
use ministark_gpu::GpuField;
use sha2::Sha256;

const BENCHMARK_TREE_DEPTH: [usize; 4] = [14, 15, 16, 17];

fn build_merkle_tree_bench<F: GpuField + Field, D: Digest + Send + Sync>(
    c: &mut Criterion,
    name: &str,
) {
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group(name);
    group.sample_size(10);

    for d in BENCHMARK_TREE_DEPTH {
        let n = 1 << d;
        let leaves: Vec<F> = (0..n).map(|_| F::rand(&mut rng)).collect();
        let matrix = Matrix::new(vec![
            leaves.to_vec_in(GpuAllocator),
            leaves.to_vec_in(GpuAllocator),
            leaves.to_vec_in(GpuAllocator),
        ]);

        group.bench_with_input(BenchmarkId::new("from_matrix", n), &n, |b, _| {
            b.iter(|| MatrixMerkleTreeImpl::<D>::from_matrix(&matrix))
        });
    }
}

fn build_merkle_tree_benches(c: &mut Criterion) {
    build_merkle_tree_bench::<Fp, Sha256>(c, "Sha256");
}

criterion_group!(benches, build_merkle_tree_benches);
criterion_main!(benches);
