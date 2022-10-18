#![feature(allocator_api)]

use ark_ff::FftField;
use ark_ff::UniformRand;
use ark_ff_optimized::fp64::Fp;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::GpuFft;
use fast_poly::plan::GpuIfft;

const BENCHMARK_INPUT_SIZES: [usize; 4] = [2048, 4096, 32768, 262144];

fn bench_fft(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group("Fast Fourier Transform");
    group.sample_size(10);

    for n in BENCHMARK_INPUT_SIZES {
        let vals = (0..n).map(|_| Fp::rand(&mut rng)).collect::<Vec<Fp>>();
        let domain = Radix2EvaluationDomain::<Fp>::new(n).unwrap();
        let coset = domain.get_coset(Fp::GENERATOR).unwrap();

        group.bench_with_input(BenchmarkId::new("FFT", n), &n, |b, _| {
            let mut coeffs = vals.to_vec_in(PageAlignedAllocator);
            b.iter(|| {
                let mut fft = GpuFft::from(domain);
                fft.encode(&mut coeffs);
                fft.execute();
            })
        });

        group.bench_with_input(BenchmarkId::new("FFT (coset)", n), &n, |b, _| {
            let mut coeffs = vals.to_vec_in(PageAlignedAllocator);
            b.iter(|| {
                let mut fft = GpuFft::from(coset);
                fft.encode(&mut coeffs);
                fft.execute();
            })
        });

        group.bench_with_input(BenchmarkId::new("IFFT", n), &n, |b, _| {
            let mut evals = vals.to_vec_in(PageAlignedAllocator);
            b.iter(|| {
                let mut fft = GpuIfft::from(domain);
                fft.encode(&mut evals);
                fft.execute();
            })
        });

        group.bench_with_input(BenchmarkId::new("IFFT (coset)", n), &n, |b, _| {
            let mut evals = vals.to_vec_in(PageAlignedAllocator);
            b.iter(|| {
                let mut fft = GpuIfft::from(coset);
                fft.encode(&mut evals);
                fft.execute();
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_fft);
criterion_main!(benches);
