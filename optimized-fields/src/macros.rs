// My main motivation with the benchmarks here was to compare the performance of arkworks
// generic field implementations with the specialized field implementations in this repo. I
// originally wanted to re-use the macros defined in the "ark-algebra-bench-templates" crate
// but Criterion's facilities for analyzing benchmarks in the aggregate didn't enable for
// comparison. I tried to use https://github.com/BurntSushi/critcmp to compare benchmarks but
// the reports where hard to read and looked bad. This is a more of less a copy of the
// benchmarks in the "ark-algebra-bench-templates" crate and parametrizes the benchmarks in
// such a way that a nice comparison table is generated in the report.
// Demo: http://andrewmilson.com/optimized-fields/criterion/report/index.html
#[macro_export]
macro_rules! field_compare {
    (prime; $test_name:expr; $mod_name:ident; $( $field:ident ),+) => {
        mod $mod_name {
            use super::*;
            use ark_ff::{Field, PrimeField, UniformRand, BigInteger};
            use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

            fn bench_compare(c: &mut Criterion) {
                const SAMPLES: usize = 1000;
                let mut group = c.benchmark_group($test_name);
                $(
                    let field_name = stringify!($field);
                    let mut rng = ark_std::test_rng();
                    let field_elements_left = (0..SAMPLES)
                        .map(|_| <$field>::rand(&mut rng))
                        .collect::<Vec<_>>();
                    let field_elements_right = (0..SAMPLES)
                        .map(|_| <$field>::rand(&mut rng))
                        .collect::<Vec<_>>();

                    // common arithmetic
                    let description = "Addition";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            field_elements_left[i] + field_elements_right[i]
                        })
                    });
                    let description = "Subtraction";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            field_elements_left[i] - field_elements_right[i]
                        })
                    });
                    let description = "Negation";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            -field_elements_left[i]
                        })
                    });
                    let description = "Double";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            field_elements_left[i].double()
                        })
                    });
                    let description = "Multiplication";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            field_elements_left[i] * field_elements_right[i]
                        })
                    });
                    let description = "Square";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            field_elements_left[i].square()
                        })
                    });
                    let description = "Inverse";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            field_elements_left[i].inverse().unwrap()
                        })
                    });
                    let description = "Sum of products of size 2";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            let j = (i + 1) % SAMPLES;
                            <$field>::sum_of_products(
                                &[field_elements_left[i], field_elements_right[j]],
                                &[field_elements_left[j], field_elements_right[i]],
                            )
                        })
                    });
                    let description = "Naive sum of products of size 2";
                    group.bench_with_input(BenchmarkId::new(field_name, description), field_name, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            let j = (i + 1) % SAMPLES;
                            field_elements_left[i] * field_elements_left[j]
                                + field_elements_right[j] * field_elements_right[i]
                        })
                    });

                    // common serialization
                    let f = field_elements_left;
                    let mut bytes = Vec::with_capacity(1000);
                    let f_compressed = f
                        .iter()
                        .map(|a| {
                            let mut bytes = Vec::with_capacity(1000);
                            a.serialize_compressed(&mut bytes).unwrap();
                            bytes
                        })
                        .collect::<Vec<_>>();
                    let f_uncompressed = f
                        .iter()
                        .map(|a| {
                            let mut bytes = Vec::with_capacity(1000);
                            a.serialize_uncompressed(&mut bytes).unwrap();
                            bytes
                        })
                        .collect::<Vec<_>>();
                    let description = "Serialize Compressed";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            bytes.clear();
                            f[i].serialize_compressed(&mut bytes).unwrap()
                        })
                    });
                    let description = "Serialize Uncompressed";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            bytes.clear();
                            f[i].serialize_uncompressed(&mut bytes).unwrap()
                        })
                    });
                    let description = "Deserialize Compressed";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            bytes.clear();
                            <$field>::deserialize_compressed(f_compressed[i].as_slice()).unwrap()
                        })
                    });
                    let description = "Deserialize Compressed Unchecked";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            bytes.clear();
                            <$field>::deserialize_compressed_unchecked(f_compressed[i].as_slice()).unwrap()
                        })
                    });
                    let description = "Deserialize Uncompressed";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            bytes.clear();
                            <$field>::deserialize_uncompressed(f_uncompressed[i].as_slice()).unwrap()
                        })
                    });
                    let description = "Deserialize Uncompressed Unchecked";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            bytes.clear();
                            <$field>::deserialize_uncompressed_unchecked(f_uncompressed[i].as_slice()).unwrap()
                        })
                    });

                    // sqrt
                    let qrs = f.iter().map(|s| s.square()).collect::<Vec<_>>();
                    let description = "Square Root for QR";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            qrs[i].sqrt().unwrap()
                        })
                    });
                    let description = "Legendre for QR";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            qrs[i].legendre()
                        })
                    });

                    // conversions
                    let bigints = f.iter().map(|f| f.into_bigint()).collect::<Vec<_>>();
                    let description = "From BigInt";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            <$field>::from_bigint(bigints[i])
                        })
                    });
                    let description = "Into BigInt";
                    group.bench_with_input(BenchmarkId::new(field_name, description), description, |b, _| {
                        let mut i = 0;
                        b.iter(|| {
                            i = (i + 1) % SAMPLES;
                            f[i].into_bigint()
                        })
                    });
                )*
                group.finish();
            }

            criterion::criterion_group!(benches, bench_compare);
        }
    };
}
