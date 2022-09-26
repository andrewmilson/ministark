use ark_algebra_bench_templates::*;
use ark_ff_optimized::field_compare;
use ark_ff_optimized::fp64::Fp as Specialized;
use criterion::criterion_main;

#[derive(ark_ff::MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct FpParams;
pub type Generic = ark_ff::Fp128<ark_ff::MontBackend<FpParams, 2>>;

field_compare!(prime; "Fp=18446744069414584321"; fp18446744069414584321; Generic, Specialized);
criterion_main!(fp18446744069414584321::benches);
