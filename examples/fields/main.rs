use ark_ff::BigInt;
use ark_ff::Field;
use ark_ff::One;
use ark_ff_optimized::fp64::Fp;

mod fq3;
use fq3::WrappedFq3 as Fq3;
// mod Fp3;
// use fp3::Fq3;

fn main() {
    println!(
        "{:?}",
        (Fp::one() + Fp::one()).pow([25769803773, 25769803769, 4294967293])
    );

    println!(
        "{:?}",
        (Fp::one() + Fp::one()).pow([25769803773, 25769803769, 4294967293])
    )

    // 16140901060737761281

    // let a = Fq3::new(Fq::one(), Fq::one(), Fq::one());
    // let b = Fq::one();

    // let domain = Radix2EvaluationDomain::<fp3>::new(8).unwrap();
    // let evals = domain.fft(&[
    //     fp3::one(),
    //     fp3::one(),
    //     fp3::one(),
    //     fp3::one(),
    //     fp3::one(),
    //     fp3::one(),
    //     fp3::one(),
    //     fp3::one(),
    // ]);

    // println!("Test: {:?}", evals);
}
