use crate::ProofOptions;
use ark_ff_optimized::fp64::Fp;
use mini_stark::Prover;

pub struct FibProver {
    options: ProofOptions,
}

impl Prover for FibProver {
    type BaseField = Fp;
    // type Air = ;
}

enum Processor {}

fn main() {}
