#![feature(int_roundings)]

use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ministark_gpu::fields::p3618502788666131213697322783095070105623107215331596699973092056135872020481::ark::Fp;
use ark_ff::One;
use ministark::air::AirConfig;
use ministark::constraints::Constraint;
use ministark::utils::FieldVariant;
use crate::rescue::Rescue;

mod rescue;

#[derive(Clone, Copy, CanonicalSerialize, CanonicalDeserialize)]
struct RescueInfo {
    input: [Fp; 2],
    output: [Fp; 2],
}

struct RescueAirConfig;

impl AirConfig for RescueAirConfig {
    const NUM_BASE_COLUMNS: usize = 0;
    type Fp = Fp;
    type Fq = Fp;
    type PublicInputs = ();

    fn constraints(_trace_len: usize) -> Vec<Constraint<FieldVariant<Self::Fp, Self::Fq>>> {
        todo!()
    }
}

fn main() {
    let state_width = 4; /* =m */
    // TODO: this may not be accurate. Generate with Algorithm 7
    let rounds = 14; /* =N */
    let security_level = 256;
    let capacity = 2;
    let digest_size = 2;

    let mut hasher = Rescue::new(state_width, capacity, rounds, security_level, digest_size);
    hasher.update(Fp::one());
    hasher.update(Fp::one());
    let _output = hasher.finish();

    todo!("Rescue example is a WIP");
}
