#![feature(test, allocator_api, const_try, int_log)]

use ark_ff::FftField;
use ark_ff_optimized::fp64;

pub mod allocator;
pub mod plan;
pub mod stage;
pub mod twiddles;
pub mod utils;

/// Represents a NTT direction
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum NttDirection {
    Forward,
    Inverse,
}

// GPU implementation of the field exists in metal/
pub trait GpuField: FftField {
    // Used to select which GPU kernel to call.
    fn field_name() -> String;
}

impl GpuField for fp64::Fp {
    fn field_name() -> String {
        "fp18446744069414584321".to_owned()
    }
}
