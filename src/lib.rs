#![warn(clippy::all, clippy::pedantic, clippy::cargo, clippy::nursery)]
#![allow(
    incomplete_features,
    clippy::must_use_candidate,
    clippy::unused_async,
    clippy::return_self_not_must_use,
    clippy::enum_glob_use,
    clippy::multiple_crate_versions,
    clippy::module_name_repetitions,
    // TODO: re-enable after writing docs
    clippy::missing_errors_doc,
    // TODO: re-enable after writing docs
    clippy::missing_panics_doc
)]
#![feature(
    allocator_api,
    let_chains,
    array_windows,
    array_chunks,
    iter_partition_in_place,
    slice_flatten,
    slice_as_chunks,
    async_fn_in_trait,
    const_trait_impl,
    cow_is_borrowed,
    exclusive_range_pattern,
    vec_into_raw_parts,
    return_position_impl_trait_in_trait,
    iter_collect_into
)]

// TODO: make some of these modules private
#[macro_use]
pub mod macros;
pub mod air;
pub mod challenges;
pub mod channel;
pub mod composer;
pub mod constraints;
pub mod debug;
pub mod eval_cpu;
pub mod eval_gpu;
pub mod expression;
pub mod fri;
pub mod hash;
pub mod hints;
pub mod matrix;
pub mod merkle;
pub mod proof;
pub mod prover;
pub mod random;
pub mod stark;
pub mod trace;
pub mod utils;
pub mod verifier;

#[macro_use]
extern crate alloc;
pub use air::Air;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use ark_poly::domain::DomainCoeff;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Mul;
use core::ops::MulAssign;
use core::ops::Sub;
use core::ops::SubAssign;
use fri::FriOptions;
pub use matrix::Matrix;
use ministark_gpu::GpuAdd;
use ministark_gpu::GpuFftField;
use ministark_gpu::GpuField;
use ministark_gpu::GpuFrom;
use ministark_gpu::GpuMul;
pub use proof::Proof;
pub use trace::Trace;

// TODO: include ability to specify:
// - base field
// - extension field
// - hashing function
#[derive(Debug, Clone, Copy, CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq)]
pub struct ProofOptions {
    pub num_queries: u8,
    pub lde_blowup_factor: u8,
    pub grinding_factor: u8,
    pub fri_folding_factor: u8,
    pub fri_max_remainder_coeffs: u8,
}

impl ProofOptions {
    pub const MIN_NUM_QUERIES: u8 = 1;
    pub const MAX_NUM_QUERIES: u8 = 128;
    pub const MIN_BLOWUP_FACTOR: u8 = 1;
    pub const MAX_BLOWUP_FACTOR: u8 = 128;
    pub const MAX_GRINDING_FACTOR: u8 = 50;

    pub const fn new(
        num_queries: u8,
        lde_blowup_factor: u8,
        grinding_factor: u8,
        fri_folding_factor: u8,
        fri_max_remainder_coeffs: u8,
    ) -> Self {
        assert!(num_queries >= Self::MIN_NUM_QUERIES);
        assert!(num_queries <= Self::MAX_NUM_QUERIES);
        assert!(lde_blowup_factor.is_power_of_two());
        assert!(lde_blowup_factor >= Self::MIN_BLOWUP_FACTOR);
        assert!(lde_blowup_factor <= Self::MAX_BLOWUP_FACTOR);
        assert!(grinding_factor <= Self::MAX_GRINDING_FACTOR);
        Self {
            num_queries,
            lde_blowup_factor,
            grinding_factor,
            fri_folding_factor,
            fri_max_remainder_coeffs,
        }
    }

    pub fn into_fri_options(self) -> FriOptions {
        // TODO: move fri params into struct
        FriOptions::new(
            self.lde_blowup_factor.into(),
            self.fri_folding_factor.into(),
            self.fri_max_remainder_coeffs.into(),
        )
    }
}

pub trait StarkExtensionOf<Fp: GpuFftField + FftField>:
    GpuField<FftField = Fp>
    + Field<BasePrimeField = Fp>
    + DomainCoeff<Fp>
    + GpuMul<Fp>
    + GpuAdd<Fp>
    + GpuFrom<Fp>
    + From<Fp>
    + MulAssign<Fp>
    + AddAssign<Fp>
    + SubAssign<Fp>
    + for<'a> MulAssign<&'a Fp>
    + for<'a> AddAssign<&'a Fp>
    + for<'a> SubAssign<&'a Fp>
    + Mul<Fp, Output = Self>
    + Add<Fp, Output = Self>
    + Sub<Fp, Output = Self>
    + for<'a> Mul<&'a Fp, Output = Self>
    + for<'a> Add<&'a Fp, Output = Self>
    + for<'a> Sub<&'a Fp, Output = Self>
{
}

impl<T, F> StarkExtensionOf<F> for T
where
    F: GpuFftField + FftField,
    T: GpuField<FftField = F>
        + Field<BasePrimeField = F>
        + DomainCoeff<F>
        + GpuMul<F>
        + GpuAdd<F>
        + GpuFrom<F>
        + MulAssign<F>
        + AddAssign<F>
        + SubAssign<F>
        + for<'a> MulAssign<&'a F>
        + for<'a> AddAssign<&'a F>
        + for<'a> SubAssign<&'a F>
        + Mul<F, Output = Self>
        + Add<F, Output = Self>
        + Sub<F, Output = Self>
        + for<'a> Mul<&'a F, Output = Self>
        + for<'a> Add<&'a F, Output = Self>
        + for<'a> Sub<&'a F, Output = Self>
        + From<F>,
{
}
