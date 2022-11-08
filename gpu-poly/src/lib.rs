#![feature(test, allocator_api, const_try, int_log, int_roundings)]

use allocator::PageAlignedAllocator;
use ark_ff::Field;
use ark_poly::domain::DomainCoeff;
use std::ops::MulAssign;

#[macro_use]
pub mod macros;
pub mod allocator;
pub mod fields;
pub mod plan;
pub mod prelude;
pub mod stage;
pub mod utils;

/// A trait to be implemented if the field can be used for FFTs on the GPU.
pub trait GpuFftField: GpuField<FftField = Self> + ark_ff::FftField {}

/// A trait to be implemented if `Self * Rhs` can be done on the GPU
pub trait GpuMulAssign<Rhs>: MulAssign<Rhs> {}

pub trait GpuDomainCoeff<F: GpuFftField>: DomainCoeff<F> + GpuMulAssign<F> {}

impl<T, F> GpuDomainCoeff<F> for T
where
    F: GpuFftField,
    T: DomainCoeff<F> + GpuMulAssign<F>,
{
}

// A trait for fields that have a GPU implementation
pub trait GpuField:
    Field + GpuMulAssign<Self> + GpuDomainCoeff<Self::FftField> + From<Self::FftField>
{
    type FftField: GpuFftField;

    // Used to select which GPU kernel to call.
    fn field_name() -> String;
}

/// Shared vec between GPU and CPU.
/// Requirement is that the vec's memory is page aligned.
pub type GpuVec<T> = Vec<T, PageAlignedAllocator>;
