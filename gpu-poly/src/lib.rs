#![feature(test, allocator_api, const_try, int_log, int_roundings)]

use allocator::PageAlignedAllocator;
use ark_ff::Field;
use ark_poly::domain::DomainCoeff;
use std::ops::Mul;

#[macro_use]
pub mod macros;
pub mod allocator;
pub mod fields;
pub mod plan;
pub mod prelude;
pub mod stage;
pub mod utils;

/// A trait to be implemented if the field can be used for FFTs on the GPU.
pub trait GpuFftField: GpuField + ark_ff::FftField {}

/// A trait to be implemented if `Self * Rhs` can be done on the GPU
pub trait GpuMul<Rhs, Output>: Mul<Rhs, Output = Output> {}

pub trait GpuDomainCoeff<F: GpuFftField>: DomainCoeff<F> + GpuMul<F, Self> {}

impl<T, F> GpuDomainCoeff<F> for T
where
    F: GpuFftField,
    T: DomainCoeff<F> + GpuMul<F, T>,
{
}

// A trait for fields that have a GPU implementation
pub trait GpuField: Field + GpuDomainCoeff<Self::FftField> + From<Self::FftField> {
    type FftField: GpuFftField;

    // Used to select which GPU kernel to call.
    fn field_name() -> String;
}

/// Shared vec between GPU and CPU.
/// Requirement is that the vec's memory is page aligned.
pub type GpuVec<T> = Vec<T, PageAlignedAllocator>;
