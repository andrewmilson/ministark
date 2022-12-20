#![feature(test, allocator_api, const_try, int_roundings)]
#![no_std]

#[macro_use]
extern crate alloc;
use alloc::string::String;
use alloc::vec::Vec;
use allocator::PageAlignedAllocator;
use ark_ff::Field;
use ark_poly::domain::DomainCoeff;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Mul;
use core::ops::MulAssign;

#[macro_use]
pub mod macros;
pub mod allocator;
pub mod fields;
pub mod plan;
pub mod prelude;
pub mod stage;
pub mod utils;

/// A marker trait to be implemented if the field can be used for FFTs on the
/// GPU.
pub trait GpuFftField: GpuField<FftField = Self> + ark_ff::FftField {}

/// A marker trait to be implemented if `Self *= Rhs` can be done on the GPU
pub trait GpuMul<Rhs>: MulAssign<Rhs> + Mul<Rhs, Output = Self> {}

/// A marker trait to be implemented if `Self += Rhs` can be done on the GPU
pub trait GpuAdd<Rhs>: AddAssign<Rhs> + Add<Rhs, Output = Self> {}

// A marker trait for fields that have a GPU implementation
// TODO: consider removing DomainCoeff
pub trait GpuField:
    Field + DomainCoeff<Self::FftField> + GpuMul<Self> + GpuAdd<Self> + From<Self::FftField>
{
    type FftField: GpuFftField;

    // Used to select which GPU kernel to call.
    fn field_name() -> String;
}

/// Shared vec between GPU and CPU.
/// Requirement is that the vec's memory is page aligned.
pub type GpuVec<T> = Vec<T, PageAlignedAllocator>;
