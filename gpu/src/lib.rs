#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(any(feature = "arkworks", feature = "winterfell")))]
compile_error!("Either feature \"arkworks\" or \"winterfell\" must be enabled for this crate.");

extern crate alloc;

#[macro_use]
pub mod macros;
pub mod fields;
pub mod plan;
pub mod prelude;
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub mod stage;
pub mod utils;

#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub use metal;

/// A trait to be implemented if the field can be used for FFTs on the GPU.
pub trait GpuFftField: GpuField<FftField = Self> {}

/// A marker trait to be implemented if `Self * Rhs` can be done on the GPU
pub trait GpuMul<Rhs> {}

/// A marker trait to be implemented if `Self + Rhs` can be done on the GPU
pub trait GpuAdd<Rhs> {}

/// A marker trait indicating `Self` can be converted to `Rhs` on the GPU
pub trait GpuFrom<Rhs> {}

/// A marker trait for fields that have a GPU implementation
pub trait GpuField: GpuMul<Self> + GpuAdd<Self> + GpuMul<Self::FftField> + Sized {
    type FftField: GpuFftField;

    // Used to select which GPU kernel to call.
    fn field_name() -> alloc::string::String;
}
