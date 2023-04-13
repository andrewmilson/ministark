#![feature(allocator_api)]
#![cfg_attr(not(feature = "std"), no_std)]

#[macro_use]
extern crate alloc;

#[macro_use]
pub mod macros;
pub mod fields;
pub mod plan;
pub mod prelude;
pub mod stage;
pub mod utils;

/// A trait to be implemented if the field can be used for FFTs on the GPU.
pub trait GpuFftField: GpuField<FftField = Self> {}

/// A marker trait to be implemented if `Self * Rhs` can be done on the GPU
pub trait GpuMul<Rhs> {}

/// A marker trait to be implemented if `Self + Rhs` can be done on the GPU
pub trait GpuAdd<Rhs> {}

// A marker trait for fields that have a GPU implementation
pub trait GpuField: GpuMul<Self> + GpuAdd<Self> + GpuMul<Self::FftField> + Sized {
    type FftField: GpuFftField;

    // Used to select which GPU kernel to call.
    fn field_name() -> alloc::string::String;
}

pub use gpu_vec::GpuVec;

mod gpu_vec {
    use ark_std::alloc::Global;
    use core::alloc::AllocError;
    use core::alloc::Allocator;
    use core::alloc::Layout;
    use core::ptr::NonNull;
    use once_cell::sync::Lazy;

    /// Shared vec between GPU and CPU.
    /// Requirement is that the vec's memory is page aligned.
    #[cfg(apple_silicon)]
    pub type GpuVec<T> = alloc::vec::Vec<T, GpuAllocator>;

    /// Allocator with page aligned allocations on Apple Silicon.
    /// Uses global allocator on all other platforms.
    pub struct GpuAllocator;

    unsafe impl Allocator for GpuAllocator {
        fn allocate(&self, layout: Layout) -> Result<NonNull<[u8]>, AllocError> {
            #[cfg(apple_silicon)]
            return page_aligned_allocator::PageAlignedAllocator.allocate(layout);
            #[cfg(not(apple_silicon))]
            return Global.allocate(layout);
        }

        unsafe fn deallocate(&self, ptr: NonNull<u8>, layout: Layout) {
            #[cfg(apple_silicon)]
            return page_aligned_allocator::PageAlignedAllocator.deallocate(ptr, layout);
            #[cfg(not(apple_silicon))]
            return Global.deallocate(ptr, layout);
        }
    }

    #[cfg(apple_silicon)]
    mod page_aligned_allocator {
        use super::*;

        pub static PAGE_SIZE: Lazy<usize> =
            Lazy::new(|| unsafe { libc::sysconf(libc::_SC_PAGESIZE).try_into().unwrap() });

        pub struct PageAlignedAllocator;

        unsafe impl Allocator for PageAlignedAllocator {
            fn allocate(&self, layout: Layout) -> Result<NonNull<[u8]>, AllocError> {
                Global.allocate(layout.align_to(*PAGE_SIZE).unwrap().pad_to_align())
            }

            unsafe fn deallocate(&self, ptr: NonNull<u8>, layout: Layout) {
                Global.deallocate(ptr, layout.align_to(*PAGE_SIZE).unwrap().pad_to_align())
            }
        }
    }
}
