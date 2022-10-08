use libc::sysconf;
use libc::_SC_PAGESIZE;
use once_cell::sync::Lazy;
use std::alloc::AllocError;
use std::alloc::Allocator;
use std::alloc::Global;
use std::alloc::Layout;
use std::ptr::NonNull;

pub static PAGE_SIZE: Lazy<usize> =
    Lazy::new(|| unsafe { sysconf(_SC_PAGESIZE).try_into().unwrap() });

pub struct PageAlignedAllocator;

unsafe impl Allocator for PageAlignedAllocator {
    fn allocate(&self, layout: Layout) -> Result<NonNull<[u8]>, AllocError> {
        Global.allocate(layout.align_to(*PAGE_SIZE).unwrap().pad_to_align())
    }

    unsafe fn deallocate(&self, ptr: NonNull<u8>, layout: Layout) {
        Global.deallocate(ptr, layout.align_to(*PAGE_SIZE).unwrap().pad_to_align())
    }
}
