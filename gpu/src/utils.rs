#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub fn bit_reverse_index(n: usize, i: usize) -> usize {
    assert!(n.is_power_of_two());
    i.reverse_bits() >> (usize::BITS - n.ilog2())
}

/// Fills a slice with twiddle factors
/// TODO: Generate of the GPU <https://kieber-emmons.medium.com/9e60b974d62> or cache
/// inverse twiddles are normalized by `1 / n`.
#[cfg(feature = "arkworks")]
pub fn fill_twiddles<F: ark_ff::FftField>(dst: &mut [F], root: F) {
    #[cfg(not(feature = "parallel"))]
    let chunk_size = dst.len();
    #[cfg(feature = "parallel")]
    let chunk_size = core::cmp::max(
        dst.len() / rayon::current_num_threads().next_power_of_two(),
        1024,
    );

    ark_std::cfg_chunks_mut!(dst, chunk_size)
        .enumerate()
        .for_each(|(chunk_offset, chunk)| {
            chunk[0] = root.pow([(chunk_offset * chunk_size) as u64]);
            for i in 1..chunk.len() {
                chunk[i] = chunk[i - 1] * root;
            }
        });
}

pub fn bit_reverse_serial<T>(v: &mut [T]) {
    assert!(v.len().is_power_of_two());
    let n = v.len();
    for i in 0..n {
        let j = bit_reverse_index(n, i);
        if j > i {
            v.swap(i, j);
        }
    }
}

#[cfg(not(feature = "parallel"))]
pub fn bit_reverse<T>(v: &mut [T]) {
    bit_reverse_serial(v);
}

/// From winterfell STARK library
#[cfg(feature = "parallel")]
pub fn bit_reverse<T: Send>(v: &mut [T]) {
    const PARALLEL_THRESHOLD: usize = 1 << 17;
    if v.len() < PARALLEL_THRESHOLD {
        return bit_reverse_serial(v);
    }
    assert!(v.len().is_power_of_two());
    let n = v.len();
    let num_batches = rayon::current_num_threads().next_power_of_two();
    assert!(num_batches <= PARALLEL_THRESHOLD);
    let batch_size = n / num_batches;
    rayon::scope(|s| {
        for batch_idx in 0..num_batches {
            // create another mutable reference to the slice of values to use in a new
            // thread; this is OK because we never write the same positions in
            // the slice from different threads
            let values = unsafe { &mut *(&mut v[..] as *mut [T]) };
            s.spawn(move |_| {
                let batch_start = batch_idx * batch_size;
                let batch_end = batch_start + batch_size;
                for i in batch_start..batch_end {
                    let j = bit_reverse_index(n, i);
                    if j > i {
                        values.swap(i, j);
                    }
                }
            });
        }
    });
}

// Copies a cpu buffer to a gpu buffer
// Never use on unified memory architechture devices (M1, M2 etc.)
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn copy_to_private_buffer<T: Sized>(
    command_queue: &metal::CommandQueue,
    v: &[T],
) -> metal::Buffer {
    let device = command_queue.device();
    let byte_len = core::mem::size_of_val(v) as metal::NSUInteger;
    let ptr = v.as_ptr() as *mut core::ffi::c_void;
    let cpu_buffer =
        device.new_buffer_with_data(ptr, byte_len, metal::MTLResourceOptions::StorageModeManaged);
    let size = cpu_buffer.length();
    let private_buffer = device.new_buffer(size, metal::MTLResourceOptions::StorageModePrivate);
    let command_buffer = command_queue.new_command_buffer();
    let blit_command_encoder = command_buffer.new_blit_command_encoder();
    blit_command_encoder.copy_from_buffer(&cpu_buffer, 0, &private_buffer, 0, size);
    blit_command_encoder.end_encoding();
    command_buffer.commit();
    command_buffer.wait_until_completed();
    private_buffer
}

/// WARNING: keep the original data around or it will be freed.
// TODO: see buffer_mut_no_copy comments
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn buffer_no_copy<T: Sized>(device: &metal::DeviceRef, v: &[T]) -> metal::Buffer {
    assert!(is_page_aligned(v));
    let byte_len = core::mem::size_of_val(v);
    device.new_buffer_with_bytes_no_copy(
        v.as_ptr() as *mut core::ffi::c_void,
        byte_len.try_into().unwrap(),
        metal::MTLResourceOptions::StorageModeShared,
        None,
    )
}

/// WARNING: keep the original data around or it will be freed.
// TODO: This method previously passed a vec instead of a slice to make sure capacity was aligned to
// the page size (as per doc requirements). Seems to work in practice (on M1 at least) if only the
// pointer is aligned. Passing a slice if handy because passing a vec with any allocator requires
// nightly allocator_api feature. https://developer.apple.com/documentation/metal/mtldevice/1433382-makebuffer
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn buffer_mut_no_copy<T: Sized>(device: &metal::DeviceRef, v: &mut [T]) -> metal::Buffer {
    assert!(is_page_aligned(v));
    // TODO: once allocator_api stabilized check capacity is aligned to page size
    // this current implementation may be brittle.
    let byte_len = core::mem::size_of_val(v);
    device.new_buffer_with_bytes_no_copy(
        v.as_mut_ptr() as *mut core::ffi::c_void,
        byte_len.try_into().unwrap(),
        metal::MTLResourceOptions::StorageModeShared,
        None,
    )
}

// adapted form arkworks
/// Multiply the `i`-th element of `coeffs` with `g^i`.
#[cfg(all(target_arch = "aarch64", target_os = "macos", feature = "arkworks"))]
pub(crate) fn distribute_powers<F: crate::GpuField + ark_ff::Field>(coeffs: &mut [F], g: F) {
    let n = coeffs.len();
    #[cfg(not(feature = "parallel"))]
    let chunk_size = n;
    #[cfg(feature = "parallel")]
    let chunk_size = core::cmp::max(n / rayon::current_num_threads().next_power_of_two(), 1024);

    ark_std::cfg_chunks_mut!(coeffs, chunk_size)
        .enumerate()
        .for_each(|(chunk_offset, chunk)| {
            let offset = g.pow([(chunk_offset * chunk_size) as u64]);
            let mut pow = offset;
            chunk.iter_mut().for_each(|coeff| {
                *coeff *= pow;
                pow *= &g
            })
        });
}

/// Returns the max FFT size each threadgroup can compute
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn threadgroup_fft_size<F: crate::GpuField>(
    max_threadgroup_mem_length: usize,
    max_threads_per_threadgroup: usize,
) -> usize {
    const MIN_THREADGROUP_FFT_SIZE: usize = 1024;

    let field_size = core::mem::size_of::<F>();
    assert!(field_size * MIN_THREADGROUP_FFT_SIZE <= max_threadgroup_mem_length);
    assert!(max_threads_per_threadgroup.is_power_of_two());
    assert!(max_threads_per_threadgroup >= MIN_THREADGROUP_FFT_SIZE / 2);

    let mut fft_size = MIN_THREADGROUP_FFT_SIZE;
    // 1. don't exceed the maximum allowed threadgroup memory
    while fft_size * 2 * field_size <= max_threadgroup_mem_length {
        fft_size *= 2;
    }

    // 2. each thread operates on 2 values so can't exceed 2 * max_threads_per_tg
    core::cmp::min(fft_size, 2 * max_threads_per_threadgroup)
}

// Converts a reference to a void pointer
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub(crate) fn void_ptr<T>(v: &T) -> *const core::ffi::c_void {
    v as *const T as *const core::ffi::c_void
}

#[repr(C, align(16384))]
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
struct Page([u8; 16384]);

/// Checks a slice is page aligned on
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn is_page_aligned<T>(v: &[T]) -> bool {
    v.as_ptr().align_offset(core::mem::align_of::<Page>()) == 0
}

/// Returns a page aligned vector of the specified length with un-initialized
/// memory. This is for Apple Silicon targets only. Apple silicon supports
/// shared memory between CPU and GPU. The data resides in system memory and is
/// visible and modifiable by both the CPU and the GPU if it's page aligned.
///
/// # Safety
/// Using values from the returned vector before initializing them will lead to
/// undefined behavior.
// #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
#[allow(clippy::uninit_vec)]
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub unsafe fn page_aligned_uninit_vector<T>(length: usize) -> alloc::vec::Vec<T> {
    #[repr(C, align(16384))]
    struct Page([u8; 16384]);
    let item_size = core::mem::size_of::<T>();
    let page_size = core::mem::size_of::<Page>();
    // assert_eq!(page_size % item_size, 0, "item size must divide page size");
    let num_pages = item_size * length / page_size + 1;
    let mut aligned: alloc::vec::Vec<Page> = alloc::vec::Vec::with_capacity(num_pages);
    let ptr = aligned.as_mut_ptr();
    let capacity = num_pages * page_size / item_size;
    core::mem::forget(aligned);
    alloc::vec::Vec::from_raw_parts(ptr as *mut T, length, capacity)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_reverse_works() {
        let mut buf = alloc::vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        bit_reverse(&mut buf);

        assert_eq!(
            alloc::vec![0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15],
            buf
        );
    }

    #[test]
    #[should_panic]
    fn bit_reversal_fails_for_non_power_of_two() {
        let mut buf = alloc::vec![0; 6];

        bit_reverse(&mut buf);
    }
}
