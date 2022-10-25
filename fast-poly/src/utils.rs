use crate::allocator::PAGE_SIZE;
use crate::GpuVec;
use ark_ff::FftField;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

fn bit_reverse_index(n: usize, i: usize) -> usize {
    assert!(n.is_power_of_two());
    i.reverse_bits() >> (usize::BITS - n.ilog2())
}

/// Fills a slice with twiddle factors
/// TODO: Generate of the GPU https://kieber-emmons.medium.com/9e60b974d62 or cache
/// [Inverse](FftDirection::Inverse) twiddles are normalized by `1 / n`.
pub fn fill_twiddles<F: FftField>(dst: &mut [F], root: F) {
    #[cfg(not(feature = "parallel"))]
    let chunk_size = dst.len();
    #[cfg(feature = "parallel")]
    let chunk_size = std::cmp::max(
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

#[cfg(not(feature = "parallel"))]
pub fn bit_reverse<T>(v: &mut [T]) {
    assert!(v.len().is_power_of_two());
    let n = v.len();
    for i in 0..n {
        let j = bit_reverse_index(n, i);
        if j > i {
            v.swap(i, j);
        }
    }
}

/// From winterfell STARK library
#[cfg(feature = "parallel")]
pub fn bit_reverse<T: Send>(v: &mut [T]) {
    assert!(v.len().is_power_of_two());
    let n = v.len();
    let num_batches = rayon::current_num_threads().next_power_of_two();
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
pub fn copy_to_private_buffer<T: Sized>(
    command_queue: &metal::CommandQueue,
    v: &GpuVec<T>,
) -> metal::Buffer {
    let device = command_queue.device();
    let shared_buffer = buffer_no_copy(device, v);
    let size = shared_buffer.length();
    let private_buffer = device.new_buffer(size, metal::MTLResourceOptions::StorageModePrivate);
    let command_buffer = command_queue.new_command_buffer();
    let blit_command_encoder = command_buffer.new_blit_command_encoder();
    blit_command_encoder.copy_from_buffer(&shared_buffer, 0, &private_buffer, 0, size);
    blit_command_encoder.end_encoding();
    command_buffer.commit();
    command_buffer.wait_until_completed();
    private_buffer
}

/// WARNING: keep the original data around or it will be freed.
pub fn buffer_no_copy<T: Sized>(device: &metal::DeviceRef, v: &GpuVec<T>) -> metal::Buffer {
    let byte_len = round_up_to_multiple(v.len() * std::mem::size_of::<T>(), *PAGE_SIZE);
    device.new_buffer_with_bytes_no_copy(
        v.as_ptr() as *mut std::ffi::c_void,
        byte_len.try_into().unwrap(),
        metal::MTLResourceOptions::StorageModeShared,
        None,
    )
}

/// WARNING: keep the original data around or it will be freed.
pub fn buffer_mut_no_copy<T: Sized>(device: &metal::DeviceRef, v: &mut GpuVec<T>) -> metal::Buffer {
    let byte_len = round_up_to_multiple(v.len() * std::mem::size_of::<T>(), *PAGE_SIZE);
    device.new_buffer_with_bytes_no_copy(
        v.as_mut_ptr() as *mut std::ffi::c_void,
        byte_len.try_into().unwrap(),
        metal::MTLResourceOptions::StorageModeShared,
        None,
    )
}

fn round_up_to_multiple(n: usize, multiple: usize) -> usize {
    if n % multiple == 0 {
        n
    } else {
        n.next_multiple_of(multiple)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_reverse_works() {
        let mut buf = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        bit_reverse(&mut buf);

        assert_eq!(
            vec![0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15],
            buf
        );
    }

    #[test]
    #[should_panic]
    fn bit_reversal_fails_for_non_power_of_two() {
        let mut buf = vec![0; 6];

        bit_reverse(&mut buf);
    }
}
