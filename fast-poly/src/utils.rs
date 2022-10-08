use crate::allocator::PageAlignedAllocator;
use crate::allocator::PAGE_SIZE;
use metal::CommandQueue;

fn bit_reverse_index(n: usize, i: usize) -> usize {
    assert!(n.is_power_of_two());
    i.reverse_bits() >> (usize::BITS - n.ilog2())
}

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

pub fn copy_to_private_buffer<T: Sized>(
    command_queue: &metal::CommandQueue,
    v: &Vec<T, PageAlignedAllocator>,
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
pub fn buffer_no_copy<T: Sized>(
    device: &metal::DeviceRef,
    v: &Vec<T, PageAlignedAllocator>,
) -> metal::Buffer {
    let byte_len = round_up_to_multiple(v.len() * std::mem::size_of::<T>(), *PAGE_SIZE);
    device.new_buffer_with_bytes_no_copy(
        v.as_ptr() as *mut std::ffi::c_void,
        byte_len.try_into().unwrap(),
        metal::MTLResourceOptions::StorageModeShared,
        None,
    )
}

/// WARNING: keep the original data around or it will be freed.
pub fn buffer_mut_no_copy<T: Sized>(
    device: &metal::DeviceRef,
    v: &mut Vec<T, PageAlignedAllocator>,
) -> metal::Buffer {
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
