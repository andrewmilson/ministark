use crate::allocator::PageAlignedAllocator;

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

pub fn buffer_no_copy<T: Sized>(
    device: &metal::DeviceRef,
    v: &mut Vec<T, PageAlignedAllocator>,
) -> metal::Buffer {
    device.new_buffer_with_bytes_no_copy(
        v.as_mut_ptr() as *mut std::ffi::c_void,
        (v.len() * std::mem::size_of::<T>()).try_into().unwrap(),
        metal::MTLResourceOptions::StorageModeShared,
        None,
    )
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
