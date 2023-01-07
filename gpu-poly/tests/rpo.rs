#![feature(allocator_api)]

use ark_ff::One;
use gpu_poly::fields::p18446744069414584321::Fp;
use gpu_poly::plan::GpuRpo128;
use gpu_poly::prelude::*;
use std::time::Instant;

#[test]
fn gpu_rpo() {
    let n = 8388608;
    let col = vec![Fp::one(); n];
    let col0 = col.to_vec_in(PageAlignedAllocator);
    let col1 = col.to_vec_in(PageAlignedAllocator);
    let col2 = col.to_vec_in(PageAlignedAllocator);
    let col3 = col.to_vec_in(PageAlignedAllocator);
    let col4 = col.to_vec_in(PageAlignedAllocator);
    let col5 = col.to_vec_in(PageAlignedAllocator);
    let col6 = col.to_vec_in(PageAlignedAllocator);
    let col7 = col.to_vec_in(PageAlignedAllocator);
    let input_size = 104;
    let requires_padding = input_size % 8 == 0;
    let mut rpo = GpuRpo128::new(n, requires_padding);

    let now = Instant::now();
    for i in 0..input_size {
        let col = match i % 8 {
            0 => &col0,
            1 => &col1,
            2 => &col2,
            3 => &col3,
            4 => &col4,
            5 => &col5,
            6 => &col6,
            7 => &col7,
            _ => unreachable!(),
        };

        rpo.update(col);
    }
    println!("Encode in {:?}", now.elapsed());

    let hashes = pollster::block_on(rpo.finish());
    println!("Hash in {:?}", now.elapsed());
    println!(
        "{}, {}, {}, {}",
        hashes[1][0], hashes[1][1], hashes[1][2], hashes[1][3]
    );
    println!("{:?}", hashes[0]);
    println!("{:?}", hashes[1]);
    println!("{:?}", hashes[2]);
    println!("{:?}", hashes[3]);
    println!("{:?}", hashes[4]);
    println!("{:?}", hashes[5]);
    println!("{:?}", hashes[6]);
    println!("{:?}", hashes[7]);
    println!("{:?}", hashes[8]);
    println!("{:?}", hashes[9]);
}
