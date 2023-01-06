#![feature(allocator_api)]

use ark_ff::One;
use ark_ff::UniformRand;
use gpu_poly::fields::p18446744069414584321::Fp;
use gpu_poly::plan::GpuRpo;
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
    let mut rpo = GpuRpo::new(n);

    let now = Instant::now();
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    rpo.update(&col0, &col1, &col2, &col3, &col4, &col5, &col6, &col7);
    println!("Encode in {:?}", now.elapsed());

    let hashes = rpo.finish();
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
