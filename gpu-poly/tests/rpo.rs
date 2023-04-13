#![feature(allocator_api, array_windows)]

use ark_ff::One;
use gpu_poly::fields::p18446744069414584321::Fp;
use gpu_poly::plan::gen_rpo_merkle_tree;
use gpu_poly::plan::GpuRpo256ColumnMajor;
use gpu_poly::plan::GpuRpo256RowMajor;
use gpu_poly::prelude::*;
use std::time::Instant;

#[test]
fn gpu_rpo_from_columns() {
    let n = 8388608;
    let col = vec![Fp::one(); n];
    let col0 = col.to_vec_in(GpuAllocator);
    let col1 = col.to_vec_in(GpuAllocator);
    let col2 = col.to_vec_in(GpuAllocator);
    let col3 = col.to_vec_in(GpuAllocator);
    let col4 = col.to_vec_in(GpuAllocator);
    let col5 = col.to_vec_in(GpuAllocator);
    let col6 = col.to_vec_in(GpuAllocator);
    let col7 = col.to_vec_in(GpuAllocator);
    let input_size = 8;
    let requires_padding = input_size % 8 != 0;
    let mut rpo = GpuRpo256ColumnMajor::new(n, requires_padding);

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

    let now = Instant::now();
    let hashes = pollster::block_on(rpo.finish());
    println!("Run in {:?}", now.elapsed());
    println!("Hashes: {:?}", hashes[0]);
    println!("Hashes1: {:?}", hashes[1]);
    println!("Hashes2: {:?}", hashes[1024]);

    let now = Instant::now();
    let merkle_tree = pollster::block_on(gen_rpo_merkle_tree(&hashes));
    println!("Root: {:?}", merkle_tree[0]);
    println!("Root1: {:?}", merkle_tree[1]);
    println!("Root1: {:?}", merkle_tree[2]);
    println!("Root1: {:?}", merkle_tree[3]);
    println!("Root1: {:?}", merkle_tree[merkle_tree.len() / 2 + 1]);
    println!("Merkle tree in {:?}", now.elapsed());
}

#[test]
fn gpu_rpo_from_rows() {
    let n = 8388608;
    let rows = vec![[Fp::one(); 8]; n].to_vec_in(GpuAllocator);
    let requires_padding = false;
    let mut rpo = GpuRpo256RowMajor::new(n, requires_padding);

    let now = Instant::now();
    rpo.update(&rows);
    println!("Encode in {:?}", now.elapsed());

    let now = Instant::now();
    let hashes = pollster::block_on(rpo.finish());
    println!("Run in {:?}", now.elapsed());
    println!("Hashes (row): {:?}", hashes[0]);
    println!("Hashes1 (row): {:?}", hashes[1]);
    println!("Hashes2 (row): {:?}", hashes[1024]);
    hashes
        .array_windows()
        .enumerate()
        .for_each(|(i, [a, b])| assert_eq!(a, b, "mismatch at {i}"));

    let now = Instant::now();
    let merkle_tree = pollster::block_on(gen_rpo_merkle_tree(&hashes));
    println!("Root (row): {:?}", merkle_tree[0]);
    println!("Root1 (row): {:?}", merkle_tree[1]);
    println!("Root1 (row): {:?}", merkle_tree[2]);
    println!("Root1 (row): {:?}", merkle_tree[3]);
    println!("Root1 (row): {:?}", merkle_tree[merkle_tree.len() / 2 + 1]);
    println!("Merkle tree in {:?}", now.elapsed());
}
