#![feature(test, allocator_api)]

extern crate test;
use ark_ff::Field;
use ark_ff_optimized::fp64::Fp;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::Planner;
use fast_poly::plan::PLANNER;
use objc::rc::autoreleasepool;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use test::Bencher;

fn gen_pcg_input<F: Field>(n: usize) -> Vec<F, PageAlignedAllocator> {
    let mut rng = Pcg64::seed_from_u64(n.try_into().unwrap());
    let mut res = Vec::new_in(PageAlignedAllocator);
    res.resize(n, F::zero());
    res.iter_mut().for_each(|v| *v = F::rand(&mut rng));
    res
}

#[bench]
fn fft_2048_vals(b: &mut Bencher) {
    autoreleasepool(|| {
        let n = 2048;
        let mut buf = gen_pcg_input::<Fp>(n);
        let mut fft = PLANNER.plan_fft(n);

        b.iter(|| fft.process(&mut buf));
    });
}
