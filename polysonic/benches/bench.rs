#![feature(test, allocator_api)]

extern crate test;
use objc::rc::autoreleasepool;
use polysonic::allocator::PageAlignedAllocator;
use polysonic::fields::fp_u128::BaseFelt;
use polysonic::fields::Felt;
use polysonic::plan::NttPlanner;
use polysonic::NttOrdering;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use test::Bencher;

fn gen_pcg_input<E: Felt>(n: usize) -> Vec<E, PageAlignedAllocator> {
    let mut rng = Pcg64::seed_from_u64(n.try_into().unwrap());
    let mut res = Vec::new_in(PageAlignedAllocator);
    res.resize(n, E::zero());
    res.iter_mut().for_each(|v| *v = E::from(rng.gen::<u128>()));
    res
}

#[bench]
fn ntt_2048_vals(b: &mut Bencher) {
    autoreleasepool(|| {
        let n = 2048;
        let mut buf = gen_pcg_input::<BaseFelt>(n);
        let planner = NttPlanner::default();
        let mut ntt = planner.plan_ntt_forward(n, NttOrdering::Natural);

        b.iter(|| ntt.process(&mut buf));
    });
}
