#![feature(allocator_api, int_log)]

use objc::rc::autoreleasepool;
use polysonic::allocator::PageAlignedAllocator;
use polysonic::fields::fp_u128::BaseFelt;
use polysonic::fields::Felt;
use polysonic::fields::StarkFelt;
use polysonic::plan::NttPlanner;
use polysonic::utils::bit_reverse;
use polysonic::NttOrdering;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use std::time::Instant;

// Number theory transform
fn ntt_control<E: StarkFelt>(values: &[E]) -> Vec<E> {
    if values.len() <= 1 {
        values.to_vec()
    } else {
        let twiddle = E::get_root_of_unity(values.len().log2());
        let half = values.len() / 2;
        let odd_values = values
            .iter()
            .skip(1)
            .step_by(2)
            .copied()
            .collect::<Vec<E>>();
        let even_values = values.iter().step_by(2).copied().collect::<Vec<E>>();
        let odds = ntt_control(&odd_values);
        let evens = ntt_control(&even_values);
        (0..values.len())
            .map(|i| evens[i % half] + twiddle.pow(i.try_into().ok().unwrap()) * odds[i % half])
            .collect()
    }
}

fn gen_pcg_input<E: Felt>(n: usize) -> Vec<E, PageAlignedAllocator> {
    let mut rng = Pcg64::seed_from_u64(42); //n.try_into().unwrap());
    let mut res = Vec::new_in(PageAlignedAllocator);
    res.resize(n, E::zero());
    res.iter_mut().for_each(|v| *v = E::from(rng.gen::<u128>()));
    res
}

#[test]
fn ntt_2048_vals() {
    autoreleasepool(|| {
        let n = 2048;
        let mut input = gen_pcg_input::<BaseFelt>(n);
        let mut expected = ntt_control(&input);
        bit_reverse(&mut expected);
        let planner = NttPlanner::default();

        let mut ntt = planner.plan_ntt_forward(n, NttOrdering::Natural);
        ntt.process(&mut input);

        for (i, (a, b)) in input.iter().zip(&expected).enumerate() {
            assert_eq!(a, b, "mismatch at index {i}");
        }
    });
}

#[test]
fn ntt_524288_vals() {
    autoreleasepool(|| {
        let n = 524288;
        let mut input = gen_pcg_input::<BaseFelt>(n);
        let mut expected = ntt_control(&input);
        bit_reverse(&mut expected);
        let planner = NttPlanner::default();

        let mut ntt = planner.plan_ntt_forward(n, NttOrdering::Natural);
        ntt.process(&mut input);

        for (i, (a, b)) in input.iter().zip(&expected).enumerate() {
            assert_eq!(a, b, "mismatch at index {i}");
        }
    });
}
