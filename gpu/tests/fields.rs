#![cfg(all(target_arch = "aarch64", target_os = "macos"))]
#![feature(allocator_api)]

use ark_ff::UniformRand;
use ministark_gpu::stage::MulPowStage;
use ministark_gpu::utils::buffer_mut_no_copy;
use ministark_gpu::utils::buffer_no_copy;

pub mod p18446744069414584321 {
    use super::*;
    use ark_ff::Field;
    use ministark_gpu::fields::p18446744069414584321::ark::Fp;
    use ministark_gpu::fields::p18446744069414584321::ark::Fq3;
    use ministark_gpu::prelude::get_planner;
    use ministark_gpu::utils::page_aligned_uninit_vector;

    #[test]
    fn mul_pow_fp() {
        let n = 2048;
        let mut rng = &mut ark_std::test_rng();
        let mut a = unsafe { page_aligned_uninit_vector(n) };
        a.fill_with(|| Fp::rand(&mut rng));
        let mut b = unsafe { page_aligned_uninit_vector(n) };
        b.fill_with(|| Fp::rand(&mut rng));
        let expected = a
            .iter()
            .copied()
            .zip(&b)
            .map(|(mut a, b)| {
                a *= b;
                a
            })
            .collect::<Vec<Fp>>();
        let planner = get_planner();
        let command_queue = &planner.command_queue;
        let mut a_buffer = buffer_mut_no_copy(command_queue.device(), &mut a);
        let b_buffer = buffer_no_copy(command_queue.device(), &b);
        let command_buffer = command_queue.new_command_buffer();

        let multiplier = MulPowStage::<Fp>::new(&planner.library, n);
        multiplier.encode(command_buffer, &mut a_buffer, &b_buffer, 1, 0);
        command_buffer.commit();
        command_buffer.wait_until_completed();

        for (i, (expected, actual)) in expected.into_iter().zip(a).enumerate() {
            assert_eq!(expected, actual, "mismatch at index {i}");
        }
    }

    #[test]
    fn mul_pow_fq3_by_fp() {
        let n = 2048;
        let mut rng = &mut ark_std::test_rng();
        let mut a = unsafe { page_aligned_uninit_vector(n) };
        a.fill_with(|| Fq3::rand(&mut rng));
        let mut b = unsafe { page_aligned_uninit_vector(n) };
        b.fill_with(|| Fp::rand(&mut rng));
        let expected = a
            .iter()
            .copied()
            .zip(&b)
            .map(|(mut a, b)| {
                a *= b;
                a
            })
            .collect::<Vec<Fq3>>();
        let planner = get_planner();
        let command_queue = &planner.command_queue;
        let mut a_buffer = buffer_mut_no_copy(command_queue.device(), &mut a);
        let b_buffer = buffer_no_copy(command_queue.device(), &b);
        let command_buffer = command_queue.new_command_buffer();

        let multiplier = MulPowStage::<Fq3, Fp>::new(&planner.library, n);
        multiplier.encode(command_buffer, &mut a_buffer, &b_buffer, 1, 0);
        command_buffer.commit();
        command_buffer.wait_until_completed();

        for (i, (expected, actual)) in expected.into_iter().zip(a).enumerate() {
            assert_eq!(expected, actual, "mismatch at index {i}");
        }
    }

    #[test]
    fn mul_pow_fq3() {
        use ark_ff::One;
        println!("{:?}", Fq3::one());

        let n = 2048;
        let mut rng = &mut ark_std::test_rng();
        let mut a = unsafe { page_aligned_uninit_vector(n) };
        a.fill_with(|| Fq3::rand(&mut rng));
        let mut b = unsafe { page_aligned_uninit_vector(n) };
        b.fill_with(|| Fq3::rand(&mut rng));
        let expected = a
            .iter()
            .copied()
            .zip(&b)
            .map(|(mut a, b)| {
                a *= b.square() * b;
                a
            })
            .collect::<Vec<Fq3>>();
        let planner = get_planner();
        let command_queue = &planner.command_queue;
        let mut a_buffer = buffer_mut_no_copy(command_queue.device(), &mut a);
        let b_buffer = buffer_no_copy(command_queue.device(), &b);
        let command_buffer = command_queue.new_command_buffer();

        let multiplier = MulPowStage::<Fq3>::new(&planner.library, n);
        multiplier.encode(command_buffer, &mut a_buffer, &b_buffer, 3, 0);
        command_buffer.commit();
        command_buffer.wait_until_completed();

        for (i, (expected, actual)) in expected.into_iter().zip(a).enumerate() {
            assert_eq!(expected, actual, "mismatch at index {i}");
        }
    }
}
