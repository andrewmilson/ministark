#![feature(int_log)]

use legacy_algebra::fp_u128::BaseFelt;
use legacy_algebra::StarkFelt;
use num_traits::One;
use num_traits::Zero;
use std::fs;
use std::io::Write;
// use zk_wasm::felt::prime_u128::BaseFelt;
// use zk_wasm::felt::StarkFelt;

fn main() -> std::io::Result<()> {
    for variant in [2048, 16384] {
        let mut elements = vec![BaseFelt::zero(); variant];
        for i in 0..variant {
            elements[i] = BaseFelt::new((i + 1) as u128);
        }

        let element_bytes = elements
            .iter()
            .flat_map(|e| e.0.to_ne_bytes())
            .collect::<Vec<u8>>();

        let root_of_unity = BaseFelt::get_root_of_unity(variant.ilog2());
        let mut acc = BaseFelt::one();
        let mut twiddles = vec![BaseFelt::zero(); variant];
        for i in 0..variant {
            acc *= root_of_unity;
            twiddles[i] = acc;
        }
        let twiddle_bytes = twiddles
            .iter()
            .flat_map(|t| t.0.to_ne_bytes())
            .collect::<Vec<u8>>();

        let mut values_file = fs::File::create(format!(
            "/tmp/ntt_montgomery_felt_u128_x{}_twiddles.bin",
            variant
        ))?;
        // Write a slice of bytes to the file
        values_file.write_all(&twiddle_bytes)?;

        let mut values_file = fs::File::create(format!(
            "/tmp/ntt_montgomery_felt_u128_x{}_input_vals.bin",
            variant
        ))?;
        // Write a slice of bytes to the file
        values_file.write_all(&element_bytes)?;

        let mut outputs_file = fs::File::create(format!(
            "/tmp/ntt_montgomery_felt_u128_x{}_output_vals.bin",
            variant
        ))?;
        // Write a slice of bytes to the file
        outputs_file.write_all(&element_bytes)?;
    }

    // for 2048

    Ok(())
}
