use super::AffinePoint;
use super::WeierstrassCurveConfig;
use crate::fields::fp_u256::U256;

const ODD_PRIMES: [u32; 53] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
];

/// Returns the size of the bound that the number of points on an elliptic curve
/// can be found within
const fn hasse_bound_size(field_size: U256) -> U256 {
    // .sqrt() is floored so add 1 to overestimate the size of the bound
    // 4 * (floor(sqrt(q)) + 1)
    field_size
        .sqrt()
        .add(U256::ONE)
        .mul(U256 { high: 0, low: 4 })
}

/// Returns the smallest list of odd primes that that multiply to at least n
fn inverse_primorial(n: U256) -> &'static [u32] {
    let mut accumulator = U256::ONE;
    let mut i = 0;
    while accumulator.lt(n) {
        accumulator = accumulator.mul(U256 {
            high: 0,
            low: ODD_PRIMES[i] as u128,
        });
        i += 1;
    }
    &ODD_PRIMES[0..i]
}

fn schoofs_algorithm() {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::fp_u256::BaseFelt;
    use crate::fields::fp_u256::U256;
    use crate::fields::PrimeFelt;

    #[test]
    fn produces_small_primes() {
        let n = hasse_bound_size(BaseFelt::MODULUS.into());

        let primes = inverse_primorial(n);

        // assert_eq!(i, 10);
        // println!("{}, {}", ODD_PRIMES[i], accumulator);
        assert_eq!(primes.len(), 26);
    }
}
