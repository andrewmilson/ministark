//! Optimized 128-bit STARK-friendly prime field with modulus `1 + 407 * 2^119`.
//!
//! The inspiration for this field came from the [Anatomy of a STARK]
//! (https://aszepieniec.github.io/stark-anatomy/basic-tools)
//! tutorial. The tutorial uses this field in Python (where integers have
//! unbounded size). This implementation has been optimized by storing field
//! elements in Montgomery form. Arithmetic in Montgomery form is significantly
//! faster due to avoiding expensive division operations. The calculation of
//! inverses has also been optimized by using the binary extended GCD algorithm.
//!
//! On my computer this optimized implementation resulted in a ~100x speedup
//! over my, first, naive implementation. For more information about Montgomery
//! arithmetic see https://en.wikipedia.org/wiki/Montgomery_modular_multiplication.
extern crate test;

use super::Felt;
use super::PrimeFelt;
use super::StarkFelt;
use num_traits::One;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;
use std::fmt::Display;
use std::hash::Hash;
use std::iter::Product;
use std::iter::Sum;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;

type PositiveInteger = u128;

/// Field modulus
const N: u128 = 1 + 407 * (1 << 119);

/// Square of auxiliary modulus `R` for montgomery reduction `R_SQUARED ≡ R^2
/// mod N`
///
/// Auxiliary modulus is `R = 1<<128` which exceeds u128::MAX which is why it
/// isn't listed as a variable. The square `R^2` is useful for converting field
/// elements into montgomery form.
///
/// Value was calculated using Python
///
/// ```python
/// R = 2**128
/// N = 2 + 407 * 2**119
/// R_SQUARED = (R * R) % N
/// print(R_SQUARED)
/// ```
const R_SQUARED: u128 = 227239200783092534449076146062029718070;

/// Integer `N'` in `[0, R − 1]` such that `N * N' ≡ −1 mod R`.
///
/// Initial calculated using Python (integers are unbounded in size in Python)
///
/// ```python
/// # Reference:
/// # https://www.geeksforgeeks.org/python-program-for-basic-and-extended-euclidean-algorithms-2/
/// def gcdExtended(a, b):
///     # Base Case
///     if a == 0:
///         return b, 0, 1
///     gcd,x1,y1 = gcdExtended(b%a, a)
///     # Update x and y using results of recursive call
///     x = y1 - (b // a) * x1
///     y = x1
///     return gcd, x, y
///
/// R = 2**128
/// N = 1 + 407 * 2**119
/// g, R_prime, N_prime = gcdExtended(R, N)
/// print(R - N_prime)
/// ```
const N_PRIME: u128 = 270497897142230380135924736767050121215;

/// Represents a base field element for the prime field with modulus `1 + 407 *
/// 2^119`
///
/// Values are stored internally in Montgomery form.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct BaseFelt(pub u128);

impl BaseFelt {
    pub const fn new(value: u128) -> BaseFelt {
        // Convert to Montgomery form
        BaseFelt(mul(value, R_SQUARED))
    }
}

impl Felt for BaseFelt {
    type PositiveInteger = PositiveInteger;

    const ELEMENT_BYTES: usize = core::mem::size_of::<u128>();

    const FIELD_ORDER_BITS: u32 = 128;

    fn inverse(&self) -> Option<Self> {
        Some(BaseFelt(modular_inverse(self.0)?))
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        match self.inverse() {
            Some(inverse) => {
                *self = inverse;
                Some(self)
            }
            None => None,
        }
    }

    fn double(&self) -> Self {
        if self.0 <= Self::MODULUS / 2 {
            BaseFelt(self.0 + self.0)
        } else {
            BaseFelt(self.0 - (Self::MODULUS - self.0))
        }
    }

    fn double_in_place(&mut self) -> &mut Self {
        *self = if self.0 <= Self::MODULUS / 2 {
            BaseFelt(self.0 + self.0)
        } else {
            BaseFelt(self.0 - (Self::MODULUS - self.0))
        };
        self
    }

    /// Converts internal value out of Montgomery form.
    fn as_integer(&self) -> Self::PositiveInteger {
        mul(1, self.0)
    }

    // TODO: find out if difference in performance if borrowed or owned self.
    fn pow(self, power: Self::PositiveInteger) -> Self {
        let mut res = Self::one();

        if power.is_zero() {
            return Self::one();
        } else if self.is_zero() {
            return Self::zero();
        }

        let mut power = power;
        let mut accumulator = self;

        while !power.is_zero() {
            if (power & Self::PositiveInteger::one()).is_one() {
                res *= accumulator;
            }
            power >>= Self::PositiveInteger::one();
            accumulator.square_in_place();
        }

        res
    }
}

impl Display for BaseFelt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_integer())
    }
}

impl From<u8> for BaseFelt {
    fn from(item: u8) -> Self {
        BaseFelt::new(item.into())
    }
}

impl From<u16> for BaseFelt {
    fn from(item: u16) -> Self {
        BaseFelt::new(item.into())
    }
}

impl From<u32> for BaseFelt {
    fn from(item: u32) -> Self {
        BaseFelt::new(item.into())
    }
}

impl From<u64> for BaseFelt {
    fn from(item: u64) -> Self {
        BaseFelt::new(item.into())
    }
}

impl From<u128> for BaseFelt {
    fn from(item: u128) -> Self {
        BaseFelt::new(item)
    }
}

impl From<usize> for BaseFelt {
    fn from(item: usize) -> Self {
        BaseFelt::new(item.try_into().unwrap())
    }
}

impl Sum for BaseFelt {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|a, b| a + b).unwrap_or_else(Self::zero)
    }
}

impl<'a> Sum<&'a Self> for BaseFelt {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.copied()
            .reduce(|a, b| a + b)
            .unwrap_or_else(Self::zero)
    }
}

impl Product for BaseFelt {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|a, b| a * b).unwrap_or_else(Self::one)
    }
}

impl<'a> Product<&'a Self> for BaseFelt {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.copied().reduce(|a, b| a * b).unwrap_or_else(Self::one)
    }
}

impl StarkFelt for BaseFelt {
    const GENERATOR: Self = BaseFelt::new(85408008396924667383611388730472331217);
    const TWO_ADICITY: u32 = 119;
    const TWO_ADIC_ROOT_OF_UNITY: Self = BaseFelt::new(85408008396924667383611388730472331217);

    /// Returns a root of unity of order `2^n`.
    fn get_root_of_unity(n: u32) -> Self {
        assert_ne!(n, 0, "n must be greater than 0");
        assert!(
            n <= Self::TWO_ADICITY,
            "n must be less than {}",
            Self::TWO_ADICITY
        );
        let power = Self::PositiveInteger::one() << (Self::TWO_ADICITY - n);
        Self::TWO_ADIC_ROOT_OF_UNITY.pow(power)
    }
}

impl PrimeFelt for BaseFelt {
    const MODULUS: u128 = N;
}

impl One for BaseFelt {
    #[inline]
    fn one() -> Self {
        BaseFelt(mul(1, R_SQUARED))
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.0 == Self::one().0
    }
}

impl Zero for BaseFelt {
    #[inline]
    fn zero() -> Self {
        BaseFelt(0)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl Neg for BaseFelt {
    type Output = Self;

    fn neg(self) -> Self::Output {
        BaseFelt(Self::MODULUS - self.0)
    }
}

impl AddAssign for BaseFelt {
    fn add_assign(&mut self, rhs: Self) {
        *self = BaseFelt(add(self.0, rhs.0))
    }
}

impl<'a> AddAssign<&'a BaseFelt> for BaseFelt {
    fn add_assign(&mut self, rhs: &Self) {
        *self = BaseFelt(add(self.0, rhs.0))
    }
}

impl Add for BaseFelt {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        BaseFelt(add(self.0, rhs.0))
    }
}

impl<'a> Add<&'a BaseFelt> for BaseFelt {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        BaseFelt(add(self.0, rhs.0))
    }
}

impl SubAssign for BaseFelt {
    fn sub_assign(&mut self, rhs: Self) {
        *self = BaseFelt(sub(self.0, rhs.0))
    }
}

impl<'a> SubAssign<&'a BaseFelt> for BaseFelt {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = BaseFelt(sub(self.0, rhs.0))
    }
}

impl Sub for BaseFelt {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        BaseFelt(sub(self.0, rhs.0))
    }
}

impl<'a> Sub<&'a BaseFelt> for BaseFelt {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        BaseFelt(sub(self.0, rhs.0))
    }
}

impl MulAssign for BaseFelt {
    fn mul_assign(&mut self, rhs: Self) {
        *self = BaseFelt(mul(self.0, rhs.0))
    }
}

impl<'a> MulAssign<&'a BaseFelt> for BaseFelt {
    fn mul_assign(&mut self, rhs: &Self) {
        *self = BaseFelt(mul(self.0, rhs.0))
    }
}

impl<'a> Mul<&'a BaseFelt> for BaseFelt {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        BaseFelt(mul(self.0, rhs.0))
    }
}

impl Mul for BaseFelt {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        BaseFelt(mul(self.0, rhs.0))
    }
}

impl DivAssign for BaseFelt {
    fn div_assign(&mut self, rhs: Self) {
        *self = BaseFelt(div(self.0, rhs.0))
    }
}

impl<'a> DivAssign<&'a BaseFelt> for BaseFelt {
    fn div_assign(&mut self, rhs: &Self) {
        *self = BaseFelt(div(self.0, rhs.0))
    }
}

impl<'a> Div<&'a BaseFelt> for BaseFelt {
    type Output = Self;

    fn div(self, rhs: &Self) -> Self::Output {
        BaseFelt(div(self.0, rhs.0))
    }
}

impl Div for BaseFelt {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        BaseFelt(div(self.0, rhs.0))
    }
}

/// Returns the multiplicative inverse in Montgomery form
///
/// `b` should be in Montgomery form. Uses binary extended GCD algorithm: https://ieeexplore.ieee.org/document/403725.
///
/// Returns None if `GCD(b, N) != 1` i.e. when `b = N` or `b = 0`.
#[inline(always)]
fn modular_inverse(b: u128) -> Option<u128> {
    let a = N;

    let mut u = a;
    let mut v = b;
    let mut r = 0;
    let mut s = 1;
    let mut k = 0;

    // First phase
    while v > 0 {
        if u & 1 == 0 {
            u /= 2;
            s *= 2;
        } else if v & 1 == 0 {
            v /= 2;
            r *= 2;
        } else if u > v {
            u = (u - v) / 2;
            r += s;
            s *= 2;
        } else {
            v = (v - u) / 2;
            s += r;

            // Can cause overflow. re-use add to reduce by the MODULUS
            r = add(r, r);
        }

        k += 1;
    }

    // GCD needs to be 1
    if u != 1 {
        return None;
    }

    if r >= a {
        r -= a;
    }

    // Second phase
    for _ in 0..(k - u128::BITS) {
        if r & 1 == 0 {
            r /= 2;
        } else {
            // Non overflowing (r + a) / 2
            r = r / 2 + a / 2 + (1 & a & r);
        }
    }

    Some(mul(a - r, R_SQUARED))
}

/// Computes `lhs * rhs mod M`
///
/// `lhs` and `rhs` are assumed to be in montgomery form. Reference:
/// - https://en.wikipedia.org/wiki/Montgomery_modular_multiplication (the REDC
///   algorithm)
/// - https://www.youtube.com/watch?v=2UmQDKcelBQ
#[inline(always)]
const fn mul(lhs: u128, rhs: u128) -> u128 {
    let lhs_low = lhs as u64 as u128;
    let lhs_high = lhs >> 64;
    let rhs_low = rhs as u64 as u128;
    let rhs_high = rhs >> 64;

    let partial_t_high = lhs_high * rhs_high;
    let partial_t_mid_a = lhs_high * rhs_low;
    let partial_t_mid_a_low = partial_t_mid_a as u64 as u128;
    let partial_t_mid_a_high = partial_t_mid_a >> 64;
    let partial_t_mid_b = rhs_high * lhs_low;
    let partial_t_mid_b_low = partial_t_mid_b as u64 as u128;
    let partial_t_mid_b_high = partial_t_mid_b >> 64;
    let partial_t_low = lhs_low * rhs_low;

    let tmp = partial_t_mid_a_low + partial_t_mid_b_low + (partial_t_low >> 64);
    let carry = tmp >> 64;
    let t_low = (tmp << 64) + partial_t_low as u64 as u128;
    let t_high = partial_t_high + partial_t_mid_a_high + partial_t_mid_b_high + carry;

    // Compute `m = T * N' mod R`
    let m = u128::wrapping_mul(t_low, N_PRIME);

    // Compute `t = (T + m * N) / R`
    let n_low = N as u64 as u128;
    let n_high = N >> 64;
    let m_low = m as u64 as u128;
    let m_high = m >> 64;

    let partial_mn_high = m_high * n_high;
    let partial_mn_mid_a = m_high * n_low;
    let partial_mn_mid_a_low = partial_mn_mid_a as u64 as u128;
    let partial_mn_mid_a_high = partial_mn_mid_a >> 64;
    let partial_mn_mid_b = n_high * m_low;
    let partial_mn_mid_b_low = partial_mn_mid_b as u64 as u128;
    let partial_mn_mid_b_high = partial_mn_mid_b >> 64;
    let partial_mn_low = m_low * n_low;

    let tmp = partial_mn_mid_a_low + partial_mn_mid_b_low + (partial_mn_low >> 64);
    let carry = tmp >> 64;
    let mn_low = (tmp << 64) + partial_mn_low as u64 as u128;
    let mn_high = partial_mn_high + partial_mn_mid_a_high + partial_mn_mid_b_high + carry;

    let (_, overflow) = u128::overflowing_add(mn_low, t_low);
    let (t, overflows_r) = u128::overflowing_add(t_high + overflow as u128, mn_high);

    let overflows_r = overflows_r as u128;
    let overflows_modulus = (t >= N) as u128;

    // TODO: overflows_r * 0u128.wrapping_sub(N) might need to be within the
    // overflows_modulus check above.
    t + overflows_r * 0u128.wrapping_sub(N) - overflows_modulus * N
}

/// Computes `lhs + rhs mod N`
#[inline(always)]
const fn add(lhs: u128, rhs: u128) -> u128 {
    let (addition, overflows_r) = u128::overflowing_add(lhs, rhs);
    let overflows_r = overflows_r as u128;
    let exceeds_modulus = (addition >= N) as u128;
    addition - exceeds_modulus * N + overflows_r * 0u128.wrapping_sub(N)
}

/// Computes `lhs - rhs mod N`
#[inline(always)]
const fn sub(lhs: u128, rhs: u128) -> u128 {
    add(lhs, BaseFelt::MODULUS - rhs)
}

/// Computes `lhs / rhs mod N`
#[inline(always)]
fn div(lhs: u128, rhs: u128) -> u128 {
    mul(lhs, modular_inverse(rhs).unwrap())
}

#[cfg(test)]
mod tests {
    use super::super::batch_inverse;
    use super::*;
    use rand::Rng;
    use rand::SeedableRng;
    use rand_pcg::Pcg64;
    use test::Bencher;

    #[test]
    fn add_adds() {
        let a = BaseFelt::new(5);
        let b = a + a;
        assert_eq!(b.as_integer(), 10);
    }

    #[test]
    fn test_batch_inverses() {
        let values = vec![BaseFelt::new(1), BaseFelt::new(2), BaseFelt::new(3)];

        let batch_inversed = batch_inverse(&values);

        assert_eq!(
            batch_inversed[0].unwrap().as_integer(),
            values[0].inverse().unwrap().as_integer()
        );
        assert_eq!(
            batch_inversed[1].unwrap().as_integer(),
            values[1].inverse().unwrap().as_integer()
        );
        assert_eq!(
            batch_inversed[2].unwrap().as_integer(),
            values[2].inverse().unwrap().as_integer()
        );
    }

    #[test]
    fn large_numbers_can_be_added() {
        // Large compared to u128::MAX
        let a = BaseFelt::new(BaseFelt::MODULUS - 2);

        assert_eq!((a + a).as_integer(), BaseFelt::MODULUS - 4);
    }

    #[test]
    fn multiplication_multiplies_large_numbers() {
        // Large compared to u128::MAX
        let a = BaseFelt::new(BaseFelt::MODULUS - 2);

        assert_eq!((a * a).as_integer(), 4);
    }

    #[test]
    fn multiplication_multiplies_odd_numbers() {
        // Large compared to u128::MAX
        let a = BaseFelt::new(BaseFelt::MODULUS - 1);
        let b = BaseFelt::new(BaseFelt::MODULUS - 2);

        assert_eq!((a * b).as_integer(), 2);
    }

    #[test]
    fn subtraction_subtracts() {
        let a = BaseFelt::new(1);
        let b = BaseFelt::new(2);

        assert_eq!((a - b).as_integer(), BaseFelt::MODULUS - 1);
    }

    #[test]
    fn division_divides() {
        let two = BaseFelt::new(2);

        assert_eq!((two / two).as_integer(), 1);
    }

    #[bench]
    fn bench_inverse_1000_items(b: &mut Bencher) {
        let mut rng = Pcg64::seed_from_u64(42);
        let items = (0..1000)
            .map(|_| BaseFelt::new(rng.gen()))
            .collect::<Vec<BaseFelt>>();

        b.iter(|| items.iter().map(|item| item.inverse()).collect::<Vec<_>>())
    }

    #[bench]
    fn bench_batch_inverse_1000_items(b: &mut Bencher) {
        let mut rng = Pcg64::seed_from_u64(42);
        let items = (0..1000)
            .map(|_| BaseFelt::new(rng.gen()))
            .collect::<Vec<BaseFelt>>();

        b.iter(|| batch_inverse(&items))
    }

    #[bench]
    fn bench_sum_20000_items(b: &mut Bencher) {
        let mut rng = Pcg64::seed_from_u64(42);
        let items = (0..20000)
            .map(|_| BaseFelt::new(rng.gen()))
            .collect::<Vec<BaseFelt>>();

        b.iter(|| items.iter().sum::<BaseFelt>())
    }

    #[bench]
    fn bench_product_1000_items(b: &mut Bencher) {
        let mut rng = Pcg64::seed_from_u64(42);
        let items = (0..1000)
            .map(|_| BaseFelt::new(rng.gen()))
            .collect::<Vec<BaseFelt>>();

        b.iter(|| items.iter().product::<BaseFelt>());

        println!("{}", items.iter().product::<BaseFelt>().as_integer());
    }
}
