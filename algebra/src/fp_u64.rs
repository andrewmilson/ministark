//! An implementation of a 64-bit STARK-friendly prime field with modulus `2^64
//! - 2^32 + 1` Code taken and adapted from the winterfell STARK library.
extern crate test;

use super::ExtensibleField;
use super::Felt;
use super::PrimeFelt;
use super::StarkFelt;
use super::Univariate;
use num_bigint::BigUint;
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

/// Field modulus `2^64 - 2^32 + 1`
const N: u64 = 0xFFFFFFFF00000001;

/// Square of auxiliary modulus `R` for montgomery reduction `R_SQUARED ≡ R^2
/// mod N`
///
/// Auxiliary modulus is `R = 1<<64` which exceeds u64::MAX which is why it
/// isn't listed as a variable. The square `R^2` is useful for converting field
/// elements into montgomery form.
///
/// Value was calculated using Python
///
/// ```python
/// R = 2**64
/// N = 2**64 - 2**32 + 1
/// R_SQUARED = (R * R) % N
/// print(R_SQUARED)
/// ```
const R_SQUARED: u64 = 18446744065119617025;

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
/// R = 2**64
/// N = 2**64 - 2**32 + 1
/// g, R_prime, N_prime = gcdExtended(R, N)
/// print(N_prime)
/// ```
// const N_PRIME: u64 = 4294967297;

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct BaseFelt(pub u64);

impl BaseFelt {
    pub const fn new(value: u64) -> BaseFelt {
        // Convert to Montgomery form
        BaseFelt(mul(value, R_SQUARED))
    }
}

impl Felt for BaseFelt {
    fn inverse(&self) -> Option<Self> {
        Some(BaseFelt(modular_inverse(self.0)))
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

    // Computes the identity in a prime field
    fn frobenius(&mut self) {}
}

impl Display for BaseFelt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.into_bigint())
    }
}

impl PartialEq for BaseFelt {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        equals(self.0, other.0) == 0xFFFFFFFFFFFFFFFF
    }
}

impl Eq for BaseFelt {}

impl Hash for BaseFelt {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.into_bigint().hash(state);
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
        BaseFelt::new(item)
    }
}

impl From<u128> for BaseFelt {
    fn from(item: u128) -> Self {
        BaseFelt::new((item % N as u128) as u64)
    }
}

impl From<usize> for BaseFelt {
    fn from(item: usize) -> Self {
        BaseFelt::new(item.try_into().unwrap())
    }
}

impl From<BigUint> for BaseFelt {
    fn from(item: BigUint) -> Self {
        todo!()
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
    /// Generates entire multiplicitive group
    const GENERATOR: Self = BaseFelt::new(7);

    const TWO_ADICITY: u32 = 32;

    /// Generates multiplicitive group of order `2^32`
    const TWO_ADIC_ROOT_OF_UNITY: Self = BaseFelt::new(1753635133440165772);

    /// Returns a root of unity of order `2^n`.
    fn get_root_of_unity(n: u32) -> Self {
        assert_ne!(n, 0, "n must be greater than 0");
        assert!(
            n <= Self::TWO_ADICITY,
            "n must be less than {}",
            Self::TWO_ADICITY
        );
        let power = 1u64 << (Self::TWO_ADICITY - n);
        Self::TWO_ADIC_ROOT_OF_UNITY.pow(&[power])
    }
}

impl PrimeFelt for BaseFelt {
    type BigInt = u64;
    const MODULUS: u64 = N;

    fn into_bigint(self) -> Self::BigInt {
        mul(1, self.0)
    }
}

impl Into<BigUint> for BaseFelt {
    fn into(self) -> BigUint {
        self.into_bigint().into()
    }
}

impl Into<u64> for BaseFelt {
    fn into(self) -> u64 {
        self.into_bigint()
    }
}

impl One for BaseFelt {
    #[inline]
    fn one() -> Self {
        BaseFelt(mul(1, R_SQUARED))
    }

    #[inline]
    fn is_one(&self) -> bool {
        equals(self.0, Self::one().0) == 0xFFFFFFFFFFFFFFFF
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
/// Returns 0 if `a` is 0
#[inline(always)]
const fn modular_inverse(a: u64) -> u64 {
    // compute base^(N - 2) using 72 multiplications
    // M - 2 = 0b1111111111111111111111111111111011111111111111111111111111111111
    let t2 = mul(mul(a, a), a);
    let t3 = mul(mul(t2, t2), a);
    let t6 = exp_acc::<3>(t3, t3);
    let t12 = exp_acc::<6>(t6, t6);
    let t24 = exp_acc::<12>(t12, t12);
    let t30 = exp_acc::<6>(t24, t6);
    let t31 = mul(mul(t30, t30), a);
    let t63 = exp_acc::<32>(t31, t31);
    mul(mul(t63, t63), a)
}

/// Squares `base` N times and multiplies the result by the tail value.
#[inline(always)]
const fn exp_acc<const N: usize>(base: u64, tail: u64) -> u64 {
    let mut result = base;
    let mut i = 0;
    while i < N {
        result = mul(result, result);
        i += 1;
    }
    mul(result, tail)
}

/// Computes `lhs * rhs mod M`
///
/// `lhs` and `rhs` are assumed to be in montgomery form.
#[inline(always)]
const fn mul(lhs: u64, rhs: u64) -> u64 {
    let x = lhs as u128 * rhs as u128;
    let xl = x as u64;
    let xh = (x >> 64) as u64;
    let (a, e) = xl.overflowing_add(xl << 32);
    let b = a.wrapping_sub(a >> 32).wrapping_sub(e as u64);
    let (r, c) = xh.overflowing_sub(b);
    r.wrapping_sub(0u32.wrapping_sub(c as u32) as u64)
}

/// Test of equality between two BaseField elements; return value is
/// 0xFFFFFFFFFFFFFFFF if the two values are equal, or 0 otherwise.
#[inline(always)]
pub fn equals(lhs: u64, rhs: u64) -> u64 {
    let t = lhs ^ rhs;
    !((((t | t.wrapping_neg()) as i64) >> 63) as u64)
}

/// Computes `lhs + rhs mod N`
#[inline(always)]
const fn add(lhs: u64, rhs: u64) -> u64 {
    let (x1, c1) = lhs.overflowing_sub(N - rhs);
    let adj = 0u32.wrapping_sub(c1 as u32);
    x1.wrapping_sub(adj as u64)
}

/// Computes `lhs - rhs mod N`
#[inline(always)]
const fn sub(lhs: u64, rhs: u64) -> u64 {
    let (x1, c1) = lhs.overflowing_sub(rhs);
    let adj = 0u32.wrapping_sub(c1 as u32);
    x1.wrapping_sub(adj as u64)
}

/// Computes `lhs / rhs mod N`
#[inline(always)]
fn div(lhs: u64, rhs: u64) -> u64 {
    mul(lhs, modular_inverse(rhs))
}

/// Irreducible polynomial for degree 4 extension field
fn irreducible_poly() -> Univariate<BaseFelt> {
    // X^3 - X + 1
    Univariate::new(vec![
        BaseFelt::one(),
        BaseFelt::zero() - BaseFelt::one(),
        BaseFelt::zero(),
        BaseFelt::one(),
    ])
}

impl ExtensibleField<4> for BaseFelt {
    fn is_zero(a: [Self; 4]) -> bool {
        a[0].is_zero() && a[1].is_zero() && a[2].is_zero() && a[3].is_zero()
    }

    fn is_one(a: [Self; 4]) -> bool {
        a[0].is_one() && a[1].is_zero() && a[2].is_zero() && a[3].is_zero()
    }

    fn mul(a: [Self; 4], b: [Self; 4]) -> [Self; 4] {
        let a = Univariate::new(a.to_vec());
        let b = Univariate::new(b.to_vec());
        let res = (a * b) % irreducible_poly();
        [
            res.coefficients[0],
            res.coefficients[1],
            res.coefficients[2],
            res.coefficients[3],
        ]
    }

    fn add(a: [Self; 4], b: [Self; 4]) -> [Self; 4] {
        let mut res = [BaseFelt(0); 4];
        for i in 0..4 {
            res[i] = a[i] + b[i];
        }
        res
    }

    fn negate(a: [Self; 4]) -> [Self; 4] {
        a.map(|v| -v)
    }

    fn inverse(_a: [Self; 4]) -> [Self; 4] {
        todo!()
    }
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
        assert_eq!(b.into_bigint(), 10);
    }

    #[test]
    fn large_numbers_can_be_added() {
        // Large compared to u128::MAX
        let a = BaseFelt::new(BaseFelt::MODULUS - 2);

        assert_eq!((a + a).into_bigint(), BaseFelt::MODULUS - 4);
    }

    #[test]
    fn multiplication_multiplies_large_numbers() {
        // Large compared to u128::MAX
        let a = BaseFelt::new(BaseFelt::MODULUS - 2);

        assert_eq!((a * a).into_bigint(), 4);
    }

    #[test]
    fn multiplication_multiplies_odd_numbers() {
        // Large compared to u128::MAX
        let a = BaseFelt::new(BaseFelt::MODULUS - 1);
        let b = BaseFelt::new(BaseFelt::MODULUS - 2);

        assert_eq!((a * b).into_bigint(), 2);
    }

    #[test]
    fn subtraction_subtracts() {
        let a = BaseFelt::new(1);
        let b = BaseFelt::new(2);

        assert_eq!((a - b).into_bigint(), BaseFelt::MODULUS - 1);
    }

    #[test]
    fn division_divides() {
        let two = BaseFelt::new(2);

        assert_eq!((two / two).into_bigint(), 1);
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

        println!("{}", items.iter().product::<BaseFelt>().into_bigint());
    }
}
