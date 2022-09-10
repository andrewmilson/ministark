#![feature(test, bigint_helper_methods, const_bigint_helper_methods)]

use bigint::BigInteger;
use core::iter::Sum;
use num_bigint::BigUint;
use num_traits::One;
use num_traits::Zero;
use rand::distributions::Standard;
use rand::prelude::Distribution;
use rand::Rng;
use serde::Deserialize;
use serde::Serialize;
use std::fmt::Debug;
use std::fmt::Display;
use std::hash::Hash;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;

mod multivariate;
pub use multivariate::Multivariate;

mod univariate;
pub use univariate::Univariate;

mod bigint;

pub mod fp_u1;
pub mod fp_u128;
pub mod fp_u256;
pub mod fp_u64;

pub trait Felt:
    Copy
    + UniformRand
    + Clone
    + Debug
    + Display
    + Zero
    + One
    + Eq
    + Hash
    + Neg<Output = Self>
    + Sized
    + Serialize
    + for<'a> Deserialize<'a>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<Self>
    + DivAssign<Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + for<'a> DivAssign<&'a Self>
    + Sum<Self>
    + for<'a> Sum<&'a Self>
    + From<usize>
    + From<u128>
    + From<u64>
    + From<u32>
    + From<u16>
    + From<u8>
{
    type BaseFelt: Felt;

    /// Returns `self * self`.
    #[must_use]
    fn square(&self) -> Self {
        *self * self
    }

    /// Squares `self` in place.
    fn square_in_place(&mut self) -> &mut Self {
        *self *= *self;
        self
    }

    /// Returns `self + self`.
    #[must_use]
    fn double(&self) -> Self {
        *self + self
    }

    /// Doubles `self` in place.
    fn double_in_place(&mut self) -> &mut Self {
        *self += *self;
        self
    }

    /// Computes the multiplicative inverse of `self` if it exists.
    ///
    /// It should exist if `self` is non-zero.
    #[must_use]
    fn inverse(&self) -> Option<Self>;

    /// Sets, and returns, `self` to its inverse if it exists. Returns `None`
    /// otherwise.
    fn inverse_in_place(&mut self) -> Option<&mut Self>;

    // TODO: find out if difference in performance if borrowed or owned self.
    fn pow<S: AsRef<[u64]>>(self, exp: S) -> Self {
        let mut res = Self::one();

        for i in BitIterator::new_be(exp) {
            res.square_in_place();

            if i {
                res *= self;
            }
        }

        res
    }

    /// Exponentiates this element by a power of the base prime modulus
    fn frobenius(&mut self);
}

/// Field element that facilitates efficient construction a stark proofs.
///
/// The efficiency comes from the fact the underlying field has a subgroup of
/// order `2^n`. This is required to compute the fast fourier transform (FFT) -
/// also called the number theory transform (NTT) in the case of finite fields.
///
/// Given `n` number of points to interpolate, fast fourier based interpolation
/// runtime is `O(n log n)` as opposed to naive implementations of `O(n^2)`.
pub trait StarkFelt: Felt {
    /// A multiplicative generator of the entire field except 0.
    const GENERATOR: Self;

    // Two adicity would be `n` given a subgroup of order `2^n`.
    const TWO_ADICITY: u32;

    /// A root of unity of the subgroup of order `2^n`.
    const TWO_ADIC_ROOT_OF_UNITY: Self;

    /// Returns a root of unity of order `2^n`.
    fn get_root_of_unity(n: u32) -> Self;
}

/// Prime field element.
pub trait PrimeFelt:
    Felt
    + From<<Self as PrimeFelt>::BigInt>
    + Into<<Self as PrimeFelt>::BigInt>
    + From<BigUint>
    + Into<BigUint>
{
    /// A `BigInteger` type that can represent elements of this field.
    type BigInt: BigInteger;

    /// Prime modulus of the field.
    const MODULUS: Self::BigInt;

    /// Returns a canonical integer representation of this field element.
    fn into_bigint(self) -> Self::BigInt;
}

/// Montgomery batch inversion.
///
/// Implementation reference:
/// - https://books.google.com.au/books?id=kGu4lTznRdgC&pg=PA54 (Sec. 5.3)
/// - https://vitalik.ca/general/2018/07/21/starks_part_3.html
pub fn batch_inverse<E: Felt>(values: &[E]) -> Vec<Option<E>> {
    if values.is_empty() {
        return vec![];
    }

    // compute running multiple of values
    // a, ab, ..., abc...xyz
    let mut accumulator = E::one();
    let mut partials = vec![Some(accumulator)];
    for value in values.iter() {
        if value.is_zero() {
            partials.push(None);
        } else {
            accumulator *= value;
            partials.push(Some(accumulator));
        }
    }

    // With Montgomery method we only need to calculate one inverse
    // 1/abc...xyz
    accumulator.inverse_in_place();

    // Calculate output values (and update the inverse as we go):
    //   - 1/z = abc...xy * 1/abx...xyz
    //   - 1/y = abc...x * 1/abx...xy
    //   - ...
    //   - 1/b = a * 1/ab
    //   - 1/a = 1 * 1/a
    let mut output = vec![None; values.len()];
    for (i, value) in values.iter().enumerate().rev() {
        if !value.is_zero() {
            if let Some(partial) = partials[i] {
                output[i] = Some(partial * accumulator);
                accumulator *= value;
            }
        }
    }

    output
}

/// Degree `N` extension of a prime field
///
/// The irreducible polynomial over the extension field is implicitly defined.
pub trait ExtensibleFelt<const N: usize>: PrimeFelt {
    fn is_zero(a: [Self; N]) -> bool;

    fn is_one(a: [Self; N]) -> bool;

    /// Returns a product of `a` and `b`
    fn mul(a: [Self; N], b: [Self; N]) -> [Self; N];

    /// Returns the addition of `a` and `b`
    fn add(a: [Self; N], b: [Self; N]) -> [Self; N];

    /// Returns the subtraction of `a` and `b`
    fn sub(a: [Self; N], b: [Self; N]) -> [Self; N] {
        ExtensibleFelt::add(a, ExtensibleFelt::negate(b))
    }

    /// Returns the negation of `a`
    fn negate(a: [Self; N]) -> [Self; N];

    // Returns the inverse of `a`
    fn inverse(a: [Self; N]) -> [Self; N];

    // Return the division result of `a / b`
    fn divide(a: [Self; N], b: [Self; N]) -> [Self; N] {
        assert!(!ExtensibleFelt::is_zero(b), "divide by zero");
        ExtensibleFelt::mul(a, ExtensibleFelt::inverse(b))
    }
}

struct BitIterator<S> {
    s: S,
    pos: usize,
    end: usize,
}

impl<S: AsRef<[u64]>> BitIterator<S> {
    // Iterator in big-endian order
    fn new_be(slice: S) -> BitIterator<S> {
        BitIterator {
            pos: slice.as_ref().len() * 64,
            s: slice,
            end: 0,
        }
    }

    // // Iterator in little-endian order
    // fn new_le(slice: S) -> BitIterator<S> {
    //     BitIterator {
    //         end: slice.as_ref().len() * 64,
    //         s: slice,
    //         pos: 0,
    //     }
    // }
}

impl<S: AsRef<[u64]>> Iterator for BitIterator<S> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.pos == self.end {
            None
        } else {
            if self.pos < self.end {
                self.pos += 1;
            } else {
                self.pos -= 1;
            }
            let part = self.pos / 64;
            let bit = self.pos - part * 64;
            Some(self.s.as_ref()[part] & (1 << bit) > 0)
        }
    }
}

/// Trait specifying a field is an extension of another
pub trait ExtensionOf<E: Felt>: Felt + From<E> {
    fn mul_base(self, other: E) -> Self;
}

/// A field is always an extension of itself.
impl<E: Felt> ExtensionOf<E> for E {
    #[inline(always)]
    fn mul_base(self, other: E) -> Self {
        self * other
    }
}

pub trait UniformRand: Sized {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

impl<T> UniformRand for T
where
    Standard: Distribution<T>,
{
    #[inline]
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        rng.sample(Standard)
    }
}
