#![feature(test, bigint_helper_methods, const_bigint_helper_methods)]

use core::iter::Sum;
use num_traits::One;
use num_traits::Zero;
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
pub(crate) use multivariate::Multivariate;

mod univariate;
pub(crate) use univariate::Univariate;

pub mod fp_u1;
pub mod fp_u128;
pub mod fp_u256;
pub mod fp_u64;

pub trait Felt:
    Copy
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
    /// A positive integer big enough to describe a field modulus for
    /// `Self::BaseField` with no loss of precision.
    /// TODO: remove `From<u128>`
    type PositiveInteger: num_traits::NumAssign + TryFrom<usize> + From<u32> + Display;

    /// Bytes needed to store the field element.
    const ELEMENT_BYTES: usize;

    /// Number of bits needed to represent the field's order.
    const FIELD_ORDER_BITS: u32;

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
    fn pow(self, power: Self::PositiveInteger) -> Self;

    /// Exponentiates this element by a power of the base prime modulus
    fn frobenius(&mut self);

    /// Returns a canonical integer representation of this field element.
    fn as_integer(&self) -> Self::PositiveInteger;
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
pub trait PrimeFelt: Felt {
    /// Prime modulus of the field.
    const MODULUS: Self::PositiveInteger;
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
pub trait ExtensibleField<const N: usize>: PrimeFelt {
    fn is_zero(a: [Self; N]) -> bool;

    fn is_one(a: [Self; N]) -> bool;

    /// Returns a product of `a` and `b`
    fn mul(a: [Self; N], b: [Self; N]) -> [Self; N];

    /// Returns the addition of `a` and `b`
    fn add(a: [Self; N], b: [Self; N]) -> [Self; N];

    /// Returns the subtraction of `a` and `b`
    fn sub(a: [Self; N], b: [Self; N]) -> [Self; N] {
        ExtensibleField::add(a, ExtensibleField::negate(b))
    }

    /// Returns the negation of `a`
    fn negate(a: [Self; N]) -> [Self; N];

    // Returns the inverse of `a`
    fn inverse(a: [Self; N]) -> [Self; N];

    // Return the division result of `a / b`
    fn divide(a: [Self; N], b: [Self; N]) -> [Self; N] {
        assert!(!ExtensibleField::is_zero(b), "divide by zero");
        ExtensibleField::mul(a, ExtensibleField::inverse(b))
    }
}
