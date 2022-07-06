use core::iter::Sum;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::{
    Add, AddAssign, BitAnd, Div, DivAssign, Mul, MulAssign, Neg, Shl, ShrAssign, Sub, SubAssign,
};

pub trait FieldElement:
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
    + From<u128>
    + From<u64>
    + From<u32>
    + From<u16>
    + From<u8>
{
    /// A type defining positive integers big enough to describe a field modulus for
    /// `Self::BaseField` with no loss of precision.
    type PositiveInteger: Debug
        + Copy
        + PartialEq
        + PartialOrd
        + ShrAssign
        + Shl<u32, Output = Self::PositiveInteger>
        + BitAnd<Output = Self::PositiveInteger>
        + One
        + Zero
        + From<u32>
        + From<u64>
        + From<u128>;

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

    /// Sets, and returns, `self` to its inverse if it exists. Returns `None` otherwise.
    fn inverse_in_place(&mut self) -> Option<&mut Self>;

    /// Montgomery batch inversion.
    ///
    /// Implementation reference:
    /// - https://books.google.com.au/books?id=kGu4lTznRdgC&pg=PA54&lpg=PA54&dq=montgomery+batch+inversion&source=bl&ots=tPJcPPOrCe&sig=Z3p_6YYwYloRU-f1K-nnv2D8lGw&hl=en&sa=X&redir_esc=y#v=onepage&q=montgomery%20batch%20inversion&f=false
    /// - https://vitalik.ca/general/2018/07/21/starks_part_3.html
    fn batch_inverse(values: &[Self]) -> Vec<Option<Self>> {
        if values.is_empty() {
            return vec![];
        }

        // compute running multiple of values
        // a, ab, ..., abc...xyz
        let mut accumulator = Self::one();
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

    /// Returns a canonical integer representation of this field element.
    fn as_integer(&self) -> Self::PositiveInteger;
}

/// Field element that facilitates efficient construction a stark proofs.
///
/// The efficiency comes from the fact the underlying field has a subgroup of order `2^n`.
/// This is required to compute fast fourier transform (FFT) also called the number theory transform (NTT)
/// in the case of finite fields.
///
/// Given `n` number of points to interpolate, fast fourier based interpolation runtime is `O(n log n)`
/// as opposed to naive implementations of `O(n^2)`.
pub trait StarkElement: FieldElement {
    /// A multiplicative generator of the entire field except 0.
    const GENERATOR: Self;

    // Two adicity would be `n` given a subgroup of order `2^n`.
    const TWO_ADICITY: u32;

    /// A root of unity of the subgroup of order `2^n`.
    const TWO_ADIC_ROOT_OF_UNITY: Self;

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

pub trait PrimeField: FieldElement {
    /// Prime modulus of the field.
    const MODULUS: Self::PositiveInteger;
}
