//! TODO
extern crate test;

use super::Felt;
use super::PrimeFelt;
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
const N: u128 = 2;

/// Represents a base field element for the prime field with modulus `2`
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct BaseFelt(bool);

impl BaseFelt {
    pub const fn new(value: usize) -> BaseFelt {
        // Convert to Montgomery form
        BaseFelt(value & 1 == 1)
    }
}

impl Felt for BaseFelt {
    type PositiveInteger = PositiveInteger;

    const ELEMENT_BYTES: usize = core::mem::size_of::<bool>();

    const FIELD_ORDER_BITS: u32 = 1;

    fn inverse(&self) -> Option<Self> {
        match self.0 {
            true => Some(*self),
            false => None,
        }
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
        *self
    }

    fn double_in_place(&mut self) -> &mut Self {
        *self = self.double();
        self
    }

    fn as_integer(&self) -> Self::PositiveInteger {
        self.0.into()
    }

    fn pow(self, power: Self::PositiveInteger) -> Self {
        if power.is_zero() {
            Self::one()
        } else {
            self
        }
    }

    // Computes the identity in a prime field
    fn frobenius(&mut self) {}
}

impl Display for BaseFelt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_integer())
    }
}

impl From<u8> for BaseFelt {
    fn from(item: u8) -> Self {
        BaseFelt::new(item as usize)
    }
}

impl From<u16> for BaseFelt {
    fn from(item: u16) -> Self {
        BaseFelt::new(item as usize)
    }
}

impl From<u32> for BaseFelt {
    fn from(item: u32) -> Self {
        BaseFelt::new(item as usize)
    }
}

impl From<u64> for BaseFelt {
    fn from(item: u64) -> Self {
        BaseFelt::new(item as usize)
    }
}

impl From<u128> for BaseFelt {
    fn from(item: u128) -> Self {
        BaseFelt::new(item as usize)
    }
}

impl From<usize> for BaseFelt {
    fn from(item: usize) -> Self {
        BaseFelt::new(item)
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

impl PrimeFelt for BaseFelt {
    const MODULUS: PositiveInteger = N;
}

impl One for BaseFelt {
    #[inline]
    fn one() -> Self {
        BaseFelt(true)
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.0 == Self::one().0
    }
}

impl Zero for BaseFelt {
    #[inline]
    fn zero() -> Self {
        BaseFelt(false)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == Self::zero().0
    }
}

impl Neg for BaseFelt {
    type Output = Self;

    fn neg(self) -> Self::Output {
        BaseFelt(!self.0)
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

/// Computes `lhs * rhs mod M`
#[inline(always)]
const fn mul(lhs: bool, rhs: bool) -> bool {
    lhs & rhs
}

/// Computes `lhs + rhs mod N`
#[inline(always)]
const fn add(lhs: bool, rhs: bool) -> bool {
    lhs ^ rhs
}

/// Computes `lhs - rhs mod N`
#[inline(always)]
const fn sub(lhs: bool, rhs: bool) -> bool {
    add(lhs, rhs)
}

/// Computes `lhs / rhs mod N`
#[inline(always)]
fn div(lhs: bool, rhs: bool) -> bool {
    mul(lhs, rhs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add() {
        let a = BaseFelt::new(1);

        assert_eq!((a + a).as_integer(), 0);
        assert_eq!((BaseFelt::zero() + a).as_integer(), 1);
    }

    #[test]
    fn subtraction() {
        let a = BaseFelt::new(1);

        assert_eq!((a - a).as_integer(), 0);
        assert_eq!((BaseFelt::zero() - a).as_integer(), 1);
    }

    #[test]
    fn division_divides() {
        assert_eq!((BaseFelt::one() / BaseFelt::one()).as_integer(), 1);
    }
}
