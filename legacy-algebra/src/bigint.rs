use core::fmt::Debug;
use num_bigint::BigUint;
use std::fmt::Display;

// Taken from the arkwork library
/// This defines a `BigInteger`, a smart wrapper around a
/// sequence of `u64` limbs, least-significant limb first.
// TODO: get rid of this trait once we can use associated constants in const
// generics.
pub trait BigInteger:
    Copy
    + Clone
    + Debug
    + Default
    + Display
    + Eq
    + Ord
    + Sized
    + 'static
    + From<u64>
    + From<u32>
    + From<u16>
    + From<u8>
    + Into<BigUint>
{
}

impl BigInteger for u128 {}
impl BigInteger for u64 {}
