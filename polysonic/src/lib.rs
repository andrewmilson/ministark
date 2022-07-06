#![feature(test, int_log, allocator_api)]

pub mod allocator;
pub mod fields;
pub mod plan;
pub mod stage;
pub mod twiddles;
pub mod utils;

use allocator::PageAlignedAllocator;
use core::num;
use fields::fp_u128::BaseFelt;
use fields::Felt;
use fields::PrimeFelt;
use fields::StarkFelt;
use std::marker::PhantomData;

/// Represents a NTT direction
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum NttDirection {
    Forward,
    Inverse,
}

/// A trait that allows NTT algorithms to report whether they compute forward
/// NTTs or inverse NTTs
pub trait Direction {
    /// Returns [Forward](NttDirection::Forward) if this instance computes
    /// forward NTT, or [Inverse](NttDirection::Inverse) for inverse NTTs
    fn ntt_direction(&self) -> NttDirection;
}

#[derive(PartialEq, Eq)]
pub enum NttOrdering {
    BitReversed,
    Natural,
}

pub trait Ordering {
    fn ntt_ordering(&self) -> NttOrdering;
}
