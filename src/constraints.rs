use crate::expression::Expr;
use crate::utils;
use alloc::collections::BTreeSet;
use ark_ff::One;
use ark_ff::Zero;
use core::iter::Product;
use core::iter::Sum;
use core::ops::Add;
use core::ops::Deref;
use core::ops::DerefMut;
use core::ops::Div;
use core::ops::Mul;
use core::ops::Neg;
use core::ops::Sub;
use num_traits::Pow;
use std::fmt::Debug;
use std::hash::Hash;

// TODO: should really remove copy as this type might change in the future
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum AlgebraicItem<T: 'static> {
    X,
    Constant(T),
    Challenge(usize),
    Periodic(PeriodicColumn<'static, T>),
    Hint(usize),
    Trace(/* =column */ usize, /* =offset */ isize),
}

impl<T> AlgebraicItem<T> {
    // Returns an upper bound on the item's degree in `x`
    const fn degree(&self, trace_degree: usize) -> Degree {
        use AlgebraicItem::*;
        match &self {
            // TODO: handle implications of a zero?
            Constant(_) | Challenge(_) | Hint(_) => Degree(0, 0),
            Trace(_, _) => Degree(trace_degree, 0),
            Periodic(col) => col.degree(trace_degree),
            X => Degree(1, 0),
        }
    }
}

impl<T: Zero> Sum<Self> for Expr<AlgebraicItem<T>> {
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        let zero = AlgebraicItem::Constant(T::zero()).into();
        iter.next()
            .map_or(zero, |expr| iter.fold(expr, |a, b| a + b))
    }
}

impl<T: One> Product<Self> for Expr<AlgebraicItem<T>> {
    fn product<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        let one = AlgebraicItem::Constant(T::one()).into();
        iter.next()
            .map_or(one, |expr| iter.fold(expr, |a, b| a * b))
    }
}

impl<T> Mul<Self> for AlgebraicItem<T> {
    type Output = Expr<Self>;

    fn mul(self, rhs: Self) -> Self::Output {
        Expr::from(self) * Expr::from(rhs)
    }
}

impl<T> Div<Self> for AlgebraicItem<T> {
    type Output = Expr<Self>;

    fn div(self, rhs: Self) -> Self::Output {
        Expr::from(self) / Expr::from(rhs)
    }
}

impl<T> Add<Self> for AlgebraicItem<T> {
    type Output = Expr<Self>;

    fn add(self, rhs: Self) -> Self::Output {
        Expr::from(self) + Expr::from(rhs)
    }
}

impl<T> Sub<Self> for AlgebraicItem<T> {
    type Output = Expr<Self>;

    fn sub(self, rhs: Self) -> Self::Output {
        Expr::from(self) - Expr::from(rhs)
    }
}

impl<T> Pow<usize> for AlgebraicItem<T> {
    type Output = Expr<Self>;

    fn pow(self, rhs: usize) -> Self::Output {
        Expr::from(self).pow(rhs)
    }
}

forward_ref_binop!(impl< T: Clone > Mul, mul for AlgebraicItem<T>, AlgebraicItem<T>);
forward_ref_binop!(impl< T: Clone > Div, div for AlgebraicItem<T>, AlgebraicItem<T>);
forward_ref_binop!(impl< T: Clone > Add, add for AlgebraicItem<T>, AlgebraicItem<T>);
forward_ref_binop!(impl< T: Clone > Sub, sub for AlgebraicItem<T>, AlgebraicItem<T>);

/// A periodic column that repeats itself every `interval_size` many rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PeriodicColumn<'a, T> {
    coeffs: &'a [T],
    interval_size: usize,
}

impl<'a, T> PeriodicColumn<'a, T> {
    /// # Panics
    /// Panics if the number of coefficients or the
    /// interval size is not a power of two.
    pub const fn new(coeffs: &'a [T], interval_size: usize) -> Self {
        assert!(coeffs.len().is_power_of_two());
        assert!(interval_size.is_power_of_two());
        assert!(coeffs.len() <= interval_size);
        Self {
            coeffs,
            interval_size,
        }
    }

    pub const fn interval_size(&self) -> usize {
        self.interval_size
    }

    pub const fn coeffs(&self) -> &'a [T] {
        self.coeffs
    }

    /// Returns an upper bound on the preiodic column's degree in `x`
    const fn degree(&self, trace_degree: usize) -> Degree {
        let trace_len = trace_degree + 1;
        assert!(trace_len.is_power_of_two());
        let poly_degree = self.coeffs.len() - 1;
        let num_intervals = trace_len / self.interval_size;
        Degree(poly_degree * num_intervals, 0)
    }
}

#[derive(Clone)]
pub struct Constraint<T: 'static>(Expr<AlgebraicItem<T>>);

impl<T> Constraint<T> {
    pub const fn new(expression: Expr<AlgebraicItem<T>>) -> Self {
        Self(expression)
    }

    /// Calculates an upper bound on the degree in X.
    /// Output is of the form `(numerator_degree, denominator_degree)`
    pub fn degree(&self, trace_degree: usize) -> (usize, usize) {
        let Degree(numerator_degree, denominator_degree) =
            self.0.eval(&mut |leaf| leaf.degree(trace_degree));
        (numerator_degree, denominator_degree)
    }

    /// Returns the power-of-2 degree blowup observed by evaluating constraints
    /// over the trace polynomials.
    pub fn blowup_factor(&self, trace_len: usize) -> usize {
        let trace_degree = trace_len - 1;
        let (numerator_degree, denominator_degree) = self.degree(trace_degree);
        blowup_factor(numerator_degree, denominator_degree, trace_degree)
    }

    /// Returns the evaluation result if the numerator is 0 when the denominator
    /// is 0 otherwise returns None. This can be used as a heuristic check by
    /// the prover to ensure they have a valid execution trace.
    // Adapted from OpenZKP
    pub fn check(&self, f: &mut impl FnMut(&AlgebraicItem<T>) -> T) -> Option<T>
    where
        T: Zero
            + Neg<Output = T>
            + Add<Output = T>
            + Mul<Output = T>
            + Div<Output = T>
            + Pow<usize, Output = T>,
    {
        pub struct CheckedEval<T: Zero>(Option<T>);

        impl<T: Zero + Neg<Output = T>> Neg for CheckedEval<T> {
            type Output = Self;

            fn neg(self) -> Self::Output {
                Self(self.0.map(|v| -v))
            }
        }

        impl<T: Zero + Add<Output = T>> Add for CheckedEval<T> {
            type Output = Self;

            fn add(self, rhs: Self) -> Self::Output {
                let Self(a) = self;
                let Self(b) = rhs;
                Self(match (a, b) {
                    (Some(a), Some(b)) => Some(a + b),
                    _ => None,
                })
            }
        }

        impl<T: Zero + Div<Output = T>> Div for CheckedEval<T> {
            type Output = Self;

            fn div(self, rhs: Self) -> Self::Output {
                let Self(a) = self;
                let Self(b) = rhs;
                Self(match (a, b) {
                    (Some(a), Some(b)) => {
                        if b.is_zero() && a.is_zero() {
                            Some(T::zero())
                        } else if b.is_zero() {
                            None
                        } else {
                            Some(a / b)
                        }
                    }
                    (Some(a), None) | (None, Some(a)) => a.is_zero().then_some(T::zero()),
                    _ => None,
                })
            }
        }

        impl<T: Zero + Mul<Output = T>> Mul for CheckedEval<T> {
            type Output = Self;

            fn mul(self, rhs: Self) -> Self::Output {
                let Self(a) = self;
                let Self(b) = rhs;
                Self(match (a, b) {
                    (Some(a), Some(b)) => Some(a * b),
                    (Some(x), None) | (None, Some(x)) => x.is_zero().then_some(x),
                    (None, None) => None,
                })
            }
        }

        impl<T: Zero + Pow<usize, Output = T>> Pow<usize> for CheckedEval<T> {
            type Output = Self;

            fn pow(self, rhs: usize) -> Self::Output {
                Self(self.0.map(|v| v.pow(rhs)))
            }
        }

        self.0.eval(&mut |leaf| CheckedEval(Some(f(leaf)))).0
    }

    // Adapted from https://github.com/0xProject/OpenZKP
    pub fn trace_arguments(&self) -> BTreeSet<(usize, isize)> {
        let mut arguments = BTreeSet::new();
        self.traverse(&mut |node| {
            if let &Expr::Leaf(AlgebraicItem::Trace(i, j)) = node {
                arguments.insert((i, j));
            }
        });
        arguments
    }
}

impl<T> From<Expr<AlgebraicItem<T>>> for Constraint<T> {
    fn from(value: Expr<AlgebraicItem<T>>) -> Self {
        Self::new(value)
    }
}

impl<T: 'static> Deref for Constraint<T> {
    type Target = Expr<AlgebraicItem<T>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for Constraint<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum CompositionItem<T: 'static> {
    Item(AlgebraicItem<T>),
    CompositionCoeff(usize),
}

impl<T> CompositionItem<T> {
    // Returns the item's corresponding degree
    const fn degree(&self, trace_degree: usize) -> Degree {
        match &self {
            Self::Item(item) => item.degree(trace_degree),
            Self::CompositionCoeff(_) => Degree(0, 0),
        }
    }
}

impl<T: Zero> Sum<Self> for Expr<CompositionItem<T>> {
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        let zero = Self::Leaf(CompositionItem::Item(AlgebraicItem::Constant(T::zero())));
        iter.next().map_or(zero, |acc| iter.fold(acc, |a, b| a + b))
    }
}

pub struct CompositionConstraint<T: 'static>(Expr<CompositionItem<T>>);

impl<T: Clone + Copy + Zero + Ord + Hash> CompositionConstraint<T> {
    pub const fn new(expression: Expr<CompositionItem<T>>) -> Self {
        Self(expression)
    }

    /// Calculates an upper bound on the degree in X.
    /// Output is of the form `(numerator_degree, denominator_degree)`
    pub fn degree(&self, trace_degree: usize) -> (usize, usize) {
        let Degree(numerator_degree, denominator_degree) =
            self.0.eval(&mut |leaf| leaf.degree(trace_degree));
        (numerator_degree, denominator_degree)
    }

    /// Returns the power-of-2 degree blowup observed by evaluating constraints
    /// over the trace polynomials.
    pub fn blowup_factor(&self, trace_len: usize) -> usize {
        let trace_degree = trace_len - 1;
        let (numerator_degree, denominator_degree) = self.degree(trace_degree);
        blowup_factor(numerator_degree, denominator_degree, trace_degree)
    }
}

impl<T> Deref for CompositionConstraint<T> {
    type Target = Expr<CompositionItem<T>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Returns the power-of-2 degree blowup observed by evaluating constraints
/// over the trace polynomials.
const fn blowup_factor(
    numerator_degree: usize,
    denominator_degree: usize,
    trace_degree: usize,
) -> usize {
    let degree = numerator_degree.saturating_sub(denominator_degree);
    utils::ceil_power_of_two(degree) / trace_degree
}

pub trait Hint {
    fn index(&self) -> usize;

    fn hint<T>(&self) -> Expr<AlgebraicItem<T>> {
        AlgebraicItem::Hint(self.index()).into()
    }
}

impl Hint for usize {
    fn index(&self) -> usize {
        *self
    }
}

pub trait VerifierChallenge {
    /// Get the challenge index
    fn index(&self) -> usize;

    /// Symbolic representation of a challenge
    // TODO: terrible name. Needs refactoring
    fn challenge<T>(&self) -> Expr<AlgebraicItem<T>> {
        AlgebraicItem::Challenge(self.index()).into()
    }
}

impl VerifierChallenge for usize {
    fn index(&self) -> usize {
        *self
    }
}

/// An interface for types that can symbolically represent a column of an
/// execution trace
pub trait ExecutionTraceColumn {
    /// Returns the execution trace column index
    fn index(&self) -> usize;

    // Create a constraint element for the current cycle
    fn curr<T>(&self) -> Expr<AlgebraicItem<T>> {
        self.offset(0)
    }

    // Create a constraint element for the next cycle
    fn next<T>(&self) -> Expr<AlgebraicItem<T>> {
        self.offset(1)
    }

    fn offset<T>(&self, offset: isize) -> Expr<AlgebraicItem<T>> {
        AlgebraicItem::Trace(self.index(), offset).into()
    }
}

impl ExecutionTraceColumn for usize {
    fn index(&self) -> usize {
        *self
    }
}

/// Degree of the form `(numerator_degree, denominator_degree)`
struct Degree(pub usize, pub usize);

impl Neg for Degree {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self
    }
}

impl Add for Degree {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let Self(an, ad) = self;
        let Self(bn, bd) = rhs;
        Self((an + bd).max(bn + ad), ad + bd)
    }
}

impl Div for Degree {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let Self(an, ad) = self;
        let Self(bn, bd) = rhs;
        Self(an + bd, ad + bn)
    }
}

impl Mul for Degree {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let Self(an, ad) = self;
        let Self(bn, bd) = rhs;
        Self(an + bn, ad + bd)
    }
}

impl Pow<usize> for Degree {
    type Output = Self;

    fn pow(self, rhs: usize) -> Self::Output {
        let Self(n, d) = self;
        Self(n * rhs, d * rhs)
    }
}
