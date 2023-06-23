use crate::expression::Expr;
use crate::utils;
use alloc::collections::BTreeSet;
use ark_ff::One;
use ark_ff::Zero;
use core::borrow::Borrow;
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
use std::hash::Hash;
use std::time::Instant;

// TODO: should really remove copy as this type might change in the future
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum AlgebraicItem<T> {
    X,
    Constant(T),
    Challenge(usize),
    Hint(usize),
    Trace(/* =column */ usize, /* =offset */ isize),
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

#[derive(Clone)]
pub struct Constraint<T>(pub Expr<AlgebraicItem<T>>);

impl<T> Constraint<T> {
    pub const fn new(expression: Expr<AlgebraicItem<T>>) -> Self {
        Self(expression)
    }

    /// Calculates an upper bound on the degree in X.
    /// Output is of the form `(numerator_degree, denominator_degree)`
    pub fn degree(&self, trace_degree: usize) -> (usize, usize) {
        /// Degree of the form `(numerator_degree, denominator_degree)`
        pub struct Degree(pub usize, pub usize);

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

        use AlgebraicItem::*;
        let Degree(numerator_degree, denominator_degree) = self.0.eval(
            // map the leaf nodes to their corresponding degree
            &mut |leaf| match leaf {
                // TODO: handle implications of a zero?
                Constant(_) | Challenge(_) | Hint(_) => Degree(0, 0),
                Trace(_, _) => Degree(trace_degree, 0),
                X => Degree(1, 0),
            },
        );

        (numerator_degree, denominator_degree)
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

impl<T> Deref for Constraint<T> {
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
pub enum CompositionItem<T> {
    Item(AlgebraicItem<T>),
    CompositionCoeff(usize),
}

impl<T: Zero> Sum<Self> for Expr<CompositionItem<T>> {
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        let zero = Self::Leaf(CompositionItem::Item(AlgebraicItem::Constant(T::zero())));
        iter.next().map_or(zero, |acc| iter.fold(acc, |a, b| a + b))
    }
}

#[derive(Clone)]
pub struct CompositionConstraint<T> {
    ce_blowup_factor: usize,
    expr: Expr<CompositionItem<T>>,
}

impl<T: Clone + Copy + Zero + Ord + Hash> CompositionConstraint<T> {
    /// Combines multiple constraints into a single constraint (the composition
    /// constraint). Constraints are composed with verifiers randomness.
    /// This verifier randomness is expressed symbolically.
    /// <https://medium.com/starkware/starkdex-deep-dive-the-stark-core-engine-497942d0f0ab>
    pub fn new(constraints: &[Constraint<T>], trace_len: usize) -> Self {
        let ce_blowup_factor = ce_blowup_factor(constraints, trace_len);
        let composition_degree = trace_len * ce_blowup_factor - 1;
        let trace_degree = trace_len - 1;
        let x = Expr::Leaf(CompositionItem::Item(AlgebraicItem::X));
        let mut composition_coeff = (0..).map(|i| Expr::Leaf(CompositionItem::CompositionCoeff(i)));
        let expr = constraints
            .iter()
            .map(|constraint| {
                let (numerator_degree, denominator_degree) = constraint.degree(trace_degree);
                let evaluation_degree = numerator_degree - denominator_degree;
                assert!(evaluation_degree <= composition_degree);
                let degree_adjustment = composition_degree - evaluation_degree;
                // TODO: if degree_adjustment is 0 then we only need one challenge
                let constraint = constraint.map_leaves(&mut |&leaf| CompositionItem::Item(leaf));
                let alpha = composition_coeff.next().unwrap();
                let beta = composition_coeff.next().unwrap();
                &constraint * (x.clone().pow(degree_adjustment) * alpha + beta)
            })
            .sum::<Expr<CompositionItem<T>>>();
        let now = Instant::now();
        let expr = expr.reuse_shared_nodes();
        println!("Reuse took: {:?}", now.elapsed());
        Self {
            ce_blowup_factor,
            expr,
        }
    }

    /// Returns the power-of-2 degree blowup observed by evaluating constraints
    /// over the trace polynomials.
    pub const fn ce_blowup_factor(&self) -> usize {
        self.ce_blowup_factor
    }
}

impl<T> Deref for CompositionConstraint<T> {
    type Target = Expr<CompositionItem<T>>;

    fn deref(&self) -> &Self::Target {
        &self.expr
    }
}

/// Constraint evaluation blowup factor
fn ce_blowup_factor<T>(constraints: &[Constraint<T>], trace_len: usize) -> usize {
    assert!(trace_len.is_power_of_two());
    let trace_degree = trace_len - 1;
    utils::ceil_power_of_two(
        constraints
            .iter()
            .map(|constraint| {
                let (numerator_degree, denominator_degree) =
                    constraint.borrow().degree(trace_degree);
                numerator_degree.saturating_sub(denominator_degree)
            })
            .max()
            // TODO: ceil_power_of_two might not be correct here. check the math
            .map_or(0, |degree| utils::ceil_power_of_two(degree) / trace_degree),
    )
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
