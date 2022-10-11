//! Implementation is adapted from the multivariable polynomial in arkworks.

use crate::Matrix;
use ark_ff::One;
use ark_ff::Zero;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::plan::PLANNER;
use fast_poly::stage::MulPowStage;
use fast_poly::utils::buffer_mut_no_copy;
use fast_poly::utils::buffer_no_copy;
use fast_poly::GpuField;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

/// A constraint element can represent several things:
/// - a column in the current cycle
/// - a column in the next cycle
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Element {
    Curr(usize),
    Next(usize),
    Challenge(usize),
}

impl<F: GpuField> From<Element> for Constraint<F> {
    fn from(element: Element) -> Self {
        Constraint::new(vec![Term::new(
            F::one(),
            Variables::new(vec![(element, 1)]),
        )])
    }
}

pub trait Challenge {
    /// Get the challenge index
    fn index(&self) -> usize;

    /// Symbolic representation of a challenge
    // TODO: terrible name. Needs refactoring
    fn get_challenge<F: GpuField>(&self) -> Constraint<F> {
        Constraint::from(Element::Challenge(self.index()))
    }
}

impl Challenge for usize {
    fn index(&self) -> usize {
        *self
    }
}

/// An interface for types that can symbolically represent a column of an
/// execution trace
pub trait Column {
    /// Get the execution trace column index
    fn index(&self) -> usize;

    // Create a constraint element for the current cycle
    fn curr<F: GpuField>(&self) -> Constraint<F> {
        Constraint::from(Element::Curr(self.index()))
    }

    // Create a constraint element for the next cycle
    fn next<F: GpuField>(&self) -> Constraint<F> {
        Constraint::from(Element::Next(self.index()))
    }
}

impl Column for usize {
    fn index(&self) -> usize {
        *self
    }
}

/// Represents the group of variables within a constraint polynomial term.
/// Each variable is of the form `(element, power)`.
#[derive(Clone, PartialEq, Eq, Default)]
struct Variables(Vec<(Element, usize)>);

impl Variables {
    /// Create a new group of variables
    fn new(mut variables: Vec<(Element, usize)>) -> Self {
        variables.retain(|(_, pow)| *pow != 0);
        variables.sort();
        Variables(Self::combine(&variables))
    }

    /// Returns the combined degree of all variables
    fn degree(&self) -> usize {
        self.0
            .iter()
            .filter(|element| match element.0 {
                // challenges are just delayed constants so they don't contribute to the degree
                Element::Challenge(_) => false,
                Element::Curr(_) | Element::Next(_) => true,
            })
            .fold(0, |sum, element| sum + element.1)
    }

    /// Sums the powers of any duplicate variables.
    /// Assumes variables are sorted by element.
    fn combine(variables: &[(Element, usize)]) -> Vec<(Element, usize)> {
        let mut combined_variables: Vec<(Element, usize)> = Vec::new();
        for (curr_element, curr_exponent) in variables {
            if let Some((prev_element, prev_exponent)) = combined_variables.last_mut() {
                if prev_element == curr_element {
                    *prev_exponent += curr_exponent;
                    continue;
                }
            }

            combined_variables.push((*curr_element, *curr_exponent));
        }
        combined_variables
    }

    fn get_elements(&self) -> Vec<Element> {
        self.0.iter().map(|element| element.0).collect()
    }
}

impl core::fmt::Debug for Variables {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> Result<(), core::fmt::Error> {
        for variable in self.0.iter() {
            let power = variable.1;
            let (symbol, index) = match variable.0 {
                Element::Curr(index) => ("x", index),
                Element::Next(index) => ("x'", index),
                Element::Challenge(index) => ("c", index),
            };

            if power.is_one() {
                write!(f, " * {symbol}_{index}")?;
            } else {
                write!(f, " * {symbol}_{index}^{power}")?;
            }
        }
        Ok(())
    }
}

impl PartialOrd for Variables {
    /// Sort by total degree. If total degree is equal then ordering
    /// is given by exponent weight in lower-numbered variables
    /// ie. `e1 > e2`, `e1^2 > e1 * e2`, etc.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.degree() != other.degree() {
            Some(self.degree().cmp(&other.degree()))
        } else {
            // Iterate through all variables and return the corresponding ordering
            // if they differ in variable numbering or power
            for (curr, other) in self.0.iter().zip(&other.0) {
                if other.0 == curr.0 {
                    if curr.1 != other.1 {
                        return Some((curr.1).cmp(&other.1));
                    }
                } else {
                    return Some((other.0).cmp(&curr.0));
                }
            }
            Some(Ordering::Equal)
        }
    }
}

impl Ord for Variables {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

/// Represents a constraint polynomial term.
/// A term is of the form `(coefficient, variables)`.
#[derive(Clone, PartialEq, Eq)]
struct Term<F>(F, Variables);

impl<F: GpuField> Term<F> {
    fn new(coefficient: F, variables: Variables) -> Self {
        Term(coefficient, variables)
    }

    fn evaluate_challenges(&self, challenges: &[F]) -> Self {
        let mut new_coefficient = self.0;
        let mut new_variables = Vec::new();
        // TODO: could turn variables into an itterator
        for variable in &(self.1).0 {
            match variable {
                (Element::Challenge(index), power) => {
                    new_coefficient *= challenges[*index].pow([*power as u64])
                }
                other => new_variables.push(*other),
            }
        }
        Term(new_coefficient, Variables(new_variables))
    }

    fn degree(&self) -> usize {
        self.1.degree()
    }
}

impl<'a, 'b, F: GpuField> Mul<&'a Term<F>> for &'b Term<F> {
    type Output = Term<F>;

    fn mul(self, rhs: &'a Term<F>) -> Self::Output {
        let vars = Variables::new(vec![(self.1).0.clone(), (rhs.1).0.clone()].concat());
        let coeff = self.0 * rhs.0;
        Term::new(coeff, vars)
    }
}

/// A multivariate constraint polynomial
#[derive(Clone)]
pub struct Constraint<F>(Vec<Term<F>>);

impl<F: GpuField> Constraint<F> {
    pub fn get_elements(&self) -> Vec<Element> {
        let mut indices = self
            .0
            .iter()
            .flat_map(|term| term.1.get_elements())
            .collect::<Vec<Element>>();
        indices.sort();
        indices.dedup();
        indices
    }

    fn new(mut terms: Vec<Term<F>>) -> Self {
        terms.sort_by(|a, b| a.1.cmp(&b.1));

        let mut constraint = Constraint(Self::combine_terms(terms));
        constraint.remove_zeros();
        constraint
    }

    pub fn degree(&self) -> usize {
        self.0.iter().map(|term| term.degree()).max().unwrap_or(0)
    }

    fn remove_zeros(&mut self) {
        self.0.retain(|Term(coeff, _)| !coeff.is_zero());
    }

    fn combine_terms(terms: Vec<Term<F>>) -> Vec<Term<F>> {
        let mut combined_terms = Vec::new();
        for Term(curr_coeff, curr_vars) in terms {
            if let Some(Term(prev_coeff, prev_vars)) = combined_terms.last_mut() {
                if prev_vars == &curr_vars {
                    *prev_coeff += curr_coeff;
                    continue;
                }
            }

            combined_terms.push(Term(curr_coeff, curr_vars));
        }
        combined_terms
    }

    fn evaluate_challenges(&self, challenges: &[F]) -> Self {
        Constraint::new(
            self.0
                .iter()
                .map(|term| term.evaluate_challenges(challenges))
                .collect(),
        )
    }

    pub fn evaluate(&self, challenges: &[F], current_row: &[F], next_row: &[F]) -> F {
        let mut result = F::zero();
        for Term(coeff, vars) in self.0.iter() {
            let mut scratch = *coeff;
            for &(element, power) in &vars.0 {
                let val = match element {
                    Element::Curr(index) => current_row[index],
                    Element::Next(index) => next_row[index],
                    Element::Challenge(index) => challenges[index],
                };
                scratch *= val.pow([power as u64]);
            }
            result += scratch;
        }
        result
    }

    pub fn evaluate_symbolic(
        &self,
        challenges: &[F],
        trace_step: usize,
        lde_matrix: &Matrix<F>,
    ) -> Matrix<F> {
        let n = lde_matrix.num_rows();
        let constraint_without_challenges = self.evaluate_challenges(challenges).0;
        assert!(!constraint_without_challenges.is_empty());
        let library = &PLANNER.library;
        let command_queue = &PLANNER.command_queue;
        let command_buffer = command_queue.new_command_buffer();
        let mut term_evaluations = Vec::new();
        let curr_multiplier = MulPowStage::<F>::new(library, n, 0);
        let next_multiplier = MulPowStage::<F>::new(library, n, trace_step);
        for term in constraint_without_challenges {
            let mut scratch = Vec::with_capacity_in(n, PageAlignedAllocator);
            scratch.resize(n, term.0);
            let mut scratch_buffer = buffer_mut_no_copy(command_queue.device(), &mut scratch);
            for (element, power) in &(term.1).0 {
                let (col_index, multiplier) = match element {
                    Element::Curr(col_index) => (col_index, &curr_multiplier),
                    Element::Next(col_index) => (col_index, &next_multiplier),
                    _ => unreachable!(),
                };
                let column_buffer = buffer_no_copy(command_queue.device(), &lde_matrix[*col_index]);
                multiplier.encode(command_buffer, &mut scratch_buffer, &column_buffer, *power);
            }
            term_evaluations.push(scratch);
        }
        command_buffer.commit();
        command_buffer.wait_until_completed();
        Matrix::new(term_evaluations).sum_columns()
    }
}

impl<F: GpuField> core::fmt::Debug for Constraint<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> Result<(), core::fmt::Error> {
        // for (coeff, term) in self.0.iter().filter(|(c, _)| !c.is_zero()) {
        for Term(coeff, variables) in self.0.iter() {
            if variables.0.is_empty() {
                write!(f, "\n{}", coeff)?;
            } else {
                write!(f, "\n{}{:?}", coeff, variables)?;
            }
        }
        Ok(())
    }
}

impl<F: GpuField> From<F> for Constraint<F> {
    fn from(element: F) -> Self {
        Constraint(vec![Term::new(element, Variables::default())])
    }
}

impl<F: GpuField> Zero for Constraint<F> {
    /// Returns the zero polynomial.
    fn zero() -> Self {
        Self(vec![])
    }

    /// Checks if the given polynomial is zero.
    fn is_zero(&self) -> bool {
        self.0.is_empty() || self.0.iter().all(|term| term.0.is_zero())
    }
}

impl<F: GpuField> Mul<&Constraint<F>> for &Constraint<F> {
    type Output = Constraint<F>;

    fn mul(self, rhs: &Constraint<F>) -> Self::Output {
        // Perform a naive n^2 multiplication of `lhs` by `rhs`.
        if self.is_zero() || rhs.is_zero() {
            Self::Output::zero()
        } else {
            let mut result_terms = Vec::new();
            for lhs_term in self.0.iter() {
                for rhs_term in rhs.0.iter() {
                    result_terms.push(lhs_term * rhs_term);
                }
            }
            // Constraint::new(result_terms)
            Constraint::new(result_terms)
        }
    }
}

impl<F: GpuField> Mul<Constraint<F>> for Constraint<F> {
    type Output = Constraint<F>;

    #[inline]
    fn mul(self, other: Constraint<F>) -> Constraint<F> {
        Mul::mul(&self, &other)
    }
}

impl<F: GpuField> Add<&Constraint<F>> for &Constraint<F> {
    type Output = Constraint<F>;

    fn add(self, rhs: &Constraint<F>) -> Self::Output {
        let mut result_terms = Vec::new();
        let mut lhs_iter = self.0.iter().peekable();
        let mut rhs_iter = rhs.0.iter().peekable();
        // Iterate over the constraints in ascending order and combine any common terms.
        loop {
            // Peek at iterators to decide which to take from
            let which = match (lhs_iter.peek(), rhs_iter.peek()) {
                (Some(cur), Some(other)) => Some((cur.1).cmp(&other.1)),
                (Some(_), None) => Some(Ordering::Less),
                (None, Some(_)) => Some(Ordering::Greater),
                (None, None) => None,
            };
            // Push the smallest element to the `result` coefficient vec
            let smallest = match which {
                Some(Ordering::Less) => lhs_iter.next().unwrap().clone(),
                Some(Ordering::Equal) => {
                    let rhs = rhs_iter.next().unwrap();
                    let lhs = lhs_iter.next().unwrap();
                    let new_coeff = rhs.0 + lhs.0;
                    // TODO: could move this to a remove_zeros call
                    if new_coeff.is_zero() {
                        continue;
                    }
                    Term(rhs.0 + lhs.0, rhs.1.clone())
                }
                Some(Ordering::Greater) => rhs_iter.next().unwrap().clone(),
                None => break,
            };
            result_terms.push(smallest);
        }
        Constraint(result_terms)
    }
}

impl<F: GpuField> Add<Constraint<F>> for Constraint<F> {
    type Output = Constraint<F>;

    #[inline]
    fn add(self, rhs: Constraint<F>) -> Self::Output {
        &self + &rhs
    }
}

impl<F: GpuField> Neg for Constraint<F> {
    type Output = Constraint<F>;

    fn neg(mut self) -> Self::Output {
        for Term(coeff, _) in &mut self.0 {
            *coeff = -*coeff;
        }
        self
    }
}

impl<F: GpuField> Neg for &Constraint<F> {
    type Output = Constraint<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        self.clone().neg()
    }
}

impl<F: GpuField> Sub<&Constraint<F>> for &Constraint<F> {
    type Output = Constraint<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: &Constraint<F>) -> Self::Output {
        self + rhs.neg()
    }
}

impl<F: GpuField> Sub<Constraint<F>> for Constraint<F> {
    type Output = Constraint<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Constraint<F>) -> Self::Output {
        self + rhs.neg()
    }
}

impl<F: GpuField> Mul<F> for Constraint<F> {
    type Output = Constraint<F>;

    fn mul(mut self, rhs: F) -> Self::Output {
        for Term(coeff, _) in &mut self.0 {
            *coeff *= rhs;
        }
        self
    }
}

impl<F: GpuField> Mul<&F> for &Constraint<F> {
    type Output = Constraint<F>;

    fn mul(self, rhs: &F) -> Self::Output {
        self.clone() * *rhs
    }
}

impl<F: GpuField> Add<F> for Constraint<F> {
    type Output = Constraint<F>;

    fn add(self, rhs: F) -> Self::Output {
        self + Constraint::from(rhs)
    }
}

impl<F: GpuField> Add<&F> for &Constraint<F> {
    type Output = Constraint<F>;

    fn add(self, rhs: &F) -> Self::Output {
        self + Constraint::from(*rhs)
    }
}

impl<F: GpuField> Sub<F> for Constraint<F> {
    type Output = Constraint<F>;

    fn sub(self, rhs: F) -> Self::Output {
        self + Constraint::from(-rhs)
    }
}

impl<F: GpuField> Sub<&F> for &Constraint<F> {
    type Output = Constraint<F>;

    fn sub(self, rhs: &F) -> Self::Output {
        self + Constraint::from(-*rhs)
    }
}

// Adapted from the `forward_ref_binop!` macro in the Rust standard library.
// Implements "&T op U", "T op &U" based on "&T op &U".
macro_rules! forward_ref_binop {
    (impl<F: GpuField> $imp:ident, $method:ident for $t:ty, $u:ty) => {
        impl<'a, F: GpuField> $imp<$u> for &'a $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, other: $u) -> <&'a $t as $imp<&'a $u>>::Output {
                $imp::$method(self, &other)
            }
        }

        impl<'a, F: GpuField> $imp<&'a $u> for $t {
            type Output = <&'a $t as $imp<&'a $u>>::Output;

            #[inline]
            fn $method(self, other: &$u) -> <&'a $t as $imp<&'a $u>>::Output {
                $imp::$method(&self, other)
            }
        }
    };
}

forward_ref_binop!(impl<F: GpuField> Mul, mul for Constraint<F>, Constraint<F>);
forward_ref_binop!(impl<F: GpuField> Add, add for Constraint<F>, Constraint<F>);
forward_ref_binop!(impl<F: GpuField> Sub, sub for Constraint<F>, Constraint<F>);
forward_ref_binop!(impl<F: GpuField> Mul, mul for Constraint<F>, F);
forward_ref_binop!(impl<F: GpuField> Add, add for Constraint<F>, F);
forward_ref_binop!(impl<F: GpuField> Sub, sub for Constraint<F>, F);

pub fn are_eq<F: GpuField>(
    a: impl Borrow<Constraint<F>>,
    b: impl Borrow<Constraint<F>>,
) -> Constraint<F> {
    a.borrow() - b.borrow()
}

/// Returns zero only when a == zero.
pub fn is_zero<F: GpuField, C: Borrow<Constraint<F>>>(a: C) -> C {
    a
}

/// Returns zero only when a == one.
pub fn is_one<F: GpuField>(a: impl Borrow<Constraint<F>>) -> Constraint<F> {
    a.borrow() - F::one()
}

/// Returns zero only when a = zero || a == one.
pub fn is_binary<F: GpuField>(a: impl Borrow<Constraint<F>>) -> Constraint<F> {
    a.borrow() * a.borrow() - a.borrow()
}
