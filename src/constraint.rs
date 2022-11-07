//! Implementation is adapted from the multivariable polynomial in arkworks.

use ark_ff::One;
use ark_ff::Zero;
use gpu_poly::prelude::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;

/// A constraint element can represent several things:
/// - a column in the current cycle
/// - a column in the next cycle
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Element {
    Curr(usize),
    Next(usize),
    Challenge(usize),
}

impl Element {
    pub fn pow<F: GpuField>(&self, exponent: usize) -> Constraint<F> {
        Constraint::new(vec![Term(
            F::one(),
            Variables::new(vec![(*self, exponent)]),
        )])
    }
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
    /// Returns the execution trace column index
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

// #[proc_macro_derive(Column, attributes(after))]
// pub fn derive_constraint_column(_item: TokenStream) -> TokenStream {
//     "fn answer() -> u32 { 42 }".parse().unwrap()
// }

/// Represents the group of variables within a constraint polynomial term.
/// Each variable is of the form `(element, power)`.
#[derive(Clone, PartialEq, Eq, Default)]
pub(crate) struct Variables(pub(crate) Vec<(Element, usize)>);

impl Variables {
    /// Create a new group of variables
    pub(crate) fn new(mut variables: Vec<(Element, usize)>) -> Self {
        variables.retain(|(_, pow)| *pow != 0);
        variables.sort();
        Variables(Self::combine(&variables))
    }

    /// Returns the combined degree of all variables
    fn degree(&self, include_challenges: bool) -> usize {
        self.0
            .iter()
            .filter(|element| match element.0 {
                // challenges are symbolic constants so they don't necessarily
                // contribute to the degree.
                Element::Challenge(_) => include_challenges,
                _ => true,
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
            write!(f, " * ")?;
            let power = variable.1;
            match variable.0 {
                Element::Curr(index) => write!(f, "x_{}", index)?,
                Element::Next(index) => write!(f, "x'_{}", index)?,
                Element::Challenge(index) => write!(f, "c_{}", index)?,
            };
            if !power.is_one() {
                write!(f, "^{power}")?;
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
        if self.degree(true) != other.degree(true) {
            Some(self.degree(true).cmp(&other.degree(true)))
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
pub struct Term<F>(pub(crate) F, pub(crate) Variables);

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
        self.1.degree(false)
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
pub struct Constraint<F>(pub(crate) Vec<Term<F>>);

impl<F: GpuField> Constraint<F> {
    pub fn pow(&self, mut exp: usize) -> Self {
        let mut res = Constraint::from(F::one());
        let mut acc = self.clone();
        while exp > 0 {
            if exp & 1 == 1 {
                res = &res * &acc;
            }
            acc = &acc * &acc;
            exp >>= 1;
        }
        res
    }

    pub fn get_elements(&self) -> Vec<Element> {
        let mut indices: Vec<Element> = self
            .0
            .iter()
            .flat_map(|term| term.1.get_elements())
            .collect();
        indices.sort();
        indices.dedup();
        indices
    }

    pub fn new(mut terms: Vec<Term<F>>) -> Self {
        terms.sort_by(|a, b| a.1.cmp(&b.1));

        let mut constraint = Constraint(Self::combine_terms(terms));
        constraint.remove_zeros();
        constraint
    }

    /// Substitutes an element with a constraint
    /// E.g. substituting `x = y^2 + y` into `3 * xy` gives `3 * y^3 + 3 * y^2`
    pub fn substitute(&mut self, element: Element, substitution: &Constraint<F>) {
        *self = Self::new(
            self.0
                .iter()
                .cloned()
                .flat_map(|term| {
                    let mut substitution_power = 0;

                    for variable in &(term.1).0 {
                        if variable.0 == element {
                            substitution_power = variable.1;
                            break;
                        }
                    }

                    if substitution_power != 0 {
                        let new_term = Term(
                            term.0,
                            Variables::new(
                                (term.1).0.into_iter().filter(|v| v.0 != element).collect(),
                            ),
                        );
                        (Constraint::new(vec![new_term]) * substitution.pow(substitution_power)).0
                    } else {
                        vec![term]
                    }
                })
                .collect(),
        );
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

    pub fn evaluate_challenges(&self, challenges: &[F]) -> Self {
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
        self - Constraint::from(rhs)
    }
}

impl<F: GpuField> Sub<&F> for &Constraint<F> {
    type Output = Constraint<F>;

    fn sub(self, rhs: &F) -> Self::Output {
        self - Constraint::from(*rhs)
    }
}

forward_ref_binop!(impl<F: GpuField> Mul, mul for Constraint<F>, Constraint<F>);
forward_ref_binop!(impl<F: GpuField> Add, add for Constraint<F>, Constraint<F>);
forward_ref_binop!(impl<F: GpuField> Sub, sub for Constraint<F>, Constraint<F>);
forward_ref_binop!(impl<F: GpuField> Mul, mul for Constraint<F>, F);
forward_ref_binop!(impl<F: GpuField> Add, add for Constraint<F>, F);
forward_ref_binop!(impl<F: GpuField> Sub, sub for Constraint<F>, F);

impl<F: GpuField> MulAssign<&Constraint<F>> for Constraint<F> {
    fn mul_assign(&mut self, other: &Constraint<F>) {
        *self = &*self * other
    }
}

impl<F: GpuField> MulAssign<&F> for Constraint<F> {
    fn mul_assign(&mut self, rhs: &F) {
        // self.clone() * *rhs
        *self = &*self * rhs
    }
}

impl<F: GpuField> AddAssign<&Constraint<F>> for Constraint<F> {
    fn add_assign(&mut self, other: &Constraint<F>) {
        *self = &*self + other
    }
}

impl<F: GpuField> AddAssign<&F> for Constraint<F> {
    fn add_assign(&mut self, rhs: &F) {
        *self = &*self + rhs
    }
}

impl<F: GpuField> SubAssign<&Constraint<F>> for Constraint<F> {
    fn sub_assign(&mut self, other: &Constraint<F>) {
        *self = &*self - other
    }
}

impl<F: GpuField> SubAssign<&F> for Constraint<F> {
    fn sub_assign(&mut self, rhs: &F) {
        *self = &*self - rhs
    }
}

// Adapted from the `forward_ref_op_assign!` macro in the Rust standard library.
// implements "T op= U", based on "T op= &U"
macro_rules! forward_ref_op_assign {
    (impl < F: GpuField > $imp:ident, $method:ident for $t:ty, $u:ty) => {
        impl<F: GpuField> $imp<$u> for $t {
            #[inline]
            fn $method(&mut self, other: $u) {
                $imp::$method(self, &other);
            }
        }
    };
}

forward_ref_op_assign!(impl<F: GpuField> AddAssign, add_assign for Constraint<F>, Constraint<F>);
forward_ref_op_assign!(impl<F: GpuField> SubAssign, sub_assign for Constraint<F>, Constraint<F>);
forward_ref_op_assign!(impl<F: GpuField> MulAssign, mul_assign for Constraint<F>, Constraint<F>);
forward_ref_op_assign!(impl<F: GpuField> AddAssign, add_assign for Constraint<F>, F);
forward_ref_op_assign!(impl<F: GpuField> SubAssign, sub_assign for Constraint<F>, F);
forward_ref_op_assign!(impl<F: GpuField> MulAssign, mul_assign for Constraint<F>, F);

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
