use ark_ff::Field;
use serde::Deserialize;
use serde::Serialize;
use std::fmt;
use std::iter::zip;
use std::ops::Add;
use std::ops::BitXor;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Rem;
use std::ops::Sub;

/// Univariate polynomial
#[derive(Clone, Serialize, Deserialize)]
pub struct Univariate<F> {
    pub coefficients: Vec<F>,
}

impl<F: Field> fmt::Display for Univariate<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let formatted_coefficients: Vec<String> = self
            .coefficients
            .iter()
            .enumerate()
            .filter(|(_, coefficient)| !coefficient.is_zero())
            .map(|(i, coefficient)| {
                if i == 0 {
                    coefficient.to_string()
                } else if i == 1 {
                    format!("{}x", coefficient)
                } else {
                    format!("{}x^{}", coefficient, i)
                }
            })
            .collect();

        write!(f, "{}", formatted_coefficients.join(" + "))
    }
}

impl<F: Field> Univariate<F> {
    pub fn degree(&self) -> isize {
        let mut all_zero = true;
        let mut max_index = 0;

        for (idx, coefficient) in self.coefficients.iter().enumerate() {
            if !coefficient.is_zero() {
                all_zero = false;
                max_index = idx;
            }
        }

        if all_zero {
            -1
        } else {
            max_index.try_into().unwrap()
        }
    }

    pub fn new(coefficients: Vec<F>) -> Self {
        Univariate { coefficients }
    }

    pub fn x() -> Self {
        Self::new(vec![F::zero(), F::one()])
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }

    fn leading_coefficient(&self) -> Option<F> {
        let degree = self.degree();
        if degree >= 0 {
            Some(self.coefficients[degree.unsigned_abs()])
        } else {
            None
        }
    }

    /// Returns quotient and remainder.
    ///
    /// Errors if division cannot be performed.
    fn divide(numerator: Self, denominator: Self) -> Result<(Self, Self), &'static str> {
        if denominator.degree() == -1 {
            return Err("cannot divide by the empty polynomial");
        }

        if numerator.degree() < denominator.degree() {
            return Ok((Univariate::new(vec![]), numerator));
        }

        let mut remainder = numerator.clone();
        let mut quotient_coefficients = vec![
            F::zero();
            (numerator.degree() - denominator.degree() + 1)
                .try_into()
                .unwrap()
        ];

        for _ in 0..quotient_coefficients.len() {
            if remainder.degree() < denominator.degree() {
                break;
            }

            let coefficient = remainder.leading_coefficient().unwrap()
                / denominator.leading_coefficient().unwrap();
            let shift: usize = (remainder.degree() - denominator.degree())
                .try_into()
                .unwrap();

            let mut subtractee_coefficients = vec![F::zero(); shift];
            subtractee_coefficients.push(coefficient);
            let subtractee = Univariate::new(subtractee_coefficients) * denominator.clone();

            quotient_coefficients[shift] = coefficient;
            remainder = remainder - subtractee;
        }

        let quotient = Univariate::new(quotient_coefficients);

        Ok((quotient, remainder))
    }

    pub fn evaluate(&self, point: F) -> F {
        let mut x_pow_i = F::one();
        let mut value = F::zero();
        for coefficient in self.coefficients.iter() {
            value += *coefficient * x_pow_i;
            x_pow_i *= point;
        }
        value
    }

    pub fn evaluate_domain(&self, domain: &[F]) -> Vec<F> {
        domain
            .iter()
            .map(|coefficient| self.evaluate(*coefficient))
            .collect()
    }

    pub fn interpolate(domain: &[F], values: &[F]) -> Self {
        assert_eq!(
            domain.len(),
            values.len(),
            "number of elements in domain does not match number of values -- cannot interpolate"
        );
        // Generate master numerator polynomial: (x - domain[0]) * (x - domain[1]) *
        // ....
        let root = Self::zerofier_domain(domain);

        // Generate the numerator for each item in the domain
        let numerators = domain
            .iter()
            .copied()
            // root / (x - domain[i])
            .map(|d| root.clone() / (Univariate::new(vec![-d, F::one()])))
            .collect::<Vec<Univariate<F>>>();

        // Generate denominators by evaluating numerator polys at each x
        let mut inverse_denominators = numerators
            .iter()
            .zip(domain.iter().copied())
            .map(|(numerator, d)| numerator.evaluate(d))
            .collect::<Vec<F>>();
        ark_ff::batch_inversion(&mut inverse_denominators);

        // Generate output polynomial
        let mut output_coefficients = vec![F::zero(); values.len()];

        for ((y, numerator), inverse_denominator) in values
            .iter()
            .copied()
            .zip(numerators)
            .zip(inverse_denominators)
        {
            let y_slice = y * inverse_denominator;
            for (j, coefficient) in numerator.coefficients.into_iter().enumerate() {
                output_coefficients[j] += coefficient * y_slice;
            }
        }

        Univariate::new(output_coefficients)
    }

    pub fn zerofier_domain(domain: &[F]) -> Self {
        let x = Self::x();
        let mut accumulator = Univariate::new(vec![F::one()]);
        for element in domain.iter() {
            accumulator = accumulator * (x.clone() - Univariate::new(vec![*element]));
        }
        accumulator
    }

    pub fn scale(&self, factor: F) -> Self {
        let mut accumulator = F::one();
        let mut res = self.coefficients.clone();
        for coefficient in &mut res {
            *coefficient *= accumulator;
            accumulator *= factor;
        }
        Self::new(res)
    }

    pub fn test_colinearity(points: Vec<(F, F)>) -> bool {
        let (domain, values): (Vec<F>, Vec<F>) = points.into_iter().unzip();
        let polynomial = Self::interpolate(&domain, &values);
        polynomial.degree() == 1
    }
}

impl<F: Field> Div for Univariate<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let (quotient, remainder) = Univariate::divide(self, rhs).unwrap();
        assert!(
            remainder.is_zero(),
            "cannot perform polynomial division because remainder is not zero"
        );
        quotient
    }
}

impl<F: Field> Rem for Univariate<F> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        let (_, remainder) = Univariate::divide(self, rhs).unwrap();
        remainder
    }
}

impl<F: Field> Neg for Univariate<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Univariate::new(
            self.coefficients
                .into_iter()
                .map(|coefficient| -coefficient)
                .collect(),
        )
    }
}

impl<F: Field> Add for Univariate<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.degree() == -1 {
            return rhs;
        }

        if rhs.degree() == -1 {
            return self;
        }

        let mut coefficients =
            vec![F::zero(); usize::max(self.coefficients.len(), rhs.coefficients.len())];

        for (idx, coefficient) in self.coefficients.into_iter().enumerate() {
            coefficients[idx] = coefficient;
        }

        for (idx, coefficient) in rhs.coefficients.into_iter().enumerate() {
            coefficients[idx] += coefficient;
        }

        Univariate::new(coefficients)
    }
}

impl<F: Field> Add for &Univariate<F> {
    type Output = Univariate<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.degree() == -1 {
            return rhs.clone();
        }

        if rhs.degree() == -1 {
            return self.clone();
        }

        let mut coefficients =
            vec![F::zero(); usize::max(self.coefficients.len(), rhs.coefficients.len())];

        for (i, coefficient) in coefficients.iter_mut().enumerate() {
            *coefficient = match (self.coefficients.get(i), rhs.coefficients.get(i)) {
                (Some(a), Some(b)) => *a + b,
                (Some(a), None) | (None, Some(a)) => *a,
                (None, None) => F::zero(),
            }
        }

        Univariate::new(coefficients)
    }
}

impl<F: Field> Sub for Univariate<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<F: Field> Mul for Univariate<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() || rhs.coefficients.is_empty() {
            return Univariate::new(vec![]);
        }

        let mut buffer = vec![F::zero(); self.coefficients.len() + rhs.coefficients.len() - 1];

        for (i, coefficient_a) in self.coefficients.into_iter().enumerate() {
            if coefficient_a.is_zero() {
                continue;
            }

            for (j, coefficient_b) in rhs.coefficients.iter().copied().enumerate() {
                buffer[i + j] += coefficient_a * coefficient_b;
            }
        }

        Univariate::new(buffer)
    }
}

impl<F: Field> Mul for &Univariate<F> {
    type Output = Univariate<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() || rhs.coefficients.is_empty() {
            return Univariate::new(vec![]);
        }

        let mut buffer = vec![F::zero(); self.coefficients.len() + rhs.coefficients.len() - 1];

        for (i, coefficient_a) in self.coefficients.iter().copied().enumerate() {
            if coefficient_a.is_zero() {
                continue;
            }

            for (j, coefficient_b) in rhs.coefficients.iter().copied().enumerate() {
                buffer[i + j] += coefficient_a * coefficient_b;
            }
        }

        Univariate::new(buffer)
    }
}

impl<F: Field> BitXor<u128> for Univariate<F> {
    type Output = Self;

    fn bitxor(self, exponent: u128) -> Self::Output {
        if self.is_zero() {
            return Univariate::new(vec![]);
        }

        if exponent == 0 {
            return Univariate::new(vec![F::one()]);
        }

        let mut accumulator = Univariate::new(vec![F::one()]);
        for i in (0..128).rev() {
            accumulator = accumulator.clone() * accumulator;
            let bit = 1u128 << i;
            if bit & exponent != 0 {
                accumulator = accumulator * self.clone();
            }
        }

        accumulator
    }
}

impl<F: Field> PartialEq for Univariate<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.degree() != other.degree() {
            false
        } else if self.degree() == -1 {
            true
        } else {
            zip(self.coefficients.iter(), other.coefficients.iter()).all(|(a, b)| a == b)
        }
    }
}
