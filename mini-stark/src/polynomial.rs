use fast_poly::fields::batch_inverse;
use fast_poly::fields::Felt;
use serde::Deserialize;
use serde::Serialize;
use std::char;
use std::cmp::Ordering;
use std::fmt;
use std::iter::once;
use std::iter::repeat;
use std::iter::zip;
use std::ops::Add;
use std::ops::BitXor;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Rem;
use std::ops::Sub;

#[derive(Clone, Serialize, Deserialize)]
pub struct Polynomial<E> {
    pub coefficients: Vec<E>,
}

impl<E: Felt> fmt::Display for Polynomial<E> {
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

// impl Deref for Polynomial {
//     type Target = Vec<Felt>;

//     fn deref(&self) -> &Self::Target {
//         &self.coefficients
//     }
// }

impl<E: Felt> Polynomial<E> {
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

    pub fn new(coefficients: Vec<E>) -> Self {
        Polynomial { coefficients }
    }

    pub fn x() -> Self {
        Self::new(vec![E::zero(), E::one()])
    }

    pub fn is_zero(&self) -> bool {
        self.degree() == -1
    }

    fn leading_coefficient(&self) -> Option<E> {
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
            return Ok((Polynomial::new(vec![]), numerator));
        }

        let mut remainder = numerator.clone();
        let mut quotient_coefficients = vec![
            E::zero();
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

            let mut subtractee_coefficients = vec![E::zero(); shift];
            subtractee_coefficients.push(coefficient);
            let subtractee = Polynomial::new(subtractee_coefficients) * denominator.clone();

            quotient_coefficients[shift] = coefficient;
            remainder = remainder - subtractee;
        }

        let quotient = Polynomial::new(quotient_coefficients);

        Ok((quotient, remainder))
    }

    pub fn evaluate(&self, point: E) -> E {
        let mut x_pow_i = E::one();
        let mut value = E::zero();
        for coefficient in self.coefficients.iter() {
            value += *coefficient * x_pow_i;
            x_pow_i *= point;
        }
        value
    }

    pub fn evaluate_domain(&self, domain: &[E]) -> Vec<E> {
        domain
            .iter()
            .map(|coefficient| self.evaluate(*coefficient))
            .collect()
    }

    pub fn interpolate(domain: &[E], values: &[E]) -> Self {
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
            .map(|d| root.clone() / (Polynomial::new(vec![-d, E::one()])))
            .collect::<Vec<Polynomial<E>>>();

        // Generate denominators by evaluating numerator polys at each x
        let inverse_denominators = batch_inverse(
            &numerators
                .iter()
                .zip(domain.iter().copied())
                .map(|(numerator, d)| numerator.evaluate(d))
                .collect::<Vec<E>>(),
        );

        // Generate output polynomial
        let mut output_coefficients = vec![E::zero(); values.len()];

        for ((y, numerator), inverse_denominator) in values
            .iter()
            .copied()
            .zip(numerators)
            .zip(inverse_denominators)
        {
            let y_slice = y * inverse_denominator.unwrap();
            for (j, coefficient) in numerator.coefficients.into_iter().enumerate() {
                output_coefficients[j] += coefficient * y_slice;
            }
        }

        Polynomial::new(output_coefficients)

        //     .map(|((value, numerator), inverse_denominator)| {
        //         Polynomial::new(vec![value])
        //             * numerator
        //             * Polynomial::new(vec![inverse_denominator])
        //     })
        //     .sum()
        //     .collect()

        // assert!(domain.len() > 0, "cannot interpolate between zero points");
        // let field = domain[0].field;
        // let x = Polynomial::x(field);
        // let mut accumulator = Polynomial::new(vec![]);

        // for i in 0..domain.len() {
        //     let mut prod = Polynomial::new(vec![values[i]]);
        //     for j in 0..domain.len() {
        //         if j == i {
        //             continue;
        //         }
        //         prod = &prod
        //             * &(x.clone() - Polynomial::new(vec![domain[j]]))
        //             * Polynomial::new(vec![(domain[i] -
        //               domain[j]).inverse()]);
        //     }
        //     accumulator = &accumulator + &prod;
        // }
        // accumulator
    }

    pub fn zerofier_domain(domain: &[E]) -> Self {
        let x = Self::x();
        let mut accumulator = Polynomial::new(vec![E::one()]);
        for element in domain.iter() {
            accumulator = accumulator * (x.clone() - Polynomial::new(vec![*element]));
        }
        accumulator
    }

    pub fn scale(&self, factor: E) -> Self {
        Polynomial::new(
            self.coefficients
                .iter()
                .copied()
                .enumerate()
                .map(|(idx, coefficient)| coefficient * factor.pow((idx as u128).into()))
                .collect(),
        )
    }

    pub fn test_colinearity(points: Vec<(E, E)>) -> bool {
        let (domain, values): (Vec<E>, Vec<E>) = points.into_iter().unzip();
        let polynomial = Self::interpolate(&domain, &values);
        polynomial.degree() == 1
    }
}

impl<E: Felt> Div for Polynomial<E> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let (quotient, remainder) = Polynomial::divide(self, rhs).unwrap();
        assert!(
            remainder.is_zero(),
            "cannot perform polynomial division because remainder is not zero"
        );
        quotient
    }
}

impl<E: Felt> Rem for Polynomial<E> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        let (_, remainder) = Polynomial::divide(self, rhs).unwrap();
        remainder
    }
}

impl<E: Felt> Neg for Polynomial<E> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Polynomial::new(
            self.coefficients
                .into_iter()
                .map(|coefficient| -coefficient)
                .collect(),
        )
    }
}

impl<E: Felt> Add for Polynomial<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.degree() == -1 {
            return rhs;
        }

        if rhs.degree() == -1 {
            return self;
        }

        let mut coefficients = vec![E::zero(); self.coefficients.len().max(rhs.coefficients.len())];

        for (idx, coefficient) in self.coefficients.into_iter().enumerate() {
            coefficients[idx] = coefficient;
        }

        for (idx, coefficient) in rhs.coefficients.into_iter().enumerate() {
            coefficients[idx] += coefficient;
        }

        Polynomial::new(coefficients)
    }
}

impl<E: Felt> Add for &Polynomial<E> {
    type Output = Polynomial<E>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.degree() == -1 {
            return rhs.clone();
        }

        if rhs.degree() == -1 {
            return self.clone();
        }

        let mut coefficients = vec![E::zero(); self.coefficients.len().max(rhs.coefficients.len())];

        for (i, coefficient) in coefficients.iter_mut().enumerate() {
            *coefficient = match (self.coefficients.get(i), rhs.coefficients.get(i)) {
                (Some(a), Some(b)) => *a + b,
                (Some(a), None) | (None, Some(a)) => *a,
                (None, None) => E::zero(),
            }
        }

        Polynomial::new(coefficients)
    }
}

impl<E: Felt> Sub for Polynomial<E> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<E: Felt> Mul for Polynomial<E> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() || rhs.coefficients.is_empty() {
            return Polynomial::new(vec![]);
        }

        let mut buffer = vec![E::zero(); self.coefficients.len() + rhs.coefficients.len() - 1];

        for (i, coefficient_a) in self.coefficients.into_iter().enumerate() {
            if coefficient_a.is_zero() {
                continue;
            }

            for (j, coefficient_b) in rhs.coefficients.iter().copied().enumerate() {
                buffer[i + j] += coefficient_a * coefficient_b;
            }
        }

        Polynomial::new(buffer)
    }
}

impl<E: Felt> Mul for &Polynomial<E> {
    type Output = Polynomial<E>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.coefficients.is_empty() || rhs.coefficients.is_empty() {
            return Polynomial::new(vec![]);
        }

        let mut buffer = vec![E::zero(); self.coefficients.len() + rhs.coefficients.len() - 1];

        for (i, coefficient_a) in self.coefficients.iter().copied().enumerate() {
            if coefficient_a.is_zero() {
                continue;
            }

            for (j, coefficient_b) in rhs.coefficients.iter().copied().enumerate() {
                buffer[i + j] += coefficient_a * coefficient_b;
            }
        }

        Polynomial::new(buffer)
    }
}

impl<E: Felt> BitXor<u128> for Polynomial<E> {
    type Output = Self;

    fn bitxor(self, exponent: u128) -> Self::Output {
        if self.is_zero() {
            return Polynomial::new(vec![]);
        }

        if exponent == 0 {
            return Polynomial::new(vec![E::one()]);
        }

        let mut accumulator = Polynomial::new(vec![E::one()]);
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

impl<E: Felt> PartialEq for Polynomial<E> {
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

#[derive(Clone)]
pub struct MultivariatePolynomial<E: Felt> {
    pub powers: Vec<Vec<u128>>,
    pub coefficients: Vec<E>,
}

impl<E: Felt> MultivariatePolynomial<E> {
    fn search_for_power<'a>(seek: &'a [u128]) -> impl FnMut(&'a Vec<u128>) -> Ordering {
        move |probe: &Vec<u128>| {
            assert_eq!(seek.len(), probe.len());

            for (a, b) in probe.iter().zip(seek.iter()) {
                let ordering = a.cmp(b);

                if ordering != Ordering::Equal {
                    return ordering;
                }
            }

            Ordering::Equal
        }
    }

    pub fn degree(&self) -> isize {
        let all_zero = !self
            .coefficients
            .iter()
            .any(|coefficient| !coefficient.is_zero());

        if all_zero {
            return -1;
        }

        self.powers
            .iter()
            .map(|pad| pad.iter().copied().sum::<u128>() as isize)
            .max()
            .unwrap_or(0)
    }

    pub fn zero() -> Self {
        Self {
            powers: vec![],
            coefficients: vec![],
        }
    }

    pub fn one() -> Self {
        Self {
            powers: vec![vec![0]],
            coefficients: vec![E::one()],
        }
    }

    pub fn is_zero(&self) -> bool {
        for coefficient in self.coefficients.iter() {
            if !coefficient.is_zero() {
                return false;
            }
        }

        true
    }

    pub fn constant(element: E) -> Self {
        Self {
            powers: vec![vec![0]],
            coefficients: vec![element],
        }
    }

    pub fn variables(num_variables: usize) -> Vec<Self> {
        (0..num_variables)
            .map(|i| Self {
                powers: vec![repeat(0)
                    .take(i)
                    .chain(once(1))
                    .chain(repeat(0).take(num_variables - i - 1))
                    .collect()],
                coefficients: vec![E::one()],
            })
            .collect()
    }

    pub fn lift(polynomial: Polynomial<E>, lift_index: usize) -> Self {
        if polynomial.is_zero() {
            return Self::zero();
        }

        let variables = Self::variables(lift_index + 1);
        let x = variables.last().unwrap();
        let mut accumulator = Self::zero();

        for (i, coefficient) in polynomial.coefficients.into_iter().enumerate() {
            accumulator =
                accumulator + Self::constant(coefficient) * (x.clone() ^ i.try_into().unwrap());
        }

        accumulator
    }

    pub fn evaluate(&self, point: &[E]) -> E {
        let mut accumulator = E::zero();
        for (pad, coefficient) in self.powers.iter().zip(self.coefficients.iter()) {
            let mut product = *coefficient;
            for (i, power) in pad.iter().enumerate() {
                product *= point[i].pow((*power).into());
            }
            accumulator += product;
        }
        accumulator
    }

    pub fn evaluate_symbolic(&self, point: &[Polynomial<E>]) -> Polynomial<E> {
        let mut accumulator = Polynomial::new(vec![]);
        for (pad, coefficient) in self.powers.iter().zip(self.coefficients.iter()) {
            let mut product = Polynomial::new(vec![*coefficient]);
            for (i, power) in pad.iter().enumerate() {
                product = product * (point[i].clone() ^ *power);
            }
            accumulator = accumulator + product;
        }
        accumulator
    }
}

impl<E: Felt> Add for MultivariatePolynomial<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut powers: Vec<Vec<u128>> = vec![];
        let mut coefficients: Vec<E> = vec![];

        let num_variables = self
            .powers
            .iter()
            .chain(rhs.powers.iter())
            .map(|pad| pad.len())
            .max()
            .unwrap_or(0);

        for (pad, coefficient) in self.powers.into_iter().zip(self.coefficients.into_iter()) {
            let num_new_variables = num_variables - pad.len();
            let pad = pad
                .into_iter()
                .chain(repeat(0).take(num_new_variables))
                .collect();
            powers.push(pad);
            coefficients.push(coefficient)
        }

        for (pad, coefficient) in rhs.powers.into_iter().zip(rhs.coefficients.into_iter()) {
            let num_new_variables = num_variables - pad.len();
            let pad = pad
                .into_iter()
                .chain(repeat(0).take(num_new_variables))
                .collect::<Vec<_>>();

            let search_result = powers.binary_search_by(Self::search_for_power(&pad));

            match search_result {
                Ok(i) => coefficients[i] += coefficient,
                Err(i) => {
                    powers.insert(i, pad);
                    coefficients.insert(i, coefficient);
                }
            }
        }

        MultivariatePolynomial {
            powers,
            coefficients,
        }
    }
}

impl<E: Felt> Sub for MultivariatePolynomial<E> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<E: Felt> Neg for MultivariatePolynomial<E> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        MultivariatePolynomial {
            powers: self.powers,
            coefficients: self
                .coefficients
                .into_iter()
                .map(|coefficient| -coefficient)
                .collect(),
        }
    }
}

impl<E: Felt> Mul for MultivariatePolynomial<E> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut powers: Vec<Vec<u128>> = vec![];
        let mut coefficients: Vec<E> = vec![];

        let num_variables = self
            .powers
            .iter()
            .chain(rhs.powers.iter())
            .map(|pad| pad.len())
            .max()
            .unwrap_or(0);

        for (lhs_pad, lhs_coefficient) in self.powers.into_iter().zip(self.coefficients.into_iter())
        {
            for (rhs_pad, rhs_coefficient) in
                rhs.powers.iter().zip(rhs.coefficients.iter().copied())
            {
                let mut pad = repeat(0).take(num_variables).collect::<Vec<_>>();

                for (i, power) in lhs_pad.iter().enumerate().chain(rhs_pad.iter().enumerate()) {
                    pad[i] += power;
                }

                let search_result = powers.binary_search_by(Self::search_for_power(&pad));

                match search_result {
                    Ok(i) => coefficients[i] += lhs_coefficient * rhs_coefficient,
                    Err(i) => {
                        powers.insert(i, pad);
                        coefficients.insert(i, lhs_coefficient * rhs_coefficient);
                    }
                }
            }
        }

        Self {
            powers,
            coefficients,
        }
    }
}

impl<E: Felt> BitXor<u128> for MultivariatePolynomial<E> {
    type Output = Self;

    fn bitxor(self, exponent: u128) -> Self::Output {
        if self.is_zero() {
            return Self::zero();
        }

        let mut accumulator = Self::one();

        for i in (0..128).rev() {
            let bit = 1u128 << i;
            accumulator = accumulator.clone() * accumulator;

            if exponent & bit != 0 {
                accumulator = accumulator * self.clone();
            }
        }

        accumulator
    }
}

impl<E: Felt> fmt::Display for MultivariatePolynomial<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let formatted_coefficients: Vec<String> = self
            .powers
            .iter()
            .zip(self.coefficients.iter())
            .filter(|(_, coefficient)| !coefficient.is_zero())
            .map(|(pad, coefficient)| {
                format!(
                    "{}{}",
                    if coefficient.is_one() && pad.iter().sum::<u128>() != 0 {
                        String::new()
                    } else {
                        coefficient.to_string()
                    },
                    pad.iter()
                        .enumerate()
                        .map(|(idx, power)| {
                            if *power == 0 {
                                return String::new();
                            }

                            let mut variable_name = char::REPLACEMENT_CHARACTER;

                            if idx == 0 {
                                variable_name = 'x';
                            } else if idx == 1 {
                                variable_name = 'y';
                            } else if idx == 2 {
                                variable_name = 'z';
                            } else if idx == 3 {
                                variable_name = 'w';
                            } else if idx == 4 {
                                variable_name = 'a';
                            } else if idx == 5 {
                                variable_name = 'b';
                            } else if idx == 6 {
                                variable_name = 'c';
                            }

                            if *power == 1 {
                                variable_name.to_string()
                            } else {
                                format!("{}^{}", variable_name, power)
                            }
                        })
                        .collect::<Vec<_>>()
                        .join("")
                )
            })
            .collect();

        write!(f, "{}", formatted_coefficients.join(" + "))
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn coliniarity() {
//         let field_elements = vec![
//             (
//                 Felt::new(227292946527062975142212960787649844118,
// Field::main()),
// Felt::new(198529779817252628145467949537111849632, Field::main()),
//             ),
//             (
//                 Felt::new(43204950615167404993711775979400277099,
// Field::main()),
// Felt::new(39518110489409123482753928896191159935, Field::main()),
//             ),
//             (
//                 Felt::new(49439774929837278997840997361743450321,
// Field::main()),
// Felt::new(203026432109815629425002792571619313110, Field::main()),
//             ),
//         ];

//         assert!(
//             Polynomial::test_colinearity(field_elements),
//             "points should be colinear"
//         );
//     }

//     #[test]
//     fn symbolic_evaluate() {
//         let field = Field::main();
//         let poly_1 =
//             MultivariatePolynomial::lift(Polynomial::new(vec![field.zero(),
// field.one()]), 0);         let poly_2 =
//             MultivariatePolynomial::lift(Polynomial::new(vec![field.zero(),
// field.one()]), 1);         let poly_3 =
// MultivariatePolynomial::lift(Polynomial::new(vec![field.one()]), 0);

//         let combined = poly_1 + poly_2 + poly_3;

//         //poly_1 + poly_2
//         println!("{}", combined);

//         println!("{:?}", combined.powers.clone());

//         println!(
//             "{}",
//             combined.evaluate_symbolic(&vec![
//                 Polynomial::x(field),
//                 Polynomial::new(vec![field.one()]),
//             ])
//         )
//     }

//     #[test]
//     fn test_evaluate() {
//         let field = Field::main();
//         let variables = MultivariatePolynomial::variables(4, field);
//         let zero = field.zero();
//         let one = field.one();
//         let two = Felt::new(2, field);
//         let five = Felt::new(5, field);

//         let mpoly1 = MultivariatePolynomial::constant(one) *
// variables[0].clone()             + MultivariatePolynomial::constant(two) *
// variables[1].clone()             + MultivariatePolynomial::constant(five) *
// (variables[2].clone() ^ 3);         let mpoly2 =
//             MultivariatePolynomial::constant(one) * variables[0].clone() *
// variables[3].clone()                 + MultivariatePolynomial::constant(five)
// * (variables[3].clone() ^ 3)                 +
// MultivariatePolynomial::constant(five);         let mpoly3 = mpoly1.clone() *
// mpoly2.clone();

//         let point = vec![zero, five, five, two];

//         let eval1 = mpoly1.evaluate(&point);
//         let eval2 = mpoly2.evaluate(&point);
//         let eval3 = mpoly3.evaluate(&point);

//         assert_eq!(
//             eval1 * eval2,
//             eval3,
//             "multivariate polynomial multiplication does not commute with
// evaluation",         );
//         assert_eq!(
//             eval1 + eval2,
//             (mpoly1 + mpoly2).evaluate(&point),
//             "multivariate polynomial addition does not commute with
// evaluation",         );

//         println!("eval3: {}", eval3.value);
//         println!("multivariate evaluate test success \\o/");
//     }

//     #[test]
//     fn test_lift() {
//         let field = Field::main();
//         let variables = MultivariatePolynomial::variables(4, field);
//         let zero = field.zero();
//         let one = field.one();
//         let two = Felt::new(2, field);
//         let five = Felt::new(5, field);

//         let upoly = Polynomial::interpolate(&vec![zero, one, two], &vec![two,
// five, five]);         let mpoly = MultivariatePolynomial::lift(upoly.clone(),
// 3);

//         println!("{}", upoly);
//         println!("{}", mpoly);

//         assert!(
//             upoly.evaluate(five) == mpoly.evaluate(&vec![zero, zero, zero,
// five]),             "lifting univariate to multivariate failed",
//         );

//         println!("lifting univariate to multivariate polynomial success
// \\o/");     }

//     #[test]
//     fn test_interpolate() {
//         let field = Field::main();
//         let domain = vec![
//             Felt::new(3, field),
//             Felt::new(4, field),
//             Felt::new(5, field),
//         ];
//         let values = vec![
//             Felt::new(1, field),
//             Felt::new(2, field),
//             Felt::new(100, field),
//         ];

//         let interpolant = Polynomial::interpolate(&domain, &values);

//         assert_eq!(interpolant.evaluate(domain[0]), values[0]);
//         assert_eq!(interpolant.evaluate(domain[1]), values[1]);
//         assert_eq!(interpolant.evaluate(domain[2]), values[2]);
//     }
// }
