use super::Felt;
use super::Univariate;
use num_traits::NumAssign;
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
use std::ops::Sub;

/// Multivariate polynomial
#[derive(Clone)]
pub struct Multivariate<E> {
    pub powers: Vec<Vec<u128>>,
    pub coefficients: Vec<E>,
}

impl<E: Felt> Multivariate<E> {
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

    pub fn lift(polynomial: Univariate<E>, lift_index: usize) -> Self {
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

    pub fn evaluate(&self, _point: &[E]) -> E {
        // let mut accumulator = E::zero();
        // for (pad, coefficient) in self.powers.iter().zip(self.coefficients.iter()) {
        //     let mut product = *coefficient;
        //     for (i, power) in pad.iter().enumerate() {
        //         product *= point[i].pow((*power).into());
        //     }
        //     accumulator += product;
        // }
        // accumulator
        todo!()
    }

    pub fn evaluate_symbolic(&self, point: &[Univariate<E>]) -> Univariate<E> {
        let mut accumulator = Univariate::new(vec![]);
        for (pad, coefficient) in self.powers.iter().zip(self.coefficients.iter()) {
            let mut product = Univariate::new(vec![*coefficient]);
            for (i, power) in pad.iter().enumerate() {
                // TODO: change all `a ^ b` to `a.pow(b)`
                product = product * (point[i].clone() ^ *power);
            }
            accumulator = accumulator + product;
        }
        accumulator
    }

    pub fn symbolic_degree_bound(self, max_degrees: &[usize]) -> usize {
        if self.degree() == 0 {
            return 0;
        }
        let mut total_degree_bound = 0;
        for (pad, coefficient) in self.powers.iter().zip(self.coefficients.iter()) {
            if coefficient.is_zero() {
                continue;
            }
            let mut term_degree_bound = 0;
            for (&exponent, &max_degree) in pad.iter().zip(max_degrees) {
                term_degree_bound += exponent as usize * max_degree;
            }
            total_degree_bound = usize::max(total_degree_bound, term_degree_bound);
        }
        total_degree_bound as usize
    }
}

impl<E: Felt> From<E> for Multivariate<E> {
    fn from(v: E) -> Self {
        Multivariate::constant(v)
    }
}

impl<E: Felt> Add for Multivariate<E> {
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

        Multivariate {
            powers,
            coefficients,
        }
    }
}

impl<E: Felt> Add<E> for Multivariate<E> {
    type Output = Self;

    fn add(self, rhs: E) -> Self::Output {
        self + Self::constant(rhs)
    }
}

impl<E: Felt> Sub for Multivariate<E> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<E: Felt> Sub<E> for Multivariate<E> {
    type Output = Self;

    fn sub(self, rhs: E) -> Self::Output {
        self + Self::constant(-rhs)
    }
}

impl<E: Felt> Neg for Multivariate<E> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Multivariate {
            powers: self.powers,
            coefficients: self
                .coefficients
                .into_iter()
                .map(|coefficient| -coefficient)
                .collect(),
        }
    }
}

impl<E: Felt> Mul for Multivariate<E> {
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

impl<E: Felt> Mul<E> for Multivariate<E> {
    type Output = Self;

    fn mul(self, rhs: E) -> Self::Output {
        self * Multivariate::constant(rhs)
    }
}

impl<E: Felt> BitXor<u128> for Multivariate<E> {
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

impl<E: Felt> fmt::Display for Multivariate<E> {
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
