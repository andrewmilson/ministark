use crate::constraints::AlgebraicItem;
use crate::constraints::PeriodicColumn;
use crate::expression::Expr;
use crate::utils::FieldVariant;
use crate::utils::GpuAllocator;
use crate::Matrix;
use crate::StarkExtensionOf;
use alloc::borrow::Cow;
use alloc::boxed::Box;
use alloc::vec::Vec;
use ark_ff::batch_inversion;
use ark_ff::FftField;
use ark_ff::Field;
use ark_poly::domain::DomainCoeff;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_std::cfg_chunks_mut;
use core::iter::zip;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Div;
use core::ops::Mul;
use core::ops::MulAssign;
use core::ops::Neg;
use ministark_gpu::GpuFftField;
use ministark_gpu::GpuField;
use num_traits::Pow;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::BTreeMap;

#[allow(clippy::too_many_arguments)]
pub fn eval<Fp: GpuFftField<FftField = Fp> + FftField, Fq: StarkExtensionOf<Fp>>(
    expr: &Expr<AlgebraicItem<FieldVariant<Fp, Fq>>>,
    challenges: &[Fq],
    hints: &[Fq],
    lde_step: usize,
    domain_offset: Fp,
    x_lde: &[Fp],
    base_trace_lde_cols: &[&[Fp]],
    extension_trace_lde_cols: Option<&[&[Fq]]>,
) -> Matrix<Fq> {
    let n = x_lde.len();
    let mut result = Vec::with_capacity_in(n, GpuAllocator);
    result.resize(n, Fq::zero());
    match n {
        1..512 => eval_impl::<Fp, Fq, 1>(
            expr,
            challenges,
            hints,
            lde_step,
            domain_offset,
            x_lde,
            base_trace_lde_cols,
            extension_trace_lde_cols,
            &mut result,
        ),
        512.. => eval_impl::<Fp, Fq, 512>(
            expr,
            challenges,
            hints,
            lde_step,
            domain_offset,
            x_lde,
            base_trace_lde_cols,
            extension_trace_lde_cols,
            &mut result,
        ),
        0 => {}
        _ => unreachable!(),
    }
    Matrix::new(vec![result])
}

#[allow(clippy::too_many_arguments)]
fn eval_impl<
    Fp: GpuFftField<FftField = Fp> + FftField,
    Fq: StarkExtensionOf<Fp>,
    const CHUNK_SIZE: usize,
>(
    expr: &Expr<AlgebraicItem<FieldVariant<Fp, Fq>>>,
    challenges: &[Fq],
    hints: &[Fq],
    lde_step: usize,
    domain_offset: Fp,
    x_lde: &[Fp],
    base_trace_lde_cols: &[&[Fp]],
    extension_trace_lde_cols: Option<&[&[Fq]]>,
    result: &mut [Fq],
) {
    use AlgebraicItem::*;
    let n = result.len();
    #[allow(clippy::cast_possible_wrap)]
    let step = lde_step as isize;
    let trace_len = n / lde_step;
    let num_base_columns = base_trace_lde_cols.len();
    let num_extension_columns = extension_trace_lde_cols.map_or(0, <[_]>::len);
    let base_column_range = 0..num_base_columns;
    let extension_column_range = num_base_columns..num_base_columns + num_extension_columns;
    let periodic_column_evals_map =
        build_periodic_column_evals_map(expr, domain_offset, trace_len, lde_step, CHUNK_SIZE);
    cfg_chunks_mut!(result, CHUNK_SIZE)
        .enumerate()
        .for_each(|(i, chunk)| {
            let chunk_offset = CHUNK_SIZE * i;
            let chunk_res: [Fq; CHUNK_SIZE] = expr
                .graph_eval(&mut |leaf| match *leaf {
                    X => EvalItem::Evals(Box::new(FieldVariant::Fp(extract_lde_chunk(
                        x_lde,
                        chunk_offset,
                    )))),
                    Constant(v) => EvalItem::Constant(v),
                    Challenge(i) => EvalItem::Constant(FieldVariant::Fq(challenges[i])),
                    Hint(i) => EvalItem::Constant(FieldVariant::Fq(hints[i])),
                    Trace(col_idx, row_offset) => {
                        let shift = step * row_offset;
                        let chunk_offset = isize::try_from(chunk_offset).unwrap();
                        #[allow(clippy::cast_possible_wrap)]
                        let position = (chunk_offset + shift).rem_euclid(n as isize) as usize;
                        if base_column_range.contains(&col_idx) {
                            let column = &base_trace_lde_cols[col_idx];
                            EvalItem::Evals(Box::new(FieldVariant::Fp(extract_lde_chunk(
                                column, position,
                            ))))
                        } else if extension_column_range.contains(&col_idx) {
                            let extension_column_idx = col_idx - num_base_columns;
                            let column = &extension_trace_lde_cols.unwrap()[extension_column_idx];
                            EvalItem::Evals(Box::new(FieldVariant::Fq(extract_lde_chunk(
                                column, position,
                            ))))
                        } else {
                            panic!("invalid column {col_idx}")
                        }
                    }
                    Periodic(col) => {
                        let lde = periodic_column_evals_map.get(&col).unwrap();
                        match lde {
                            FieldVariant::Fp(lde) => EvalItem::Evals(Box::new(FieldVariant::Fp(
                                extract_lde_chunk(lde, chunk_offset),
                            ))),
                            FieldVariant::Fq(lde) => EvalItem::Evals(Box::new(FieldVariant::Fq(
                                extract_lde_chunk(lde, chunk_offset),
                            ))),
                        }
                    }
                })
                .into_fq_array();
            chunk.copy_from_slice(&chunk_res);
        });
}

/// Extracts a chunk of evaluations from a low-degree-extension
#[inline]
pub fn extract_lde_chunk<F: Field, const CHUNK_SIZE: usize>(
    lde: &'_ [F],
    offset: usize,
) -> Cow<'_, [F; CHUNK_SIZE]> {
    let n = lde.len();
    let lde_offset = offset % n;
    if n >= lde_offset + CHUNK_SIZE {
        Cow::Borrowed(
            (&lde[lde_offset..lde_offset + CHUNK_SIZE])
                .try_into()
                .unwrap(),
        )
    } else {
        let prefix = &lde[lde_offset..];
        let suffix = &lde[0..lde_offset + CHUNK_SIZE - n];
        let mut chunk = [F::ZERO; CHUNK_SIZE];
        chunk[0..prefix.len()].copy_from_slice(prefix);
        chunk[prefix.len()..].copy_from_slice(suffix);
        Cow::Owned(chunk)
    }
}

/// Build a map from periodic column to evaluations
#[allow(clippy::type_complexity)]
pub fn build_periodic_column_evals_map<
    Fp: GpuFftField<FftField = Fp> + FftField,
    Fq: StarkExtensionOf<Fp>,
>(
    expr: &Expr<AlgebraicItem<FieldVariant<Fp, Fq>>>,
    domain_offset: Fp,
    trace_len: usize,
    blowup_factor: usize,
    min_domain_size: usize,
) -> BTreeMap<PeriodicColumn<'static, FieldVariant<Fp, Fq>>, FieldVariant<Vec<Fp>, Vec<Fq>>> {
    let mut res = BTreeMap::new();
    expr.traverse(&mut |node| {
        if let &Expr::Leaf(AlgebraicItem::Periodic(col)) = node {
            let interval_size = col.interval_size();
            let coeffs = col.coeffs();
            let is_fp = |&v| match v {
                FieldVariant::Fp(_) => true,
                FieldVariant::Fq(_) => false,
            };

            let lde = if coeffs.iter().all(is_fp) {
                let coeffs: Vec<Fp> = coeffs
                    .iter()
                    .map(|v| match v {
                        FieldVariant::Fp(v) => *v,
                        FieldVariant::Fq(_) => unreachable!(),
                    })
                    .collect();
                let col = PeriodicColumn::new(&coeffs, interval_size);
                let lde = eval_periodic_column(
                    domain_offset,
                    trace_len,
                    blowup_factor,
                    col,
                    min_domain_size,
                );
                FieldVariant::Fp(lde)
            } else {
                let coeffs: Vec<Fq> = coeffs.iter().map(FieldVariant::as_fq).collect();
                let col = PeriodicColumn::new(&coeffs, interval_size);
                let lde = eval_periodic_column(
                    domain_offset,
                    trace_len,
                    blowup_factor,
                    col,
                    min_domain_size,
                );
                FieldVariant::Fq(lde)
            };

            res.insert(col, lde);
        }
    });
    res
}

/// Generates a preiodic low degree extension of a periodic column of values
pub fn eval_periodic_column<F: GpuField + Field + DomainCoeff<F::FftField>>(
    domain_offset: F::FftField,
    trace_len: usize,
    blowup_factor: usize,
    periodic_column: PeriodicColumn<'_, F>,
    min_len: usize,
) -> Vec<F>
where
    F::FftField: FftField,
{
    let interval_size = periodic_column.interval_size();
    let domain_size = interval_size * blowup_factor;
    let domain_offset = domain_offset.pow([(trace_len / interval_size) as u64]);
    let domain = Radix2EvaluationDomain::new_coset(domain_size, domain_offset).unwrap();
    let mut evals = domain.fft(periodic_column.coeffs());
    let mut i = 0;
    while evals.len() < min_len {
        evals.push(evals[i]);
        i += 1;
    }
    evals
}

#[derive(Clone)]
#[allow(clippy::type_complexity)]
enum EvalItem<'a, Fp: Field, Fq: Field, const N: usize> {
    Constant(FieldVariant<Fp, Fq>),
    Evals(Box<FieldVariant<Cow<'a, [Fp; N]>, Cow<'a, [Fq; N]>>>),
}

impl<'a, Fp: Field, Fq: Field + From<Fp>, const N: usize> EvalItem<'a, Fp, Fq, N> {
    pub fn into_fq_array(self) -> [Fq; N] {
        match self {
            // TODO: think if constant case should be considered. maybe panic?
            Self::Constant(v) => match v {
                FieldVariant::Fp(v) => [v.into(); N],
                FieldVariant::Fq(v) => [v; N],
            },
            Self::Evals(evals) => match *evals {
                FieldVariant::Fp(evals) => evals.into_owned().map(|v| Fq::from(v)),
                FieldVariant::Fq(evals) => evals.into_owned(),
            },
        }
    }

    fn inverse(self) -> Self {
        match self {
            Self::Constant(v) => Self::Constant(v.inverse().unwrap()),
            Self::Evals(evals) => match *evals {
                FieldVariant::Fp(mut evals) => {
                    batch_inversion(evals.to_mut());
                    Self::Evals(Box::new(FieldVariant::Fp(evals)))
                }
                FieldVariant::Fq(mut evals) => {
                    batch_inversion(evals.to_mut());
                    Self::Evals(Box::new(FieldVariant::Fq(evals)))
                }
            },
        }
    }
}

impl<
        'a,
        Fp: Field,
        Fq: Field + From<Fp> + Add<Fp, Output = Fq> + AddAssign<Fp>,
        const N: usize,
    > Add for EvalItem<'a, Fp, Fq, N>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        use EvalItem::*;
        match (self, rhs) {
            (Constant(lhs), Constant(rhs)) => Constant(lhs + rhs),
            (Evals(evals), Constant(c)) | (Constant(c), Evals(evals)) => match (*evals, c) {
                (FieldVariant::Fp(mut evals), FieldVariant::Fp(c)) => {
                    evals.to_mut().iter_mut().for_each(|v| *v += c);
                    Evals(Box::new(FieldVariant::Fp(evals)))
                }
                (FieldVariant::Fp(evals), FieldVariant::Fq(c)) => {
                    let evals = Cow::Owned(evals.into_owned().map(|v| Fq::from(v) + c));
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
                (FieldVariant::Fq(mut evals), FieldVariant::Fp(c)) => {
                    evals.to_mut().iter_mut().for_each(|v| *v += c);
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
                (FieldVariant::Fq(mut evals), FieldVariant::Fq(c)) => {
                    evals.to_mut().iter_mut().for_each(|v| *v += c);
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
            },
            (Evals(lhs), Evals(rhs)) => match (*lhs, *rhs) {
                (FieldVariant::Fq(mut lhs), FieldVariant::Fp(rhs))
                | (FieldVariant::Fp(rhs), FieldVariant::Fq(mut lhs)) => {
                    zip(lhs.to_mut(), &*rhs).for_each(|(lhs, rhs)| *lhs += *rhs);
                    Evals(Box::new(FieldVariant::Fq(lhs)))
                }
                (FieldVariant::Fq(lhs), FieldVariant::Fq(rhs)) => {
                    let (mut lhs, rhs) = if lhs.is_owned() {
                        (lhs, rhs)
                    } else {
                        (rhs, lhs)
                    };
                    zip(lhs.to_mut(), &*rhs).for_each(|(lhs, rhs)| *lhs += *rhs);
                    Evals(Box::new(FieldVariant::Fq(lhs)))
                }
                (FieldVariant::Fp(lhs), FieldVariant::Fp(rhs)) => {
                    let (mut lhs, rhs) = if lhs.is_owned() {
                        (lhs, rhs)
                    } else {
                        (rhs, lhs)
                    };
                    zip(lhs.to_mut(), &*rhs).for_each(|(lhs, rhs)| *lhs += *rhs);
                    Evals(Box::new(FieldVariant::Fp(lhs)))
                }
            },
        }
    }
}

impl<
        'a,
        Fp: Field,
        Fq: Field + From<Fp> + Mul<Fp, Output = Fq> + MulAssign<Fp>,
        const N: usize,
    > Mul for EvalItem<'a, Fp, Fq, N>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        use EvalItem::*;
        // Don't mult if one side is `1`
        // let (lhs, rhs) = match (self, rhs) {
        //     (Constant(c), other) | (other, Constant(c)) => {
        //         if c.is_one() {
        //             return other;
        //         }
        //         (Constant(c), other)
        //     }
        //     (lhs, rhs) => (lhs, rhs),
        // };

        match (self, rhs) {
            (Constant(lhs), Constant(rhs)) => Constant(lhs * rhs),
            (Evals(evals), Constant(c)) | (Constant(c), Evals(evals)) => match (*evals, c) {
                (FieldVariant::Fp(mut evals), FieldVariant::Fp(c)) => {
                    for v in evals.to_mut() {
                        *v *= c;
                    }
                    Evals(Box::new(FieldVariant::Fp(evals)))
                }
                (FieldVariant::Fp(evals), FieldVariant::Fq(c)) => {
                    let evals = Cow::Owned(evals.into_owned().map(|v| Fq::from(v) * c));
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
                (FieldVariant::Fq(mut evals), FieldVariant::Fp(c)) => {
                    evals.to_mut().iter_mut().for_each(|v| *v *= c);
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
                (FieldVariant::Fq(mut evals), FieldVariant::Fq(c)) => {
                    evals.to_mut().iter_mut().for_each(|v| *v *= c);
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
            },
            (Evals(lhs), Evals(rhs)) => match (*lhs, *rhs) {
                (FieldVariant::Fq(mut lhs), FieldVariant::Fp(rhs))
                | (FieldVariant::Fp(rhs), FieldVariant::Fq(mut lhs)) => {
                    zip(lhs.to_mut(), &*rhs).for_each(|(lhs, rhs)| *lhs *= *rhs);
                    Evals(Box::new(FieldVariant::Fq(lhs)))
                }
                (FieldVariant::Fq(lhs), FieldVariant::Fq(rhs)) => {
                    let (mut lhs, rhs) = if lhs.is_owned() {
                        (lhs, rhs)
                    } else {
                        (rhs, lhs)
                    };
                    zip(lhs.to_mut(), &*rhs).for_each(|(lhs, rhs)| *lhs *= *rhs);
                    Evals(Box::new(FieldVariant::Fq(lhs)))
                }
                (FieldVariant::Fp(lhs), FieldVariant::Fp(rhs)) => {
                    let (mut lhs, rhs) = if lhs.is_owned() {
                        (lhs, rhs)
                    } else {
                        (rhs, lhs)
                    };
                    zip(lhs.to_mut(), &*rhs).for_each(|(lhs, rhs)| *lhs *= *rhs);
                    Evals(Box::new(FieldVariant::Fp(lhs)))
                }
            },
        }
    }
}

impl<
        'a,
        Fp: Field,
        Fq: Field + From<Fp> + Mul<Fp, Output = Fq> + MulAssign<Fp>,
        const N: usize,
    > Div for EvalItem<'a, Fp, Fq, N>
{
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse()
    }
}

// TODO: add type safety to all the GPU shaders
// could create TypedBuffer which is metal::Buffer wrapper with phantom data
impl<'a, Fp: Field, Fq: Field, const N: usize> Pow<usize> for EvalItem<'a, Fp, Fq, N> {
    type Output = Self;

    fn pow(self, exp: usize) -> Self::Output {
        use EvalItem::*;
        match self {
            Constant(v) => Constant(v.pow(exp)),
            Evals(evals) => match *evals {
                FieldVariant::Fp(mut evals) => {
                    evals
                        .to_mut()
                        .iter_mut()
                        .for_each(|v| *v = v.pow([exp as u64]));
                    Evals(Box::new(FieldVariant::Fp(evals)))
                }
                FieldVariant::Fq(mut evals) => {
                    evals
                        .to_mut()
                        .iter_mut()
                        .for_each(|v| *v = v.pow([exp as u64]));
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
            },
        }
    }
}

impl<'a, Fp: Field, Fq: Field, const N: usize> Neg for EvalItem<'a, Fp, Fq, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        use EvalItem::*;
        match self {
            Constant(v) => Constant(-v),
            Evals(evals) => match *evals {
                FieldVariant::Fp(mut evals) => {
                    evals.to_mut().iter_mut().for_each(|v| *v = -*v);
                    Evals(Box::new(FieldVariant::Fp(evals)))
                }
                FieldVariant::Fq(mut evals) => {
                    evals.to_mut().iter_mut().for_each(|v| *v = -*v);
                    Evals(Box::new(FieldVariant::Fq(evals)))
                }
            },
        }
    }
}
