use crate::constraints::AlgebraicItem;
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
use ark_std::cfg_chunks_mut;
use core::iter::zip;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Div;
use core::ops::Mul;
use core::ops::MulAssign;
use core::ops::Neg;
use ministark_gpu::GpuFftField;
use num_traits::Pow;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub fn eval<Fp: GpuFftField<FftField = Fp> + FftField, Fq: StarkExtensionOf<Fp>>(
    expr: &Expr<AlgebraicItem<FieldVariant<Fp, Fq>>>,
    challenges: &[Fq],
    hints: &[Fq],
    lde_step: usize,
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
    x_lde: &[Fp],
    base_trace_lde_cols: &[&[Fp]],
    extension_trace_lde_cols: Option<&[&[Fq]]>,
    result: &mut [Fq],
) {
    use AlgebraicItem::*;
    let n = result.len();
    #[allow(clippy::cast_possible_wrap)]
    let step = lde_step as isize;
    let num_base_columns = base_trace_lde_cols.len();
    let num_extension_columns = extension_trace_lde_cols.map_or(0, <[_]>::len);
    let base_column_range = 0..num_base_columns;
    let extension_column_range = num_base_columns..num_base_columns + num_extension_columns;
    cfg_chunks_mut!(result, CHUNK_SIZE)
        .enumerate()
        .for_each(|(i, chunk)| {
            let chunk_offset = CHUNK_SIZE * i;
            let chunk_res = expr
                .graph_eval(&mut |leaf| match *leaf {
                    X => {
                        if n >= chunk_offset + CHUNK_SIZE {
                            EvalItem::Evals(Box::new(FieldVariant::Fp(Cow::Borrowed(
                                (&x_lde[chunk_offset..chunk_offset + CHUNK_SIZE])
                                    .try_into()
                                    .unwrap(),
                            ))))
                        } else {
                            let prefix = &x_lde[chunk_offset..];
                            let suffix = &x_lde[0..chunk_offset + CHUNK_SIZE - n];
                            let mut chunk = [Fp::zero(); CHUNK_SIZE];
                            chunk[0..prefix.len()].copy_from_slice(prefix);
                            chunk[prefix.len()..].copy_from_slice(suffix);
                            EvalItem::Evals(Box::new(FieldVariant::Fp(Cow::Owned(chunk))))
                        }
                    }
                    Constant(v) => EvalItem::Constant(v),
                    Challenge(i) => EvalItem::Constant(FieldVariant::Fq(challenges[i])),
                    Hint(i) => EvalItem::Constant(FieldVariant::Fq(hints[i])),
                    Trace(col_idx, row_offset) => {
                        let shift = step * row_offset;
                        #[allow(clippy::cast_possible_wrap)]
                        let position =
                            (chunk_offset as isize + shift).rem_euclid(n as isize) as usize;
                        #[allow(clippy::cast_possible_wrap)]
                        if base_column_range.contains(&col_idx) {
                            let column = &base_trace_lde_cols[col_idx];
                            if n >= position + CHUNK_SIZE {
                                EvalItem::Evals(Box::new(FieldVariant::Fp(Cow::Borrowed(
                                    (&column[position..position + CHUNK_SIZE])
                                        .try_into()
                                        .unwrap(),
                                ))))
                            } else {
                                let prefix = &column[position..];
                                let suffix = &column[0..position + CHUNK_SIZE - n];
                                let mut chunk = [Fp::zero(); CHUNK_SIZE];
                                chunk[0..prefix.len()].copy_from_slice(prefix);
                                chunk[prefix.len()..].copy_from_slice(suffix);
                                EvalItem::Evals(Box::new(FieldVariant::Fp(Cow::Owned(chunk))))
                            }
                        } else if extension_column_range.contains(&col_idx) {
                            let extension_column_offset = col_idx - num_base_columns;
                            let column =
                                &extension_trace_lde_cols.unwrap()[extension_column_offset];
                            if n >= position + CHUNK_SIZE {
                                EvalItem::Evals(Box::new(FieldVariant::Fq(Cow::Borrowed(
                                    (&column[position..position + CHUNK_SIZE])
                                        .try_into()
                                        .unwrap(),
                                ))))
                            } else {
                                let prefix = &column[position..];
                                let suffix = &column[0..position + CHUNK_SIZE - n];
                                let mut chunk = [Fq::zero(); CHUNK_SIZE];
                                chunk[0..prefix.len()].copy_from_slice(prefix);
                                chunk[prefix.len()..].copy_from_slice(suffix);
                                EvalItem::Evals(Box::new(FieldVariant::Fq(Cow::Owned(chunk))))
                            }
                        } else {
                            panic!("invalid column {col_idx}")
                        }
                    }
                })
                .into_fq_array();
            chunk.copy_from_slice(&chunk_res);
        });
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
