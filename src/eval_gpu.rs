#![cfg(feature = "gpu")]

use crate::constraints::AlgebraicItem;
use crate::expression::Expr;
use crate::utils::FieldType;
use crate::utils::FieldVariant;
use crate::utils::GpuAllocator;
use crate::utils::GpuVec;
use crate::Matrix;
use crate::StarkExtensionOf;
use alloc::collections::BTreeMap;
use alloc::rc::Rc;
use alloc::rc::Weak;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::Field;
use core::cell::RefCell;
use core::cmp::Ordering;
use core::ops::Add;
use core::ops::Div;
use core::ops::Mul;
use core::ops::Neg;
use ministark_gpu::metal;
use ministark_gpu::prelude::*;
use ministark_gpu::stage::AddAssignConstStage;
use ministark_gpu::stage::AddIntoConstStage;
use ministark_gpu::stage::AddIntoStage;
use ministark_gpu::stage::ConvertIntoStage;
use ministark_gpu::stage::ExpInPlaceStage;
use ministark_gpu::stage::ExpIntoStage;
use ministark_gpu::stage::InverseInPlaceStage;
use ministark_gpu::stage::InverseIntoStage;
use ministark_gpu::stage::MulAssignConstStage;
use ministark_gpu::stage::MulAssignStage;
use ministark_gpu::stage::MulIntoConstStage;
use ministark_gpu::stage::MulIntoStage;
use ministark_gpu::stage::NegInPlaceStage;
use ministark_gpu::stage::NegIntoStage;
use ministark_gpu::utils::buffer_no_copy;
use ministark_gpu::GpuAdd;
use ministark_gpu::GpuFftField;
use ministark_gpu::GpuFrom;
use ministark_gpu::GpuMul;
use num_traits::Pow;

pub fn eval<Fp: GpuFftField<FftField = Fp> + FftField, Fq: StarkExtensionOf<Fp>>(
    expr: &Expr<AlgebraicItem<FieldVariant<Fp, Fq>>>,
    challenges: &[Fq],
    hints: &[Fq],
    lde_step: usize,
    x_lde: GpuVec<Fp>,
    base_trace_lde: &Matrix<Fp>,
    extension_trace_lde: Option<&Matrix<Fq>>,
) -> Matrix<Fq> {
    use AlgebraicItem::*;
    let library = &get_planner().library;
    let command_queue = &get_planner().command_queue;
    let device = command_queue.device();
    let step = lde_step as isize;
    let lde_size = x_lde.len();
    let mut x_lde = Some(x_lde);
    let lde_calculator = GpuLdeCalculator::new(library, lde_size);
    let lde_cache = RefCell::new(LdeCache::<Fp, Fq>::new(lde_size));
    let command_buffer = command_queue.new_command_buffer();

    let mut trace_ldes = Vec::new();
    let mut trace_ldes_map = BTreeMap::new();

    // TODO: this is really bad but can be changed when GPU constraint
    // evaluation is refactored.
    for lde in base_trace_lde.clone().0 {
        let gpu_buffer = buffer_no_copy(device, &lde);
        trace_ldes.push(Some(FieldVariant::Fp(Lde(lde, gpu_buffer))));
    }

    for lde in extension_trace_lde.cloned().into_iter().flatten() {
        let gpu_buffer = buffer_no_copy(device, &lde);
        trace_ldes.push(Some(FieldVariant::Fq(Lde(lde, gpu_buffer))));
    }

    let res = expr.graph_eval(&mut |leaf| match leaf {
        &Constant(v) => {
            EvaluationItem::new_constant(&lde_calculator, &lde_cache, command_buffer, v)
        }
        &Challenge(i) => EvaluationItem::new_constant(
            &lde_calculator,
            &lde_cache,
            command_buffer,
            FieldVariant::Fq(challenges[i]),
        ),
        &Hint(i) => EvaluationItem::new_constant(
            &lde_calculator,
            &lde_cache,
            command_buffer,
            FieldVariant::Fq(hints[i]),
        ),
        &Trace(i, j) => {
            // // Clippy's suggestion "Using Option::map_or_else" doesn't work because
            // // `.insert` requires mutable access to seen (which is already borrowed).
            #[allow(clippy::option_if_let_else)]
            let lde = if let Some(lde) = trace_ldes_map.get(&i) {
                Weak::upgrade(lde).unwrap()
            } else {
                let lde = lde_cache
                    .borrow_mut()
                    .add(Option::take(&mut trace_ldes[i]).unwrap());
                trace_ldes_map.insert(i, Rc::downgrade(&lde));
                lde
            };
            EvaluationItem::new_lde(&lde_calculator, &lde_cache, command_buffer, lde, j * step)
        }
        &Periodic(_col) => {
            todo!()
        }
        X => {
            // generate an LDE for the only X (we called reuse_shared_nodes)
            let mut x_lde = Option::take(&mut x_lde).unwrap();
            let buffer = buffer_mut_no_copy(device, &mut x_lde);
            let lde = lde_cache
                .borrow_mut()
                .add(FieldVariant::Fp(Lde(x_lde, buffer)));
            EvaluationItem::new_lde(&lde_calculator, &lde_cache, command_buffer, lde, 0)
        }
    });
    let item = res.item.clone();
    command_buffer.commit();
    command_buffer.wait_until_completed();
    drop(res);
    drop(lde_cache);
    Matrix::new(vec![item.into_fq_vec()])
}

/// Holds GPU shaders for performing arithmetic on LDEs
struct GpuLdeCalculator<Fp, Fq> {
    lde_size: usize,
    mul_into_const_fp: MulIntoConstStage<Fp>,
    mul_into_const_fq: MulIntoConstStage<Fq>,
    mul_into_const_fq_fp: MulIntoConstStage<Fq, Fp>,
    mul_assign_const_fp: MulAssignConstStage<Fp>,
    mul_assign_const_fq: MulAssignConstStage<Fq>,
    mul_assign_const_fq_fp: MulAssignConstStage<Fq, Fp>,
    mul_assign_fp: MulAssignStage<Fp>,
    mul_assign_fq: MulAssignStage<Fq>,
    mul_assign_fq_fp: MulAssignStage<Fq, Fp>,
    mul_into_fp: MulIntoStage<Fp>,
    mul_into_fq: MulIntoStage<Fq>,
    mul_into_fq_fp: MulIntoStage<Fq, Fp>,
    add_assign_fp: AddAssignStage<Fp>,
    add_assign_fq: AddAssignStage<Fq>,
    add_assign_fq_fp: AddAssignStage<Fq, Fp>,
    add_into_fp: AddIntoStage<Fp>,
    add_into_fq: AddIntoStage<Fq>,
    add_into_fq_fp: AddIntoStage<Fq, Fp>,
    add_into_const_fp: AddIntoConstStage<Fp>,
    add_into_const_fq: AddIntoConstStage<Fq>,
    add_into_const_fq_fp: AddIntoConstStage<Fq, Fp>,
    add_assign_const_fp: AddAssignConstStage<Fp>,
    add_assign_const_fq: AddAssignConstStage<Fq>,
    add_assign_const_fq_fp: AddAssignConstStage<Fq, Fp>,
    // TODO: this is problematic if Fp==Fq
    convert_fp_into_fq: ConvertIntoStage<Fq, Fp>,
    inverse_in_place_fp: InverseInPlaceStage<Fp>,
    inverse_into_fp: InverseIntoStage<Fp>,
    neg_in_place_fp: NegInPlaceStage<Fp>,
    neg_in_place_fq: NegInPlaceStage<Fq>,
    neg_into_fp: NegIntoStage<Fp>,
    neg_into_fq: NegIntoStage<Fq>,
    exp_in_place_fp: ExpInPlaceStage<Fp>,
    exp_in_place_fq: ExpInPlaceStage<Fq>,
    exp_into_fp: ExpIntoStage<Fp>,
    exp_into_fq: ExpIntoStage<Fq>,
}

impl<Fp: Field + GpuField, Fq: Field + GpuField + GpuMul<Fp> + GpuAdd<Fp> + GpuFrom<Fp>>
    GpuLdeCalculator<Fp, Fq>
{
    pub fn new(library: &metal::LibraryRef, lde_size: usize) -> Self {
        Self {
            lde_size,
            mul_into_const_fp: MulIntoConstStage::new(library, lde_size),
            mul_into_const_fq: MulIntoConstStage::new(library, lde_size),
            mul_into_const_fq_fp: MulIntoConstStage::new(library, lde_size),
            mul_assign_const_fp: MulAssignConstStage::new(library, lde_size),
            mul_assign_const_fq: MulAssignConstStage::new(library, lde_size),
            mul_assign_const_fq_fp: MulAssignConstStage::new(library, lde_size),
            mul_assign_fp: MulAssignStage::new(library, lde_size),
            mul_assign_fq: MulAssignStage::new(library, lde_size),
            mul_assign_fq_fp: MulAssignStage::new(library, lde_size),
            mul_into_fp: MulIntoStage::new(library, lde_size),
            mul_into_fq: MulIntoStage::new(library, lde_size),
            mul_into_fq_fp: MulIntoStage::new(library, lde_size),
            add_assign_fp: AddAssignStage::new(library, lde_size),
            add_assign_fq: AddAssignStage::new(library, lde_size),
            add_assign_fq_fp: AddAssignStage::new(library, lde_size),
            add_into_fp: AddIntoStage::new(library, lde_size),
            add_into_fq: AddIntoStage::new(library, lde_size),
            add_into_fq_fp: AddIntoStage::new(library, lde_size),
            add_into_const_fp: AddIntoConstStage::new(library, lde_size),
            add_into_const_fq: AddIntoConstStage::new(library, lde_size),
            add_into_const_fq_fp: AddIntoConstStage::new(library, lde_size),
            add_assign_const_fp: AddAssignConstStage::new(library, lde_size),
            add_assign_const_fq: AddAssignConstStage::new(library, lde_size),
            add_assign_const_fq_fp: AddAssignConstStage::new(library, lde_size),
            convert_fp_into_fq: ConvertIntoStage::new(library, lde_size),
            inverse_in_place_fp: InverseInPlaceStage::new(library, lde_size),
            inverse_into_fp: InverseIntoStage::new(library, lde_size),
            neg_in_place_fp: NegInPlaceStage::new(library, lde_size),
            neg_in_place_fq: NegInPlaceStage::new(library, lde_size),
            neg_into_fp: NegIntoStage::new(library, lde_size),
            neg_into_fq: NegIntoStage::new(library, lde_size),
            exp_in_place_fp: ExpInPlaceStage::new(library, lde_size),
            exp_in_place_fq: ExpInPlaceStage::new(library, lde_size),
            exp_into_fp: ExpIntoStage::new(library, lde_size),
            exp_into_fq: ExpIntoStage::new(library, lde_size),
        }
    }
}

/// Rust vector and its corresponding GPU buffer
pub(crate) struct Lde<F>(pub GpuVec<F>, pub metal::Buffer);

#[derive(Clone)]
enum EvaluationVariant<Fp: Field, Fq: Field> {
    Constant(FieldVariant<Fp, Fq>),
    Lde(Rc<FieldVariant<Lde<Fp>, Lde<Fq>>>, /* =offset */ isize),
}

impl<Fp: Field, Fq: Field + From<Fp>> EvaluationVariant<Fp, Fq> {
    pub fn into_fq_vec(self) -> GpuVec<Fq> {
        use EvaluationVariant::*;
        match self {
            Constant(_) => {
                // let mut x_lde = Vec::new_in(GpuAllocator);
                // x_lde.resize(lde_size, v.as_fq());
                // x_lde
                panic!()
            }
            Lde(lde, offset) => {
                let vec = match Rc::try_unwrap(lde) {
                    Err(lde) => match lde.as_ref() {
                        FieldVariant::Fp(lde) => FieldVariant::Fp(lde.0.to_vec_in(GpuAllocator)),
                        FieldVariant::Fq(lde) => FieldVariant::Fq(lde.0.to_vec_in(GpuAllocator)),
                    },
                    Ok(lde) => match lde {
                        FieldVariant::Fp(lde) => FieldVariant::Fp(lde.0),
                        FieldVariant::Fq(lde) => FieldVariant::Fq(lde.0),
                    },
                };

                let mut res = match vec {
                    FieldVariant::Fq(lde) => lde,
                    FieldVariant::Fp(lde) => {
                        let mut res = Vec::new_in(GpuAllocator);
                        lde.into_iter().for_each(|v| res.push(Fq::from(v)));
                        res
                    }
                };

                match offset.cmp(&0) {
                    Ordering::Greater => res.rotate_left(offset.unsigned_abs()),
                    Ordering::Less => res.rotate_right(offset.unsigned_abs()),
                    Ordering::Equal => {}
                }

                res
            }
        }
    }
}

#[derive(Clone)]
struct EvaluationItem<'a, Fp: Field, Fq: Field> {
    calculator: &'a GpuLdeCalculator<Fp, Fq>,
    cache: &'a RefCell<LdeCache<Fp, Fq>>,
    command_buffer: &'a metal::CommandBufferRef,
    item: EvaluationVariant<Fp, Fq>,
}

impl<'a, Fp: Field + GpuField, Fq: Field + GpuField> EvaluationItem<'a, Fp, Fq> {
    pub const fn new_constant(
        calculator: &'a GpuLdeCalculator<Fp, Fq>,
        cache: &'a RefCell<LdeCache<Fp, Fq>>,
        command_buffer: &'a metal::CommandBufferRef,
        v: FieldVariant<Fp, Fq>,
    ) -> Self {
        EvaluationItem {
            calculator,
            cache,
            command_buffer,
            item: EvaluationVariant::Constant(v),
        }
    }

    fn new_lde(
        calculator: &'a GpuLdeCalculator<Fp, Fq>,
        cache: &'a RefCell<LdeCache<Fp, Fq>>,
        command_buffer: &'a metal::CommandBufferRef,
        item: Rc<FieldVariant<Lde<Fp>, Lde<Fq>>>,
        offset: isize,
    ) -> Self {
        EvaluationItem {
            calculator,
            cache,
            command_buffer,
            item: EvaluationVariant::Lde(item, offset),
        }
    }

    fn inverse(self) -> Self {
        let calculator = self.calculator;
        let command_buffer = self.command_buffer;
        let cache = self.cache;

        let GpuLdeCalculator {
            inverse_in_place_fp,
            inverse_into_fp,
            ..
        } = calculator;

        let item = match self.item {
            EvaluationVariant::Constant(v) => EvaluationVariant::Constant(v.inverse().unwrap()),
            EvaluationVariant::Lde(lde, offset) => {
                let lde_ref_count = Rc::strong_count(&lde);
                match lde.as_ref() {
                    FieldVariant::Fp(a) => {
                        if lde_ref_count <= 2 {
                            inverse_in_place_fp.encode(command_buffer, &a.1);
                            EvaluationVariant::Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(b) = dst.as_ref() else {
                                panic!()
                            };
                            inverse_into_fp.encode(command_buffer, &b.1, &a.1);
                            EvaluationVariant::Lde(dst, offset)
                        }
                    }
                    FieldVariant::Fq(_) => todo!(),
                }
            }
        };

        EvaluationItem {
            calculator,
            cache,
            command_buffer,
            item,
        }
    }
    // fn get_lde_ref
}

impl<
        'a,
        Fp: Field + GpuField,
        Fq: Field + GpuField + Add<Fp, Output = Fq> + GpuFrom<Fp> + GpuAdd<Fp>,
    > Add for EvaluationItem<'a, Fp, Fq>
{
    type Output = Self;

    #[allow(clippy::too_many_lines)]
    fn add(self, rhs: Self) -> Self::Output {
        use EvaluationVariant::*;
        let calculator = self.calculator;
        let command_buffer = self.command_buffer;
        let cache = self.cache;

        let GpuLdeCalculator {
            lde_size,
            add_assign_const_fp,
            convert_fp_into_fq,
            add_into_const_fp,
            add_into_const_fq,
            add_into_const_fq_fp,
            add_assign_const_fq,
            add_assign_const_fq_fp,
            add_assign_fp,
            add_into_fp,
            add_assign_fq,
            add_into_fq,
            add_assign_fq_fp,
            add_into_fq_fp,
            ..
        } = calculator;

        let item = match (self.item, rhs.item) {
            (Constant(a), Constant(b)) => Constant(a + b),
            (Lde(lde, offset), Constant(b)) | (Constant(b), Lde(lde, offset)) => {
                let lde_ref_count = Rc::strong_count(&lde);
                match (lde.as_ref(), b) {
                    (FieldVariant::Fp(a), FieldVariant::Fp(b)) => {
                        if lde_ref_count <= 2 {
                            add_assign_const_fp.encode(command_buffer, &a.1, &b);
                            Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_const_fp.encode(command_buffer, &c.1, &a.1, b);
                            Lde(dst, offset)
                        }
                    }
                    (FieldVariant::Fp(a), FieldVariant::Fq(b)) => {
                        let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                        let FieldVariant::Fq(c) = dst.as_ref() else {
                            panic!()
                        };
                        convert_fp_into_fq.encode(command_buffer, &c.1, &a.1);
                        add_assign_const_fq.encode(command_buffer, &c.1, &b);
                        Lde(dst, offset)
                    }
                    (FieldVariant::Fq(a), FieldVariant::Fp(b)) => {
                        if lde_ref_count <= 2 {
                            add_assign_const_fq_fp.encode(command_buffer, &a.1, &b);
                            Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_const_fq_fp.encode(command_buffer, &c.1, &a.1, b);
                            Lde(dst, offset)
                        }
                    }
                    (FieldVariant::Fq(a), FieldVariant::Fq(b)) => {
                        if lde_ref_count <= 2 {
                            add_assign_const_fq.encode(command_buffer, &a.1, &b);
                            Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_const_fq.encode(command_buffer, &c.1, &a.1, b);
                            Lde(dst, offset)
                        }
                    }
                }
            }
            (Lde(lhs, lhs_offset), Lde(rhs, rhs_offset)) => {
                let lhs_ref_count = Rc::strong_count(&lhs);
                let rhs_ref_count = Rc::strong_count(&rhs);
                let lhs_offset = lhs_offset % *lde_size as isize;
                let rhs_offset = rhs_offset % *lde_size as isize;
                match (lhs.as_ref(), rhs.as_ref()) {
                    (FieldVariant::Fq(a), FieldVariant::Fq(b)) => {
                        if lhs_ref_count <= 2 {
                            let offset_diff = rhs_offset - lhs_offset;
                            add_assign_fq.encode(command_buffer, &a.1, &b.1, offset_diff);
                            // TODO: simplify offset to lhs_offset
                            Lde(lhs, rhs_offset - offset_diff)
                        } else if rhs_ref_count <= 2 {
                            let offset_diff = lhs_offset - rhs_offset;
                            add_assign_fq.encode(command_buffer, &b.1, &a.1, offset_diff);
                            Lde(rhs, lhs_offset - offset_diff)
                        } else {
                            let offset_diff = rhs_offset - lhs_offset;
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_fq.encode(command_buffer, &c.1, &a.1, &b.1, offset_diff);
                            Lde(dst, rhs_offset - offset_diff)
                        }
                    }
                    (FieldVariant::Fp(a), FieldVariant::Fp(b)) => {
                        if lhs_ref_count <= 2 {
                            let offset_diff = rhs_offset - lhs_offset;
                            add_assign_fp.encode(command_buffer, &a.1, &b.1, offset_diff);
                            Lde(lhs, rhs_offset - offset_diff)
                        } else if rhs_ref_count <= 2 {
                            let offset_diff = lhs_offset - rhs_offset;
                            add_assign_fp.encode(command_buffer, &b.1, &a.1, offset_diff);
                            Lde(rhs, lhs_offset - offset_diff)
                        } else {
                            let offset_diff = rhs_offset - lhs_offset;
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_fp.encode(command_buffer, &c.1, &a.1, &b.1, offset_diff);
                            Lde(dst, rhs_offset - offset_diff)
                        }
                    }
                    (FieldVariant::Fp(a), FieldVariant::Fq(b)) => {
                        let offset_diff = lhs_offset - rhs_offset;
                        if rhs_ref_count <= 2 {
                            add_assign_fq_fp.encode(command_buffer, &b.1, &a.1, offset_diff);
                            Lde(rhs, lhs_offset - offset_diff)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_fq_fp.encode(command_buffer, &c.1, &b.1, &a.1, offset_diff);
                            Lde(dst, lhs_offset - offset_diff)
                        }
                    }
                    (FieldVariant::Fq(a), FieldVariant::Fp(b)) => {
                        let offset_diff = rhs_offset - lhs_offset;
                        if lhs_ref_count <= 2 {
                            add_assign_fq_fp.encode(command_buffer, &a.1, &b.1, offset_diff);
                            Lde(lhs, rhs_offset - offset_diff)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            add_into_fq_fp.encode(command_buffer, &c.1, &a.1, &b.1, offset_diff);
                            Lde(dst, rhs_offset - offset_diff)
                        }
                    }
                }
            }
        };

        Self {
            calculator,
            cache,
            command_buffer,
            item,
        }
    }
}

impl<
        'a,
        Fp: Field + GpuField,
        Fq: Field + GpuField + GpuFrom<Fp> + GpuMul<Fp> + Mul<Fp, Output = Fq>,
    > Mul for EvaluationItem<'a, Fp, Fq>
{
    type Output = Self;

    #[allow(clippy::too_many_lines)]
    fn mul(self, rhs: Self) -> Self::Output {
        use EvaluationVariant::*;
        let calculator = self.calculator;
        let command_buffer = self.command_buffer;
        let cache = self.cache;
        let GpuLdeCalculator {
            lde_size,
            mul_assign_const_fp,
            convert_fp_into_fq,
            mul_into_const_fp,
            mul_into_const_fq,
            mul_into_const_fq_fp,
            mul_assign_const_fq,
            mul_assign_const_fq_fp,
            mul_assign_fp,
            mul_into_fp,
            mul_assign_fq,
            mul_into_fq,
            mul_assign_fq_fp,
            mul_into_fq_fp,
            ..
        } = calculator;

        let item = match (self.item, rhs.item) {
            (Constant(a), Constant(b)) => Constant(a * b),
            (Lde(lde, offset), Constant(b)) | (Constant(b), Lde(lde, offset)) => {
                let lde_ref_count = Rc::strong_count(&lde);
                match (lde.as_ref(), b) {
                    (FieldVariant::Fp(a), FieldVariant::Fp(b)) => {
                        if lde_ref_count <= 2 {
                            mul_assign_const_fp.encode(command_buffer, &a.1, b);
                            Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_const_fp.encode(command_buffer, &c.1, &a.1, &b);
                            Lde(dst, offset)
                        }
                    }
                    (FieldVariant::Fp(a), FieldVariant::Fq(b)) => {
                        let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                        let FieldVariant::Fq(c) = dst.as_ref() else {
                            panic!()
                        };
                        convert_fp_into_fq.encode(command_buffer, &c.1, &a.1);
                        mul_assign_const_fq.encode(command_buffer, &c.1, b);
                        Lde(dst, offset)
                    }
                    (FieldVariant::Fq(a), FieldVariant::Fp(b)) => {
                        if lde_ref_count <= 2 {
                            mul_assign_const_fq_fp.encode(command_buffer, &a.1, b);
                            Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_const_fq_fp.encode(command_buffer, &c.1, &a.1, &b);
                            Lde(dst, offset)
                        }
                    }
                    (FieldVariant::Fq(a), FieldVariant::Fq(b)) => {
                        if lde_ref_count <= 2 {
                            mul_assign_const_fq.encode(command_buffer, &a.1, b);
                            Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_const_fq.encode(command_buffer, &c.1, &a.1, &b);
                            Lde(dst, offset)
                        }
                    }
                }
            }
            (Lde(lhs, lhs_offset), Lde(rhs, rhs_offset)) => {
                let lhs_ref_count = Rc::strong_count(&lhs);
                let rhs_ref_count = Rc::strong_count(&rhs);
                let lhs_offset = lhs_offset % *lde_size as isize;
                let rhs_offset = rhs_offset % *lde_size as isize;
                match (lhs.as_ref(), rhs.as_ref()) {
                    (FieldVariant::Fq(a), FieldVariant::Fq(b)) => {
                        if lhs_ref_count <= 2 {
                            let offset_diff = rhs_offset - lhs_offset;
                            mul_assign_fq.encode(command_buffer, &a.1, &b.1, offset_diff);
                            // TODO: simplify offset to lhs_offset
                            Lde(lhs, rhs_offset - offset_diff)
                        } else if rhs_ref_count <= 2 {
                            let offset_diff = lhs_offset - rhs_offset;
                            mul_assign_fq.encode(command_buffer, &b.1, &a.1, offset_diff);
                            Lde(rhs, lhs_offset - offset_diff)
                        } else {
                            let offset_diff = rhs_offset - lhs_offset;
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_fq.encode(command_buffer, &c.1, &a.1, &b.1, offset_diff);
                            Lde(dst, rhs_offset - offset_diff)
                        }
                    }
                    (FieldVariant::Fp(a), FieldVariant::Fp(b)) => {
                        if lhs_ref_count <= 2 {
                            let offset_diff = rhs_offset - lhs_offset;
                            mul_assign_fp.encode(command_buffer, &a.1, &b.1, offset_diff);
                            Lde(lhs, rhs_offset - offset_diff)
                        } else if rhs_ref_count <= 2 {
                            let offset_diff = lhs_offset - rhs_offset;
                            mul_assign_fp.encode(command_buffer, &b.1, &a.1, offset_diff);
                            Lde(rhs, lhs_offset - offset_diff)
                        } else {
                            let offset_diff = rhs_offset - lhs_offset;
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_fp.encode(command_buffer, &c.1, &a.1, &b.1, offset_diff);
                            Lde(dst, rhs_offset - offset_diff)
                        }
                    }
                    (FieldVariant::Fp(a), FieldVariant::Fq(b)) => {
                        let offset_diff = lhs_offset - rhs_offset;
                        if rhs_ref_count <= 2 {
                            mul_assign_fq_fp.encode(command_buffer, &b.1, &a.1, offset_diff);
                            Lde(rhs, lhs_offset - offset_diff)
                        } else {
                            // TODO::
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_fq_fp.encode(command_buffer, &c.1, &b.1, &a.1, offset_diff);
                            Lde(dst, lhs_offset - offset_diff)
                        }
                    }
                    (FieldVariant::Fq(a), FieldVariant::Fp(b)) => {
                        let offset_diff = rhs_offset - lhs_offset;
                        if lhs_ref_count <= 2 {
                            mul_assign_fq_fp.encode(command_buffer, &a.1, &b.1, offset_diff);
                            Lde(lhs, rhs_offset - offset_diff)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(c) = dst.as_ref() else {
                                panic!()
                            };
                            mul_into_fq_fp.encode(command_buffer, &c.1, &a.1, &b.1, offset_diff);
                            Lde(dst, rhs_offset - offset_diff)
                        }
                    }
                }
            }
        };

        Self {
            calculator,
            cache,
            command_buffer,
            item,
        }
    }
}

impl<
        'a,
        Fp: Field + GpuField,
        Fq: Field + GpuField + GpuFrom<Fp> + GpuMul<Fp> + Mul<Fp, Output = Fq>,
    > Div for EvaluationItem<'a, Fp, Fq>
{
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse()
    }
}

// TODO: add type safety to all the GPU shaders
// could create TypedBuffer which is metal::Buffer wrapper with phantom data
impl<'a, Fp: Field + GpuField, Fq: Field + GpuField> Pow<usize> for EvaluationItem<'a, Fp, Fq> {
    type Output = Self;

    fn pow(self, rhs: usize) -> Self::Output {
        let calculator = self.calculator;
        let command_buffer = self.command_buffer;
        let cache = self.cache;

        let GpuLdeCalculator {
            exp_in_place_fp,
            exp_in_place_fq,
            exp_into_fp,
            exp_into_fq,
            ..
        } = calculator;

        let item = match self.item {
            EvaluationVariant::Constant(v) => EvaluationVariant::Constant(v.pow(rhs)),
            EvaluationVariant::Lde(lde, offset) => {
                let lde_ref_count = Rc::strong_count(&lde);
                match lde.as_ref() {
                    FieldVariant::Fp(a) => {
                        if lde_ref_count <= 2 {
                            exp_in_place_fp.encode(command_buffer, &a.1, rhs);
                            EvaluationVariant::Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(b) = dst.as_ref() else {
                                panic!()
                            };
                            exp_into_fp.encode(command_buffer, &b.1, &a.1, rhs);
                            EvaluationVariant::Lde(dst, offset)
                        }
                    }
                    FieldVariant::Fq(a) => {
                        if lde_ref_count <= 2 {
                            exp_in_place_fq.encode(command_buffer, &a.1, rhs);
                            EvaluationVariant::Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(b) = dst.as_ref() else {
                                panic!()
                            };
                            exp_into_fq.encode(command_buffer, &b.1, &a.1, rhs);
                            EvaluationVariant::Lde(dst, offset)
                        }
                    }
                }
            }
        };

        Self {
            calculator,
            cache,
            command_buffer,
            item,
        }
    }
}

impl<'a, Fp: Field + GpuField, Fq: Field + GpuField> Neg for EvaluationItem<'a, Fp, Fq> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let calculator = self.calculator;
        let command_buffer = self.command_buffer;
        let cache = self.cache;

        let GpuLdeCalculator {
            neg_in_place_fp,
            neg_in_place_fq,
            neg_into_fp,
            neg_into_fq,
            ..
        } = calculator;

        let item = match self.item {
            EvaluationVariant::Constant(v) => EvaluationVariant::Constant(-v),
            EvaluationVariant::Lde(lde, offset) => {
                let lde_ref_count = Rc::strong_count(&lde);
                match lde.as_ref() {
                    FieldVariant::Fp(a) => {
                        if lde_ref_count <= 2 {
                            neg_in_place_fp.encode(command_buffer, &a.1);
                            EvaluationVariant::Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fp);
                            let FieldVariant::Fp(b) = dst.as_ref() else {
                                panic!()
                            };
                            neg_into_fp.encode(command_buffer, &b.1, &a.1);
                            EvaluationVariant::Lde(dst, offset)
                        }
                    }
                    FieldVariant::Fq(a) => {
                        if lde_ref_count <= 2 {
                            neg_in_place_fq.encode(command_buffer, &a.1);
                            EvaluationVariant::Lde(lde, offset)
                        } else {
                            let dst = cache.borrow_mut().get_buffer(FieldType::Fq);
                            let FieldVariant::Fq(b) = dst.as_ref() else {
                                panic!()
                            };
                            neg_into_fq.encode(command_buffer, &b.1, &a.1);
                            EvaluationVariant::Lde(dst, offset)
                        }
                    }
                }
            }
        };

        Self {
            calculator,
            cache,
            command_buffer,
            item,
        }
    }
}

pub struct LdeCache<Fp, Fq> {
    // TODO: make a type for vec and gpu buffer
    lde_size: usize,
    buffers: Vec<Rc<FieldVariant<Lde<Fp>, Lde<Fq>>>>,
}

impl<Fp: GpuField, Fq: GpuField> LdeCache<Fp, Fq> {
    const fn new(lde_size: usize) -> Self {
        Self {
            lde_size,
            buffers: Vec::new(),
        }
    }

    fn add(&mut self, lde: FieldVariant<Lde<Fp>, Lde<Fq>>) -> Rc<FieldVariant<Lde<Fp>, Lde<Fq>>> {
        let res = Rc::new(lde);
        self.buffers.push(Rc::clone(&res));
        res
    }

    fn get_buffer(&mut self, ty: FieldType) -> Rc<FieldVariant<Lde<Fp>, Lde<Fq>>> {
        let command_queue = &get_planner().command_queue;
        let device = command_queue.device();
        // TODO: make O(1)
        self.buffers
            .iter()
            .find_map(|lde| {
                if Rc::strong_count(lde) == 1 {
                    match (lde.as_ref(), ty) {
                        (FieldVariant::Fp(..), FieldType::Fp)
                        | (FieldVariant::Fq(..), FieldType::Fq) => Some(Rc::clone(lde)),
                        _ => None,
                    }
                } else {
                    None
                }
            })
            .unwrap_or_else(|| {
                // if a buffer can't be found in the pool allocate new memory
                let n = self.lde_size;
                let buffer = match ty {
                    FieldType::Fp => {
                        let mut buffer = GpuVec::<Fp>::with_capacity_in(n, GpuAllocator);
                        // ok because all buffers are treated as uninitialized
                        unsafe { buffer.set_len(n) }
                        let gpu_buffer = buffer_no_copy(device, &buffer);
                        FieldVariant::Fp(Lde(buffer, gpu_buffer))
                    }
                    FieldType::Fq => {
                        let mut buffer = GpuVec::<Fq>::with_capacity_in(n, GpuAllocator);
                        // ok because all buffers are treated as uninitialized
                        unsafe { buffer.set_len(n) }
                        let gpu_buffer = buffer_no_copy(device, &buffer);
                        FieldVariant::Fq(Lde(buffer, gpu_buffer))
                    }
                };

                let res = Rc::new(buffer);
                self.buffers.push(Rc::clone(&res));
                res
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::ExecutionTraceColumn;
    use crate::constraints::AlgebraicItem;
    use crate::expression::Expr;
    use crate::utils::FieldVariant;
    use crate::utils::tests::gen_fib_matrix;
    use ark_ff::One;
    use ark_ff::Zero;
    use ark_poly::EvaluationDomain;
    use ark_poly::Radix2EvaluationDomain;
    use ministark_gpu::fields::p18446744069414584321::ark::Fp;
    use ministark_gpu::fields::p3618502788666131213697322783095070105623107215331596699973092056135872020481::ark::Fp as Fp252;
    use ministark_gpu::fields::p18446744069414584321::ark::Fq3;

    #[test]
    fn evaluate_x_lde() {
        use AlgebraicItem::*;
        let lde_blowup_factor = 4;
        let trace_len = 2048;
        let five = Fp::from(5u32);
        let expr = (X.pow(3) / X - X + Constant(FieldVariant::Fp(five))).pow(21) / X;
        let lde_step = lde_blowup_factor;
        let offset = Fp::GENERATOR;
        let lde_domain = Radix2EvaluationDomain::new_coset(trace_len * lde_step, offset).unwrap();
        let x_lde = lde_domain
            .elements()
            .collect::<Vec<_>>()
            .to_vec_in(GpuAllocator);

        let result = eval::<Fp, Fp>(
            &expr.reuse_shared_nodes(),
            &[],
            &[],
            lde_step,
            x_lde,
            &Matrix(vec![]),
            None,
        );

        for (i, (v, x)) in result.0[0].iter().zip(lde_domain.elements()).enumerate() {
            assert_eq!(*v, (x.pow([2]) - x + five).pow([21]) / x, "mismatch at {i}");
        }
    }

    #[test]
    fn evaluate_x_lde_with_fp_and_fq() {
        let lde_blowup_factor = 4;
        let trace_len = 2048;
        let five = Fp::from(5u32);
        let extension_element = Fq3::from_base_prime_field_elems(&[five, five, five]).unwrap();
        let expr = 1.curr() * (0.curr() * AlgebraicItem::X) / (0.next().pow(5) - AlgebraicItem::X)
            * (1.next() * AlgebraicItem::Constant(FieldVariant::Fq(extension_element)))
            / AlgebraicItem::X;
        let lde_step = lde_blowup_factor;
        let offset = Fp::GENERATOR;
        let lde_domain = Radix2EvaluationDomain::new_coset(trace_len * lde_step, offset).unwrap();
        let x_lde = lde_domain
            .elements()
            .collect::<Vec<_>>()
            .to_vec_in(GpuAllocator);
        let base_lde = Matrix::new(vec![gen_random_col::<Fp>(lde_domain.size())]);
        let extension_lde = Matrix::new(vec![gen_random_col::<Fq3>(lde_domain.size())]);

        let result = eval::<Fp, Fq3>(
            &expr.reuse_shared_nodes(),
            &[],
            &[],
            lde_step,
            x_lde,
            &base_lde,
            Some(&extension_lde),
        );

        for i in 0..lde_domain.size() - lde_step {
            let expected = expr
                .eval(&mut |leaf| match leaf {
                    AlgebraicItem::X => FieldVariant::Fp(lde_domain.element(i)),
                    &AlgebraicItem::Constant(c) => c,
                    AlgebraicItem::Challenge(_) | AlgebraicItem::Hint(_) => unreachable!(),
                    &AlgebraicItem::Trace(c, o) => match c {
                        0 => FieldVariant::Fp(base_lde.0[0][i + o.unsigned_abs() * lde_step]),
                        1 => FieldVariant::Fq(extension_lde.0[0][i + o.unsigned_abs() * lde_step]),
                        _ => unreachable!(),
                    },
                })
                .as_fq();
            let actual = result.0[0][i];
            assert_eq!(expected, actual, "mismatch at {i}");
        }
    }

    #[test]
    fn evaluate_x_inverse_lde() {
        let lde_blowup_factor = 4;
        let trace_len = 2048;
        let expr = AlgebraicItem::Constant(FieldVariant::Fp(Fp::one())) / AlgebraicItem::X;
        let lde_step = lde_blowup_factor;
        let offset = Fp::GENERATOR;
        let lde_domain = Radix2EvaluationDomain::new_coset(trace_len * lde_step, offset).unwrap();
        let x_lde = lde_domain
            .elements()
            .collect::<Vec<_>>()
            .to_vec_in(GpuAllocator);

        let result = eval::<Fp, Fp>(
            &expr.reuse_shared_nodes(),
            &[],
            &[],
            lde_step,
            x_lde,
            &Matrix(vec![]),
            None,
        );

        for (i, (v, x)) in result.0[0].iter().zip(lde_domain.elements()).enumerate() {
            assert_eq!(*v, Fp::one() / x, "mismatch at {i}");
        }
    }

    #[test]
    fn evaluate_trace_lde() {
        let lde_blowup_factor = 1u8;
        let trace_len = 2048;
        let n = trace_len * lde_blowup_factor as usize;
        let one = AlgebraicItem::Constant(FieldVariant::Fp(Fp::one()));
        let expr = 0.next() - 1.curr() - 0.curr() + one;
        let trace = gen_fib_matrix(n);
        let lde_step = lde_blowup_factor as usize;
        let offset = Fp::GENERATOR;
        let lde_domain = Radix2EvaluationDomain::new_coset(trace_len * lde_step, offset).unwrap();
        let x_lde = lde_domain
            .elements()
            .collect::<Vec<_>>()
            .to_vec_in(GpuAllocator);

        let result = eval::<Fp, Fp>(
            &expr.reuse_shared_nodes(),
            &[],
            &[],
            lde_step,
            x_lde,
            &trace,
            None,
        );

        for v in &result.0[0][0..result.num_rows() - 1] {
            assert_eq!(*v, Fp::one());
        }
    }

    #[test]
    fn evaluate_constant_lde() {
        let lde_blowup_factor = 1u8;
        let trace_len = 2048;
        let n = trace_len * lde_blowup_factor as usize;
        let one = AlgebraicItem::Constant(FieldVariant::Fp(Fp252::one()));
        let expr = Expr::from(one) - 0.curr();
        let trace = Matrix::new(vec![vec![Fp252::one(); n].to_vec_in(GpuAllocator)]);
        let lde_step = lde_blowup_factor as usize;
        let offset = Fp252::GENERATOR;
        let lde_domain = Radix2EvaluationDomain::new_coset(trace_len * lde_step, offset).unwrap();
        let x_lde = lde_domain
            .elements()
            .collect::<Vec<_>>()
            .to_vec_in(GpuAllocator);

        let result = eval::<Fp252, Fp252>(
            &expr.reuse_shared_nodes(),
            &[],
            &[],
            lde_step,
            x_lde,
            &trace,
            None,
        );

        for v in &result.0[0][0..result.num_rows() - 1] {
            assert_eq!(*v, Fp252::zero());
        }
    }

    fn gen_random_col<F: Field>(n: usize) -> GpuVec<F> {
        let mut rng = ark_std::test_rng();
        (0..n)
            .map(|_| F::rand(&mut rng))
            .collect::<Vec<F>>()
            .to_vec_in(GpuAllocator)
    }
}
