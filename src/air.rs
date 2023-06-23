use crate::challenges::Challenges;
use crate::composer::DeepCompositionCoeffs;
use crate::constraints::AlgebraicItem;
use crate::constraints::CompositionConstraint;
use crate::constraints::CompositionItem;
use crate::constraints::Constraint;
use crate::expression::Expr;
use crate::hints::Hints;
use crate::random::PublicCoin;
use crate::utils::FieldVariant;
use crate::utils::GpuVec;
use crate::Matrix;
use crate::ProofOptions;
use crate::StarkExtensionOf;
use alloc::collections::BTreeSet;
use alloc::vec::Vec;
use ark_ff::FftField;
use ark_ff::UniformRand;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use core::ops::Range;
use digest::Digest;
use ministark_gpu::GpuFftField;

pub trait AirConfig: Send + Sync + Sized + 'static {
    const NUM_BASE_COLUMNS: usize;
    const NUM_EXTENSION_COLUMNS: usize = 0;
    type Fp: GpuFftField<FftField = Self::Fp> + FftField;
    type Fq: StarkExtensionOf<Self::Fp>;
    type PublicInputs: CanonicalSerialize + CanonicalDeserialize + Clone;

    fn constraints(trace_len: usize) -> Vec<Constraint<FieldVariant<Self::Fp, Self::Fq>>>;

    fn gen_hints(
        _trace_len: usize,
        _public_inputs: &Self::PublicInputs,
        _challenges: &Challenges<Self::Fq>,
    ) -> Hints<Self::Fq> {
        Hints::default()
    }

    fn domain_offset() -> Self::Fp {
        Self::Fp::GENERATOR
    }

    // TODO: maybe move this into a constraint evaluator
    #[allow(clippy::too_many_arguments)]
    fn eval_constraint(
        composition_constraint: &CompositionConstraint<FieldVariant<Self::Fp, Self::Fq>>,
        challenges: &[Self::Fq],
        hints: &[Self::Fq],
        composition_constraint_coeffs: &[Self::Fq],
        lde_step: usize,
        x_lde: GpuVec<Self::Fp>,
        base_trace_lde: &Matrix<Self::Fp>,
        extension_trace_lde: Option<&Matrix<Self::Fq>>,
    ) -> Matrix<Self::Fq> {
        let eval_expr = composition_constraint.map_leaves(&mut |leaf| match leaf {
            CompositionItem::Item(item) => *item,
            CompositionItem::CompositionCoeff(i) => {
                AlgebraicItem::Constant(FieldVariant::Fq(composition_constraint_coeffs[*i]))
            }
        });
        // TODO: add back in
        // .reuse_shared_nodes();
        // TODO: GPU constraint eval is currently slower than CPU
        // #[cfg(feature = "gpu")]
        // return crate::eval_gpu::eval::<Self::Fp, Self::Fq>(
        //     &eval_expr,
        //     challenges,
        //     hints,
        //     lde_step,
        //     x_lde,
        //     base_trace_lde,
        //     extension_trace_lde,
        // );
        // #[cfg(not(feature = "gpu"))]
        return crate::eval_cpu::eval::<Self::Fp, Self::Fq>(
            &eval_expr,
            challenges,
            hints,
            lde_step,
            &x_lde,
            base_trace_lde,
            extension_trace_lde,
        );
    }
}

pub fn trace_domain<A: AirConfig>(trace_len: usize) -> Radix2EvaluationDomain<A::Fp> {
    Radix2EvaluationDomain::new(trace_len).unwrap()
}

pub struct Air<C: AirConfig> {
    constraints: Vec<Constraint<FieldVariant<C::Fp, C::Fq>>>,
    composition_constraint: CompositionConstraint<FieldVariant<C::Fp, C::Fq>>,
    trace_len: usize,
    options: ProofOptions,
    public_inputs: C::PublicInputs,
}

impl<C: AirConfig> Air<C> {
    pub fn new(trace_len: usize, public_inputs: C::PublicInputs, options: ProofOptions) -> Self {
        let constraints = C::constraints(trace_len);
        let composition_constraint = CompositionConstraint::new(&constraints, trace_len);
        assert!(composition_constraint.ce_blowup_factor() <= options.lde_blowup_factor.into());

        Self {
            constraints,
            composition_constraint,
            trace_len,
            options,
            public_inputs,
        }
    }

    pub const fn trace_len(&self) -> usize {
        self.trace_len
    }

    pub const fn options(&self) -> ProofOptions {
        self.options
    }

    pub const fn public_inputs(&self) -> &C::PublicInputs {
        &self.public_inputs
    }

    pub const fn ce_blowup_factor(&self) -> usize {
        self.composition_constraint.ce_blowup_factor()
    }

    /// Returns a degree that all constraint polynomials must be normalized to.
    pub const fn composition_degree(&self) -> usize {
        let ce_domain_size = self.trace_len * self.ce_blowup_factor();
        ce_domain_size - 1
    }

    pub fn gen_challenges(&self, public_coin: &mut PublicCoin<impl Digest>) -> Challenges<C::Fq> {
        let mut num_challenges = 0;
        for constraint in &self.constraints {
            constraint.traverse(&mut |node| {
                if let Expr::Leaf(AlgebraicItem::Challenge(i)) = node {
                    num_challenges = core::cmp::max(num_challenges, *i + 1);
                }
            });
        }

        if num_challenges == 0 {
            Challenges::default()
        } else {
            let mut rng = public_coin.draw_rng();
            Challenges::new(&mut rng, num_challenges)
        }
    }

    pub fn gen_hints(&self, challenges: &Challenges<C::Fq>) -> Hints<C::Fq> {
        C::gen_hints(self.trace_len(), self.public_inputs(), challenges)
    }

    pub fn gen_composition_constraint_coeffs(
        &self,
        public_coin: &mut PublicCoin<impl Digest>,
    ) -> Vec<C::Fq> {
        let mut rng = public_coin.draw_rng();
        let mut num_coeffs = 0;
        self.composition_constraint.traverse(&mut |node| {
            if let Expr::Leaf(CompositionItem::CompositionCoeff(i)) = node {
                num_coeffs = num_coeffs.max(i + 1);
            }
        });
        (0..num_coeffs).map(|_| C::Fq::rand(&mut rng)).collect()
    }

    // TODO: make this generic
    /// Output is of the form `(trace_coeffs, composition_coeffs,
    /// degree_adjustment_coeffs)`
    pub fn gen_deep_composition_coeffs(
        &self,
        public_coin: &mut PublicCoin<impl Digest>,
    ) -> DeepCompositionCoeffs<C::Fq> {
        let mut rng = public_coin.draw_rng();

        // execution trace coeffs
        let mut execution_trace_coeffs = Vec::new();
        for _ in self.trace_arguments() {
            execution_trace_coeffs.push(C::Fq::rand(&mut rng));
        }

        // composition trace coeffs
        let num_composition_trace_cols = self.ce_blowup_factor();
        let mut composition_trace_coeffs = Vec::new();
        for _ in 0..num_composition_trace_cols {
            composition_trace_coeffs.push(C::Fq::rand(&mut rng));
        }

        DeepCompositionCoeffs {
            execution_trace: execution_trace_coeffs,
            composition_trace: composition_trace_coeffs,
            degree: (C::Fq::rand(&mut rng), C::Fq::rand(&mut rng)),
        }
    }

    pub fn trace_domain(&self) -> Radix2EvaluationDomain<C::Fp> {
        trace_domain::<C>(self.trace_len)
    }

    /// Low degree extension domain
    pub fn lde_domain(&self) -> Radix2EvaluationDomain<C::Fp> {
        let offset = C::domain_offset();
        let trace_len = self.trace_len();
        let lde_blowup_factor = self.lde_blowup_factor();
        Radix2EvaluationDomain::new_coset(trace_len * lde_blowup_factor, offset).unwrap()
    }

    /// Low degree extension domain
    #[inline]
    pub const fn lde_blowup_factor(&self) -> usize {
        self.options.lde_blowup_factor as usize
    }

    pub const fn composition_constraint(
        &self,
    ) -> &CompositionConstraint<FieldVariant<C::Fp, C::Fq>> {
        &self.composition_constraint
    }

    pub fn trace_arguments(&self) -> BTreeSet<(usize, isize)> {
        self.constraints
            .iter()
            .map(Constraint::trace_arguments)
            .fold(BTreeSet::new(), |a, b| &a | &b)
    }

    // #[cfg(all(feature = "std", debug_assertions))]
    // fn validate_constraints(
    //     &self,
    //     challenges: &Challenges<C::Fq>,
    //     hints: &Hints<C::Fq>,
    //     base_trace: &crate::Matrix<C::Fp>,
    //     extension_trace: Option<&crate::Matrix<C::Fq>>,
    // ) {
    //     use AlgebraicItem::*;
    //     use Expr::*;

    //     let num_execution_trace_columns = C::NUM_BASE_COLUMNS +
    // C::NUM_EXTENSION_COLUMNS;     let mut col_indicies = vec![false;
    // num_execution_trace_columns];     let mut challenge_indicies =
    // vec![false; challenges.len()];     let mut hint_indicies = vec![false;
    // hints.len()];

    //     for constraint in &self.constraints {
    //         constraint.traverse(&mut |node| match node {
    //             Leaf(Challenge(i)) => challenge_indicies[*i] = true,
    //             Leaf(Trace(i, _)) => col_indicies[*i] = true,
    //             Leaf(Hint(i)) => hint_indicies[*i] = true,
    //             _ => {}
    //         })
    //     }

    //     for (index, exists) in col_indicies.into_iter().enumerate() {
    //         if !exists {
    //             // TODO: make assertion
    //             println!("WARN: no constraints for execution trace column
    // {index}");         }
    //     }

    //     for (index, exists) in challenge_indicies.into_iter().enumerate() {
    //         if !exists {
    //             // TODO: make assertion
    //             println!("WARN: challenge at index {index} never used");
    //         }
    //     }

    //     for (index, exists) in hint_indicies.into_iter().enumerate() {
    //         if !exists {
    //             // TODO: make assertion
    //             println!("WARN: hint at index {index} never used");
    //         }
    //     }

    //     let trace_domain = self.trace_domain();
    //     let base_column_range = Self::base_column_range();
    //     let extension_column_range = Self::extension_column_range();

    //     // helper function to get a value from the execution trace
    //     let get_trace_value = |row: usize, col: usize, offset: isize| {
    //         let pos = (row as isize + offset).rem_euclid(trace_domain.size() as
    // isize) as usize;         if base_column_range.contains(&col) {
    //             FieldVariant::Fp(base_trace.0[col][pos])
    //         } else if extension_column_range.contains(&col) {
    //             let col = col - C::NUM_BASE_COLUMNS;
    //             FieldVariant::Fq(extension_trace.unwrap().0[col][pos])
    //         } else {
    //             unreachable!("requested column {col} does not exist")
    //         }
    //     };

    //     for (c_idx, constraint) in self.constraints().into_iter().enumerate() {
    //         for (row, x) in trace_domain.elements().enumerate() {
    //             let is_valid = constraint
    //                 .check(&mut |leaf| match leaf {
    //                     X => FieldVariant::Fp(x),
    //                     &Hint(i) => FieldVariant::Fq(hints[i]),
    //                     &Challenge(i) => FieldVariant::Fq(challenges[i]),
    //                     &Trace(col, offset) => get_trace_value(row, col, offset),
    //                     &Constant(c) => c,
    //                 })
    //                 .is_some();

    //             if !is_valid {
    //                 let mut vals = vec![format!("x = {x}")];
    //                 constraint.traverse(&mut |node| match *node {
    //                     // get a description of each leaf node
    //                     Leaf(Trace(col, offset)) => vals.push(format!(
    //                         "Trace(col={col:0>3}, offset={offset:0>3}) = {}",
    //                         get_trace_value(row, col, offset)
    //                     )),
    //                     Leaf(Challenge(i)) => {
    //                         vals.push(format!("Challenge({i}) = {}",
    // challenges[i]))                     }
    //                     Leaf(Hint(i)) => vals.push(format!("Hint({i}) = {}",
    // hints[i])),                     // skip tree nodes
    //                     _ => (),
    //                 });

    //                 vals.sort();
    //                 vals.dedup();

    //                 // TODO: display constraint? eprintln!("Constraint
    // is:\n{constraint}\n");                 #[cfg(feature = "std")]
    //                 eprint!("Constraint {c_idx} does not evaluate to a low degree
    // polynomial. ");                 #[cfg(feature = "std")]
    //                 eprintln!("Divide by zero occurs at row {row}.\n");
    //                 #[cfg(feature = "std")]
    //                 eprintln!("Expression values:\n{}", vals.join("\n"));
    //                 panic!();
    //             }
    //         }
    //     }
    // }

    pub const fn base_column_range() -> Range<usize> {
        0..C::NUM_BASE_COLUMNS
    }

    pub const fn extension_column_range() -> Range<usize> {
        C::NUM_BASE_COLUMNS..C::NUM_BASE_COLUMNS + C::NUM_EXTENSION_COLUMNS
    }
}
