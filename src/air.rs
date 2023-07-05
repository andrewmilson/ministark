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

pub struct Air<AC: AirConfig> {
    constraints: Vec<Constraint<FieldVariant<AC::Fp, AC::Fq>>>,
    composition_constraint: CompositionConstraint<FieldVariant<AC::Fp, AC::Fq>>,
    trace_len: usize,
    options: ProofOptions,
    public_inputs: AC::PublicInputs,
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
}
