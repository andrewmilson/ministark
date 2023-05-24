use crate::air::AirConfig;
use crate::utils::divide_out_point_into;
use crate::utils::horner_evaluate;
use crate::utils::GpuAllocator;
use crate::utils::GpuVec;
use crate::Air;
use crate::Matrix;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct DeepPolyComposer<'a, A: AirConfig> {
    z: A::Fq,
    air: &'a Air<A>,
    base_trace_polys: &'a Matrix<A::Fp>,
    extension_trace_polys: Option<&'a Matrix<A::Fq>>,
    composition_trace_polys: Matrix<A::Fq>,
}

impl<'a, A: AirConfig> DeepPolyComposer<'a, A> {
    pub const fn new(
        air: &'a Air<A>,
        z: A::Fq,
        base_trace_polys: &'a Matrix<A::Fp>,
        extension_trace_polys: Option<&'a Matrix<A::Fq>>,
        composition_trace_polys: Matrix<A::Fq>,
    ) -> Self {
        Self {
            z,
            air,
            base_trace_polys,
            extension_trace_polys,
            composition_trace_polys,
        }
    }

    /// Output is of the form `(execution_trace_evals, composition_trace_evals)`
    pub fn get_ood_evals(&mut self) -> (Vec<A::Fq>, Vec<A::Fq>) {
        let Self {
            z,
            air,
            base_trace_polys,
            extension_trace_polys,
            composition_trace_polys,
            ..
        } = self;

        let trace_domain = air.trace_domain();
        let g = trace_domain.group_gen();
        let g_inv = trace_domain.group_gen_inv();

        // generate ood evaluations for the execution trace polynomials
        let base_column_range = Air::<A>::base_column_range();
        let extension_column_range = Air::<A>::extension_column_range();
        let execution_trace_evals = ark_std::cfg_into_iter!(air.trace_arguments())
            .map(|(column, offset)| {
                let generator = if offset >= 0 { g } else { g_inv };
                let offset = offset.unsigned_abs() as u64;
                let x = *z * generator.pow([offset]);
                if base_column_range.contains(&column) {
                    let coeffs = &base_trace_polys[column];
                    horner_evaluate(coeffs, &x)
                } else if extension_column_range.contains(&column) {
                    let coeffs = &extension_trace_polys.unwrap()[column - A::NUM_BASE_COLUMNS];
                    horner_evaluate(coeffs, &x)
                } else {
                    panic!(
                        "column is {column} but there are only {} columns",
                        A::NUM_BASE_COLUMNS + A::NUM_EXTENSION_COLUMNS
                    )
                }
            })
            .collect();

        // generate ood evaluations for the composition trace polynomials
        let z_n = self.z.pow([composition_trace_polys.num_cols() as u64]);
        let composition_trace_evals = ark_std::cfg_iter!(composition_trace_polys)
            .map(|column| horner_evaluate(column, &z_n))
            .collect();

        (execution_trace_evals, composition_trace_evals)
    }

    pub fn into_deep_poly(self, composition_coeffs: DeepCompositionCoeffs<A::Fq>) -> Matrix<A::Fq> {
        let Self {
            z,
            air,
            base_trace_polys,
            extension_trace_polys,
            composition_trace_polys,
            ..
        } = self;

        let DeepCompositionCoeffs {
            execution_trace: execution_trace_alphas,
            composition_trace: composition_trace_alphas,
            degree: (degree_alpha, degree_beta),
        } = composition_coeffs;

        let trace_domain = air.trace_domain();
        let g = trace_domain.group_gen();
        let g_inv = trace_domain.group_gen_inv();

        // divide out OOD point from composition trace polys
        let z_n = self.z.pow([composition_trace_polys.num_cols() as u64]);
        let composition_trace_quotients = Matrix::new(
            ark_std::cfg_into_iter!(composition_trace_polys.0)
                .zip(composition_trace_alphas)
                .map(|(coeffs, alpha)| {
                    let mut res = Vec::new_in(GpuAllocator);
                    res.resize(trace_domain.size(), A::Fq::zero());
                    divide_out_point_into(&mut res, &coeffs, &z_n, &alpha);
                    res
                })
                .collect(),
        );

        // divide out OOD points from execution trace polys
        let base_column_range = Air::<A>::base_column_range();
        let extension_column_range = Air::<A>::extension_column_range();
        // NOTE: ark_std::cfg_into_iter! doesn't work with
        // .zip() on BTreeSet but works with Vec.
        #[allow(clippy::needless_collect)]
        let trace_arguments = air.trace_arguments().into_iter().collect::<Vec<_>>();
        let execution_trace_quotients = Matrix::new(
            ark_std::cfg_into_iter!(trace_arguments)
                .zip(execution_trace_alphas)
                .map(|((col, offset), alpha)| {
                    let mut res = Vec::new_in(GpuAllocator);
                    res.resize(trace_domain.size(), A::Fq::zero());
                    let generator = if offset >= 0 { g } else { g_inv };
                    let offset = offset.unsigned_abs() as u64;
                    let x = z * generator.pow([offset]);
                    if base_column_range.contains(&col) {
                        let coeffs = &base_trace_polys[col];
                        divide_out_point_into(&mut res, coeffs, &x, &alpha);
                    } else if extension_column_range.contains(&col) {
                        let coeffs = &extension_trace_polys.unwrap()[col - A::NUM_BASE_COLUMNS];
                        divide_out_point_into(&mut res, coeffs, &x, &alpha);
                    } else {
                        panic!(
                            "column is {col} but there are only {} columns",
                            A::NUM_BASE_COLUMNS + A::NUM_EXTENSION_COLUMNS
                        )
                    }
                    res
                })
                .collect(),
        );

        let quotients = Matrix::join(vec![execution_trace_quotients, composition_trace_quotients]);
        let mut combined_coeffs = GpuVec::try_from(quotients.sum_columns()).unwrap();

        // Adjust the degree
        // P(x) * (alpha + x * beta)
        let mut last = A::Fq::zero();
        for coeff in &mut combined_coeffs {
            let tmp = *coeff;
            *coeff *= degree_alpha;
            *coeff += last * degree_beta;
            last = tmp;
        }

        Matrix::new(vec![combined_coeffs])
    }
}

pub struct DeepCompositionCoeffs<F> {
    /// Execution trace poly coefficients
    pub execution_trace: Vec<F>,
    /// Composition trace poly coefficients
    pub composition_trace: Vec<F>,
    /// Degree adjustment coefficients
    pub degree: (F, F),
}
