use crate::air::AirConfig;
use crate::utils::divide_out_point_into;
use crate::utils::divide_out_points_into;
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
use std::iter::zip;

pub struct DeepPolyComposer<'a, A: AirConfig> {
    z: A::Fq,
    air: &'a Air<A>,
    base_trace_polys: Matrix<A::Fp>,
    extension_trace_polys: Option<Matrix<A::Fq>>,
    composition_trace_polys: Matrix<A::Fq>,
}

impl<'a, A: AirConfig> DeepPolyComposer<'a, A> {
    pub const fn new(
        air: &'a Air<A>,
        z: A::Fq,
        base_trace_polys: Matrix<A::Fp>,
        extension_trace_polys: Option<Matrix<A::Fq>>,
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
        } = self;

        let trace_domain = air.trace_domain();
        let g = trace_domain.group_gen();
        let g_inv = trace_domain.group_gen_inv();

        let num_columns = A::NUM_BASE_COLUMNS + A::NUM_EXTENSION_COLUMNS;
        let base_column_range = 0..A::NUM_BASE_COLUMNS;
        let extension_column_range = A::NUM_BASE_COLUMNS..num_columns;

        // generate ood evaluations for the execution trace polynomials
        let execution_trace_evals = ark_std::cfg_into_iter!(air.trace_arguments())
            .map(|(col_idx, offset)| {
                let generator = if offset >= 0 { g } else { g_inv };
                let offset = offset.unsigned_abs() as u64;
                let x = *z * generator.pow([offset]);
                if base_column_range.contains(&col_idx) {
                    let coeffs = &base_trace_polys[col_idx];
                    horner_evaluate(coeffs, &x)
                } else if extension_column_range.contains(&col_idx) {
                    let coeffs =
                        &extension_trace_polys.as_deref().unwrap()[col_idx - A::NUM_BASE_COLUMNS];
                    horner_evaluate(coeffs, &x)
                } else {
                    panic!("column is {col_idx} but there are only {num_columns} columns")
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

    // <https://medium.com/starkware/starkdex-deep-dive-the-stark-core-engine-497942d0f0ab>
    pub fn into_deep_poly(self, composition_coeffs: DeepCompositionCoeffs<A::Fq>) -> Matrix<A::Fq> {
        let Self {
            z,
            air,
            base_trace_polys,
            extension_trace_polys,
            composition_trace_polys,
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
        let composition_trace_quotients = ark_std::cfg_into_iter!(composition_trace_polys.0)
            .zip(composition_trace_alphas)
            .map(|(mut coeffs, alpha)| {
                divide_out_point_into(&mut coeffs, &z_n, &alpha);
                coeffs
            });

        let num_columns = A::NUM_BASE_COLUMNS + A::NUM_EXTENSION_COLUMNS;
        let base_column_range = 0..A::NUM_BASE_COLUMNS;
        let extension_column_range = A::NUM_BASE_COLUMNS..num_columns;
        let trace_arguments = air.trace_arguments();
        let execution_trace_xs_and_alphas = |col_idx| {
            let mut xs = Vec::new();
            let mut alphas = Vec::new();
            for (&(col, offset), &alpha) in zip(&trace_arguments, &execution_trace_alphas) {
                if col == col_idx {
                    let generator = if offset >= 0 { g } else { g_inv };
                    let offset = offset.unsigned_abs() as u64;
                    let x = z * generator.pow([offset]);
                    xs.push(x);
                    alphas.push(alpha);
                }
            }
            (xs, alphas)
        };

        let base_trace_quotients = ark_std::cfg_into_iter!(base_trace_polys.0)
            .zip(base_column_range)
            .map(|(coeffs, col_idx)| {
                let (xs, alphas) = execution_trace_xs_and_alphas(col_idx);
                // TODO: inefficient when Fp != Fq
                let mut coeffs = coeffs
                    .into_iter()
                    .map(A::Fq::from)
                    .collect::<Vec<_>>()
                    .to_vec_in(GpuAllocator);
                divide_out_points_into(&mut coeffs, &xs, &alphas);
                coeffs
            });

        let extension_trace_quotients =
            ark_std::cfg_into_iter!(extension_trace_polys.map_or(vec![], |t| t.0))
                .zip(extension_column_range)
                .map(|(mut coeffs, col_idx)| {
                    let (xs, alphas) = execution_trace_xs_and_alphas(col_idx);
                    divide_out_points_into(&mut coeffs, &xs, &alphas);
                    coeffs
                });

        let quotients = Matrix::new(
            composition_trace_quotients
                .chain(base_trace_quotients)
                .chain(extension_trace_quotients)
                .collect(),
        );
        let mut combined_coeffs = GpuVec::try_from(quotients.sum_columns()).unwrap();

        let chunk_size = 1 << 16;
        if degree_beta.is_zero() {
            // P(x) * alpha
            ark_std::cfg_chunks_mut!(combined_coeffs, chunk_size).for_each(|coeff_chunk| {
                for coeff in coeff_chunk {
                    *coeff *= degree_alpha;
                }
            });
        } else {
            // Adjust the degree
            // P(x) * (alpha + x * beta)
            let mut last = A::Fq::zero();
            for coeff in &mut combined_coeffs {
                let tmp = *coeff;
                *coeff *= degree_alpha;
                *coeff += last * degree_beta;
                last = tmp;
            }
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
