use crate::Air;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use fast_poly::plan::Fft;
use fast_poly::GpuField;

// pub struct StarkDomain<F> {
//     trace_fft: Fft<F>,
//     ce_fft: Fft<F>,
//     lde_fft: Fft<F>,
//     ce_domain_size: usize,
//     ce_to_lde_blowup_size: usize,
// }

// impl<F: GpuField> StarkDomain<F> {
//     fn new<A: Air<Fp = F>>(air: &A) -> Self {
//         let trace_domain = air.trace_domain();
//         let ce_domain = air.ce_domain();
//         let
//         StarkDomain {}
//     }
// }
