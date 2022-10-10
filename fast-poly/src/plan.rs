use crate::allocator::PageAlignedAllocator;
use crate::stage::BitReverseGpuStage;
use crate::stage::FftGpuStage;
use crate::stage::MulPowStage;
use crate::stage::ScaleAndNormalizeGpuStage;
use crate::stage::Variant;
use crate::twiddles::fill_twiddles;
use crate::utils::bit_reverse;
use crate::utils::buffer_mut_no_copy;
use crate::utils::copy_to_private_buffer;
use crate::FftDirection;
use crate::GpuField;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use once_cell::sync::Lazy;
use std::sync::Arc;

const LIBRARY_DATA: &[u8] = include_bytes!("metal/fft.metallib");

pub struct FftEncoder<'a, F: GpuField> {
    n: usize,
    command_queue: Arc<metal::CommandQueue>,
    twiddles_buffer: metal::Buffer,
    scale_and_normalize_stage: Option<ScaleAndNormalizeGpuStage<F>>,
    butterfly_stages: Vec<FftGpuStage<F>>,
    bit_reverse_stage: BitReverseGpuStage<F>,
    command_buffer: &'a metal::CommandBufferRef,
}

// https://github.com/gfx-rs/metal-rs/issues/40
unsafe impl<'a, F: GpuField> Send for FftEncoder<'a, F> {}
unsafe impl<'a, F: GpuField> Sync for FftEncoder<'a, F> {}

impl<'a, F: GpuField> FftEncoder<'a, F> {
    fn encode_butterfly_stages(&mut self, input_buffer: &mut metal::Buffer) {
        for stage in &self.butterfly_stages {
            stage.encode(self.command_buffer, input_buffer, &self.twiddles_buffer);
        }
    }

    fn encode_bit_reverse_stage(&mut self, input_buffer: &mut metal::Buffer) {
        self.bit_reverse_stage
            .encode(self.command_buffer, input_buffer);
    }

    fn encode_scale_stage(&mut self, input_buffer: &mut metal::Buffer) {
        if let Some(scale_stage) = &self.scale_and_normalize_stage {
            scale_stage.encode(self.command_buffer, input_buffer);
        }
    }

    // TODO: change to &mut
    pub fn execute(self) {
        self.command_buffer.commit();
        self.command_buffer.wait_until_completed();
    }
}

pub struct GpuFft<'a, F: GpuField>(FftEncoder<'a, F>);

impl<'a, F: GpuField> GpuFft<'a, F> {
    pub fn encode(&mut self, buffer: &mut Vec<F, PageAlignedAllocator>) {
        assert!(self.0.n >= buffer.len());
        buffer.resize(self.0.n, F::zero());
        let mut input_buffer = buffer_mut_no_copy(self.0.command_queue.device(), buffer);
        self.0.encode_scale_stage(&mut input_buffer);
        self.0.encode_butterfly_stages(&mut input_buffer);
        self.0.encode_bit_reverse_stage(&mut input_buffer);
    }

    pub fn execute(self) {
        self.0.execute()
    }
}

impl<'a, F: GpuField> From<Radix2EvaluationDomain<F>> for GpuFft<'a, F> {
    fn from(domain: Radix2EvaluationDomain<F>) -> Self {
        PLANNER.plan_fft(domain)
    }
}

pub struct GpuIfft<'a, F: GpuField>(FftEncoder<'a, F>);

impl<'a, F: GpuField> GpuIfft<'a, F> {
    pub fn encode(&mut self, input: &mut Vec<F, PageAlignedAllocator>) {
        assert_eq!(self.0.n, input.len());
        let mut input_buffer = buffer_mut_no_copy(self.0.command_queue.device(), input);
        self.0.encode_butterfly_stages(&mut input_buffer);
        self.0.encode_bit_reverse_stage(&mut input_buffer);
        self.0.encode_scale_stage(&mut input_buffer);
    }

    pub fn execute(self) {
        self.0.execute()
    }
}

impl<'a, F: GpuField> From<Radix2EvaluationDomain<F>> for GpuIfft<'a, F> {
    fn from(domain: Radix2EvaluationDomain<F>) -> Self {
        PLANNER.plan_ifft(domain)
    }
}

// struct MulPow<'a, F: GpuField> {
//     stage: MulPowStage<F>,
//     command_queue: Arc<metal::CommandQueue>,
//     command_buffer: &'a metal::CommandBufferRef,
// }

// impl<'a, F: GpuField> MulPow<'a, F> {
//     fn encode(
//         &mut self,
//         dst: &mut Vec<F, PageAlignedAllocator>,
//         src: &mut Vec<F, PageAlignedAllocator>,
//     ) {
//         assert_eq!(dst.len(), src.len());
//         let mut dst_buffer = buffer_mut_no_copy(self.command_queue.device(),
// dst);         let mut src_buffer =
// buffer_mut_no_copy(self.command_queue.device(), src);     }

//     fn execute() {}
// }

pub static PLANNER: Lazy<Planner> = Lazy::new(Planner::default);

pub struct Planner {
    pub library: metal::Library,
    pub command_queue: Arc<metal::CommandQueue>,
}

unsafe impl Send for Planner {}
unsafe impl Sync for Planner {}

impl Planner {
    pub fn new(device: &metal::DeviceRef) -> Self {
        let library = device.new_library_with_data(LIBRARY_DATA).unwrap();
        let command_queue = Arc::new(device.new_command_queue());
        Self {
            library,
            command_queue,
        }
    }

    pub fn plan_fft<F: GpuField>(&self, domain: Radix2EvaluationDomain<F>) -> GpuFft<F> {
        GpuFft(self.create_fft_encoder(FftDirection::Forward, domain))
    }

    pub fn plan_ifft<F: GpuField>(&self, domain: Radix2EvaluationDomain<F>) -> GpuIfft<F> {
        GpuIfft(self.create_fft_encoder(FftDirection::Inverse, domain))
    }

    // TODO: move to FftEncoder struct
    fn create_fft_encoder<F: GpuField>(
        &self,
        direction: FftDirection,
        domain: Radix2EvaluationDomain<F>,
    ) -> FftEncoder<F> {
        let n = domain.size();
        assert!(n >= 2048);

        // generate twiddles buffer
        let mut twiddles = Vec::with_capacity_in(n, PageAlignedAllocator);
        twiddles.resize(n / 2, F::zero());
        fill_twiddles(&mut twiddles, n, direction);
        bit_reverse(&mut twiddles);
        let twiddles_buffer = copy_to_private_buffer(&self.command_queue, &twiddles);

        // in-place FFT requires a bit reversal
        let bit_reverse_stage = BitReverseGpuStage::new(&self.library, n);

        // scale and normalise
        let scale_and_normalize_stage = if direction == FftDirection::Forward {
            if domain.offset.is_one() {
                None
            } else {
                Some(ScaleAndNormalizeGpuStage::new(
                    &self.library,
                    &self.command_queue,
                    n,
                    domain.offset,
                    F::one(),
                ))
            }
        } else {
            Some(ScaleAndNormalizeGpuStage::new(
                &self.library,
                &self.command_queue,
                n,
                domain.offset_inv,
                domain.size_inv,
            ))
        };

        // stages that involve an FFT butterfly
        let mut butterfly_stages = Vec::new();
        for stage in 0..n.ilog2() {
            let variant = match n >> stage {
                2048 => Variant::Multiple,
                _ => Variant::Single,
            };
            butterfly_stages.push(FftGpuStage::new(&self.library, n, 1 << stage, variant));
            if let Variant::Multiple = variant {
                break;
            }
        }

        FftEncoder {
            n,
            twiddles_buffer,
            scale_and_normalize_stage,
            butterfly_stages,
            bit_reverse_stage,
            command_queue: Arc::clone(&self.command_queue),
            command_buffer: self.command_queue.new_command_buffer(),
        }
    }
}

impl Default for Planner {
    fn default() -> Self {
        Planner::new(&metal::Device::system_default().expect("no device found"))
    }
}
