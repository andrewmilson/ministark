use crate::allocator::PageAlignedAllocator;
use crate::stage::BitReverseGpuStage;
use crate::stage::FftGpuStage;
use crate::stage::PolyScaleGpuStage;
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
use std::ops::Deref;
use std::ops::DerefMut;
use std::sync::Arc;
use std::time::Instant;

const LIBRARY_DATA: &[u8] = include_bytes!("metal/fft.metallib");

pub struct Fft<'a, F: GpuField> {
    direction: FftDirection,
    command_queue: Arc<metal::CommandQueue>,
    twiddles_buffer: metal::Buffer,
    poly_scale_stage: Option<PolyScaleGpuStage<F>>,
    grid_dim: metal::MTLSize,
    threadgroup_dim: metal::MTLSize,
    butterfly_stages: Vec<FftGpuStage<F>>,
    bit_reverse_stage: BitReverseGpuStage<F>,
    command_buffer: &'a metal::CommandBufferRef,
}

// https://github.com/gfx-rs/metal-rs/issues/40
unsafe impl<'a, F: GpuField> Send for Fft<'a, F> {}
unsafe impl<'a, F: GpuField> Sync for Fft<'a, F> {}

impl<'a, F: GpuField> Fft<'a, F> {
    pub fn encode(&mut self, buffer: &mut Vec<F, PageAlignedAllocator>) {
        let start = Instant::now();
        let mut input_buffer = buffer_mut_no_copy(self.command_queue.device(), buffer);
        match self.direction {
            FftDirection::Forward => self.encode_scale_stage(&mut input_buffer),
            FftDirection::Inverse => self.encode_bit_reverse_stage(&mut input_buffer),
        }
        self.encode_butterfly_stages(&mut input_buffer);
        match self.direction {
            FftDirection::Forward => self.encode_bit_reverse_stage(&mut input_buffer),
            FftDirection::Inverse => self.encode_scale_stage(&mut input_buffer),
        }
        println!("Preperation: {:?}", start.elapsed());
    }

    fn encode_butterfly_stages(&self, input_buffer: &mut metal::Buffer) {
        for stage in &self.butterfly_stages {
            stage.encode(
                self.command_buffer,
                self.grid_dim,
                self.threadgroup_dim,
                input_buffer,
                &self.twiddles_buffer,
            );
        }
    }

    fn encode_bit_reverse_stage(&self, input_buffer: &mut metal::Buffer) {
        self.bit_reverse_stage
            .encode(self.command_buffer, input_buffer);
    }

    fn encode_scale_stage(&self, input_buffer: &mut metal::Buffer) {
        if let Some(scale_stage) = &self.poly_scale_stage {
            scale_stage.encode(self.command_buffer, input_buffer);
        }
    }

    // TODO: allow mutable
    pub fn execute(self) {
        self.command_buffer.commit();
        self.command_buffer.wait_until_completed();
    }
}

impl<'a, F: GpuField> From<Radix2EvaluationDomain<F>> for Fft<'a, F> {
    fn from(domain: Radix2EvaluationDomain<F>) -> Self {
        PLANNER.plan_fft(domain)
    }
}

pub struct Ifft<'a, F: GpuField>(Fft<'a, F>);

impl<'a, F: GpuField> Ifft<'a, F> {
    pub fn execute(self) {
        self.0.execute()
    }
}

impl<'a, F: GpuField> From<Radix2EvaluationDomain<F>> for Ifft<'a, F> {
    fn from(domain: Radix2EvaluationDomain<F>) -> Self {
        PLANNER.plan_ifft(domain)
    }
}

impl<'a, F: GpuField> DerefMut for Ifft<'a, F> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<'a, F: GpuField> Deref for Ifft<'a, F> {
    type Target = Fft<'a, F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub static PLANNER: Lazy<Planner> = Lazy::new(Planner::default);

pub struct Planner {
    library: metal::Library,
    command_queue: Arc<metal::CommandQueue>,
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

    fn plan_generic_fft<F: GpuField>(
        &self,
        direction: FftDirection,
        domain: Radix2EvaluationDomain<F>,
    ) -> Fft<F> {
        let n = domain.size();
        assert!(n >= 2048);

        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new((n / 2).try_into().unwrap(), 1, 1);

        // generate twiddles buffer
        let mut twiddles = Vec::with_capacity_in(n, PageAlignedAllocator);
        twiddles.resize(n / 2, F::zero());
        fill_twiddles(&mut twiddles, n, direction);
        bit_reverse(&mut twiddles);
        let twiddles_buffer = copy_to_private_buffer(&self.command_queue, &twiddles);

        // Cooley-Tuckey fft requires a bit reversal
        let bit_reverse_stage = BitReverseGpuStage::new(&self.library, n);

        // need to scale the polynomial if the domain is a coset
        let poly_scale_factor = match direction {
            FftDirection::Forward => domain.offset,
            FftDirection::Inverse => domain.offset_inv,
        };
        let poly_scale_stage = (!poly_scale_factor.is_one())
            .then(|| PolyScaleGpuStage::new(&self.library, &self.command_queue, n, domain.offset));

        // stages that involve an FFT butterfly
        let mut butterfly_stages = Vec::new();
        for stage in 0..n.ilog2() {
            let num_boxes = match direction {
                FftDirection::Forward => 1 << (stage),
                FftDirection::Inverse => n >> (stage + 1),
            };

            println!("stage:{stage} num_boxes:{num_boxes}");

            butterfly_stages.push(FftGpuStage::new(
                &self.library,
                direction,
                n,
                num_boxes,
                Variant::Single,
            ));

            // let box_len = n / num_boxes;
            // let variant = match box_len {
            //     2048 => Variant::Multiple,
            //     _ => Variant::Single,
            // };
            // println!("{}", box_len);
            // butterfly_stages.push(FftGpuStage::new(
            //     &self.library,
            //     direction,
            //     n,
            //     num_boxes,
            //     variant,
            // ));
            // match variant {
            //     Variant::Single => num_boxes *= 2,
            //     Variant::Multiple => break,
            // }
        }

        Fft {
            direction,
            grid_dim,
            threadgroup_dim,
            twiddles_buffer,
            poly_scale_stage,
            butterfly_stages,
            bit_reverse_stage,
            command_queue: Arc::clone(&self.command_queue),
            command_buffer: self.command_queue.new_command_buffer(),
        }
    }

    pub fn plan_fft<F: GpuField>(&self, domain: Radix2EvaluationDomain<F>) -> Fft<F> {
        self.plan_generic_fft(FftDirection::Forward, domain)
    }

    pub fn plan_ifft<F: GpuField>(&self, domain: Radix2EvaluationDomain<F>) -> Ifft<F> {
        Ifft(self.plan_generic_fft(FftDirection::Inverse, domain))
    }
}

impl Default for Planner {
    fn default() -> Self {
        Planner::new(&metal::Device::system_default().expect("no device found"))
    }
}
