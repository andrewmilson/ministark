use crate::allocator::PageAlignedAllocator;
use crate::stage::BitReverseGpuStage;
use crate::stage::FftGpuStage;
use crate::stage::Variant;
use crate::twiddles::fill_twiddles;
use crate::utils::bit_reverse;
use crate::utils::buffer_mut_no_copy;
use crate::utils::buffer_no_copy;
use crate::FftDirection;
use crate::GpuField;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use once_cell::sync::Lazy;
use std::sync::Arc;
use std::time::Instant;

pub struct Fft<F: GpuField> {
    domain: Radix2EvaluationDomain<F>,
    command_queue: Arc<metal::CommandQueue>,
    direction: FftDirection,
    twiddles: Vec<F, PageAlignedAllocator>,
    grid_dim: metal::MTLSize,
    threadgroup_dim: metal::MTLSize,
    n: usize,
    stages: Vec<FftGpuStage<F>>,
    bit_reverse_stage: BitReverseGpuStage<F>,
}

// https://github.com/gfx-rs/metal-rs/issues/40
unsafe impl<F: GpuField> Send for Fft<F> {}
unsafe impl<F: GpuField> Sync for Fft<F> {}

impl<F: GpuField> Fft<F> {
    pub fn process(&self, buffer: &mut Vec<F, PageAlignedAllocator>) {
        let start = Instant::now();
        let mut input_buffer = buffer_mut_no_copy(self.command_queue.device(), buffer);
        let twiddles_buffer = buffer_no_copy(self.command_queue.device(), &self.twiddles);
        let command_buffer = self.command_queue.new_command_buffer();
        for stage in &self.stages {
            stage.encode(
                command_buffer,
                self.grid_dim,
                self.threadgroup_dim,
                &mut input_buffer,
                &twiddles_buffer,
            );
        }
        self.bit_reverse_stage
            .encode(command_buffer, &mut input_buffer);
        println!("Preperation: {:?}", start.elapsed());
        let start = Instant::now();
        command_buffer.commit();
        println!("Commitment: {:?}", start.elapsed());
        let start = Instant::now();
        command_buffer.wait_until_completed();
        println!("Completion: {:?}", start.elapsed());
    }
}

impl<F: GpuField> From<Radix2EvaluationDomain<F>> for Fft<F> {
    fn from(domain: Radix2EvaluationDomain<F>) -> Self {
        PLANNER.plan_fft(domain)
    }
}

const LIBRARY_DATA: &[u8] = include_bytes!("metal/fft.metallib");

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

    pub fn plan_fft<F: GpuField>(&self, domain: Radix2EvaluationDomain<F>) -> Fft<F> {
        let n = domain.size();
        assert!(n >= 2048);
        let direction = FftDirection::Forward;
        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new((n / 2).try_into().unwrap(), 1, 1);
        let start = Instant::now();
        let mut twiddles = Vec::with_capacity_in(n, PageAlignedAllocator);
        twiddles.resize(n / 2, F::zero());
        println!("Allocation: {:?}", start.elapsed());
        let start = Instant::now();
        fill_twiddles(&mut twiddles, n, direction);
        println!("Generation: {:?}", start.elapsed());
        let start = Instant::now();
        bit_reverse(&mut twiddles);
        println!("Reversal: {:?}", start.elapsed());
        let stages = match n {
            33554432 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 512, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 1024, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2048, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4096, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8192, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16384, Variant::Multiple),
            ],
            8388608 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 512, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 1024, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2048, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4096, Variant::Multiple),
            ],
            4194304 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 512, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 1024, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2048, Variant::Multiple),
            ],
            2097152 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 512, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 1024, Variant::Multiple),
            ],
            1048576 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 512, Variant::Multiple),
            ],
            524288 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                // FftGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 256, Variant::Multiple),
            ],
            262144 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 128, Variant::Multiple),
            ],
            131072 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 64, Variant::Multiple),
            ],
            65536 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 32, Variant::Multiple),
            ],
            32768 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 16, Variant::Multiple),
            ],
            16384 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 8, Variant::Multiple),
            ],
            8192 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 4, Variant::Multiple),
            ],
            4096 => vec![
                FftGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                FftGpuStage::new(&self.library, direction, n, 2, Variant::Multiple),
            ],
            2048 => vec![FftGpuStage::new(
                &self.library,
                direction,
                n,
                1,
                Variant::Multiple,
            )],
            // TODO: change to an error
            _ => panic!("invalid FFT size of {n}"),
        };
        let bit_reverse_stage = BitReverseGpuStage::new(&self.library, n);
        Fft {
            domain,
            command_queue: self.command_queue.clone(),
            grid_dim,
            threadgroup_dim,
            twiddles,
            direction,
            n,
            stages,
            bit_reverse_stage,
        }
    }
}

impl Default for Planner {
    fn default() -> Self {
        Planner::new(&metal::Device::system_default().expect("no device found"))
    }
}
