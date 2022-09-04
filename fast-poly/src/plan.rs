use crate::allocator::PageAlignedAllocator;
use crate::stage::NttGpuStage;
use crate::stage::Variant;
use crate::twiddles::fill_twiddles;
use crate::utils::bit_reverse;
use crate::utils::buffer_no_copy;
use crate::NttDirection;
use crate::NttOrdering;
use algebra::PrimeFelt;
use algebra::StarkFelt;
use std::time::Instant;

pub struct Ntt<'a, E> {
    command_queue: &'a metal::CommandQueueRef,
    direction: NttDirection,
    input_order: NttOrdering,
    twiddles: Vec<E, PageAlignedAllocator>,
    grid_dim: metal::MTLSize,
    threadgroup_dim: metal::MTLSize,
    n: usize,
    stages: Vec<NttGpuStage<E>>,
}

impl<'a, E: StarkFelt + PrimeFelt> Ntt<'a, E> {
    pub fn process(&mut self, buffer: &mut Vec<E, PageAlignedAllocator>) {
        let mut input_buffer = buffer_no_copy(self.command_queue.device(), buffer);
        let mut twiddles_buffer = buffer_no_copy(self.command_queue.device(), &mut self.twiddles);
        let command_buffer = self.command_queue.new_command_buffer();
        for stage in &self.stages {
            stage.encode(
                command_buffer,
                self.grid_dim,
                self.threadgroup_dim,
                &mut input_buffer,
                &mut twiddles_buffer,
            );
        }
        command_buffer.commit();
        command_buffer.wait_until_completed();
    }
}

pub struct NttPlanner {
    library: metal::Library,
    command_queue: metal::CommandQueue,
}

impl NttPlanner {
    pub fn new(device: &metal::DeviceRef) -> Self {
        let library_data = include_bytes!("metal/ntt.metallib");
        let library = device.new_library_with_data(library_data).unwrap();
        let command_queue = device.new_command_queue();
        Self {
            library,
            command_queue,
        }
    }

    pub fn plan_ntt_forward<E: StarkFelt + PrimeFelt>(
        &self,
        n: usize,
        input_order: NttOrdering,
    ) -> Ntt<E> {
        assert!(n.is_power_of_two(), "must be a power of two");
        assert!(n >= 2048);
        let direction = NttDirection::Forward;
        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new((n / 2).try_into().unwrap(), 1, 1);
        let start = Instant::now();
        let mut twiddles = Vec::with_capacity_in(n, PageAlignedAllocator);
        twiddles.resize(n / 2, E::zero());
        println!("Allocation: {:?}", start.elapsed());
        let start = Instant::now();
        fill_twiddles(&mut twiddles, n, direction);
        println!("Generation: {:?}", start.elapsed());
        let start = Instant::now();
        match input_order {
            NttOrdering::Natural => bit_reverse(&mut twiddles),
            NttOrdering::BitReversed => {}
        }
        println!("Reversal: {:?}", start.elapsed());
        let stages = match n {
            4194304 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 512, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 1024, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2048, Variant::Multiple),
            ],
            524288 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 128, Variant::Single),
                // NttGpuStage::new(&self.library, direction, n, 256, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 256, Variant::Multiple),
            ],
            262144 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 64, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 128, Variant::Multiple),
            ],
            131072 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 32, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 64, Variant::Multiple),
            ],
            65536 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 16, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 32, Variant::Multiple),
            ],
            32768 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 16, Variant::Multiple),
            ],
            16384 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 8, Variant::Multiple),
            ],
            8192 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 4, Variant::Multiple),
            ],
            4096 => vec![
                NttGpuStage::new(&self.library, direction, n, 1, Variant::Single),
                NttGpuStage::new(&self.library, direction, n, 2, Variant::Multiple),
            ],
            2048 => vec![NttGpuStage::new(
                &self.library,
                direction,
                n,
                1,
                Variant::Multiple,
            )],
            // TODO: change to an error
            _ => panic!("invalid ntt size of {n}"),
        };
        Ntt {
            command_queue: &self.command_queue,
            grid_dim,
            threadgroup_dim,
            twiddles,
            direction,
            input_order,
            n,
            stages,
        }
    }
}

impl Default for NttPlanner {
    fn default() -> Self {
        NttPlanner::new(&metal::Device::system_default().expect("no device found"))
    }
}
