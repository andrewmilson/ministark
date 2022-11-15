#![cfg(target_arch = "aarch64")]

use super::GpuField;
use crate::allocator::PageAlignedAllocator;
use crate::utils::buffer_no_copy;
use crate::utils::distribute_powers;
use crate::utils::void_ptr;
use crate::GpuMulAssign;
use crate::GpuVec;
use std::marker::PhantomData;
use std::mem::size_of;

#[derive(Clone, Copy, Debug)]
pub enum Variant {
    Multiple,
    Single,
}

/// GPU FFT kernel name as declared at the bottom of `fft.metal`
fn fft_kernel_name<F: GpuField>(variant: Variant) -> String {
    format!(
        "fft_{}_{}",
        match variant {
            Variant::Multiple => "multiple",
            Variant::Single => "single",
        },
        F::field_name()
    )
}

pub struct FftGpuStage<E> {
    variant: Variant,
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    threadgroup_fft_size: usize,
    _phantom: PhantomData<E>,
}

impl<F: GpuField> FftGpuStage<F> {
    pub fn new(
        library: &metal::LibraryRef,
        n: usize,
        num_boxes: usize,
        variant: Variant,
        threadgroup_fft_size: usize,
    ) -> FftGpuStage<F> {
        use metal::MTLDataType::UInt;
        assert!(n.is_power_of_two());
        assert!(num_boxes.is_power_of_two());
        assert!(threadgroup_fft_size.is_power_of_two());
        assert!(num_boxes < n);
        assert!((2048..=1073741824).contains(&n));

        // Create the compute pipeline
        let fft_constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        let num_boxes = num_boxes as u32;
        let tg_fft_size = threadgroup_fft_size as u32;
        fft_constants.set_constant_value_at_index(void_ptr(&n), UInt, 0);
        fft_constants.set_constant_value_at_index(void_ptr(&num_boxes), UInt, 1);
        fft_constants.set_constant_value_at_index(void_ptr(&tg_fft_size), UInt, 2);
        let func = library
            .get_function(&fft_kernel_name::<F>(variant), Some(fft_constants))
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();
        let max_threadgroup_threads = pipeline.max_total_threads_per_threadgroup();
        assert!(threadgroup_fft_size <= (max_threadgroup_threads / 2) as usize);

        // each thread operates on two values each round
        let threadgroup_dim = metal::MTLSize::new((tg_fft_size / 2).try_into().unwrap(), 1, 1);
        let grid_dim = metal::MTLSize::new((n / 2).try_into().unwrap(), 1, 1);

        FftGpuStage {
            variant,
            pipeline,
            threadgroup_dim,
            grid_dim,
            threadgroup_fft_size,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        input_buffer: &mut metal::BufferRef,
        twiddles_buffer: &metal::BufferRef,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        if let Variant::Multiple = self.variant {
            let field_size = size_of::<F>();
            let num_bytes = (self.threadgroup_fft_size * field_size).try_into().unwrap();
            command_encoder.set_threadgroup_memory_length(0, num_bytes);
        }
        command_encoder.set_buffer(0, Some(input_buffer), 0);
        command_encoder.set_buffer(1, Some(twiddles_buffer), 0);
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[input_buffer]);
        command_encoder.end_encoding()
    }
}

pub struct ScaleAndNormalizeGpuStage<F> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _scale_factors: GpuVec<F>,
    scale_factors_buffer: metal::Buffer,
    _phantom: PhantomData<F>,
}

impl<F: GpuField> ScaleAndNormalizeGpuStage<F> {
    pub fn new(
        library: &metal::LibraryRef,
        command_queue: &metal::CommandQueue,
        n: usize,
        scale_factor: F,
        norm_factor: F,
    ) -> Self {
        // Create the compute pipeline
        let kernel_name = format!("mul_assign_LHS_{}_RHS_{}", F::field_name(), F::field_name());
        let func = library.get_function(&kernel_name, None).unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let mut _scale_factors = Vec::with_capacity_in(n, PageAlignedAllocator);
        _scale_factors.resize(n, norm_factor);
        if !scale_factor.is_one() {
            distribute_powers(&mut _scale_factors, scale_factor);
        }
        let scale_factors_buffer = buffer_no_copy(command_queue.device(), &_scale_factors);

        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(n.try_into().unwrap(), 1, 1);

        ScaleAndNormalizeGpuStage {
            pipeline,
            threadgroup_dim,
            grid_dim,
            _scale_factors,
            scale_factors_buffer,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        input_buffer: &mut metal::BufferRef,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(input_buffer), 0);
        command_encoder.set_buffer(1, Some(&self.scale_factors_buffer), 0);
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[input_buffer]);
        command_encoder.end_encoding()
    }
}

/// FFT stage to perform a bit reversal of an input array in place
pub struct BitReverseGpuStage<F> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<F>,
}

impl<F: GpuField> BitReverseGpuStage<F> {
    pub fn new(library: &metal::LibraryRef, n: usize) -> Self {
        use metal::MTLDataType::UInt;
        assert!(n.is_power_of_two());
        assert!((2048..=1073741824).contains(&n));

        // Create the compute pipeline
        let fft_constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        let num_boxes = 5u32;
        fft_constants.set_constant_value_at_index(void_ptr(&n), UInt, 0);
        fft_constants.set_constant_value_at_index(void_ptr(&num_boxes), UInt, 1);
        let kernel_name = format!("bit_reverse_{}", F::field_name());
        let func = library
            .get_function(&kernel_name, Some(fft_constants))
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(n.try_into().unwrap(), 1, 1);

        BitReverseGpuStage {
            pipeline,
            threadgroup_dim,
            grid_dim,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        input_buffer: &mut metal::BufferRef,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(input_buffer), 0);
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[input_buffer]);
        command_encoder.end_encoding()
    }
}

pub struct MulPowStage<Lhs, Rhs = Lhs> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<(Lhs, Rhs)>,
}

impl<LhsF: GpuField, RhsF: GpuField> MulPowStage<LhsF, RhsF>
where
    LhsF: GpuMulAssign<RhsF>,
{
    pub fn new(library: &metal::LibraryRef, n: usize) -> Self {
        // Create the compute pipeline
        let constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        constants.set_constant_value_at_index(void_ptr(&n), metal::MTLDataType::UInt, 0);
        let kernel_name = format!(
            "mul_pow_LHS_{}_RHS_{}",
            LhsF::field_name(),
            RhsF::field_name()
        );
        let func = library.get_function(&kernel_name, Some(constants)).unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(n.try_into().unwrap(), 1, 1);

        MulPowStage {
            threadgroup_dim,
            pipeline,
            grid_dim,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        dst_buffer: &mut metal::BufferRef,
        src_buffer: &metal::BufferRef,
        power: usize,
        shift: usize,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(dst_buffer), 0);
        command_encoder.set_buffer(1, Some(src_buffer), 0);
        let power = power as u32;
        command_encoder.set_bytes(2, size_of::<u32>().try_into().unwrap(), void_ptr(&power));
        let shift = shift as u32;
        command_encoder.set_bytes(3, size_of::<u32>().try_into().unwrap(), void_ptr(&shift));
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[dst_buffer]);
        command_encoder.end_encoding();
    }
}

pub struct AddAssignStage<F> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<F>,
}

impl<F: GpuField> AddAssignStage<F> {
    pub fn new(library: &metal::LibraryRef, n: usize) -> Self {
        // Create the compute pipeline
        let func = library
            .get_function(
                &format!("add_assign_LHS_{}_RHS_{}", F::field_name(), F::field_name()),
                None,
            )
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let n = n as u32;
        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(n.try_into().unwrap(), 1, 1);

        AddAssignStage {
            threadgroup_dim,
            pipeline,
            grid_dim,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        dst_buffer: &mut metal::BufferRef,
        src_buffer: &metal::BufferRef,
    ) {
        // TODO: why is `metal::MTLDispatchType::Concurrent` slower?
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(dst_buffer), 0);
        command_encoder.set_buffer(1, Some(src_buffer), 0);
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[dst_buffer]);
        command_encoder.end_encoding()
    }
}

pub struct FillBuffStage<F> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<F>,
}

impl<F: GpuField> FillBuffStage<F> {
    pub fn new(library: &metal::LibraryRef, n: usize) -> Self {
        // Create the compute pipeline
        let func = library
            .get_function(&format!("fill_buff_{}", F::field_name()), None)
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let n = n as u32;
        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(n.try_into().unwrap(), 1, 1);

        FillBuffStage {
            threadgroup_dim,
            pipeline,
            grid_dim,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        dst_buffer: &mut metal::BufferRef,
        value: F,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(dst_buffer), 0);
        command_encoder.set_bytes(1, size_of::<F>().try_into().unwrap(), void_ptr(&value));
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[dst_buffer]);
        command_encoder.end_encoding()
    }
}

pub struct GenerateTwiddlesStage<F> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<F>,
}

impl<F: GpuField> GenerateTwiddlesStage<F> {
    pub fn new(library: &metal::LibraryRef, n: usize) -> Self {
        // Create the compute pipeline
        let constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        constants.set_constant_value_at_index(
            &n as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            0,
        );
        let func = library
            .get_function(
                &format!("generate_twiddles_{}", F::field_name()),
                Some(constants),
            )
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(n.try_into().unwrap(), 1, 1);

        GenerateTwiddlesStage {
            threadgroup_dim,
            pipeline,
            grid_dim,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        dst_buffer: &mut metal::BufferRef,
        value: F,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(dst_buffer), 0);
        command_encoder.set_bytes(1, size_of::<F>().try_into().unwrap(), void_ptr(&value));
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[dst_buffer]);
        command_encoder.end_encoding()
    }
}
