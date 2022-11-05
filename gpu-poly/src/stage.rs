#![cfg(target_arch = "aarch64")]

use super::GpuField;
use crate::allocator::PageAlignedAllocator;
use crate::utils::buffer_no_copy;
use crate::utils::distribute_powers;
use crate::GpuVec;
use std::marker::PhantomData;

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
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<E>,
}

impl<F: GpuField> FftGpuStage<F> {
    pub fn new(
        library: &metal::LibraryRef,
        n: usize,
        num_boxes: usize,
        variant: Variant,
    ) -> FftGpuStage<F> {
        assert!(n.is_power_of_two());
        assert!(num_boxes.is_power_of_two());
        assert!(num_boxes < n);
        assert!((2048..=1073741824).contains(&n));

        // Create the compute pipeline
        let fft_constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        let num_boxes = num_boxes as u32;
        fft_constants.set_constant_value_at_index(
            &n as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            0,
        );
        fft_constants.set_constant_value_at_index(
            &num_boxes as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            1,
        );
        let func = library
            .get_function(&fft_kernel_name::<F>(variant), Some(fft_constants))
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new((n / 2).try_into().unwrap(), 1, 1);

        FftGpuStage {
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
        twiddles_buffer: &metal::BufferRef,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_threadgroup_memory_length(
            0,
            (2048 * std::mem::size_of::<F>()).try_into().unwrap(),
        );
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
        let func = library
            .get_function(&format!("mul_assign_{}", F::field_name()), None)
            .unwrap();
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
        assert!(n.is_power_of_two());
        assert!((2048..=1073741824).contains(&n));

        // Create the compute pipeline
        let fft_constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        let num_boxes = 5u32;
        fft_constants.set_constant_value_at_index(
            &n as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            0,
        );
        fft_constants.set_constant_value_at_index(
            &num_boxes as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            1,
        );
        let func = library
            .get_function(
                &format!("bit_reverse_{}", F::field_name()),
                Some(fft_constants),
            )
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

pub struct MulPowStage<F> {
    pipeline: metal::ComputePipelineState,
    threadgroup_dim: metal::MTLSize,
    grid_dim: metal::MTLSize,
    _phantom: PhantomData<F>,
}

impl<F: GpuField> MulPowStage<F> {
    pub fn new(library: &metal::LibraryRef, n: usize) -> Self {
        // Create the compute pipeline
        let constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        constants.set_constant_value_at_index(
            &n as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            0,
        );
        // Create the compute pipeline
        let func = library
            .get_function(&format!("mul_pow_{}", F::field_name()), Some(constants))
            .unwrap();
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
        command_encoder.set_bytes(
            2,
            std::mem::size_of::<u32>() as u64,
            &power as *const u32 as *const std::ffi::c_void,
        );
        let shift = shift as u32;
        command_encoder.set_bytes(
            3,
            std::mem::size_of::<u32>() as u64,
            &shift as *const u32 as *const std::ffi::c_void,
        );
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[dst_buffer]);
        command_encoder.end_encoding()
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
            .get_function(&format!("add_assign_{}", F::field_name()), None)
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
        // let command_encoder = command_buffer.new_compute_command_encoder();
        let command_encoder = command_buffer
            .compute_command_encoder_with_dispatch_type(metal::MTLDispatchType::Serial);
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(dst_buffer), 0);
        command_encoder.set_bytes(
            1,
            std::mem::size_of::<F>() as u64,
            &value as *const F as *const std::ffi::c_void,
        );
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
        command_encoder.set_bytes(
            1,
            std::mem::size_of::<F>() as u64,
            &value as *const F as *const std::ffi::c_void,
        );
        command_encoder.dispatch_threads(self.grid_dim, self.threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[dst_buffer]);
        command_encoder.end_encoding()
    }
}
