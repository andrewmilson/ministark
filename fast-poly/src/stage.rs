use super::GpuField;
use crate::FftDirection;
use std::marker::PhantomData;

#[derive(Clone, Copy)]
pub enum Variant {
    Multiple,
    Single,
}

/// GPU kernel name
///
/// Kernel names declared at the bottom of `fft.metal`
fn fft_kernel_name<F: GpuField>(direction: FftDirection, variant: Variant) -> String {
    format!(
        "{}_{}_{}",
        match direction {
            FftDirection::Forward => "fft",
            FftDirection::Inverse => "ifft",
        },
        match variant {
            Variant::Multiple => "multiple",
            Variant::Single => "single",
        },
        F::field_name()
    )
}

pub struct FftGpuStage<E> {
    pipeline: metal::ComputePipelineState,
    n: u32,
    num_boxes: u32,
    direction: FftDirection,
    variant: Variant,
    _phantom: PhantomData<E>,
}

impl<E: GpuField> FftGpuStage<E> {
    pub fn new(
        library: &metal::LibraryRef,
        direction: FftDirection,
        n: usize,
        num_boxes: usize,
        variant: Variant,
    ) -> FftGpuStage<E> {
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
            .get_function(
                &fft_kernel_name::<E>(direction, variant),
                Some(fft_constants),
            )
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        FftGpuStage {
            pipeline,
            n,
            num_boxes,
            variant,
            direction,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        grid_dim: metal::MTLSize,
        threadgroup_dim: metal::MTLSize,
        input_buffer: &mut metal::BufferRef,
        twiddles_buffer: &mut metal::BufferRef,
    ) {
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_threadgroup_memory_length(
            0,
            (2048 * std::mem::size_of::<E>()).try_into().unwrap(),
        );
        command_encoder.set_buffer(0, Some(input_buffer), 0);
        command_encoder.set_buffer(1, Some(twiddles_buffer), 0);
        command_encoder.dispatch_threads(grid_dim, threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[input_buffer]);
        command_encoder.end_encoding()
    }
}

/// FFT stage to perform a bit reversal of an input array in place
pub struct BitReverseGpuStage<F> {
    pipeline: metal::ComputePipelineState,
    n: u32,
    _phantom: PhantomData<F>,
}

impl<E: GpuField> BitReverseGpuStage<E> {
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
                &format!("bit_reverse_{}", E::field_name()),
                Some(fft_constants),
            )
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        BitReverseGpuStage {
            pipeline,
            n,
            _phantom: PhantomData,
        }
    }

    pub fn encode(
        &self,
        command_buffer: &metal::CommandBufferRef,
        input_buffer: &mut metal::BufferRef,
    ) {
        let threadgroup_dim = metal::MTLSize::new(1024, 1, 1);
        let grid_dim = metal::MTLSize::new(self.n.try_into().unwrap(), 1, 1);
        let command_encoder = command_buffer.new_compute_command_encoder();
        command_encoder.set_compute_pipeline_state(&self.pipeline);
        command_encoder.set_buffer(0, Some(input_buffer), 0);
        command_encoder.dispatch_threads(grid_dim, threadgroup_dim);
        command_encoder.memory_barrier_with_resources(&[input_buffer]);
        command_encoder.end_encoding()
    }
}
