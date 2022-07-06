use crate::fields::PrimeFelt;
use crate::fields::StarkFelt;
use crate::NttDirection;
use std::marker::PhantomData;

#[derive(Clone, Copy)]
pub enum Variant {
    Multiple,
    Single,
}

/// GPU kernel name
fn ntt_kernel_name<E: PrimeFelt>(direction: NttDirection, variant: Variant) -> String {
    format!(
        "{}_{}_fp{}",
        match direction {
            NttDirection::Forward => "ntt",
            NttDirection::Inverse => "intt",
        },
        match variant {
            Variant::Multiple => "multiple",
            Variant::Single => "single",
        },
        E::MODULUS
    )
}

pub struct NttGpuStage<E> {
    pipeline: metal::ComputePipelineState,
    n: u32,
    num_boxes: u32,
    direction: NttDirection,
    variant: Variant,
    _phantom: PhantomData<E>,
}

impl<E: StarkFelt + PrimeFelt> NttGpuStage<E> {
    pub fn new(
        library: &metal::LibraryRef,
        direction: NttDirection,
        n: usize,
        num_boxes: usize,
        variant: Variant,
    ) -> NttGpuStage<E> {
        assert!(n.is_power_of_two());
        assert!(num_boxes.is_power_of_two());
        assert!(num_boxes < n);
        assert!((2048..=1073741824).contains(&n));

        // Create the compute pipeline
        let ntt_constants = metal::FunctionConstantValues::new();
        let n = n as u32;
        let num_boxes = num_boxes as u32;
        ntt_constants.set_constant_value_at_index(
            &n as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            0,
        );
        ntt_constants.set_constant_value_at_index(
            &num_boxes as *const u32 as *const std::ffi::c_void,
            metal::MTLDataType::UInt,
            1,
        );
        let func = library
            .get_function(
                &ntt_kernel_name::<E>(direction, variant),
                Some(ntt_constants),
            )
            .unwrap();
        let pipeline = library
            .device()
            .new_compute_pipeline_state_with_function(&func)
            .unwrap();

        NttGpuStage {
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
