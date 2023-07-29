#![cfg(all(target_arch = "aarch64", target_os = "macos"))]
#[cfg(feature = "arkworks")]
use crate::stage::BitReverseGpuStage;
#[cfg(feature = "arkworks")]
use crate::stage::FftGpuStage;
#[cfg(feature = "arkworks")]
use crate::stage::FftVariant;
use crate::stage::Rpo256AbsorbColumnsStage;
use crate::stage::Rpo256AbsorbRowsStage;
use crate::stage::Rpo256GenMerkleNodesFirstRowStage;
use crate::stage::Rpo256GenMerkleNodesRowStage;
#[cfg(feature = "arkworks")]
use crate::stage::ScaleAndNormalizeGpuStage;
use crate::utils::buffer_mut_no_copy;
use crate::utils::buffer_no_copy;
use crate::utils::is_page_aligned;
use crate::utils::page_aligned_uninit_vector;
use crate::GpuField;
use alloc::rc::Rc;
use alloc::vec::Vec;
#[cfg(feature = "arkworks")]
use ark_ff::One;
#[cfg(feature = "arkworks")]
use ark_poly::EvaluationDomain;
#[cfg(feature = "arkworks")]
use ark_poly::Radix2EvaluationDomain;
use metal::CommandBufferRef;
use once_cell::sync::Lazy;

const LIBRARY_DATA: &[u8] = include_bytes!("metal/shaders.metallib");

pub struct GpuRpo256ColumnMajor<'a, F: GpuField> {
    n: usize,
    _requires_padding: bool,
    stage: Rpo256AbsorbColumnsStage<F>,
    state: Vec<&'a [F]>,
    command_buffer: Option<&'a CommandBufferRef>,
}

impl<'a, F: GpuField + From<u32> + Copy> GpuRpo256ColumnMajor<'a, F> {
    pub const RATE: usize = 8;

    pub fn new(n: usize, requires_padding: bool) -> Self {
        Self {
            n,
            _requires_padding: requires_padding,
            stage: Rpo256AbsorbColumnsStage::new(&get_planner().library, n, requires_padding),
            state: Vec::new(),
            command_buffer: None,
        }
    }

    pub fn update(&mut self, col: &'a [F]) {
        assert!(is_page_aligned(col));
        self.state.push(col);
        if self.state.len() % Self::RATE == 0 {
            let command_buffer = get_planner().command_queue.new_command_buffer();
            #[cfg(debug_assertions)]
            command_buffer.set_label("rpo update columns");
            let state = &core::mem::take(&mut self.state)[0..8];
            self.stage.encode(command_buffer, state.try_into().unwrap());
            command_buffer.commit();
            self.command_buffer = Some(command_buffer);
        }
    }

    pub async fn finish(mut self) -> Vec<[F; 4]> {
        // Wait for all update stages to finish. Stages run sequentially so waiting for
        // the last stage to finish is sufficient
        if let Some(cb) = self.command_buffer {
            cb.wait_until_completed()
        } else {
            // TODO: error? "The zero-length input is not allowed."
        }

        // return if no padding is required
        // TODO: check self.requires_padding == false
        if self.state.is_empty() {
            return self.stage.digests;
        }

        // padding rule: "a single 1 element followed by as many zeros as are necessary
        // to make the input length a multiple of the rate." - https://eprint.iacr.org/2022/1577.pdf
        // TODO: check self.requires_padding == true
        let mut ones = unsafe { page_aligned_uninit_vector(self.n) };
        ones.fill(F::from(1));
        self.state.push(&ones);

        let mut zeros: Vec<F>;
        if self.state.len() != Self::RATE {
            // only access memory for zeros if needed
            zeros = unsafe { page_aligned_uninit_vector(self.n) };
            zeros.fill(F::from(0));
            while self.state.len() != 8 {
                self.state.push(&zeros);
            }
        }

        let planner = get_planner();
        let command_buffer = planner.command_queue.new_command_buffer();
        let state = &self.state[0..8];
        self.stage.encode(command_buffer, state.try_into().unwrap());
        command_buffer.commit();
        command_buffer.wait_until_completed();
        self.stage.digests
    }
}

pub struct GpuRpo256RowMajor<'a, F: GpuField> {
    _requires_padding: bool,
    stage: Rpo256AbsorbRowsStage<F>,
    command_buffer: Option<&'a CommandBufferRef>,
}

impl<'a, F: GpuField + From<u32> + Copy> GpuRpo256RowMajor<'a, F> {
    pub const RATE: usize = 8;

    pub fn new(n: usize, requires_padding: bool) -> Self {
        Self {
            _requires_padding: requires_padding,
            stage: Rpo256AbsorbRowsStage::new(&get_planner().library, n, requires_padding),
            command_buffer: None,
        }
    }

    pub fn update(&mut self, rows: &'a [[F; 8]]) {
        assert!(is_page_aligned(rows));
        let planner = get_planner();
        let command_buffer = planner.command_queue.new_command_buffer();
        #[cfg(debug_assertions)]
        command_buffer.set_label("rpo update rows");
        self.stage.encode(command_buffer, rows);
        command_buffer.commit();
        self.command_buffer = Some(command_buffer);
    }

    pub async fn finish(self) -> Vec<[F; 4]> {
        // Wait for all update stages to finish. Stages run sequentially so waiting for
        // the last stage to finish is sufficient
        if let Some(cb) = self.command_buffer {
            cb.wait_until_completed();
            self.stage.digests
        } else {
            // TODO: error? "The zero-length input is not allowed."
            panic!()
        }
    }
}

pub async fn gen_rpo_merkle_tree<F: GpuField + From<u32> + Copy>(leaves: &[[F; 4]]) -> Vec<[F; 4]> {
    assert!(is_page_aligned(leaves));
    let planner = get_planner();
    let num_leaves = leaves.len();
    let leaves_buffer = buffer_no_copy(planner.library.device(), leaves);
    let mut nodes = unsafe { page_aligned_uninit_vector(num_leaves) };
    // TODO: might be unnecessary. only zero first item?
    nodes.fill([F::from(0); 4]);
    let nodes_buffer = buffer_mut_no_copy(planner.library.device(), &mut nodes);

    let first_row_stage = Rpo256GenMerkleNodesFirstRowStage::<F>::new(&planner.library, num_leaves);
    let nth_row_stage = Rpo256GenMerkleNodesRowStage::<F>::new(&planner.library, num_leaves);

    let command_buffer = planner.command_queue.new_command_buffer();
    #[cfg(debug_assertions)]
    command_buffer.set_label("rpo merkle tree");
    first_row_stage.encode(command_buffer, &leaves_buffer, &nodes_buffer);
    for row in 2..=num_leaves.ilog2() {
        nth_row_stage.encode(command_buffer, &nodes_buffer, row);
    }
    command_buffer.commit();
    command_buffer.wait_until_completed();

    nodes
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
#[cfg(feature = "arkworks")]
enum FftDirection {
    /// FFT
    Forward,
    /// IFFT
    Inverse,
}

#[cfg(feature = "arkworks")]
pub struct FftEncoder<'a, F: GpuField + ark_ff::Field>
where
    F::FftField: ark_ff::FftField,
{
    n: usize,
    command_queue: Rc<metal::CommandQueue>,
    // twiddles_buffer references this memory
    // field exists to keep the memory around
    _twiddles: Vec<F::FftField>,
    twiddles_buffer: metal::Buffer,
    scale_and_normalize_stage: Option<ScaleAndNormalizeGpuStage<F, F::FftField>>,
    butterfly_stages: Vec<FftGpuStage<F>>,
    bit_reverse_stage: BitReverseGpuStage<F>,
    command_buffer: &'a metal::CommandBufferRef,
}

// // https://github.com/gfx-rs/metal-rs/issues/40
// unsafe impl<'a, F: GpuField> Send for FftEncoder<'a, F> {}
// unsafe impl<'a, F: GpuField> Sync for FftEncoder<'a, F> {}

#[cfg(feature = "arkworks")]
impl<'a, F: GpuField + ark_ff::Field> FftEncoder<'a, F>
where
    F::FftField: ark_ff::FftField,
{
    fn encode_butterfly_stages(&self, input_buffer: &mut metal::Buffer) {
        for stage in &self.butterfly_stages {
            stage.encode(self.command_buffer, input_buffer, &self.twiddles_buffer);
        }
    }

    fn encode_bit_reverse_stage(&self, input_buffer: &mut metal::Buffer) {
        self.bit_reverse_stage
            .encode(self.command_buffer, input_buffer);
    }

    fn encode_scale_stage(&self, input_buffer: &mut metal::Buffer) {
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

#[cfg(feature = "arkworks")]
pub struct GpuFft<'a, F: GpuField + ark_ff::Field>
where
    F::FftField: ark_ff::FftField,
{
    encoder: FftEncoder<'a, F>,
}

#[cfg(feature = "arkworks")]
impl<'a, F: GpuField + ark_ff::Field> GpuFft<'a, F>
where
    F::FftField: ark_ff::FftField,
{
    pub const MIN_SIZE: usize = 2048;

    fn new(encoder: FftEncoder<'a, F>) -> Self {
        GpuFft { encoder }
    }

    pub fn encode(&mut self, buffer: &mut [F]) {
        assert!(is_page_aligned(buffer));
        let encoder = &self.encoder;
        assert_eq!(encoder.n, buffer.len());
        let mut input_buffer =
            crate::utils::buffer_mut_no_copy(encoder.command_queue.device(), buffer);
        encoder.encode_scale_stage(&mut input_buffer);
        encoder.encode_butterfly_stages(&mut input_buffer);
        encoder.encode_bit_reverse_stage(&mut input_buffer);
    }

    pub fn execute(self) {
        self.encoder.execute()
    }
}

#[cfg(feature = "arkworks")]
impl<'a, F: GpuField + ark_ff::Field> From<Radix2EvaluationDomain<F::FftField>> for GpuFft<'a, F>
where
    F::FftField: ark_ff::FftField,
{
    fn from(domain: Radix2EvaluationDomain<F::FftField>) -> Self {
        let planner = get_planner();
        planner.plan_fft(domain)
    }
}

#[cfg(feature = "arkworks")]
pub struct GpuIfft<'a, F: GpuField + ark_ff::Field>
where
    F::FftField: ark_ff::FftField,
{
    encoder: FftEncoder<'a, F>,
}

#[cfg(feature = "arkworks")]
impl<'a, F: GpuField + ark_ff::Field> GpuIfft<'a, F>
where
    F::FftField: ark_ff::FftField,
{
    pub const MIN_SIZE: usize = 2048;

    fn new(encoder: FftEncoder<'a, F>) -> Self {
        GpuIfft { encoder }
    }

    pub fn encode(&mut self, input: &mut [F]) {
        assert!(is_page_aligned(input));
        let encoder = &self.encoder;
        assert_eq!(encoder.n, input.len());
        let mut input_buffer =
            crate::utils::buffer_mut_no_copy(encoder.command_queue.device(), input);
        encoder.encode_butterfly_stages(&mut input_buffer);
        encoder.encode_bit_reverse_stage(&mut input_buffer);
        encoder.encode_scale_stage(&mut input_buffer);
    }

    pub fn execute(self) {
        self.encoder.execute()
    }
}

#[cfg(feature = "arkworks")]
impl<'a, F: GpuField + ark_ff::Field> From<Radix2EvaluationDomain<F::FftField>> for GpuIfft<'a, F>
where
    F::FftField: ark_ff::FftField,
{
    fn from(domain: Radix2EvaluationDomain<F::FftField>) -> Self {
        let planner = get_planner();
        planner.plan_ifft(domain)
    }
}

static PLANNER: Lazy<Planner> = Lazy::new(Planner::default);

pub fn get_planner() -> &'static Planner {
    &PLANNER
}

pub struct Planner {
    pub library: metal::Library,
    pub command_queue: Rc<metal::CommandQueue>,
}

// TODO: unsafe
unsafe impl Send for Planner {}
unsafe impl Sync for Planner {}

impl Planner {
    pub fn new(device: &metal::DeviceRef) -> Self {
        let library = device.new_library_with_data(LIBRARY_DATA).unwrap();
        let command_queue = Rc::new(device.new_command_queue());
        Self {
            library,
            command_queue,
        }
    }

    #[cfg(feature = "arkworks")]
    pub fn plan_fft<F: GpuField + ark_ff::Field>(
        &self,
        domain: Radix2EvaluationDomain<F::FftField>,
    ) -> GpuFft<F>
    where
        F::FftField: ark_ff::FftField,
    {
        assert!(domain.size() >= GpuFft::<F>::MIN_SIZE);
        GpuFft::new(self.create_fft_encoder(FftDirection::Forward, domain))
    }

    #[cfg(feature = "arkworks")]
    pub fn plan_ifft<F: GpuField + ark_ff::Field>(
        &self,
        domain: Radix2EvaluationDomain<F::FftField>,
    ) -> GpuIfft<F>
    where
        F::FftField: ark_ff::FftField,
    {
        assert!(domain.size() >= GpuIfft::<F>::MIN_SIZE);
        GpuIfft::new(self.create_fft_encoder(FftDirection::Inverse, domain))
    }

    // TODO: move to FftEncoder struct
    #[cfg(feature = "arkworks")]
    fn create_fft_encoder<F: GpuField + ark_ff::Field>(
        &self,
        direction: FftDirection,
        domain: Radix2EvaluationDomain<F::FftField>,
    ) -> FftEncoder<F>
    where
        F::FftField: ark_ff::FftField,
    {
        let n = domain.size();
        let device = self.command_queue.device();

        let root = match direction {
            FftDirection::Forward => domain.group_gen,
            FftDirection::Inverse => domain.group_gen_inv,
        };

        // generate twiddles buffer
        let mut _twiddles = unsafe { page_aligned_uninit_vector(n / 2) };
        crate::utils::fill_twiddles(&mut _twiddles, root);
        crate::utils::bit_reverse(&mut _twiddles);
        let twiddles_buffer = crate::utils::buffer_no_copy(device, &_twiddles);

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
                    F::FftField::one(),
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
        let threadgroup_mem_len = device.max_threadgroup_memory_length() as usize;
        // TODO: get max_threads_per_threadgroup from metal api. Depends on pipeline
        let threadgroup_fft_size =
            crate::utils::threadgroup_fft_size::<F>(threadgroup_mem_len, 1024);
        for stage in 0..n.ilog2() {
            let variant = if n >> stage == threadgroup_fft_size {
                FftVariant::Multiple
            } else {
                FftVariant::Single
            };

            butterfly_stages.push(FftGpuStage::new(
                &self.library,
                n,
                1 << stage,
                variant,
                threadgroup_fft_size,
            ));

            if let FftVariant::Multiple = variant {
                break;
            }
        }

        FftEncoder {
            n,
            _twiddles,
            twiddles_buffer,
            scale_and_normalize_stage,
            butterfly_stages,
            bit_reverse_stage,
            command_queue: Rc::clone(&self.command_queue),
            command_buffer: self.command_queue.new_command_buffer(),
        }
    }
}

impl Default for Planner {
    fn default() -> Self {
        Planner::new(&metal::Device::system_default().expect("no device found"))
    }
}
