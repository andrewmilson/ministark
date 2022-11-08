// A suite of functions for performing decimation-in-time (DIT)
// and decimation-in-frequency (DIF) number theory transforms.
//
// Code is has been optimized and may be difficult to reason 
// about. Refer to README.md for a higher lever explanation.

#include <metal_stdlib>
#include "felt_u128.h.metal"
#include "felt_u64.h.metal"
#include "permute.h.metal"
using namespace metal;

// Number of input items being transformed
constant unsigned N [[ function_constant(0) ]];

// Number of boxes for the FFT/IFFT itteration e.g.
// ┌─────────────────────────┐ ┌─────────────────────┐ ┌───────────────────┐
// │ a[0]               a[0] │ │ a[0]           a[0] │ │ a[0]         a[0] │
// │      ╲           ╱      │ │      ╲       ╱      │ │      ❭ εϊз ❬      │
// │ a[1]  ╲         ╱  a[1] │ │ a[1]  ❭ εϊз ❬  a[1] │ │ a[1]         a[1] │
// │      ╲ ╲       ╱ ╱      │ │      ╳       ╳      │ ├───────────────────┤
// │ a[2]  ╲ ❭ εϊз ❬ ╱  a[2] │ │ a[2]  ❭ εϊз ❬  a[2] │ │ a[2]         a[2] │
// │      ╲ ╳       ╳ ╱      │ │      ╱       ╲      │ │      ❭ εϊз ❬      │
// │ a[3]  ╳ ❭ εϊз ❬ ╳  a[3] │ │ a[3]           a[3] │ │ a[3]         a[3] │
// │      ╳ ╳       ╳ ╳      │ ├─────────────────────┤ ├───────────────────┤
// │ a[4]  ╳ ❭ εϊз ❬ ╳  a[4] │ │ a[4]           a[4] │ │ a[4]         a[4] │
// │      ╱ ╳       ╳ ╲      │ │      ╲       ╱      │ │      ❭ εϊз ❬      │
// │ a[5]  ╱ ❭ εϊз ❬ ╲  a[5] │ │ a[5]  ❭ εϊз ❬  a[5] │ │ a[5]         a[5] │
// │      ╱ ╱       ╲ ╲      │ │      ╳       ╳      │ ├───────────────────┤
// │ a[6]  ╱         ╲  a[6] │ │ a[6]  ❭ εϊз ❬  a[6] │ │ a[6]         a[6] │
// │      ╱           ╲      │ │      ╱       ╲      │ │      ❭ εϊз ❬      │
// │ a[7]               a[7] │ │ a[7]           a[7] │ │ a[7]         a[7] │
// ├─────────────────────────┤ ├─────────────────────┤ ├───────────────────┤
// │       NUM_BOXES=1       │ │     NUM_BOXES=2     │ │    NUM_BOXES=4    │
// └─────────────────────────┘ └─────────────────────┘ └───────────────────┘
constant unsigned NUM_BOXES [[ function_constant(1) ]];

// Performs a single itteration of Cooley-Tuckey radix-2 decimation-in-time (DIT)
template<typename CoeffFieldT, typename TwiddleFieldT = CoeffFieldT> kernel void
FftSingle(device CoeffFieldT *vals [[ buffer(0) ]],
        constant TwiddleFieldT *twiddles [[ buffer(1) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    unsigned input_step = (N / NUM_BOXES) / 2;
    unsigned box_id = global_tid / input_step;
    unsigned target_index = box_id * input_step * 2 + (global_tid % input_step);

    TwiddleFieldT twiddle = twiddles[box_id];
    CoeffFieldT p = vals[target_index];
    CoeffFieldT tmp = vals[target_index + input_step];
    CoeffFieldT q = tmp * twiddle;

    vals[target_index] = p + q;
    vals[target_index + input_step] = p - q;
}

// Performs bit reversal.
// A useful transformation after a Cooley-Tuckey FFT to put outputs in order.
template<typename FieldT> kernel void
BitReverse(device FieldT *vals [[ buffer(0) ]],
        unsigned i [[ thread_position_in_grid ]]) {
    // ctz(N) is essentially equal to log2(N) since N is a power of two
    unsigned ri = reverse_bits(i) >> (sizeof(i) * 8 - ctz(N));

    if (i < ri) {
        // Swap positions
        FieldT tmp = vals[i];
        vals[i] = vals[ri];
        vals[ri] = tmp;
    }
}

// Performs bit reversal 
template<typename LHSFieldT, typename RHSFieldT = LHSFieldT> kernel void
MulAssign(device LHSFieldT *lhs [[ buffer(0) ]],
        constant RHSFieldT *rhs [[ buffer(1) ]],
        unsigned i [[ thread_position_in_grid ]]) {
    LHSFieldT lhs_val = lhs[i];
    RHSFieldT rhs_val = rhs[i];
    lhs[i] = lhs_val * rhs_val;
}

// Performs bit reversal 
template<typename FieldT> kernel void
AddAssign(device FieldT *lhs_vals [[ buffer(0) ]],
        constant FieldT *rhs_vals [[ buffer(1) ]],
        unsigned i [[ thread_position_in_grid ]]) {
    FieldT lhs = lhs_vals[i];
    FieldT rhs = rhs_vals[i];
    lhs_vals[i] = lhs + rhs;
}

// TODO: I want to move fft unrelated kernels into their own .metal
// lhs[i] *= rhs[i + shift] ^ exponent
template<typename LHSFieldT, typename RHSFieldT = LHSFieldT> kernel void
MulPow(device LHSFieldT *lhs [[ buffer(0) ]],
        constant RHSFieldT *rhs [[ buffer(1) ]],
        constant unsigned &exponent [[ buffer(2) ]],
        constant unsigned &shift [[ buffer(3) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    unsigned lhs_idx = global_tid;
    unsigned rhs_idx = (global_tid + shift) % N;

    LHSFieldT lhs_val = lhs[lhs_idx];
    RHSFieldT rhs_val = rhs[rhs_idx];
    lhs[lhs_idx] = lhs_val * rhs_val.pow(exponent);
}

template<typename FieldT> kernel void
FillBuff(device FieldT *dst [[ buffer(0) ]],
        constant FieldT &value [[ buffer(1) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    dst[global_tid] = value;
}

// TODO: not being used. Consider removing
template<typename FieldT> kernel void
GenerateTwiddles(device FieldT *dst [[ buffer(0) ]],
        constant FieldT &root [[ buffer(1) ]],
        unsigned i [[ thread_position_in_grid ]]) {
    unsigned ri = reverse_bits(i) >> (sizeof(i) * 8 - ctz(N));
    FieldT tmp = root;
    dst[i] = tmp.pow(ri);
}

// Performs multiple itteration stages of Cooley-Tuckey radix-2 decimation-in-time (DIT)
//
// Template param "F" represents the type of field element.
//
// TODO: Figure out poor perf reasons. Unrolls might cause instruction cache misses.
// TODO: Theoretically should be faster due to use of threadgroup memory... but it's not :(
template<typename CoeffFieldT, typename TwiddleFieldT = CoeffFieldT> kernel void
FftMultiple(device CoeffFieldT *vals [[ buffer(0) ]],
        constant TwiddleFieldT *twiddles [[ buffer(1) ]],
        threadgroup CoeffFieldT *shared_array [[ threadgroup(0) ]],
        unsigned group_id [[ threadgroup_position_in_grid ]],
        unsigned local_tid [[ thread_index_in_threadgroup ]]) {
#pragma unroll
    for (unsigned iteration_num = 0; iteration_num < (N / 1024 / NUM_BOXES); iteration_num++) {
        unsigned global_tid = local_tid + iteration_num * 1024;
        shared_array[global_tid] = vals[global_tid + group_id * (N / NUM_BOXES)];
    }

// #pragma unroll
    for (unsigned boxes = NUM_BOXES; boxes < N; boxes *= 2) {
        unsigned input_step = (N / boxes) / 2;

#pragma unroll
        for (unsigned iteration_num = 0; iteration_num < (N / 1024 / NUM_BOXES) / 2; iteration_num++) {
            unsigned global_tid = local_tid + iteration_num * 1024;
            unsigned box_id = global_tid / input_step;
            unsigned target_index = box_id * input_step * 2 + (global_tid % input_step);

            CoeffFieldT p = shared_array[target_index];
            TwiddleFieldT twiddle = twiddles[box_id + group_id * (boxes / NUM_BOXES)];
            CoeffFieldT tmp = shared_array[target_index + input_step];
            CoeffFieldT q = tmp * twiddle;

            shared_array[target_index] = p + q;
            shared_array[target_index + input_step] = p - q;
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

#pragma unroll
    for (unsigned iteration_num = 0; iteration_num < (N / 1024 / NUM_BOXES); iteration_num++) {
        // copy back to global from shared
        unsigned global_tid = local_tid + iteration_num * 1024;
        vals[global_tid + group_id * (N / NUM_BOXES)] = shared_array[global_tid];
    }
}


// ===========================================================
// FFT for Fp=270497897142230380135924736767050121217
// - 128 bit prime field
// - from Stark Anatomy series
template [[ host_name("fft_single_fp270497897142230380135924736767050121217") ]] kernel void
FftSingle<p270497897142230380135924736767050121217::Fp>(
        device p270497897142230380135924736767050121217::Fp*,
        constant p270497897142230380135924736767050121217::Fp*,
        unsigned);
template [[ host_name("fft_multiple_fp270497897142230380135924736767050121217") ]] kernel void
FftMultiple<p270497897142230380135924736767050121217::Fp>(
        device p270497897142230380135924736767050121217::Fp*,
        constant p270497897142230380135924736767050121217::Fp*,
        threadgroup p270497897142230380135924736767050121217::Fp*,
        unsigned,
        unsigned);
// ===========================================================
// FFT for Fp=18446744069414584321
// - 64 bit prime field (2^64−2^32+1 = 18446744069414584321)
// - Polygon filed (usesed by Miden and Zero)
// - Prime has many nice properties
template [[ host_name("bit_reverse_p18446744069414584321_fp") ]] kernel void
BitReverse<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        unsigned);
template [[ host_name("add_assign_LHS_p18446744069414584321_fp_RHS_p18446744069414584321_fp") ]] kernel void
AddAssign<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp*,
        unsigned);
template [[ host_name("mul_assign_LHS_p18446744069414584321_fp_RHS_p18446744069414584321_fp") ]] kernel void
MulAssign<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp*,
        unsigned);
template [[ host_name("mul_pow_LHS_p18446744069414584321_fp_RHS_p18446744069414584321_fp") ]] kernel void
MulPow<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp*,
        constant unsigned&,
        constant unsigned&,
        unsigned);
template [[ host_name("fill_buff_p18446744069414584321_fp") ]] kernel void
FillBuff<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp&,
        unsigned);
template [[ host_name("generate_twiddles_p18446744069414584321_fp") ]] kernel void
GenerateTwiddles<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp&,
        unsigned);
template [[ host_name("fft_single_p18446744069414584321_fp") ]] kernel void
FftSingle<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp*,
        unsigned);
template [[ host_name("fft_multiple_p18446744069414584321_fp") ]] kernel void
FftMultiple<p18446744069414584321::Fp>(
        device p18446744069414584321::Fp*,
        constant p18446744069414584321::Fp*,
        threadgroup p18446744069414584321::Fp*,
        unsigned,
        unsigned);
// ===========================================================
// FFT for cubic extension of Fp=18446744069414584321
template [[ host_name("bit_reverse_p18446744069414584321_fq3") ]] kernel void
BitReverse<p18446744069414584321::Fq3>(
        device p18446744069414584321::Fq3*,
        unsigned);
template [[ host_name("add_assign_LHS_p18446744069414584321_fq3_RHS_p18446744069414584321_fq3") ]] kernel void
AddAssign<p18446744069414584321::Fq3>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fq3*,
        unsigned);
template [[ host_name("mul_assign_LHS_p18446744069414584321_fq3_RHS_p18446744069414584321_fq3") ]] kernel void
MulAssign<p18446744069414584321::Fq3>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fq3*,
        unsigned);
template [[ host_name("fill_buff_p18446744069414584321_fq3") ]] kernel void
FillBuff<p18446744069414584321::Fq3>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fq3&,
        unsigned);
template [[ host_name("mul_assign_LHS_p18446744069414584321_fq3_RHS_p18446744069414584321_fp") ]] kernel void
MulAssign<p18446744069414584321::Fq3, p18446744069414584321::Fp>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fp*,
        unsigned);
template [[ host_name("mul_pow_LHS_p18446744069414584321_fq3_RHS_p18446744069414584321_fq3") ]] kernel void
MulPow<p18446744069414584321::Fq3>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fq3*,
        constant unsigned&,
        constant unsigned&,
        unsigned);
template [[ host_name("mul_pow_LHS_p18446744069414584321_fq3_RHS_p18446744069414584321_fp") ]] kernel void
MulPow<p18446744069414584321::Fq3, p18446744069414584321::Fp>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fp*,
        constant unsigned&,
        constant unsigned&,
        unsigned);
template [[ host_name("fft_single_p18446744069414584321_fq3") ]] kernel void
FftSingle<p18446744069414584321::Fq3, p18446744069414584321::Fp>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fp*,
        unsigned);
template [[ host_name("fft_multiple_p18446744069414584321_fq3") ]] kernel void
FftMultiple<p18446744069414584321::Fq3, p18446744069414584321::Fp>(
        device p18446744069414584321::Fq3*,
        constant p18446744069414584321::Fp*,
        threadgroup p18446744069414584321::Fq3*,
        unsigned,
        unsigned);
// ===========================================================
