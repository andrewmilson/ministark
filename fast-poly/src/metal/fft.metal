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
template<typename FieldT> kernel void
FftSingle(device FieldT *vals [[ buffer(0) ]],
        constant FieldT *twiddles [[ buffer(1) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    unsigned input_step = (N / NUM_BOXES) / 2;
    unsigned box_id = global_tid / input_step;
    unsigned target_index = box_id * input_step * 2 + (global_tid % input_step);

    FieldT twiddle = twiddles[box_id];
    FieldT p = vals[target_index];
    FieldT tmp = vals[target_index + input_step];
    FieldT q = tmp * twiddle;

    vals[target_index] = p + q;
    vals[target_index + input_step] = p - q;
}

// Performs bit reversal 
template<typename FieldT> kernel void
BitReverse(device FieldT *vals [[ buffer(0) ]],
        unsigned i [[ thread_position_in_grid ]]) {
    unsigned ri = reverse_bits(i) >> (32 - log2_floor(N));

    if (i < ri) {
        // Swap positions
        FieldT tmp = vals[i];
        vals[i] = vals[ri];
        vals[ri] = tmp;
    }
}

// Performs bit reversal 
template<typename FieldT> kernel void
MulAssign(device FieldT *lhs_vals [[ buffer(0) ]],
        constant FieldT *rhs_vals [[ buffer(1) ]],
        unsigned i [[ thread_position_in_grid ]]) {
    FieldT lhs = lhs_vals[i];
    FieldT rhs = rhs_vals[i];
    lhs_vals[i] = lhs * rhs;
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

// dst[i] *= src[i] ^ exponent
// preconditions: 1 <= exponent <= 16
template<typename FieldT> kernel void
MulPow(device FieldT *dst [[ buffer(0) ]],
        constant FieldT *src [[ buffer(1) ]],
        constant unsigned &exponent [[ buffer(2) ]],
        constant unsigned &shift [[ buffer(3) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    unsigned dst_idx = global_tid;
    unsigned src_idx = global_tid + shift;
    if (src_idx >= N) {
        src_idx -= N;
    }

    FieldT src_val = src[src_idx];
    FieldT acc = src_val;
    
    // TODO: optimize
    for (unsigned i = 1; i < exponent; i++) {
        acc = acc * src_val;
    }

    FieldT dst_val = dst[dst_idx];
    dst[dst_idx] = dst_val * acc;
}


// Performs multiple itteration stages of Cooley-Tuckey radix-2 decimation-in-time (DIT)
//
// Template param "F" represents the type of field element.
//
// TODO: Figure out poor perf reasons. Unrolls might cause instruction cache misses.
// TODO: Theoretically should be faster due to use of threadgroup memory... but it's not :(
template<typename FieldT> kernel void
FftMultiple(device FieldT *vals [[ buffer(0) ]],
        constant FieldT *twiddles [[ buffer(1) ]],
        threadgroup FieldT *shared_array [[ threadgroup(0) ]],
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

            FieldT p = shared_array[target_index];
            FieldT twiddle = twiddles[box_id + group_id * (boxes / NUM_BOXES)];
            FieldT tmp = shared_array[target_index + input_step];
            FieldT q = tmp * twiddle;

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
FftSingle<FP270497897142230380135924736767050121217>(
        device FP270497897142230380135924736767050121217*,
        constant FP270497897142230380135924736767050121217*,
        unsigned);
template [[ host_name("fft_multiple_fp270497897142230380135924736767050121217") ]] kernel void
FftMultiple<FP270497897142230380135924736767050121217>(
        device FP270497897142230380135924736767050121217*,
        constant FP270497897142230380135924736767050121217*,
        threadgroup FP270497897142230380135924736767050121217*,
        unsigned,
        unsigned);
// ===========================================================
// FFT for Fp=18446744069414584321
// - 64 bit prime field (2^64−2^32+1 = 18446744069414584321)
// - Polygon filed (usesed by Miden and Zero)
// - Prime has many nice properties
template [[ host_name("bit_reverse_fp18446744069414584321") ]] kernel void
BitReverse<FP18446744069414584321>(
        device FP18446744069414584321*,
        unsigned);
template [[ host_name("add_assign_fp18446744069414584321") ]] kernel void
AddAssign<FP18446744069414584321>(
        device FP18446744069414584321*,
        constant FP18446744069414584321*,
        unsigned);
template [[ host_name("mul_assign_fp18446744069414584321") ]] kernel void
MulAssign<FP18446744069414584321>(
        device FP18446744069414584321*,
        constant FP18446744069414584321*,
        unsigned);
template [[ host_name("mul_pow_fp18446744069414584321") ]] kernel void
MulPow<FP18446744069414584321>(
        device FP18446744069414584321*,
        constant FP18446744069414584321*,
        constant unsigned&,
        constant unsigned&,
        unsigned);
template [[ host_name("fft_single_fp18446744069414584321") ]] kernel void
FftSingle<FP18446744069414584321>(
        device FP18446744069414584321*,
        constant FP18446744069414584321*,
        unsigned);
template [[ host_name("fft_multiple_fp18446744069414584321") ]] kernel void
FftMultiple<FP18446744069414584321>(
        device FP18446744069414584321*,
        constant FP18446744069414584321*,
        threadgroup FP18446744069414584321*,
        unsigned,
        unsigned);
// ===========================================================
