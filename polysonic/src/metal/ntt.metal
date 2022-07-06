// A suite of functions for performing decimation-in-time (DIT)
// and decimation-in-frequency (DIF) number theory transforms.
//
// Code is has been optimized and may be difficult to reason 
// about. Refer to README.md for a higher lever explanation.

#include <metal_stdlib>
#include "felt_u128.h"
#include "permute.h"
using namespace metal;

// Number of input items being transformed
constant unsigned N [[ function_constant(0) ]];

// Number of boxes for the NTT/INTT itteration e.g.
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
template<typename E> kernel void
NttSingle(device E *vals [[ buffer(0) ]],
        constant E *twiddles [[ buffer(1) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    unsigned input_step = (N / NUM_BOXES) / 2;
    unsigned box_id = global_tid / input_step;
    unsigned target_index = box_id * input_step * 2 + (global_tid % input_step);

    E twiddle = twiddles[box_id];
    E p = vals[target_index];
    E tmp = vals[target_index + input_step];
    E q = tmp * twiddle;

    vals[target_index] = p + q;
    vals[target_index + input_step] = p - q;
}

//
constexpr unsigned inv_twiddle_idx(unsigned idx) {
    return permute_index(N, N - permute_index(N, idx)) % N;
}

// Performs a single itteration of Cooley-Tuckey radix-2 decimation-in-frequency (DIF)
//
// Inverse number theory transform. Opperates on input `vals` in bit-reversed order and
// transformed values are in natural order. Operations happen in-place.
//
// ## Example: 8 elemet array in bit-reversed order
//
// ```
// [
//     0, // 000->000
//     4, // 100->001
//     2, // 010->010
//     6, // 110->011
//     1, // 001->100
//     5, // 101->101
//     3, // 011->110
//     7  // 111->111
// ]
// ```
template<typename E> kernel void
INttSingle(device E *vals [[ buffer(0) ]],
        constant E *inv_twiddles [[ buffer(1) ]],
        unsigned global_tid [[ thread_position_in_grid ]]) {
    unsigned input_step = (N / NUM_BOXES) / 2;
    unsigned box_id = global_tid / input_step;
    unsigned target_index = 2 * box_id * input_step + (global_tid % input_step);

    E inv_twiddle = inv_twiddles[2 * box_id];
    E p = vals[target_index];
    E q = vals[target_index + input_step];
    
    vals[target_index] = p + q;
    vals[target_index + input_step] = (p - q) * inv_twiddle;
}


// Performs multiple itteration stages of Cooley-Tuckey radix-2 decimation-in-time (DIT)
//
// Template param "E" represents the type of field element.
template<typename E> kernel void
NttMultiple(device E *vals [[ buffer(0) ]],
        constant E *twiddles [[ buffer(1) ]],
        threadgroup E *shared_array [[ threadgroup(0) ]],
        unsigned group_id [[ threadgroup_position_in_grid ]],
        unsigned local_tid [[ thread_index_in_threadgroup ]]) {
#pragma unroll
    for (unsigned iteration_num = 0; iteration_num < (N / 1024 / NUM_BOXES); iteration_num++) {
        unsigned global_tid = local_tid + iteration_num * 1024;
        shared_array[global_tid] = vals[global_tid + group_id * (N / NUM_BOXES)];
    }

#pragma unroll
    for (unsigned boxes = NUM_BOXES; boxes < N; boxes *= 2) {
        unsigned input_step = (N / boxes) / 2;

#pragma unroll
        for (unsigned iteration_num = 0; iteration_num < (N / 1024 / NUM_BOXES) / 2; iteration_num++) {
            unsigned global_tid = local_tid + iteration_num * 1024;
            unsigned box_id = global_tid / input_step;
            unsigned target_index = box_id * input_step * 2 + (global_tid % input_step);

            E p = shared_array[target_index];
            E twiddle = twiddles[box_id + group_id * (boxes / NUM_BOXES)];
            E tmp = shared_array[target_index + input_step];
            E q = tmp * twiddle;

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

template [[ host_name("ntt_single_fp270497897142230380135924736767050121217") ]] kernel void
NttSingle<FP270497897142230380135924736767050121217>(
        device FP270497897142230380135924736767050121217*,
        constant FP270497897142230380135924736767050121217*,
        unsigned);
template [[ host_name("intt_single_fp270497897142230380135924736767050121217") ]] kernel void
INttSingle<FP270497897142230380135924736767050121217>(
        device FP270497897142230380135924736767050121217*,
        constant FP270497897142230380135924736767050121217*,
        unsigned);
template [[ host_name("ntt_multiple_fp270497897142230380135924736767050121217") ]] kernel void
NttMultiple<FP270497897142230380135924736767050121217>(
        device FP270497897142230380135924736767050121217*,
        constant FP270497897142230380135924736767050121217*,
        threadgroup FP270497897142230380135924736767050121217*,
        unsigned,
        unsigned);
