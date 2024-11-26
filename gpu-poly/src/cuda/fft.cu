#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <cuda.h>
#include "frarith.h"
#include <iostream>


/*==================== CUDA KERNELS ====================*/

/*
	Inner FFT loop
*/
__device__ void inplace_fft_inner(uint64_t *__restrict__ A, uint64_t *__restrict__ tw, int j, int k, int m, int n) {
	if (j + k + m / 2 < n && k < m / 2) {
		unsigned r = j + k + m / 2;
		unsigned l = j + k;
		unsigned g = k * (n / m);

		Mul(A + (r << 2), tw + (g << 2), A + (r << 2));

		uint64_t t[4]{
			A[(l << 2) + 0],
			A[(l << 2) + 1],
			A[(l << 2) + 2],
			A[(l << 2) + 3]
		};

		Add(A + (l << 2), A + (r << 2), A + (l << 2));
		Sub(t, A + (r << 2), A + (r << 2));
	}
}

/*
	Reorders array by bit-reversing the indices
*/
__global__ void bit_reverse(uint64_t *__restrict__ A, uint64_t *__restrict__ tmp, int s, size_t nthr) {
	int id = blockIdx.x * nthr + threadIdx.x;
	int n = 1 << s;
	int shifted = __brev(id) >> (32 - s);
	if (id < n && shifted < n) {
		A[shifted * 4 + 0] = tmp[id * 4 + 0];
		A[shifted * 4 + 1] = tmp[id * 4 + 1];
		A[shifted * 4 + 2] = tmp[id * 4 + 2];
		A[shifted * 4 + 3] = tmp[id * 4 + 3];
	}
}

/*
	FFT if we have enough threads
*/
__global__ void inplace_fft(uint64_t *__restrict__ A, uint64_t *__restrict__ tw, int j, int m, int n, size_t nthr) {
	int k = blockIdx.x * nthr + threadIdx.x;
	inplace_fft_inner(A, tw, j, k, m, n);
}

/*
	FFT if we don't have enough threads
*/
__global__ void inplace_fft_outer(uint64_t *__restrict__ A, uint64_t *__restrict__ tw, int m, int n, size_t nthr) {
	int j = (blockIdx.x * nthr + threadIdx.x) * m;
	for (int k = 0; k < m / 2; k++) {
		inplace_fft_inner(A, tw, j, k, m, n);
	}
}

// #define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)
// void checkLast(const char* const file, const int line)
// {
//     cudaError_t err{cudaGetLastError()};
//     if (err != cudaSuccess)
//     {
//         std::cerr << "CUDA Runtime Error at: " << file << ":" << line
//                   << std::endl;
//         std::cerr << cudaGetErrorString(err) << std::endl;
//         // We don't exit when we encounter CUDA errors in this example.
//         // std::exit(EXIT_FAILURE);
//     }
// }

extern "C" {
	void fft(uint64_t *a, uint64_t *omegas, size_t n, bool invert, int balance, size_t threads) {
		// NOTE: n is the number of ELEMENTS,
		// The total list size is n * S
		uint64_t list_size = n * S;
		
		// Create array from input vector
		uint64_t buffer_size = list_size * sizeof(uint64_t);

		// Allocate array
		uint64_t *data_array = (uint64_t *) malloc(buffer_size);
		for (int i = 0; i < list_size; i++) {
			data_array[i] = a[i];
		}

		// Copy data to FPGA using 2 arrays for bitreverse
		uint64_t *A, *tmp, *tw;
		cudaMalloc((void **) &A, buffer_size);
		cudaMalloc((void **) &tmp, buffer_size);
		cudaMalloc((void **) &tw, buffer_size); // TODO: I think this will/can be smaller than buffer_size 
		cudaMemcpy(tmp, data_array, buffer_size, cudaMemcpyHostToDevice);
		cudaMemcpy(tw, omegas, buffer_size / 2, cudaMemcpyHostToDevice);

		// Bit reverse ordering
		int s = log2(n);
		bit_reverse<<<ceil((float)n / threads), threads>>>(A, tmp, s, threads);
		
		// Synchronize 
		cudaDeviceSynchronize();
		// Iterative FFT with loop parallelism balancing
		for (int i = 1; i <= s; i++) {
			int m = 1 << i;
			if (n / m > balance) {
				inplace_fft_outer<<<ceil((float)n / m / threads), threads>>>(A, tw, m, n, threads);
			} else {
				for (int j = 0; j < n; j += m) {
					float repeats = m / 2;
					inplace_fft<<<ceil(repeats / threads), threads>>>(A, tw, j, m, n, threads);
				}
			}
		}

		// Copy back result
		uint64_t *result;
		result = (uint64_t *) malloc(buffer_size);
		cudaMemcpy(result, A, buffer_size, cudaMemcpyDeviceToHost);

		for (int i = 0; i < list_size; i++) {
			a[i] = result[i];
		}

		free(data_array);
		free(result);
		cudaFree(A);
		cudaFree(tmp);
		cudaFree(tw);

		return;
	}
}