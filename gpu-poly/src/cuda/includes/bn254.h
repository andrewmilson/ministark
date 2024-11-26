#ifndef __BN254_H__
#define __BN254_H__

// K -> Num bits to store a field element
#define K 256

// S -> Number of words to store a field element
#define S 4

#include <cuda.h>

// Constants for BN254
__device__ uint64_t MODULUS[S]{
	4891460686036598785,
	2896914383306846353,
	13281191951274694749,
	3486998266802970665
};
__device__ uint64_t R[S]{
	12436184717236109307,
	3962172157175319849,
	7381016538464732718,
	1011752739694698287
};
__device__ uint64_t UNITY[S]{
	7164790868263648668,
	11685701338293206998,
	6216421865291908056,
	1756667274303109607
};
__device__ uint64_t montConstant = 14042775128853446655;

#endif