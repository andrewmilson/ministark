#ifndef __FR_ARITH_H__
#define __FR_ARITH_H__

#define LO 1
#define HI 0
#define SUM 0
#define CARRY 1
#define DIFF 0
#define BORROW 1

#include <stdint.h>
#include <stdbool.h>
#include <cuda.h>

#include "bn254.h"

/*=============== DEVICE FUNCTIONS =================*/
/* Modulus Arithmetic */
static __device__ inline void Add64(uint64_t a, uint64_t b, uint64_t carry, uint64_t *result) {
	uint64_t sum = a + b + carry;
	carry = ((a & b) | ((a | b) & ~sum)) >> 63;

	result[SUM] = sum;
	result[CARRY] = carry;
}

static __device__ inline void Sub64(uint64_t a, uint64_t b, uint64_t borrow, uint64_t *result) {
	uint64_t diff = a - b - borrow;
	borrow = ((~a & b) | (~(a ^ b) & diff)) >> 63;

	result[DIFF] = diff;
	result[BORROW] = borrow;
}

static __device__ inline void Mul64(uint64_t a, uint64_t b, uint64_t *result) {
	const uint64_t one = 1;
	const uint64_t mask32 = (one << 32) - 1;
	uint64_t a0 = a & mask32;
	uint64_t a1 = a >> 32;
	uint64_t b0 = b & mask32;
	uint64_t b1 = b >> 32;
	
	uint64_t w0 = a0 * b0;
	uint64_t t = a1 * b0 + (w0 >> 32);
	uint64_t w1 = t & mask32;
	uint64_t w2 = t >> 32;
	w1 += a0 * b1;

	result[HI] = a1 * b1 + w2 + (w1 >> 32);
	result[LO] = a * b;
}

static __device__ inline void madd0(uint64_t a, uint64_t b, uint64_t c, uint64_t &hi) {
	uint64_t mul_result[2];
	Mul64(a, b, mul_result);

	uint64_t add_result[2];
	Add64(mul_result[LO], c, 0, add_result);

	uint64_t add_result_2[2];
	Add64(mul_result[HI], 0, add_result[CARRY], add_result_2);

	hi = add_result_2[SUM];
}

static __device__ inline void madd1(uint64_t a, uint64_t b, uint64_t c, uint64_t &hi, uint64_t &lo) {
	uint64_t mul_result[2];
	Mul64(a, b, mul_result);

	uint64_t add_result[2];
	Add64(mul_result[LO], c, 0, add_result);

	uint64_t add_result_2[2];
	Add64(mul_result[HI], 0, add_result[CARRY], add_result_2);

	hi = add_result_2[SUM];
	lo = add_result[SUM];
}

static __device__ inline void madd2(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t &hi, uint64_t &lo) {
	uint64_t mul_result[2];
	Mul64(a, b, mul_result);

	uint64_t add_result[2];
	Add64(c, d, 0, add_result);

	uint64_t add_result_2[2];
	Add64(mul_result[HI], 0, add_result[CARRY], add_result_2);
	
	uint64_t add_result_3[2];
	Add64(mul_result[LO], add_result[SUM], 0, add_result_3);

	uint64_t add_result_4[2];
	Add64(add_result_2[SUM], 0, add_result_3[CARRY], add_result_4);

	hi = add_result_4[SUM];
	lo = add_result_3[SUM];
}

static __device__ inline void madd3(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e, uint64_t &hi, uint64_t &lo) {
	uint64_t mul_result[2];
	Mul64(a, b, mul_result);

	uint64_t add_result[2];
	Add64(c, d, 0, add_result);

	uint64_t add_result_2[2];
	Add64(mul_result[HI], 0, add_result[CARRY], add_result_2);

	uint64_t add_result_3[2];
	Add64(mul_result[LO], add_result[SUM], 0, add_result_3);

	uint64_t add_result_4[2];
	Add64(add_result_2[SUM], e, add_result_3[CARRY], add_result_4);

	hi = add_result_4[SUM];
	lo = add_result_3[SUM];
}

static __device__ inline bool smallerThanModulus(uint64_t* z) {
	return (
		(z[3] < MODULUS[3] || 
		(z[3] == MODULUS[3] && 
		(z[2] < MODULUS[2] || 
		(z[2] == MODULUS[2] && 
		(z[1] < MODULUS[1] || 
		(z[1] == MODULUS[1] && 
		(z[0] < MODULUS[0])))))))
	);
}

static __device__ inline void Add(uint64_t *a, uint64_t* b, uint64_t* out) {
	uint64_t carry = 0;
	uint64_t add_container[2];

	Add64(a[0], b[0], carry, add_container);
	out[0] = add_container[SUM];
	carry = add_container[CARRY];

	Add64(a[1], b[1], carry, add_container);
	out[1] = add_container[SUM];
	carry = add_container[CARRY];

	Add64(a[2], b[2], carry, add_container);
	out[2] = add_container[SUM];
	carry = add_container[CARRY];

	Add64(a[3], b[3], carry, add_container);
	out[3] = add_container[SUM];

	// Reduce if necessary
	if (!smallerThanModulus(out)) {
		uint64_t borrow = 0;
		uint64_t sub_container[2];

		Sub64(out[0], MODULUS[0], borrow, sub_container);
		out[0] = sub_container[DIFF];
		borrow = sub_container[BORROW];

		Sub64(out[1], MODULUS[1], borrow, sub_container);
		out[1] = sub_container[DIFF];
		borrow = sub_container[BORROW];

		Sub64(out[2], MODULUS[2], borrow, sub_container);
		out[2] = sub_container[DIFF];
		borrow = sub_container[BORROW];

		Sub64(out[3], MODULUS[3], borrow, sub_container);
		out[3] = sub_container[DIFF];
	}
}

static __device__ inline void Sub(uint64_t* a, uint64_t* b, uint64_t* out) {
	uint64_t borrow = 0;
	uint64_t sub_container[2];

	Sub64(a[0], b[0], borrow, sub_container);
	out[0] = sub_container[DIFF];
	borrow = sub_container[BORROW];

	Sub64(a[1], b[1], borrow, sub_container);
	out[1] = sub_container[DIFF];
	borrow = sub_container[BORROW];

	Sub64(a[2], b[2], borrow, sub_container);
	out[2] = sub_container[DIFF];
	borrow = sub_container[BORROW];

	Sub64(a[3], b[3], borrow, sub_container);
	out[3] = sub_container[DIFF];
	borrow = sub_container[BORROW];

	if (borrow != 0) {
		uint64_t carry = 0;
		uint64_t add_container[2];

		Add64(out[0], MODULUS[0], carry, add_container);
		out[0] = add_container[SUM];
		carry = add_container[CARRY];

		Add64(out[1], MODULUS[1], carry, add_container);
		out[1] = add_container[SUM];
		carry = add_container[CARRY];

		Add64(out[2], MODULUS[2], carry, add_container);
		out[2] = add_container[SUM];
		carry = add_container[CARRY];

		Add64(out[3], MODULUS[3], carry, add_container);
		out[3] = add_container[SUM];
	}
}

static __device__ inline void Mul(uint64_t* a, uint64_t* b, uint64_t* out) {
	// Implements CIOS multiplication -- section 2.3.2 of Tolga Acar's thesis
	// https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf

	uint64_t t[4]{0, 0, 0, 0};
	uint64_t c[3]{0, 0, 0};

	{
		// round 0
		uint64_t v = a[0];
		uint64_t mul_container[2];
		Mul64(v, b[0], mul_container);
		c[1] = mul_container[HI];
		c[0] = mul_container[LO];

		uint64_t m = c[0] * montConstant;
		madd0(m, MODULUS[0], c[0], c[2]);
		madd1(v, b[1], c[1], c[1], c[0]);
		madd2(m, MODULUS[1], c[2], c[0], c[2], t[0]);
		madd1(v, b[2], c[1], c[1], c[0]);
		madd2(m, MODULUS[2], c[2], c[0], c[2], t[1]);
		madd1(v, b[3], c[1], c[1], c[0]);
		madd3(m, MODULUS[3], c[0], c[2], c[1], t[3], t[2]);
	}
	{
		// round 1
		uint64_t v = a[1];
		madd1(v, b[0], t[0], c[1], c[0]);
		uint64_t m = c[0] * montConstant;
		madd0(m, MODULUS[0], c[0], c[2]);
		madd2(v, b[1], c[1], t[1], c[1], c[0]);
		madd2(m, MODULUS[1], c[2], c[0], c[2], t[0]);
		madd2(v, b[2], c[1], t[2], c[1], c[0]);
		madd2(m, MODULUS[2], c[2], c[0], c[2], t[1]);
		madd2(v, b[3], c[1], t[3], c[1], c[0]);
		madd3(m, MODULUS[3], c[0], c[2], c[1], t[3], t[2]);
	}
	{
		// round 2
		uint64_t v = a[2];
		madd1(v, b[0], t[0], c[1], c[0]);
		uint64_t m = c[0] * montConstant;
		madd0(m, MODULUS[0], c[0], c[2]);
		madd2(v, b[1], c[1], t[1], c[1], c[0]);
		madd2(m, MODULUS[1], c[2], c[0], c[2], t[0]);
		madd2(v, b[2], c[1], t[2], c[1], c[0]);
		madd2(m, MODULUS[2], c[2], c[0], c[2], t[1]);
		madd2(v, b[3], c[1], t[3], c[1], c[0]);
		madd3(m, MODULUS[3], c[0], c[2], c[1], t[3], t[2]);
	}
	{
		// round 3
		uint64_t v = a[3];
		madd1(v, b[0], t[0], c[1], c[0]);
		uint64_t m = c[0] * montConstant;
		madd0(m, MODULUS[0], c[0], c[2]);
		madd2(v, b[1], c[1], t[1], c[1], c[0]);
		madd2(m, MODULUS[1], c[2], c[0], c[2], out[0]);
		madd2(v, b[2], c[1], t[2], c[1], c[0]);
		madd2(m, MODULUS[2], c[2], c[0], c[2], out[1]);
		madd2(v, b[3], c[1], t[3], c[1], c[0]);
		madd3(m, MODULUS[3], c[0], c[2], c[1], out[3], out[2]);
	}

	//if z >= q → z -= q
	if (!smallerThanModulus(out)) {
		uint64_t borrow = 0;
		uint64_t sub_container[2];

		Sub64(out[0], MODULUS[0], borrow, sub_container);
		out[0] = sub_container[DIFF];
		borrow = sub_container[BORROW];

		Sub64(out[1], MODULUS[1], borrow, sub_container);
		out[1] = sub_container[DIFF];
		borrow = sub_container[BORROW];

		Sub64(out[2], MODULUS[2], borrow, sub_container);
		out[2] = sub_container[DIFF];
		borrow = sub_container[BORROW];

		Sub64(out[3], MODULUS[3], borrow, sub_container);
		out[3] = sub_container[DIFF];
	}
}
#endif