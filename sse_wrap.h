/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * sse_wrap.h
 *
 * Routines to wrap Streaming SIMD Extensions (SSE) emmintrin.h
 * for an Intel x86 CPU and SIMD Everywhere (simde) for other CPUs.
 */

#ifndef SSE_WRAP_H_
#define SSE_WRAP_H_

#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
#include "simde/x86/sse2.h"
#else
#include <emmintrin.h>
#endif

#define NBYTES_PER_REG 16

#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
typedef simde__m128i SSERegI;
#define sse_adds_epi16(x, y) simde_mm_adds_epi16(x, y)
#define sse_adds_epu8(x, y) simde_mm_adds_epu8(x, y)
#define sse_cmpeq_epi16(x, y) simde_mm_cmpeq_epi16(x, y)
#define sse_cmpeq_epi8(x, y) simde_mm_cmpeq_epi8(x, y)
#define sse_cmpgt_epi16(x, y) simde_mm_cmpgt_epi16(x, y)
#define sse_cmpgt_epi8(x, y) simde_mm_cmpgt_epi8(x, y)
#define sse_cmplt_epi16(x, y) simde_mm_cmplt_epi16(x, y)
#define sse_cmplt_epu8(x, y) simde_mm_cmplt_epu8(x, y)
#define sse_extract_epi16(x, y) simde_mm_extract_epi16(x, y)
#define sse_insert_epi16(x, y, z) simde_mm_insert_epi16(x, y, z)
#define sse_load_siall(x) simde_mm_load_si128(x)
#define sse_max_epi16(x, y) simde_mm_max_epi16(x, y)
#define sse_max_epu8(x, y) simde_mm_max_epu8(x, y)
#define sse_movemask_epi8(x) simde_mm_movemask_epi8(x)
#define sse_or_siall(x, y) simde_mm_or_si128(x, y)
#define sse_setzero_siall() simde_mm_setzero_si128()
#define sse_shuffle_epi32(x, y) simde_mm_shuffle_epi32(x, y)
#define sse_shufflelo_epi16(x, y) simde_mm_shufflelo_epi16(x, y)
#define sse_slli_epi16(x, y) simde_mm_slli_epi16(x, y)
#define sse_slli_siall(x, y) simde_mm_slli_si128(x, y)
#define sse_srli_epi16(x, y) simde_mm_srli_epi16(x, y)
#define sse_srli_epu8(x, y) simde_mm_srli_epu8(x, y)
#define sse_srli_siall(x, y) simde_mm_srli_si128(x, y)
#define sse_store_siall(x, y) simde_mm_store_si128(x, y)
#define sse_subs_epi16(x, y) simde_mm_subs_epi16(x, y)
#define sse_subs_epu8(x, y) simde_mm_subs_epu8(x, y)
#define sse_xor_siall(x, y) simde_mm_xor_si128(x, y)
#define sse_set1_epi16(x) simde_mm_set1_epi16(x)

#else
typedef __m128i SSERegI;
#define sse_adds_epi16(x, y) _mm_adds_epi16(x, y)
#define sse_adds_epu8(x, y) _mm_adds_epu8(x, y)
#define sse_cmpeq_epi16(x, y) _mm_cmpeq_epi16(x, y)
#define sse_cmpeq_epi8(x, y) _mm_cmpeq_epi8(x, y)
#define sse_cmpgt_epi16(x, y) _mm_cmpgt_epi16(x, y)
#define sse_cmpgt_epi8(x, y) _mm_cmpgt_epi8(x, y)
#define sse_cmplt_epi16(x, y) _mm_cmplt_epi16(x, y)
#define sse_cmplt_epu8(x, y) _mm_cmplt_epu8(x, y)
#define sse_extract_epi16(x, y) _mm_extract_epi16(x, y)
#define sse_insert_epi16(x, y, z) _mm_insert_epi16(x, y, z)
#define sse_load_siall(x) _mm_load_si128(x)
#define sse_max_epi16(x, y) _mm_max_epi16(x, y)
#define sse_max_epu8(x, y) _mm_max_epu8(x, y)
#define sse_movemask_epi8(x) _mm_movemask_epi8(x)
#define sse_or_siall(x, y) _mm_or_si128(x, y)
#define sse_setzero_siall() _mm_setzero_si128()
#define sse_shuffle_epi32(x, y) _mm_shuffle_epi32(x, y)
#define sse_shufflelo_epi16(x, y) _mm_shufflelo_epi16(x, y)
#define sse_slli_epi16(x, y) _mm_slli_epi16(x, y)
#define sse_slli_siall(x, y) _mm_slli_si128(x, y)
#define sse_srli_epi16(x, y) _mm_srli_epi16(x, y)
#define sse_srli_epu8(x, y) _mm_srli_epu8(x, y)
#define sse_srli_siall(x, y) _mm_srli_si128(x, y)
#define sse_store_siall(x, y) _mm_store_si128(x, y)
#define sse_subs_epi16(x, y) _mm_subs_epi16(x, y)
#define sse_subs_epu8(x, y) _mm_subs_epu8(x, y)
#define sse_xor_siall(x, y) _mm_xor_si128(x, y)
#define sse_set1_epi16(x) _mm_set1_epi16(x)

#endif

/* Fill all elements in outval with inval */
/* opt version will check for special ivals that can use shortcuts */
#define sse_fill_i16(inval, outval) outval=sse_set1_epi16(inval)

#define sse_fill_i16_opt(inval, outval) { \
	if (inval==0xffff) sse_cmpeq_epi16(outval, outval); \
	else if (inval==0) sse_xor_siall(outval, outval); \
	else sse_fill_i16(inval, outval); \
}

#define sse_fill_u8(inval, outval) {\
	int invalloc = inval; \
	int dup = (invalloc << 8) | (invalloc & 0x00ff); \
	sse_fill_i16(dup, outval); \
}

#define sse_fill_u8_opt(inval, outval) {\
	if (inval==0xff) sse_cmpeq_epi16(outval, outval); \
	else if (inval==0) sse_xor_siall(outval, outval); \
	else sse_fill_u8(inval, outval); \
}

/* Set the low element with invl, all others to 0 */
#define sse_set_low_i16(inval, outval) { \
	outval = sse_setzero_siall(); \
	outval = sse_insert_epi16(outval, inval, 0); \
}

#define sse_set_low_u8(inval, outval) { \
	outval = sse_setzero_siall(); \
	outval = sse_insert_epi16(outval, inval, 0); \
}

/* compute the max val of a vector */
#define sse_max_score_i16(inval, outval) { \
		SSERegI vlmax = inval; \
		SSERegI vltmp = sse_srli_siall(vlmax, 8); \
		vlmax = sse_max_epu8(vlmax, vltmp); \
		vltmp = sse_srli_siall(vlmax, 4); \
		vlmax = sse_max_epu8(vlmax, vltmp); \
		vltmp = sse_srli_siall(vlmax, 2); \
		vlmax = sse_max_epu8(vlmax, vltmp); \
		vltmp = sse_srli_siall(vlmax, 1); \
		vlmax = sse_max_epu8(vlmax, vltmp); \
		outval = sse_extract_epi16(vlmax, 0); \
}

#endif /* SSE_WRAP_H_ */
