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

#ifdef __aarch64__
#include "simde/x86/sse2.h"
#else
#include <emmintrin.h>
#endif

#ifdef __aarch64__
typedef  __m128i simdesimde__m128i;
#define simde_mm_adds_epi16(x, y) _mm_adds_epi16(x, y)
#define simde_mm_adds_epu8(x, y) _mm_adds_epu8(x, y)
#define simde_mm_cmpeq_epi16(x, y) _mm_cmpeq_epi16(x, y)
#define simde_mm_cmpeq_epi8(x, y) _mm_cmpeq_epi8(x, y)
#define simde_mm_cmpgt_epi16(x, y) _mm_cmpgt_epi16(x, y)
#define simde_mm_cmpgt_epi8(x, y) _mm_cmpgt_epi8(x, y)
#define simde_mm_cmplt_epi16(x, y) _mm_cmplt_epi16(x, y)
#define simde_mm_cmplt_epu8(x, y) _mm_cmplt_epu8(x, y)
#define simde_mm_extract_epi16(x, y) _mm_extract_epi16(x, y)
#define simde_mm_insert_epi16(x, y, z) _mm_insert_epi16(x, y, z)
#define simde_mm_load_si128(x) _mm_load_si128(x)
#define simde_mm_max_epi16(x, y) _mm_max_epi16(x, y)
#define simde_mm_max_epu8(x, y) _mm_max_epu8(x, y)
#define simde_mm_movemask_epi8(x) _mm_movemask_epi8(x)
#define simde_mm_or_si128(x, y) _mm_or_si128(x, y)
#define simde_mm_setzero_si128() _mm_setzero_si128()
#define simde_mm_shuffle_epi32(x, y) _mm_shuffle_epi32(x, y)
#define simde_mm_shufflelo_epi16(x, y) _mm_shufflelo_epi16(x, y)
#define simde_mm_slli_epi16(x, y) _mm_slli_epi16(x, y)
#define simde_mm_slli_si128(x, y) _mm_slli_si128(x, y)
#define simde_mm_srli_epi16(x, y) _mm_srli_epi16(x, y)
#define simde_mm_srli_epu8(x, y) _mm_srli_epu8(x, y)
#define simde_mm_srli_si128(x, y) _mm_srli_si128(x, y)
#define simde_mm_store_si128(x, y) _mm_store_si128(x, y)
#define simde_mm_subs_epi16(x, y) _mm_subs_epi16(x, y)
#define simde_mm_subs_epu8(x, y) _mm_subs_epu8(x, y)
#define simde_mm_xor_si128(x, y) _mm_xor_si128(x, y)

#endif

#endif /* SSE_WRAP_H_ */
