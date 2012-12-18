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

#ifndef ENDIAN_SWAP_H
#define ENDIAN_SWAP_H

#include <stdint.h>
#include <inttypes.h>

/**
 * Return true iff the machine running this program is big-endian.
 */
static inline bool currentlyBigEndian() {
	static uint8_t endianCheck[] = {1, 0, 0, 0};
	return *((uint32_t*)endianCheck) != 1;
}

/**
 * Return copy of uint32_t argument with byte order reversed.
 */
static inline uint32_t endianSwapU32(uint32_t u) {
	uint32_t tmp = 0;
	tmp |= ((u >> 24) & (0xff <<  0));
	tmp |= ((u >>  8) & (0xff <<  8));
	tmp |= ((u <<  8) & (0xff << 16));
	tmp |= ((u << 24) & (0xff << 24));
	return tmp;
}

/**
 * Return copy of uint64_t argument with byte order reversed.
 */
static inline uint64_t endianSwapU64(uint64_t u) {
	uint64_t tmp = 0;
	tmp |= ((u >> 56) & (0xffull <<  0));
	tmp |= ((u >> 40) & (0xffull <<  8));
	tmp |= ((u >> 24) & (0xffull << 16));
	tmp |= ((u >>  8) & (0xffull << 24));
	tmp |= ((u <<  8) & (0xffull << 32));
	tmp |= ((u << 24) & (0xffull << 40));
	tmp |= ((u << 40) & (0xffull << 48));
	tmp |= ((u << 56) & (0xffull << 56));
	return tmp;
}


/**
 * Return copy of int32_t argument with byte order reversed.
 */
static inline int32_t endianSwapI32(int32_t i) {
	int32_t tmp = 0;
	tmp |= ((i >> 24) & (0xff <<  0));
	tmp |= ((i >>  8) & (0xff <<  8));
	tmp |= ((i <<  8) & (0xff << 16));
	tmp |= ((i << 24) & (0xff << 24));
	return tmp;
}

/**
 * Convert uint32_t argument to the specified endianness.  It's assumed
 * that u currently has the endianness of the current machine.
 */
static inline uint32_t endianizeU32(uint32_t u, bool toBig) {
	if(toBig == currentlyBigEndian()) {
		return u;
	}
	return endianSwapU32(u);
}

/**
 * Convert int32_t argument to the specified endianness.  It's assumed
 * that u currently has the endianness of the current machine.
 */
static inline int32_t endianizeI32(int32_t i, bool toBig) {
	if(toBig == currentlyBigEndian()) {
		return i;
	}
	return endianSwapI32(i);
}

#endif
