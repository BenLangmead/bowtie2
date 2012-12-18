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

#ifndef BITPACK_H_
#define BITPACK_H_

#include <stdint.h>
#include "assert_helpers.h"

/**
 * Routines for marshalling 2-bit values into and out of 8-bit or
 * 32-bit hosts
 */

static inline void pack_2b_in_8b(const int two, uint8_t& eight, const int off) {
	assert_lt(two, 4);
	assert_lt(off, 4);
	eight |= (two << (off*2));
}

static inline int unpack_2b_from_8b(const uint8_t eight, const int off) {
	assert_lt(off, 4);
	return ((eight >> (off*2)) & 0x3);
}

static inline void pack_2b_in_32b(const int two, uint32_t& thirty2, const int off) {
	assert_lt(two, 4);
	assert_lt(off, 16);
	thirty2 |= (two << (off*2));
}

static inline int unpack_2b_from_32b(const uint32_t thirty2, const int off) {
	assert_lt(off, 16);
	return ((thirty2 >> (off*2)) & 0x3);
}

#endif /*BITPACK_H_*/
