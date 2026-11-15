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

#ifndef MM_H_
#define MM_H_

/**
 * mm.h:
 *
 * Defines that make it easier to handle files in the two different MM
 * contexts: i.e. on Linux and Mac where MM is supported and POSIX I/O
 * functions work as expected, and on Windows where MM is not supported
 * and where there isn't POSIX I/O,
 */

#include <stddef.h>
#include <stdint.h>

#define MM_READ(file, dest, sz) fread(dest, 1, sz, file)
#define MM_IS_IO_ERR(file_hd, ret, count) is_fread_err(file_hd, ret, count)

/**
 * Alignment, in bytes, of the large index arrays read from disk into the
 * heap (BWT, SA sample, ftab, eftab, bit-packed reference).  A BWT side
 * is 64 bytes (128 for 64-bit offsets) and is read at a random place, so
 * a side that starts on a cache line touches one line (two) instead of
 * two (three).  glibc serves large new[] blocks with mmap and returns an
 * address 16 bytes into a page, which is never a line boundary.  Only
 * where the arrays start in memory changes; their bytes are the file's.
 */
static const size_t MM_INDEX_ALIGN = 64;

/**
 * Round p up to the next multiple of MM_INDEX_ALIGN bytes.  The block p
 * points into must extend MM_INDEX_ALIGN bytes past what is used.
 */
template<typename T>
static inline T* mmAlignUp(T* p) {
	const uintptr_t a = MM_INDEX_ALIGN;
	return (T*)(((uintptr_t)p + a - 1) & ~(a - 1));
}

#endif /* MM_H_ */
