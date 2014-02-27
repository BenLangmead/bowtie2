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

#ifndef WORD_IO_H_
#define WORD_IO_H_

#include <stdint.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "btypes.h"

/**
 * Write a 32/64 bit unsigned to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
template <typename T>
static inline void writeU(std::ostream& out, T x, bool toBigEndian) {
	T y = endianizeU<T>(x, toBigEndian);
	out.write((const char*)&y, sizeof(T));
}

/**
 * Write a 32/64 bit unsigned to an output stream using the native
 * endianness.
 */
template <typename T>
static inline void writeU(std::ostream& out, T x) {
	out.write((const char*)&x, sizeof(T));
}

/**
 * Write a 32/64 bit signed int to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
template <typename T>
static inline void writeI(std::ostream& out, T x, bool toBigEndian) {
	T y = endianizeI<T>(x, toBigEndian);
	out.write((const char*)&y, sizeof(T));
}

/**
 * Write a 32/64 bit unsigned to an output stream using the native
 * endianness.
 */
template <typename T>
static inline void writeI(std::ostream& out, T x) {
	out.write((const char*)&x, sizeof(T));
}

/**
 * Read a 32/64 bit unsigned from an input stream, inverting endianness
 * if necessary.
 */
//template <typename T>
//static inline T readU(std::istream& in, bool swap) {
//	T x;
//	in.read((char *)&x, OFF_SIZE);
//	assert_eq(OFF_SIZE, in.gcount());
//	if(swap) {
//		return endianSwapU(x);
//	} else {
//		return x;
//	}
//}
template <typename T>
static inline T readU(std::istream& in, bool swap) {
	T x;
	in.read((char *)&x, sizeof(T));
	assert_eq(sizeof(T), in.gcount());
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapU32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapU64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}


/**
 * Read a 32/64 bit unsigned from a FILE*, optionally inverting
 * endianness.
 */
//template <typename T>
//static inline T readU(FILE* in, bool swap) {
//	T x;
//	if(fread((void *)&x, 1, OFF_SIZE, in) != OFF_SIZE) {
//		assert(false);
//	}
//	if(swap) {
//		return endianSwapU(x);
//	} else {
//		return x;
//	}
//}
template <typename T>
static inline T readU(FILE* in, bool swap) {
	T x;
	if(fread((void *)&x, 1, sizeof(T), in) != sizeof(T)) {
		assert(false);
	}
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapU32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapU64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}


/**
 * Read a 32/64 bit signed from an input stream, inverting endianness
 * if necessary.
 */
//template <typename T>
//static inline T readI(std::istream& in, bool swap) {
//	T x;
//	in.read((char *)&x, OFF_SIZE);
//	assert_eq(OFF_SIZE, in.gcount());
//	if(swap) {
//		return endianSwapI(x);
//	} else {
//		return x;
//	}
//}
template <typename T>
static inline T readI(std::istream& in, bool swap) {
	T x;
	in.read((char *)&x, sizeof(T));
	assert_eq(sizeof(T), in.gcount());
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapI32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapI64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}


/**
 * Read a 32/64 bit unsigned from a FILE*, optionally inverting
 * endianness.
 */
//template <typename T>
//static inline T readI(FILE* in, bool swap) {
//	T x;
//	if(fread((void *)&x, 1, OFF_SIZE, in) != OFF_SIZE) {
//		assert(false);
//	}
//	if(swap) {
//		return endianSwapI(x);
//	} else {
//		return x;
//	}
//}
template <typename T>
static inline T readI(FILE* in, bool swap) {
	T x;
	if(fread((void *)&x, 1, sizeof(T), in) != sizeof(T)) {
		assert(false);
	}
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapI32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapI64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}


#endif /*WORD_IO_H_*/
