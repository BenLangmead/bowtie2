/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
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

#ifndef UTIL_H_
#define UTIL_H_

#include <stdlib.h>
#include <limits>

/**
 * C++ version char* style "itoa":
 */
template<typename T>
char* itoa10(const T& value, char* result) {
	// Check that base is valid
	char* out = result;
	T quotient = value;
	if(std::numeric_limits<T>::is_signed) {
		if(quotient <= 0) quotient = -quotient;
	}
	// Now write each digit from most to least significant
	do {
		*out = "0123456789"[quotient % 10];
		++out;
		quotient /= 10;
	} while (quotient > 0);
	// Only apply negative sign for base 10
	if(std::numeric_limits<T>::is_signed) {
		// Avoid compiler warning in cases where T is unsigned
		if (value <= 0 && value != 0) *out++ = '-';
	}
	reverse( result, out );
	*out = 0; // terminator
	return out;
}

#endif /*ndef UTIL_H_*/
