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

#ifndef STR_UTIL_H_
#define STR_UTIL_H_

#include <stdint.h>
#include <string>

/**
 * Given a string, return an int hash for it.
 */
static inline int
hash_string(const std::string& s) {
	int ret = 0;
	int a = 63689;
	int b = 378551;
	for(size_t i = 0; i < s.length(); i++) {
		ret = (ret * a) + (int)s[i];
		if(a == 0) {
			a += b;
		} else {
			a *= b;
		}
		if(a == 0) {
			a += b;
		}
	}
	return ret;
}

/*
 * Implementation of the FNV-1a hash function.
 * https://tools.etf.org/html/draft/eastlake-fnv-03
*/
static inline uint32_t hash_str(const char *str) {
	int c;
	const uint32_t FNV_PRIME = 0x010000193;
	uint32_t hash = 0x811C9Dc5;


	while ((c = *str++) != '\0')
		hash = (hash ^ c) * FNV_PRIME;
	return hash;
}
#endif /* STR_UTIL_H_ */
