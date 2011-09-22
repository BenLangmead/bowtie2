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

#ifndef BITSET_H_
#define BITSET_H_

#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string.h>
#include <stdexcept>
#include "assert_helpers.h"
#include "threading.h"

/**
 * Given a words array and a size, allocate a new, larger array, moving
 * data from the old to the new array, and set all newly-allocated
 * words to 0.  Return the new, larger array, which can be substituted
 * for the old one.  The new array is larger than the old by about 50%.
 */
static inline uint32_t*
bitsetRealloc(uint32_t& sz, uint32_t* words, const char *errmsg = NULL) {
	uint32_t oldsz = sz;
	if(sz > 0) {
		sz += (sz >> 1) + 31; // Add 50% more elements, plus a bit
		sz &= ~31;            // Make sure it's 32-aligned
	} else {
		sz = 1024; // Start off at 1024 bits to avoid many expansions
	}
	assert_gt(sz, oldsz);
	assert_eq(0, (sz & 31));
	uint32_t *newwords;
	try {
		newwords = new uint32_t[sz >> 5 /* convert to words */];
	} catch(std::bad_alloc& ba) {
		if(errmsg != NULL) {
			// Output given error message
			std::cerr << errmsg;
		}
		throw 1;
	}
	if(oldsz > 0) {
		// Move old values into new array
		memcpy(newwords, words, oldsz >> 3 /* convert to bytes */);
	}
	// Initialize all new words to 0
	memset(newwords + (oldsz >> 5 /*convert to words*/), 0,
	       (sz - oldsz) >> 3 /* convert to bytes */);
	return newwords; // return new array
}

/**
 * A simple unsynchronized bitset class.
 */
class Bitset {

public:
	Bitset(uint32_t sz, const char *errmsg = NULL) : _errmsg(errmsg) {
		uint32_t nwords = (sz >> 5)+1;
		try {
			_words = new uint32_t[nwords];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			throw 1;
		}
		assert(_words != NULL);
		memset(_words, 0, nwords * 4);
		_sz = nwords << 5;
		_cnt = 0;
	}

	Bitset(const Bitset& o) : _words(NULL) {
		this->operator=(o);
	}

	~Bitset() {
		delete[] _words;
	}

	/**
	 * Test whether the given bit is set.
	 */
	bool test(uint32_t i) const {
		bool ret = false;
		if(i < _sz) {
			ret = ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		}
		return ret;
	}

	/**
	 * Set a bit in the vector that hasn't been set before.  Assert if
	 * it has been set.
	 */
	void set(uint32_t i) {
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(uint32_t oldsz = _sz);
			expand();
			assert_gt(_sz, oldsz);
		}
		// Fast path
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
		_cnt++;
		_words[i >> 5] |= (1 << (i & 0x1f));
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	/**
	 * Set a bit in the vector that might have already been set.
	 */
	void setOver(uint32_t i) {
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(uint32_t oldsz = _sz);
			expand();
			assert_gt(_sz, oldsz);
		}
		// Fast path
		if(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0) _cnt++;
		_words[i >> 5] |= (1 << (i & 0x1f));
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	/**
	 * Unset all entries.  Don't adjust size.
	 */
	void clear() {
		for(size_t i = 0; i < ((_sz+31)>>5); i++) {
			_words[i] = 0;
		}
		_cnt = 0;
	}

	/**
	 * Return the number of set bits.
	 */
	uint32_t count() const {
		return _cnt;
	}

	/**
	 * Return true iff no bits are set.
	 */
	bool empty() const {
		return _cnt == 0;
	}

	/**
	 * Deep copy from given Bitset to this one.
	 */
	Bitset& operator=(const Bitset& o) {
		_errmsg = o._errmsg;
		_sz = o._sz;
		_cnt = o._cnt;
		if(_words != NULL) delete[] _words;
		_words = new uint32_t[(_sz+31)>>5];
		for(size_t i = 0; i < (_sz+31)>>5; i++) {
			_words[i] = o._words[i];
		}
		return *this;
	}

private:

	/**
	 * Expand the size of the _words array by 50% to accommodate more
	 * bits.
	 */
	void expand() {
		uint32_t *newwords = bitsetRealloc(_sz, _words, _errmsg);
		delete[] _words;   // delete old array
		_words = newwords; // install new array
	}

	uint32_t _cnt;       // number of set bits
	const char *_errmsg; // error message if an allocation fails
	uint32_t _sz;        // size as # of bits
	uint32_t *_words;    // storage
};

/**
 * A simple fixed-length unsynchronized bitset class.
 */
template<int LEN>
class FixedBitset {

public:
	FixedBitset() : _cnt(0), _size(0) {
		memset(_words, 0, ((LEN>>5)+1) * 4);
	}

	/**
	 * Unset all bits.
	 */
	void clear() {
		memset(_words, 0, ((LEN>>5)+1) * 4);
	}

	/**
	 * Return true iff the bit at offset i has been set.
	 */
	bool test(uint32_t i) const {
		bool ret = false;
		assert_lt(i, LEN);
		ret = ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		return ret;
	}

	/**
	 * Set the bit at offset i.  Assert if the bit was already set.
	 */
	void set(uint32_t i) {
		// Fast path
		assert_lt(i, LEN);
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
		_words[i >> 5] |= (1 << (i & 0x1f));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	/**
	 * Set the bit at offset i.  Do not assert if the bit was already
	 * set.
	 */
	void setOver(uint32_t i) {
		// Fast path
		assert_lt(i, LEN);
		_words[i >> 5] |= (1 << (i & 0x1f));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	uint32_t count() const { return _cnt; }
	uint32_t size() const  { return _size; }

	/**
	 * Return true iff this FixedBitset has the same bits set as
	 * FixedBitset 'that'.
	 */
	bool operator== (const FixedBitset<LEN>& that) const {
		for(uint32_t i = 0; i < (LEN>>5)+1; i++) {
			if(_words[i] != that._words[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true iff this FixedBitset does not have the same bits set
	 * as FixedBitset 'that'.
	 */
	bool operator!= (const FixedBitset<LEN>& that) const {
		for(uint32_t i = 0; i < (LEN>>5)+1; i++) {
			if(_words[i] != that._words[i]) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Return a string-ized version of this FixedBitset.
	 */
	std::string str() const {
		std::ostringstream oss;
		for(int i = (int)size()-1; i >= 0; i--) {
			oss << (test(i)? "1" : "0");
		}
		return oss.str();
	}

private:
	uint32_t _cnt;
	uint32_t _size;
	uint32_t _words[(LEN>>5)+1]; // storage
};

/**
 * A simple fixed-length unsynchronized bitset class.
 */
class FixedBitset2 {

public:
	FixedBitset2(uint32_t len) : len_(len), _cnt(0), _size(0) {
		_words = new uint32_t[((len_ >> 5)+1)];
		memset(_words, 0, ((len_ >> 5)+1) * 4);
	}

	~FixedBitset2() { delete[] _words; }

	/**
	 * Unset all bits.
	 */
	void clear() {
		memset(_words, 0, ((len_ >> 5)+1) * 4);
		_cnt = 0;
		_size = 0;
	}

	/**
	 * Return true iff the bit at offset i has been set.
	 */
	bool test(uint32_t i) const {
		bool ret = false;
		assert_lt(i, len_);
		ret = ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		return ret;
	}

	/**
	 * Set the bit at offset i.  Assert if the bit was already set.
	 */
	void set(uint32_t i) {
		// Fast path
		assert_lt(i, len_);
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
		_words[i >> 5] |= (1 << (i & 0x1f));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	/**
	 * Clear the bit at offset i.  Assert if the bit was not already set.
	 */
	void clear(uint32_t i) {
		// Fast path
		assert_lt(i, len_);
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
		_words[i >> 5] &= ~(1 << (i & 0x1f));
		_cnt--;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
	}

	/**
	 * Set the bit at offset i.  Do not assert if the bit was already
	 * set.
	 */
	void setOver(uint32_t i) {
		// Fast path
		assert_lt(i, len_);
		if(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0) {
			_words[i >> 5] |= (1 << (i & 0x1f));
			_cnt++;
		}
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	uint32_t count() const { return _cnt; }
	uint32_t size() const  { return _size; }

	/**
	 * Return true iff this FixedBitset has the same bits set as
	 * FixedBitset 'that'.
	 */
	bool operator== (const FixedBitset2& that) const {
		for(uint32_t i = 0; i < (len_>>5)+1; i++) {
			if(_words[i] != that._words[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true iff this FixedBitset does not have the same bits set
	 * as FixedBitset 'that'.
	 */
	bool operator!= (const FixedBitset2& that) const {
		for(uint32_t i = 0; i < (len_>>5)+1; i++) {
			if(_words[i] != that._words[i]) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Return a string-ized version of this FixedBitset.
	 */
	std::string str() const {
		std::ostringstream oss;
		for(int i = (int)size()-1; i >= 0; i--) {
			oss << (test(i)? "1" : "0");
		}
		return oss.str();
	}

private:
	const uint32_t len_;
	uint32_t _cnt;
	uint32_t _size;
	uint32_t *_words; // storage
};

#endif /* BITSET_H_ */
