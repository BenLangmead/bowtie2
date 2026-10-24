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

#ifndef RANDOM_UTIL_H_
#define RANDOM_UTIL_H_

#include <algorithm>
#include <limits>
#include "random_source.h"
#include "ds.h"
#if defined(SSE_AVX2) && defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

/**
 * Return true iff v equals one of s[0], ..., s[k-1].  The same answer as
 * comparing them one at a time; only the number of compares per step differs.
 */
static inline bool random1toNSeen(const size_t* s, size_t k, size_t v) {
	size_t i = 0;
#if defined(SSE_AVX2) && defined(__AVX2__)
	if(sizeof(size_t) == 8) {
		const __m256i nv = _mm256_set1_epi64x((long long)v);
		for(; i + 8 <= k; i += 8) {
			__m256i a = _mm256_cmpeq_epi64(_mm256_loadu_si256((const __m256i*)(s + i)), nv);
			__m256i b = _mm256_cmpeq_epi64(_mm256_loadu_si256((const __m256i*)(s + i + 4)), nv);
			a = _mm256_or_si256(a, b);
			if(!_mm256_testz_si256(a, a)) {
				return true;
			}
		}
		if(i + 4 <= k) {
			__m256i a = _mm256_cmpeq_epi64(_mm256_loadu_si256((const __m256i*)(s + i)), nv);
			if(!_mm256_testz_si256(a, a)) {
				return true;
			}
			i += 4;
		}
	}
#elif defined(__SSE2__)
	if(sizeof(size_t) == 8) {
		// SSE2 has no 64-bit compare: a 64-bit lane is equal iff both of its
		// 32-bit halves are
		const __m128i nv = _mm_set1_epi64x((long long)v);
		for(; i + 4 <= k; i += 4) {
			__m128i a = _mm_cmpeq_epi32(_mm_loadu_si128((const __m128i*)(s + i)), nv);
			__m128i b = _mm_cmpeq_epi32(_mm_loadu_si128((const __m128i*)(s + i + 2)), nv);
			a = _mm_and_si128(a, _mm_shuffle_epi32(a, _MM_SHUFFLE(2, 3, 0, 1)));
			b = _mm_and_si128(b, _mm_shuffle_epi32(b, _MM_SHUFFLE(2, 3, 0, 1)));
			if(_mm_movemask_epi8(_mm_or_si128(a, b)) != 0) {
				return true;
			}
		}
	}
#endif
	for(; i < k; i++) {
		if(s[i] == v) {
			return true;
		}
	}
	return false;
}

/**
 * Return a random integer in [1, N].  Each time it's called it samples again
 * without replacement.  done() indicates when all elements have been given
 * out.
 */
class Random1toN {

	typedef size_t T;

public:

	// A set with fewer than this many elements should kick us into swap-list
	// mode immediately.  Otherwise we start in seen-list mode and then
	// possibly proceed to swap-list mode later.
	static const size_t SWAPLIST_THRESH;
	
	// Convert seen-list to swap-list after this many entries in the seen-list.
	static const size_t CONVERSION_THRESH;

	// Convert seen-list to swap-list after this (this times n_) many entries
	// in the seen-list.
	static const float CONVERSION_FRAC;

	Random1toN(int cat = 0) :
		sz_(0), n_(0), cur_(0),
		list_(SWAPLIST_THRESH, cat), seen_(CONVERSION_THRESH, cat),
		thresh_(0) {}
	
	Random1toN(size_t n, int cat = 0) :
		sz_(0), n_(n), cur_(0),
		list_(SWAPLIST_THRESH, cat), seen_(CONVERSION_THRESH, cat),
		thresh_(0) {}

	/**
	 * Initialize the set of pseudo-randoms to be given out without replacement.
	 */
	void init(size_t n, bool withoutReplacement) {
		sz_ = n_ = n;
		converted_ = false;
		swaplist_ = n < SWAPLIST_THRESH || withoutReplacement;
		cur_ = 0;
		list_.clear();
		seen_.clear();
		thresh_ = std::max(CONVERSION_THRESH, (size_t)(CONVERSION_FRAC * n));
	}
	
	/**
	 * Reset in preparation for giving out a fresh collection of pseudo-randoms
	 * without replacement.
	 */
	void reset() {
		sz_ = n_ = cur_ = 0; swaplist_ = converted_ = false;
		list_.clear(); seen_.clear();
		thresh_ = 0;
	}

	/**
	 * Get next pseudo-random element without replacement.
	 */
	T next(RandomSource& rnd) {
		assert(!done());
		if(cur_ == 0 && !converted_) {
			// This is the first call to next()
			if(n_ == 1) {
				// Trivial case: set of 1
				cur_ = 1;
				return 0;
			}
			if(swaplist_) {
				// The set is small, so we go immediately to the random
				// swapping list
				list_.resize(n_);
				for(size_t i = 0; i < n_; i++) {
					list_[i] = (T)i;
				}
			}
		}
		if(swaplist_) {
			// Get next pseudo-random using the swap-list
			size_t r = cur_ + (rnd.nextU32() % (n_ - cur_));
			if(r != cur_) {
				std::swap(list_[cur_], list_[r]);
			}
			return list_[cur_++];
		} else {
			assert(!converted_);
			// Get next pseudo-random but reject it if it's in the seen-list
			T rn = 0;
			size_t seenSz = seen_.size();
			const T* seen = seen_.ptr();
			do {
				rn = rnd.nextU32() % (T)n_;
			} while(random1toNSeen(seen, seenSz, rn));
			// Add it to the seen-list
			seen_.push_back(rn);
			cur_++;
			assert_leq(cur_, n_);
			// Move on to using the swap-list?
			assert_gt(thresh_, 0);
			if(seen_.size() >= thresh_ && cur_ < n_) {
				// Add all elements not already in the seen list to the
				// swap-list
				assert(!seen_.empty());
				// The elements not in the seen list, in increasing order, as
				// sorting the seen list and filling the gaps gave: write 0..n-1,
				// overwrite the seen ones with a value no element has, and keep
				// the others in order
				const T gone = std::numeric_limits<T>::max();
				list_.resize(n_);
				T* l = list_.ptr();
				for(size_t j = 0; j < n_; j++) {
					l[j] = (T)j;
				}
				const T* sp = seen_.ptr(); // seenSz + 1 entries now
				for(size_t i = 0; i <= seenSz; i++) {
					l[sp[i]] = gone;
				}
				size_t cur = 0;
				for(size_t j = 0; j < n_; j++) {
					const T v = l[j];
					l[cur] = v;
					cur += (v != gone);
				}
				assert_eq(cur, n_ - cur_);
				list_.resize(cur);
				seen_.clear();
				cur_ = 0;
				n_ = list_.size();
				converted_ = true;
				swaplist_ = true;
			}
			return rn;
		}
	}
	
	/**
	 * Return true iff the generator was initialized.
	 */
	bool inited() const { return n_ > 0; }
	
	/**
	 * Set so that there are no pseudo-randoms remaining.
	 */
	void setDone() { assert(inited()); cur_ = n_; }
	
	/**
	 * Return true iff all pseudo-randoms have already been given out.
	 */
	bool done() const { return inited() && cur_ >= n_; }

	/**
	 * Return the total number of pseudo-randoms we are initialized to give
	 * out, including ones already given out.
	 */
	size_t size() const { return n_; }
	
	/**
	 * Return the number of pseudo-randoms left to give out.
	 */
	size_t left() const { return n_ - cur_; }

	/**
	 * Return the total size occupued by the Descent driver and all its
	 * constituent parts.
	 */
	size_t totalSizeBytes() const {
		return list_.totalSizeBytes() +
		       seen_.totalSizeBytes();
	}

	/**
	 * Return the total capacity of the Descent driver and all its constituent
	 * parts.
	 */
	size_t totalCapacityBytes() const {
		return list_.totalCapacityBytes() +
		       seen_.totalCapacityBytes();
	}

protected:

	size_t   sz_;        // domain to pick elts from
	size_t   n_;         // number of elements in active list
	bool     swaplist_;  // if small, use swapping
	bool     converted_; // true iff seen-list was converted to swap-list
	size_t   cur_;       // # times next() was called
	EList<T> list_;      // pseudo-random swapping list
	EList<T> seen_;      // prior to swaplist_ mode, list of
	                     // pseudo-randoms given out
	size_t   thresh_;    // conversion threshold for this instantiation, which
	                     // depends both on CONVERSION_THRESH and on n_
};

#endif
