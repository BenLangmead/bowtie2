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

#ifndef RANDOM_UTIL_H_
#define RANDOM_UTIL_H_

/**
 * Return a random integer in [1, N].  Each time it's called it samples again
 * without replacement.  done() indicates when all elements have been given
 * out.
 */
class Random1toN {

	typedef uint32_t T;

public:

	// A set with fewer than this many elements should kick us into swap-list
	// mode immediately.  Otherwise we start in seen-list mode and then
	// possibly proceed to swap-list mode later.
	static const size_t SWAPLIST_THRESH = 128;
	
	// Convert seen-list to swap-list after this many 
	static const size_t CONVERSION_THRESH = 16;

	Random1toN(int cat = 0) :
		sz_(0), n_(0), cur_(0),
		list_(SWAPLIST_THRESH, cat), seen_(CONVERSION_THRESH, cat) {}
	
	Random1toN(size_t n, int cat = 0) :
		sz_(0), n_(n), cur_(0),
		list_(SWAPLIST_THRESH, cat), seen_(CONVERSION_THRESH, cat) {}

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
	}
	
	/**
	 * Reset in preparation for giving out a fresh collection of pseudo-randoms
	 * without replacement.
	 */
	void reset() {
		sz_ = n_ = cur_ = 0; swaplist_ = converted_ = false;
		list_.clear(); seen_.clear();
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
				swap(list_[cur_], list_[r]);
			}
			return list_[cur_++];
		} else {
			assert(!converted_);
			// Get next pseudo-random but reject it if it's in the seen-list
			bool again = true;
			T rn = 0;
			size_t seenSz = seen_.size();
			while(again) {
				rn = rnd.nextU32() % (T)n_;
				again = false;
				for(size_t i = 0; i < seenSz; i++) {
					if(seen_[i] == rn) {
						again = true;
						break;
					}
				}
			}
			// Add it to the seen-list
			seen_.push_back(rn);
			cur_++;
			assert_leq(cur_, n_);
			// Move on to using the swap-list?
			if(seen_.size() >= CONVERSION_THRESH && cur_ < n_) {
				// Add all elements not already in the seen list to the
				// swap-list
				assert(!seen_.empty());
				seen_.sort();
				list_.resize(n_ - cur_);
				size_t prev = 0;
				size_t cur = 0;
				for(size_t i = 0; i <= seenSz; i++) {
					// Add all the elements between the previous element and
					// this one
					for(size_t j = prev; j < seen_[i]; j++) {
						list_[cur++] = (T)j;
					}
					prev = seen_[i]+1;
				}
				for(size_t j = prev; j < n_; j++) {
					list_[cur++] = (T)j;
				}
				assert_eq(cur, n_ - cur_);
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
	bool inited() { return n_ > 0; }
	
	/**
	 * Set so that there are no pseudo-randoms remaining.
	 */
	void setDone() { assert(inited()); cur_ = n_; }
	
	/**
	 * Return true iff all pseudo-randoms have already been given out.
	 */
	bool done() { return inited() && cur_ >= n_; }

	/**
	 * Return the total number of pseudo-randoms we are initialized to give
	 * out, including ones already given out.
	 */
	size_t size() const { return n_; }
	
	/**
	 * Return the number of pseudo-randoms left to give out.
	 */
	size_t left() const { return n_ - cur_; }

protected:

	size_t          sz_;        // domain to pick elts from
	size_t          n_;         // number of elements in active list
	bool            swaplist_;  // if small, use swapping
	bool            converted_; // true iff seen-list was converted to swap-list
	size_t          cur_;       // # times next() was called
	EList<T> list_;             // pseudo-random swapping list
	EList<T> seen_;             // prior to swaplist_ mode, list of
	                            // pseudo-randoms given out
};

#endif
