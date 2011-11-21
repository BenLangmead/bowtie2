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

public:
	static const size_t THRESH = 128;

	Random1toN(int cat = 0) : n_(0), cur_(0), list_(THRESH, cat) {}
	Random1toN(size_t n, int cat = 0) : n_(n), cur_(0), list_(THRESH, cat) {}

	void init(size_t n) {
		n_ = n;
		small_ = n < THRESH;
		cur_ = 0;
		list_.clear();
	}
	
	void reset() {
		n_ = cur_ = 0; small_ = false; list_.clear();
	}

	uint32_t next(RandomSource& rnd) {
		assert(!done());
		if(cur_ == 0) {
			if(n_ == 1) {
				cur_ = 1;
				return 0;
			}
			if(small_) {
				// Initialize list
				list_.resize(n_);
				for(size_t i = 0; i < n_; i++) {
					list_[i] = ((uint32_t)i);
				}
			}
		}
		if(small_) {
			size_t r = cur_ + (rnd.nextU32() % (n_ - cur_));
			if(r != cur_) {
				uint32_t tmp = list_[cur_];
				list_[cur_] = list_[r];
				list_[r] = tmp;
			}
			return list_[cur_++];
		} else {
			// We have a large domain, so we'll draw a random without regard to
			// previous draws
			return rnd.nextU32() % (uint32_t)n_;
		}
	}
	
	bool inited() { return n_ > 0; }
	
	void setDone() { assert(inited()); cur_ = n_; }
	
	bool done() { return inited() && cur_ >= n_; }

	size_t size() const { return n_; }
	
	size_t left() const { return n_ - cur_; }

protected:

	size_t          n_;     // domain to pick elts from
	bool            small_; // if small, use swapping
	size_t          cur_;   // # times next() was called
	EList<uint32_t> list_;  // state to remember what's been picked
};

#endif
