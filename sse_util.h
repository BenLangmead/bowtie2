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

#ifndef SSE_UTIL_H_
#define SSE_UTIL_H_

#include "assert_helpers.h"
#include "ds.h"
#include "limit.h"
#include <iostream>
#include <emmintrin.h>

class EList_m128i {
public:

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList_m128i(int cat = 0) :
		cat_(cat), last_alloc_(NULL), list_(NULL), sz_(0), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Destructor.
	 */
	~EList_m128i() { free(); }

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return number of elements allocated.
	 */
	inline size_t capacity() const { return sz_; }
	
	/**
	 * Ensure that there is sufficient capacity to expand to include
	 * 'thresh' more elements without having to expand.
	 */
	inline void ensure(size_t thresh) {
		if(list_ == NULL) lazyInit();
		expandCopy(cur_ + thresh);
	}

	/**
	 * Ensure that there is sufficient capacity to include 'newsz' elements.
	 * If there isn't enough capacity right now, expand capacity to exactly
	 * equal 'newsz'.
	 */
	inline void reserveExact(size_t newsz) {
		if(list_ == NULL) lazyInitExact(newsz);
		expandCopyExact(newsz);
	}

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }
	
	/**
	 * Return true iff list hasn't been initialized yet.
	 */
	inline bool null() const { return list_ == NULL; }

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandCopy(sz);
		}
		cur_ = sz;
	}
	
	/**
	 * Zero out contents of vector.
	 */
	void zero() {
		if(cur_ > 0) {
			memset(list_, 0, cur_ * sizeof(__m128i));
		}
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.  Do not copy the elements over.
	 */
	void resizeNoCopy(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandNoCopy(sz);
		}
		cur_ = sz;
	}

	/**
	 * If size is less than requested size, resize up to exactly sz and set
	 * cur_ to requested sz.
	 */
	void resizeExact(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInitExact(sz);
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopyExact(sz);
		cur_ = sz;
	}

	/**
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline __m128i& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline __m128i operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline __m128i& get(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	inline __m128i get(size_t i) const {
		return operator[](i);
	}

	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	__m128i *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	const __m128i *ptr() const { return list_; }

	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }

private:

	/**
	 * Initialize memory for EList.
	 */
	void lazyInit() {
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}

	/**
	 * Initialize exactly the prescribed number of elements for EList.
	 */
	void lazyInitExact(size_t sz) {
		assert_gt(sz, 0);
		assert(list_ == NULL);
		sz_ = sz;
		list_ = alloc(sz);
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	__m128i *alloc(size_t sz) {
		__m128i* last_alloc_;
		try {
			last_alloc_ = new __m128i[sz + 2];
		} catch(std::bad_alloc& e) {
			std::cerr << "Error: Out of memory allocating " << sz << " __m128i's for DP matrix: '" << e.what() << "'" << std::endl;
			throw e;
		}
		__m128i* tmp = last_alloc_;
		size_t tmpint = (size_t)tmp;
		// Align it!
		if((tmpint & 0xf) != 0) {
			tmpint += 15;
			tmpint &= (~0xf);
			tmp = reinterpret_cast<__m128i*>(tmpint);
		}
		assert_eq(0, (tmpint & 0xf)); // should be 16-byte aligned
		assert(tmp != NULL);
		gMemTally.add(cat_, sz);
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] last_alloc_;
			gMemTally.del(cat_, sz_);
			list_ = NULL;
			sz_ = cur_ = 0;
		}
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.  Size
	 * increases quadratically with number of expansions.  Copy old contents
	 * into new buffer using operator=.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandCopyExact(newsz);
	}

	/**
	 * Expand the list_ buffer until it has exactly 'newsz' elements.  Copy
	 * old contents into new buffer using operator=.
	 */
	void expandCopyExact(size_t newsz) {
		if(newsz <= sz_) return;
		__m128i* tmp = alloc(newsz);
		assert(tmp != NULL);
		size_t cur = cur_;
		if(list_ != NULL) {
 			for(size_t i = 0; i < cur_; i++) {
				// Note: operator= is used
				tmp[i] = list_[i];
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
		cur_ = cur;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Size increases quadratically with number of expansions.  Don't copy old
	 * contents into the new buffer.
	 */
	void expandNoCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandNoCopyExact(newsz);
	}

	/**
	 * Expand the list_ buffer until it has exactly 'newsz' elements.  Don't
	 * copy old contents into the new buffer.
	 */
	void expandNoCopyExact(size_t newsz) {
		assert(list_ != NULL);
		assert_gt(newsz, 0);
		free();
		__m128i* tmp = alloc(newsz);
		assert(tmp != NULL);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	int      cat_;        // memory category, for accounting purposes
	__m128i* last_alloc_; // what new[] originally returns
	__m128i *list_;       // list ptr, aligned version of what new[] returns
	size_t   sz_;         // capacity
	size_t   cur_;        // occupancy (AKA size)
};

struct  CpQuad {
	CpQuad() { reset(); }
	
	void reset() { sc[0] = sc[1] = sc[2] = sc[3] = 0; }
	
	bool operator==(const CpQuad& o) const {
		return sc[0] == o.sc[0] &&
		       sc[1] == o.sc[1] &&
			   sc[2] == o.sc[2] &&
			   sc[3] == o.sc[3];
	}

	int16_t sc[4];
};

/**
 * Encapsulates a collection of checkpoints.  Assumes the scheme is to
 * checkpoint adjacent pairs of anti-diagonals.
 */
class Checkpointer {

public:

	Checkpointer() { reset(); }
	
	/**
	 * Set the checkpointer up for a new rectangle.
	 */
	void init(
		size_t nrow,          // # of rows
		size_t ncol,          // # of columns
		size_t perpow2,       // checkpoint every 1 << perpow2 diags (& next)
		int64_t perfectScore, // what is a perfect score?  for sanity checks
		bool is8,             // 8-bit?
		bool doTri,           // triangle shaped?
		bool local,           // is alignment local?  for sanity checks
		bool debug)           // gather debug checkpoints?
	{
		assert_gt(perpow2, 0);
		nrow_ = nrow;
		ncol_ = ncol;
		perpow2_ = perpow2;
		per_ = 1 << perpow2;
		lomask_ = ~(0xffffffff << perpow2);
		perf_ = perfectScore;
		local_ = local;
		ndiag_ = (ncol + nrow - 1 + 1) / per_;
		locol_ = MAX_SIZE_T;
		hicol_ = MIN_SIZE_T;
		debug_ = debug;
		commitMap_.clear();
		firstCommit_ = true;
		size_t perword = (is8 ? 16 : 8);
		is8_ = is8;
		niter_ = ((nrow_ + perword - 1) / perword);
		if(doTri) {
			// Save a pair of anti-diagonals every per_ anti-diagonals for
			// backtrace purposes
			qdiag1s_.resize(ndiag_ * nrow_);
			qdiag2s_.resize(ndiag_ * nrow_);
		} else {
			// Save every per_ columns and rows for backtrace purposes
			qrows_.resize((nrow_ / per_) * ncol_);
			qcols_.resize((ncol_ / per_) * (niter_ << 2));
		}
		if(debug_) {
			// Save all columns for debug purposes
			qcolsD_.resize(ncol_ * (niter_ << 2));
		}
	}
	
	/**
	 * Return true iff we've been collecting debug cells.
	 */
	bool debug() const { return debug_; }
	
	/**
	 * Check whether the given score matches the saved score at row, col, hef.
	 */
	int64_t debugCell(size_t row, size_t col, int hef) const {
		assert(debug_);
		const __m128i* ptr = qcolsD_.ptr() + hef;
		// Fast forward to appropriate column
		ptr += ((col * niter_) << 2);
		size_t mod = row % niter_; // which m128i
		size_t div = row / niter_; // offset into m128i
		// Fast forward to appropriate word
		ptr += (mod << 2);
		// Extract score
		int16_t sc = (is8_ ? ((uint8_t*)ptr)[div] : ((int16_t*)ptr)[div]);
		int64_t asc = MIN_I64;
		// Convert score
		if(is8_) {
			if(local_) {
				asc = sc;
			} else {
				if(sc == 0) asc = MIN_I64;
				else asc = sc - 0xff;
			}
		} else {
			if(local_) {
				asc = sc + 0x8000;
			} else {
				if(sc != MIN_I16) asc = sc - 0x7fff;
			}
		}
		return asc;
	}
	
	/**
	 * Return true iff the given row/col is checkpointed.
	 */
	bool isCheckpointed(size_t row, size_t col) const {
		assert_leq(col, hicol_);
		assert_geq(col, locol_);
		size_t mod = (row + col) & lomask_;
		assert_lt(mod, per_);
		return mod >= per_ - 2;
	}

	/**
	 * Return the checkpointed H, E, or F score from the given cell.
	 */
	inline int64_t scoreTriangle(size_t row, size_t col, int hef) const {
		assert(isCheckpointed(row, col));
		bool diag1 = ((row + col) & lomask_) == per_ - 2;
		size_t off = (row + col) >> perpow2_;
		if(diag1) {
			if(qdiag1s_[off * nrow_ + row].sc[hef] == MIN_I16) {
				return MIN_I64;
			} else {
				return qdiag1s_[off * nrow_ + row].sc[hef];
			}
		} else {
			if(qdiag2s_[off * nrow_ + row].sc[hef] == MIN_I16) {
				return MIN_I64;
			} else {
				return qdiag2s_[off * nrow_ + row].sc[hef];
			}
		}
	}

	/**
	 * Return the checkpointed H, E, or F score from the given cell.
	 */
	inline int64_t scoreSquare(size_t row, size_t col, int hef) const {
		// Is it in a checkpointed row?  Note that checkpointed rows don't
		// necessarily have the horizontal contributions calculated, so we want
		// to use the column info in that case.
		if((row & lomask_) == lomask_ && hef != 1) {
			int64_t sc = qrows_[(row >> perpow2_) * ncol_ + col].sc[hef];
			if(sc == MIN_I16) return MIN_I64;
			return sc;
		}
		hef--;
		if(hef == -1) hef = 2;
		// It must be in a checkpointed column
		assert_eq(lomask_, (col & lomask_));
		// Fast forward to appropriate column
		const __m128i* ptr = qcols_.ptr() + hef;
		ptr += (((col >> perpow2_) * niter_) << 2);
		size_t mod = row % niter_; // which m128i
		size_t div = row / niter_; // offset into m128i
		// Fast forward to appropriate word
		ptr += (mod << 2);
		// Extract score
		int16_t sc = (is8_ ? ((uint8_t*)ptr)[div] : ((int16_t*)ptr)[div]);
		int64_t asc = MIN_I64;
		// Convert score
		if(is8_) {
			if(local_) {
				asc = sc;
			} else {
				if(sc == 0) asc = MIN_I64;
				else asc = sc - 0xff;
			}
		} else {
			if(local_) {
				asc = sc + 0x8000;
			} else {
				if(sc != MIN_I16) asc = sc - 0x7fff;
			}
		}
		return asc;
	}

	/**
	 * Given a column of filled-in cells, save the checkpointed cells in cs_.
	 */
	void commitCol(__m128i *pvH, __m128i *pvE, __m128i *pvF, size_t coli);
	
	/**
	 * Reset the state of the Checkpointer.
	 */
	void reset() {
		perpow2_ = per_ = lomask_ = nrow_ = ncol_ = 0;
		local_ = false;
		niter_ = ndiag_ = locol_ = hicol_ = 0;
		perf_ = 0;
		firstCommit_ = true;
		is8_ = debug_ = false;
	}
	
	/**
	 * Return true iff the Checkpointer has been initialized.
	 */
	bool inited() const {
		return nrow_ > 0;
	}
	
	size_t per()     const { return per_;     }
	size_t perpow2() const { return perpow2_; }
	size_t lomask()  const { return lomask_;  }
	size_t locol()   const { return locol_;   }
	size_t hicol()   const { return hicol_;   }
	size_t nrow()    const { return nrow_;    }
	size_t ncol()    const { return ncol_;    }
	
	const CpQuad* qdiag1sPtr() const { return qdiag1s_.ptr(); }
	const CpQuad* qdiag2sPtr() const { return qdiag2s_.ptr(); }

	size_t   perpow2_;   // 1 << perpow2_ - 2 is the # of uncheckpointed
	                     // anti-diags between checkpointed anti-diag pairs
	size_t   per_;       // 1 << perpow2_
	size_t   lomask_;    // mask for extracting low bits
	size_t   nrow_;      // # rows in current rectangle
	size_t   ncol_;      // # cols in current rectangle
	int64_t  perf_;      // perfect score
	bool     local_;     // local alignment?
	
	size_t   ndiag_;     // # of double-diags
	
	size_t   locol_;     // leftmost column committed
	size_t   hicol_;     // rightmost column committed

	// Map for committing scores from vector columns to checkpointed diagonals
	EList<size_t> commitMap_;
	bool          firstCommit_;
	
	EList<CpQuad> qdiag1s_; // checkpoint H/E/F values for diagonal 1
	EList<CpQuad> qdiag2s_; // checkpoint H/E/F values for diagonal 2

	EList<CpQuad> qrows_;   // checkpoint H/E/F values for rows
	
	// We store columns in this way to reduce overhead of populating them
	bool          is8_;     // true -> fill used 8-bit cells
	size_t        niter_;   // # __m128i words per column
	EList_m128i   qcols_;   // checkpoint E/F/H values for select columns
	
	bool          debug_;   // get debug checkpoints? (i.e. fill qcolsD_?)
	EList_m128i   qcolsD_;  // checkpoint E/F/H values for all columns (debug)
};

#endif
