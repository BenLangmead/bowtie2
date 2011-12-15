//
//  sse_util.h
//  bowtie2
//
//  Created by Benjamin Langmead on 11/8/11.
//  Copyright 2011 Johns Hopkins University. All rights reserved.
//

#ifndef SSE_UTIL_H_
#define SSE_UTIL_H_

#include "assert_helpers.h"
#include "ds.h"
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
			std::cerr << "Error: Out of memory allocating __m128i's for DP matrix: '" << e.what() << "'" << std::endl;
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

#endif
