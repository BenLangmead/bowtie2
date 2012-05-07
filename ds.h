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

#ifndef DS_H_
#define DS_H_

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <stdint.h>
#include <string.h>
#include "assert_helpers.h"
#include "threading.h"

/**
 * Tally how much memory is allocated to certain 
 */
class MemoryTally {

public:

	MemoryTally() : tot_(0), peak_(0) {
		memset(tots_,  0, 256 * sizeof(uint64_t));
		memset(peaks_, 0, 256 * sizeof(uint64_t));
		MUTEX_INIT(lock_);
	}
	
	/**
	 * Tally a memory allocation of size amt bytes.
	 */
	void add(int cat, uint64_t amt) {
		ThreadSafe ts(&lock_);
		tots_[cat] += amt;
		tot_ += amt;
		if(tots_[cat] > peaks_[cat]) {
			peaks_[cat] = tots_[cat];
		}
		if(tot_ > peak_) {
			peak_ = tot_;
		}
	}

	/**
	 * Tally a memory free of size amt bytes.
	 */
	void del(int cat, uint64_t amt) {
		ThreadSafe ts(&lock_);
		assert_geq(tots_[cat], amt);
		assert_geq(tot_, amt);
		tots_[cat] -= amt;
		tot_ -= amt;
	}
	
	/**
	 * Return the total amount of memory allocated.
	 */
	uint64_t total() { return tot_; }

	/**
	 * Return the total amount of memory allocated in a particular
	 * category.
	 */
	uint64_t total(int cat) { return tots_[cat]; }

	/**
	 * Return the peak amount of memory allocated.
	 */
	uint64_t peak() { return peak_; }

	/**
	 * Return the peak amount of memory allocated in a particular
	 * category.
	 */
	uint64_t peak(int cat) { return peaks_[cat]; }
	
	/**
	 * Check that memory tallies are internally consistent;
	 */
	bool repOk() const {
		uint64_t tot = 0;
		for(int i = 0; i < 256; i++) {
			assert_leq(tots_[i], peaks_[i]);
			tot += tots_[i];
		}
		assert_eq(tot, tot_);
		return true;
	}

protected:

	MUTEX_T lock_;
	uint64_t tots_[256];
	uint64_t tot_;
	uint64_t peaks_[256];
	uint64_t peak_;
};

extern MemoryTally gMemTally;

/**
 * A simple fixed-length array of type T, automatically freed in the
 * destructor.
 */
template<typename T>
class AutoArray {
public:

	AutoArray(size_t sz, int cat = 0) : cat_(cat) {
		t_ = NULL;
		t_ = new T[sz];
		gMemTally.add(cat_, sz);
		memset(t_, 0, sz * sizeof(T));
		sz_ = sz;
	}
	
	~AutoArray() {
		if(t_ != NULL) {
			delete[] t_;
			gMemTally.del(cat_, sz_);
		}
	}
	
	T& operator[](size_t sz) {
		return t_[sz];
	}
	
	const T& operator[](size_t sz) const {
		return t_[sz];
	}
	
	size_t size() const { return sz_; }

private:
	int cat_;
	T *t_;
	size_t sz_;
};

/**
 * A wrapper for a non-array pointer that associates it with a memory
 * category for tracking purposes and calls delete on it when the
 * PtrWrap is destroyed.
 */
template<typename T>
class PtrWrap {
public:

	explicit PtrWrap(
		T* p,
		bool freeable = true,
		int cat = 0) :
		cat_(cat),
		p_(NULL)
	{
		init(p, freeable);
	}

	explicit PtrWrap(int cat = 0) :
		cat_(cat),
		p_(NULL)
	{
		reset();
	}

	void reset() {
		free();
		init(NULL);
	}

	~PtrWrap() { free(); }
	
	void init(T* p, bool freeable = true) {
		assert(p_ == NULL);
		p_ = p;
		freeable_ = freeable;
		if(p != NULL && freeable_) {
			gMemTally.add(cat_, sizeof(T));
		}
	}
	
	void free() {
		if(p_ != NULL) {
			if(freeable_) {
				delete p_;
				gMemTally.del(cat_, sizeof(T));
			}
			p_ = NULL;
		}
	}
	
	inline T* get() { return p_; }
	inline const T* get() const { return p_; }

private:
	int cat_;
	T *p_;
	bool freeable_;
};

/**
 * A wrapper for an array pointer that associates it with a memory
 * category for tracking purposes and calls delete[] on it when the
 * PtrWrap is destroyed.
 */
template<typename T>
class APtrWrap {
public:

	explicit APtrWrap(
		T* p,
		size_t sz,
		bool freeable = true,
		int cat = 0) :
		cat_(cat),
		p_(NULL)
	{
		init(p, sz, freeable);
	}

	explicit APtrWrap(int cat = 0) :
		cat_(cat),
		p_(NULL)
	{
		reset();
	}
	
	void reset() {
		free();
		init(NULL, 0);
	}

	~APtrWrap() { free(); }
	
	void init(T* p, size_t sz, bool freeable = true) {
		assert(p_ == NULL);
		p_ = p;
		sz_ = sz;
		freeable_ = freeable;
		if(p != NULL && freeable_) {
			gMemTally.add(cat_, sizeof(T) * sz_);
		}
	}
	
	void free() {
		if(p_ != NULL) {
			if(freeable_) {
				delete[] p_;
				gMemTally.del(cat_, sizeof(T) * sz_);
			}
			p_ = NULL;
		}
	}
	
	inline T* get() { return p_; }
	inline const T* get() const { return p_; }

private:
	int cat_;
	T *p_;
	bool freeable_;
	size_t sz_;
};

/**
 * An EList<T> is an expandable list with these features:
 *
 *  - Payload type is a template parameter T.
 *  - Initial size can be specified at construction time, otherwise
 *    default of 128 is used.
 *  - When allocated initially or when expanding, the new[] operator is
 *    used, which in turn calls the default constructor for T.
 *  - All copies (e.g. assignment of a const T& to an EList<T> element,
 *    or during expansion) use operator=.
 *  - When the EList<T> is resized to a smaller size (or cleared, which
 *    is like resizing to size 0), the underlying containing is not
 *    reshaped.  Thus, ELists<T>s never release memory before
 *    destruction.
 *
 * And these requirements:
 *
 *  - Payload type T must have a default constructor.
 *
 * For efficiency reasons, ELists should not be declared on the stack
 * in often-called worker functions.  Best practice is to declare
 * ELists at a relatively stable layer of the stack (such that it
 * rarely bounces in and out of scope) and let the worker function use
 * it and *expand* it only as needed.  The effect is that only
 * relatively few allocations and copies will be incurred, and they'll
 * occur toward the beginning of the computation before stabilizing at
 * a "high water mark" for the remainder of the computation.
 *
 * A word about multidimensional lists.  One way to achieve a
 * multidimensional lists is to nest ELists.  This works, but it often
 * involves a lot more calls to the default constructor and to
 * operator=, especially when the outermost EList needs expanding, than
 * some of the alternatives.  One alternative is use a most specialized
 * container that still uses ELists but knows to use xfer instead of
 * operator= when T=EList.
 *
 * The 'cat_' fiends encodes a category.  This makes it possible to
 * distinguish between object subgroups in the global memory tally.
 *
 * Memory allocation is lazy.  Allocation is only triggered when the
 * user calls push_back, expand, resize, or another function that
 * increases the size of the list.  This saves memory and also makes it
 * easier to deal with nested ELists, since the default constructor
 * doesn't set anything in stone.
 */
template <typename T, int S = 128>
class EList {

public:

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList() :
		cat_(0), allocCat_(-1), list_(NULL), sz_(S), cur_(0) { }

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList(int cat) :
		cat_(cat), allocCat_(-1), list_(NULL), sz_(S), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	explicit EList(size_t isz, int cat = 0) :
		cat_(cat), allocCat_(-1), list_(NULL), sz_(isz), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Copy from another EList using operator=.
	 */
	EList(const EList<T, S>& o) :
		cat_(0), allocCat_(-1), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}

	/**
	 * Copy from another EList using operator=.
	 */
	explicit EList(const EList<T, S>& o, int cat) :
		cat_(cat), allocCat_(-1), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}

	/**
	 * Destructor.
	 */
	~EList() { free(); }

	/**
	 * Make this object into a copy of o by allocat
	 */
	EList<T, S>& operator=(const EList<T, S>& o) {
		assert_eq(cat_, o.cat());
		if(o.cur_ == 0) {
			// Nothing to copy
			cur_ = 0;
			return *this;
		}
		if(list_ == NULL) {
			// cat_ should already be set
			lazyInit();
		}
		if(sz_ < o.cur_) expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for(size_t i = 0; i < cur_; i++) {
			list_[i] = o.list_[i];
		}
		return *this;
	}
	
	/**
	 * Transfer the guts of another EList into this one without using
	 * operator=, etc.  We have to set EList o's list_ field to NULL to
	 * avoid o's destructor from deleting list_ out from under us.
	 */
	void xfer(EList<T, S>& o) {
		// What does it mean to transfer to a different-category list?
		assert_eq(cat_, o.cat());
		// Can only transfer into an empty object
		free();
		allocCat_ = cat_;
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
		o.allocCat_ = -1;
	}

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
	 * Add an element to the back and immediately initialize it via
	 * operator=.
	 */
	void push_back(const T& el) {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		list_[cur_++] = el;
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void expand() {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void fill(size_t begin, size_t end, const T& v) {
		assert_leq(begin, end);
		assert_leq(end, cur_);
		for(size_t i = begin; i < end; i++) {
			list_[i] = v;
		}
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void fill(const T& v) {
		for(size_t i = 0; i < cur_; i++) {
			list_[i] = v;
		}
	}

	/**
	 * Set all bits in specified range of elements in list array to 0.
	 */
	void fillZero(size_t begin, size_t end) {
		assert_leq(begin, end);
		memset(&list_[begin], 0, sizeof(T) * (end-begin));
	}

	/**
	 * Set all bits in the list array to 0.
	 */
	void fillZero() {
		memset(list_, 0, sizeof(T) * cur_);
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resizeNoCopy(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandNoCopy(sz);
		cur_ = sz;
	}

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
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
	}

	/**
	 * Erase range of elements starting at offset idx and going for len.
	 */
	void erase(size_t idx, size_t len) {
		assert_geq(len, 0);
		if(len == 0) {
			return;
		}
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-len; i++) {
			list_[i] = list_[i+len];
		}
		cur_ -= len;
	}

	/**
	 * Insert value 'el' at offset 'idx'
	 */
	void insert(const T& el, size_t idx) {
		if(list_ == NULL) lazyInit();
		assert_leq(idx, cur_);
		if(cur_ == sz_) expandCopy(sz_+1);
		for(size_t i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
	}

	/**
	 * Insert contents of list 'l' at offset 'idx'
	 */
	void insert(const EList<T>& l, size_t idx) {
		if(list_ == NULL) lazyInit();
		assert_lt(idx, cur_);
		if(l.cur_ == 0) return;
		if(cur_ + l.cur_ > sz_) expandCopy(cur_ + l.cur_);
		for(size_t i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--) {
			list_[i] = list_[i - l.cur_];
		}
		for(size_t i = 0; i < l.cur_; i++) {
			list_[i+idx] = l.list_[i];
		}
		cur_ += l.cur_;
	}

	/**
	 * Remove an element from the top of the stack.
	 */
	void pop_back() {
		assert_gt(cur_, 0);
		cur_--;
	}

	/**
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Get the element on the top of the stack.
	 */
	inline T& back() {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Reverse list elements.
	 */
	void reverse() {
		if(cur_ > 1) {
			size_t n = cur_ >> 1;
			for(size_t i = 0; i < n; i++) {
				T tmp = list_[i];
				list_[i] = list_[cur_ - i - 1];
				list_[cur_ - i - 1] = tmp;
			}
		}
	}

	/**
	 * Get the element on the top of the stack, const version.
	 */
	inline const T& back() const {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the frontmost element (bottom of stack).
	 */
	inline T& front() {
		assert_gt(cur_, 0);
		return list_[0];
	}

	/**
	 * Get the element on the bottom of the stack, const version.
	 */
	inline const T& front() const { return front(); }

	/**
	 * Return true iff this list and list o contain the same elements in the
	 * same order according to type T's operator==.
	 */
	bool operator==(const EList<T, S>& o) const {
		if(size() != o.size()) {
			return false;
		}
		for(size_t i = 0; i < size(); i++) {
			if(!(get(i) == o.get(i))) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true iff this list contains all of the elements in o according to
	 * type T's operator==.
	 */
	bool isSuperset(const EList<T, S>& o) const {
		if(o.size() > size()) {
			// This can't be a superset if the other set contains more elts
			return false;
		}
		// For each element in o
		for(size_t i = 0; i < o.size(); i++) {
			bool inthis = false;
			// Check if it's in this
			for(size_t j = 0; j < size(); j++) {
				if(o[i] == (*this)[j]) {
					inthis = true;
					break;
				}
			}
			if(!inthis) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline T& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const T& operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline T& get(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	inline const T& get(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	T& getSlow(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	const T& getSlow(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Sort some of the contents.
	 */
	void sortPortion(size_t begin, size_t num) {
		assert_leq(begin+num, cur_);
		if(num < 2) return;
		std::sort(list_ + begin, list_ + begin + num);
	}
	
	/**
	 * Sort contents
	 */
	void sort() {
		sortPortion(0, cur_);
	}

	/**
	 * Return true iff every element is < its successor.  Only operator< is
	 * used.
	 */
	bool sorted() const {
		for(size_t i = 1; i < cur_; i++) {
			if(!(list_[i-1] < list_[i])) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Delete element at position 'idx'; slide subsequent chars up.
	 */
	void remove(size_t idx) {
		assert_lt(idx, cur_);
		assert_gt(cur_, 0);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
	}
	
	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	T *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	const T *ptr() const { return list_; }

	/**
	 * Set the memory category for this object.
	 */
	void setCat(int cat) {
		// What does it mean to set the category after the list_ is
		// already allocated?
		assert(null());
		assert_gt(cat, 0); cat_ = cat;
	}

	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }

	/**
	 * Perform a binary search for the first element that is not less
	 * than 'el'.  Return cur_ if all elements are less than el.
	 */
	size_t bsearchLoBound(const T& el) const {
		size_t hi = cur_;
		size_t lo = 0;
		while(true) {
			if(lo == hi) {
				return lo;
			}
			size_t mid = lo + ((hi-lo)>>1);
			assert_neq(mid, hi);
			if(list_[mid] < el) {
				if(lo == mid) {
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
	}

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
	T *alloc(size_t sz) {
		T* tmp = new T[sz];
		assert(tmp != NULL);
		gMemTally.add(cat_, sz);
		allocCat_ = cat_;
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			assert_neq(-1, allocCat_);
			assert_eq(allocCat_, cat_);
			delete[] list_;
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
		T* tmp = alloc(newsz);
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
		T* tmp = alloc(newsz);
		assert(tmp != NULL);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	int cat_;      // memory category, for accounting purposes
	int allocCat_; // category at time of allocation
	T *list_;      // list pointer, returned from new[]
	size_t sz_;    // capacity
	size_t cur_;   // occupancy (AKA size)
};

/**
 * An ELList<T> is an expandable list of lists with these features:
 *
 *  - Payload type of the inner list is a template parameter T.
 *  - Initial size can be specified at construction time, otherwise
 *    default of 128 is used.
 *  - When allocated initially or when expanding, the new[] operator is
 *    used, which in turn calls the default constructor for EList<T>.
 *  - Upon expansion, instead of copies, xfer is used.
 *  - When the ELList<T> is resized to a smaller size (or cleared,
 *    which is like resizing to size 0), the underlying containing is
 *    not reshaped.  Thus, ELLists<T>s never release memory before
 *    destruction.
 *
 * And these requirements:
 *
 *  - Payload type T must have a default constructor.
 *
 */
template <typename T, int S1 = 128, int S2 = 128>
class ELList {

public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	explicit ELList(int cat = 0) :
		cat_(cat), list_(NULL), sz_(S2), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	explicit ELList(size_t isz, int cat = 0) :
		cat_(cat), list_(NULL), sz_(isz), cur_(0)
	{
		assert_gt(isz, 0);
		assert_geq(cat, 0);
	}

	/**
	 * Copy from another ELList using operator=.
	 */
	ELList(const ELList<T, S1, S2>& o) :
		cat_(0), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}

	/**
	 * Copy from another ELList using operator=.
	 */
	explicit ELList(const ELList<T, S1, S2>& o, int cat) :
		cat_(cat), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}

	/**
	 * Destructor.
	 */
	~ELList() { free(); }

	/**
	 * Make this object into a copy of o by allocating enough memory to
	 * fit the number of elements in o (note: the number of elements
	 * may be substantially less than the memory allocated in o) and
	 * using operator= to copy them over.
	 */
	ELList<T, S1, S2>& operator=(const ELList<T, S1, S2>& o) {
		assert_eq(cat_, o.cat());
		if(list_ == NULL) {
			lazyInit();
		}
		if(o.cur_ == 0) {
			cur_ = 0;
			return *this;
		}
		if(sz_ < o.cur_) expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for(size_t i = 0; i < cur_; i++) {
			// Note: using operator=, not xfer
			assert_eq(list_[i].cat(), o.list_[i].cat());
			list_[i] = o.list_[i];
		}
		return *this;
	}
	
	/**
	 * Transfer the guts of another EList into this one without using
	 * operator=, etc.  We have to set EList o's list_ field to NULL to
	 * avoid o's destructor from deleting list_ out from under us.
	 */
	void xfer(ELList<T, S1, S2>& o) {
		assert_eq(cat_, o.cat());
		list_ = o.list_; // list_ is an array of EList<T>s
		sz_   = o.sz_;
		cur_  = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }

	/**
	 * Return true iff list hasn't been initialized yet.
	 */
	inline bool null() const { return list_ == NULL; }

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void expand() {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
	}

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
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Get the element on the top of the stack.
	 */
	inline EList<T, S1>& back() {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the element on the top of the stack, const version.
	 */
	inline const EList<T, S1>& back() const {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the frontmost element (bottom of stack).
	 */
	inline EList<T, S1>& front() {
		assert_gt(cur_, 0);
		return list_[0];
	}

	/**
	 * Get the element on the bottom of the stack, const version.
	 */
	inline const EList<T, S1>& front() const { return front(); }

	/**
	 * Return a reference to the ith element.
	 */
	inline EList<T, S1>& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const EList<T, S1>& operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline EList<T, S1>& get(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	inline const EList<T, S1>& get(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	EList<T, S1>& getSlow(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	const EList<T, S1>& getSlow(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	EList<T, S1> *ptr() { return list_; }
	
	/**
	 * Set the memory category for this object and all children.
	 */
	void setCat(int cat) {
		assert_gt(cat, 0);
		cat_ = cat;
		if(cat_ != 0) {
			for(size_t i = 0; i < sz_; i++) {
				assert(list_[i].null());
				list_[i].setCat(cat_);
			}
		}
	}

	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }

protected:

	/**
	 * Initialize memory for EList.
	 */
	void lazyInit() {
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	EList<T, S1> *alloc(size_t sz) {
		assert_gt(sz, 0);
		EList<T, S1> *tmp = new EList<T, S1>[sz];
		gMemTally.add(cat_, sz);
		if(cat_ != 0) {
			for(size_t i = 0; i < sz; i++) {
				assert(tmp[i].ptr() == NULL);
				tmp[i].setCat(cat_);
			}
		}
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] list_;
			gMemTally.del(cat_, sz_);
			list_ = NULL;
		}
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Copy old contents into new buffer
	 * using operator=.
	 */
	void expandCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		EList<T, S1>* tmp = alloc(newsz);
		if(list_ != NULL) {
			for(size_t i = 0; i < cur_; i++) {
				assert_eq(cat_, tmp[i].cat());
				tmp[i].xfer(list_[i]);
				assert_eq(cat_, tmp[i].cat());
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Don't copy old contents over.
	 */
	void expandNoCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		free();
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		EList<T, S1>* tmp = alloc(newsz);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	int cat_;    // memory category, for accounting purposes
	EList<T, S1> *list_; // list pointer, returned from new[]
	size_t sz_;  // capacity
	size_t cur_; // occupancy (AKA size)

};

/**
 * An ELLList<T> is an expandable list of expandable lists with these
 * features:
 *
 *  - Payload type of the innermost list is a template parameter T.
 *  - Initial size can be specified at construction time, otherwise
 *    default of 128 is used.
 *  - When allocated initially or when expanding, the new[] operator is
 *    used, which in turn calls the default constructor for ELList<T>.
 *  - Upon expansion, instead of copies, xfer is used.
 *  - When the ELLList<T> is resized to a smaller size (or cleared,
 *    which is like resizing to size 0), the underlying containing is
 *    not reshaped.  Thus, ELLLists<T>s never release memory before
 *    destruction.
 *
 * And these requirements:
 *
 *  - Payload type T must have a default constructor.
 *
 */
template <typename T, int S1 = 128, int S2 = 128, int S3 = 128>
class ELLList {

public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	explicit ELLList(int cat = 0) :
		cat_(cat), list_(NULL), sz_(S3), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	explicit ELLList(size_t isz, int cat = 0) :
		cat_(cat), list_(NULL), sz_(isz), cur_(0)
	{
		assert_geq(cat, 0);
		assert_gt(isz, 0);
	}

	/**
	 * Copy from another ELLList using operator=.
	 */
	ELLList(const ELLList<T, S1, S2, S3>& o) :
		cat_(0), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}

	/**
	 * Copy from another ELLList using operator=.
	 */
	explicit ELLList(const ELLList<T, S1, S2, S3>& o, int cat) :
		cat_(cat), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}

	/**
	 * Destructor.
	 */
	~ELLList() { free(); }

	/**
	 * Make this object into a copy of o by allocating enough memory to
	 * fit the number of elements in o (note: the number of elements
	 * may be substantially less than the memory allocated in o) and
	 * using operator= to copy them over.
	 */
	ELLList<T, S1, S2, S3>& operator=(const ELLList<T, S1, S2, S3>& o) {
		assert_eq(cat_, o.cat());
		if(list_ == NULL) lazyInit();
		if(o.cur_ == 0) {
			cur_ = 0;
			return *this;
		}
		if(sz_ < o.cur_) expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for(size_t i = 0; i < cur_; i++) {
			// Note: using operator=, not xfer
			assert_eq(list_[i].cat(), o.list_[i].cat());
			list_[i] = o.list_[i];
		}
		return *this;
	}
	
	/**
	 * Transfer the guts of another EList into this one without using
	 * operator=, etc.  We have to set EList o's list_ field to NULL to
	 * avoid o's destructor from deleting list_ out from under us.
	 */
	void xfer(ELLList<T, S1, S2, S3>& o) {
		assert_eq(cat_, o.cat());
		list_ = o.list_; // list_ is an array of EList<T>s
		sz_   = o.sz_;
		cur_  = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }

	/**
	 * Return true iff list hasn't been initialized yet.
	 */
	inline bool null() const { return list_ == NULL; }

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void expand() {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
	}

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
		if(sz_ < sz) expandCopy(sz);
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
	 * Get the element on the top of the stack.
	 */
	inline ELList<T, S1, S2>& back() {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the element on the top of the stack, const version.
	 */
	inline const ELList<T, S1, S2>& back() const {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the frontmost element (bottom of stack).
	 */
	inline ELList<T, S1, S2>& front() {
		assert_gt(cur_, 0);
		return list_[0];
	}

	/**
	 * Get the element on the bottom of the stack, const version.
	 */
	inline const ELList<T, S1, S2>& front() const { return front(); }

	/**
	 * Return a reference to the ith element.
	 */
	inline ELList<T, S1, S2>& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const ELList<T, S1, S2>& operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline ELList<T, S1, S2>& get(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	inline const ELList<T, S1, S2>& get(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	ELList<T, S1, S2>& getSlow(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	const ELList<T, S1, S2>& getSlow(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	ELList<T, S1, S2> *ptr() { return list_; }

	/**
	 * Set the memory category for this object and all children.
	 */
	void setCat(int cat) {
		assert_gt(cat, 0);
		cat_ = cat;
		if(cat_ != 0) {
			for(size_t i = 0; i < sz_; i++) {
				assert(list_[i].null());
				list_[i].setCat(cat_);
			}
		}
	}
	
	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }

protected:

	/**
	 * Initialize memory for EList.
	 */
	void lazyInit() {
		assert(null());
		list_ = alloc(sz_);
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	ELList<T, S1, S2> *alloc(size_t sz) {
		assert_gt(sz, 0);
		ELList<T, S1, S2> *tmp = new ELList<T, S1, S2>[sz];
		gMemTally.add(cat_, sz);
		if(cat_ != 0) {
			for(size_t i = 0; i < sz; i++) {
				assert(tmp[i].ptr() == NULL);
				tmp[i].setCat(cat_);
			}
		}
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] list_;
			gMemTally.del(cat_, sz_);
			list_ = NULL;
		}
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Copy old contents into new buffer
	 * using operator=.
	 */
	void expandCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ELList<T, S1, S2>* tmp = alloc(newsz);
		if(list_ != NULL) {
			for(size_t i = 0; i < cur_; i++) {
				assert_eq(cat_, tmp[i].cat());
				tmp[i].xfer(list_[i]);
				assert_eq(cat_, tmp[i].cat());
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Don't copy old contents over.
	 */
	void expandNoCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		free();
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ELList<T, S1, S2>* tmp = alloc(newsz);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	int cat_;    // memory category, for accounting purposes
	ELList<T, S1, S2> *list_; // list pointer, returned from new[]
	size_t sz_;  // capacity
	size_t cur_; // occupancy (AKA size)

};

/**
 * Expandable set using a heap-allocated sorted array.
 *
 * Note that the copy constructor and operator= routines perform
 * shallow copies (w/ memcpy).
 */
template <typename T>
class ESet {
public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	ESet(int cat = 0) :
		cat_(cat),
		list_(NULL),
		sz_(0),
		cur_(0)
	{
		if(sz_ > 0) {
			list_ = alloc(sz_);
		}
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	ESet(size_t isz, int cat = 0) :
		cat_(cat),
		list_(NULL),
		sz_(isz),
		cur_(0)
	{
		assert_gt(isz, 0);
		if(sz_ > 0) {
			list_ = alloc(sz_);
		}
	}

	/**
	 * Copy from another ESet.
	 */
	ESet(const ESet<T>& o, int cat = 0) :
		cat_(cat), list_(NULL)
	{
		assert_eq(cat_, o.cat());
		*this = o;
	}

	/**
	 * Destructor.
	 */
	~ESet() { free(); }

	/**
	 * Copy contents of given ESet into this ESet.
	 */
	ESet& operator=(const ESet<T>& o) {
		assert_eq(cat_, o.cat());
		sz_ = o.sz_;
		cur_ = o.cur_;
		free();
		if(sz_ > 0) {
			list_ = alloc(sz_);
			memcpy(list_, o.list_, cur_ * sizeof(T));
		} else {
			list_ = NULL;
		}
		return *this;
	}

	/**
	 * Return number of elements.
	 */
	size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	bool empty() const { return cur_ == 0; }

	/**
	 * Return true iff list isn't initialized yet.
	 */
	bool null() const { return list_ == NULL; }

	/**
	 * Insert a new element into the set in sorted order.
	 */
	bool insert(const T& el) {
		size_t i = 0;
		if(cur_ == 0) {
			insert(el, 0);
			return true;
		}
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		if(i < cur_ && list_[i] == el) return false;
		insert(el, i);
		return true;
	}

	/**
	 * Return true iff this set contains 'el'.
	 */
	bool contains(const T& el) const {
		if(cur_ == 0) {
			return false;
		}
		else if(cur_ == 1) {
			return el == list_[0];
		}
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i] == el;
	}

	/**
	 * Remove element from set.
	 */
	void remove(const T& el) {
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		assert(i != cur_ && list_[i] == el);
		erase(i);
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz <= cur_) return;
		if(sz_ < sz) expandCopy(sz);
	}

	/**
	 * Clear set without deallocating (or setting) anything.
	 */
	void clear() { cur_ = 0; }

	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }
	
	/**
	 * Set the memory category for this object.
	 */
	void setCat(int cat) {
		cat_ = cat;
	}

	/**
	 * Transfer the guts of another EList into this one without using
	 * operator=, etc.  We have to set EList o's list_ field to NULL to
	 * avoid o's destructor from deleting list_ out from under us.
	 */
	void xfer(ESet<T>& o) {
		// What does it mean to transfer to a different-category list?
		assert_eq(cat_, o.cat());
		// Can only transfer into an empty object
		free();
		list_ = o.list_;
		sz_ = o.sz_;
		cur_ = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}

	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	T *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	const T *ptr() const { return list_; }

private:

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	T *alloc(size_t sz) {
		assert_gt(sz, 0);
		T *tmp = new T[sz];
		gMemTally.add(cat_, sz);
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] list_;
			gMemTally.del(cat_, sz_);
			list_ = NULL;
		}
	}

	/**
	 * Simple linear scan that returns the index of the first element
	 * of list_ that is not less than el, or cur_ if all elements are
	 * less than el.
	 */
	size_t scanLoBound(const T& el) const {
		for(size_t i = 0; i < cur_; i++) {
			if(!(list_[i] < el)) {
				// Shouldn't be equal
				return i;
			}
		}
		return cur_;
	}

	/**
	 * Perform a binary search for the first element that is not less
	 * than 'el'.  Return cur_ if all elements are less than el.
	 */
	size_t bsearchLoBound(const T& el) const {
		size_t hi = cur_;
		size_t lo = 0;
		while(true) {
			if(lo == hi) {
#ifndef NDEBUG
				if((rand() % 10) == 0) {
					assert_eq(lo, scanLoBound(el));
				}
#endif
				return lo;
			}
			size_t mid = lo + ((hi-lo)>>1);
			assert_neq(mid, hi);
			if(list_[mid] < el) {
				if(lo == mid) {
#ifndef NDEBUG
					if((rand() % 10) == 0) {
						assert_eq(hi, scanLoBound(el));
					}
#endif
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
	}

	/**
	 * Return true if sorted, assert otherwise.
	 */
	bool sorted() const {
		if(cur_ <= 1) return true;
#ifndef NDEBUG
		if((rand() % 20) == 0) {
			for(size_t i = 0; i < cur_-1; i++) {
				assert(list_[i] < list_[i+1]);
			}
		}
#endif
		return true;
	}

	/**
	 * Insert value 'el' at offset 'idx'.  It's OK to insert at cur_,
	 * which is equivalent to appending.
	 */
	void insert(const T& el, size_t idx) {
		assert_leq(idx, cur_);
		if(cur_ == sz_) {
			expandCopy(sz_+1);
			assert(sorted());
		}
		for(size_t i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
		assert(sorted());
	}

	/**
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
		assert(sorted());
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) {
			newsz *= 2;
		}
		T* tmp = alloc(newsz);
		for(size_t i = 0; i < cur_; i++) {
			tmp[i] = list_[i];
		}
		free();
		list_ = tmp;
		sz_ = newsz;
	}

	int cat_;    // memory category, for accounting purposes
	T *list_;    // list pointer, returned from new[]
	size_t sz_;  // capacity
	size_t cur_; // occupancy (AKA size)
};

template <typename T, int S = 128>
class ELSet {

public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	explicit ELSet(int cat = 0) :
		cat_(cat), list_(NULL), sz_(S), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	explicit ELSet(size_t isz, int cat = 0) :
		cat_(cat), list_(NULL), sz_(isz), cur_(0)
	{
		assert_gt(isz, 0);
		assert_geq(cat, 0);
	}

	/**
	 * Copy from another ELList using operator=.
	 */
	ELSet(const ELSet<T, S>& o) :
		cat_(0), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
	}

	/**
	 * Copy from another ELList using operator=.
	 */
	explicit ELSet(const ELSet<T, S>& o, int cat) :
		cat_(cat), list_(NULL), sz_(0), cur_(0)
	{
		*this = o;
		assert_geq(cat, 0);
	}

	/**
	 * Destructor.
	 */
	~ELSet() { free(); }

	/**
	 * Make this object into a copy of o by allocating enough memory to
	 * fit the number of elements in o (note: the number of elements
	 * may be substantially less than the memory allocated in o) and
	 * using operator= to copy them over.
	 */
	ELSet<T, S>& operator=(const ELSet<T, S>& o) {
		assert_eq(cat_, o.cat());
		if(list_ == NULL) {
			lazyInit();
		}
		if(o.cur_ == 0) {
			cur_ = 0;
			return *this;
		}
		if(sz_ < o.cur_) expandNoCopy(o.cur_ + 1);
		assert_geq(sz_, o.cur_);
		cur_ = o.cur_;
		for(size_t i = 0; i < cur_; i++) {
			// Note: using operator=, not xfer
			assert_eq(list_[i].cat(), o.list_[i].cat());
			list_[i] = o.list_[i];
		}
		return *this;
	}
	
	/**
	 * Transfer the guts of another ESet into this one without using
	 * operator=, etc.  We have to set ESet o's list_ field to NULL to
	 * avoid o's destructor from deleting list_ out from under us.
	 */
	void xfer(ELSet<T, S>& o) {
		assert_eq(cat_, o.cat());
		list_ = o.list_; // list_ is an array of ESet<T>s
		sz_   = o.sz_;
		cur_  = o.cur_;
		o.list_ = NULL;
		o.sz_ = o.cur_ = 0;
	}

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }

	/**
	 * Return true iff list hasn't been initialized yet.
	 */
	inline bool null() const { return list_ == NULL; }

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void expand() {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
	}

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
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Get the element on the top of the stack.
	 */
	inline ESet<T>& back() {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the element on the top of the stack, const version.
	 */
	inline const ESet<T>& back() const {
		assert_gt(cur_, 0);
		return list_[cur_-1];
	}

	/**
	 * Get the frontmost element (bottom of stack).
	 */
	inline ESet<T>& front() {
		assert_gt(cur_, 0);
		return list_[0];
	}

	/**
	 * Get the element on the bottom of the stack, const version.
	 */
	inline const ESet<T>& front() const { return front(); }

	/**
	 * Return a reference to the ith element.
	 */
	inline ESet<T>& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const ESet<T>& operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline ESet<T>& get(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	inline const ESet<T>& get(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	ESet<T>& getSlow(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.  This version is not
	 * inlined, which guarantees we can use it from the debugger.
	 */
	const ESet<T>& getSlow(size_t i) const {
		return operator[](i);
	}
	
	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	ESet<T> *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	const ESet<T> *ptr() const { return list_; }

	/**
	 * Set the memory category for this object and all children.
	 */
	void setCat(int cat) {
		assert_gt(cat, 0);
		cat_ = cat;
		if(cat_ != 0) {
			for(size_t i = 0; i < sz_; i++) {
				assert(list_[i].null());
				list_[i].setCat(cat_);
			}
		}
	}

	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }

protected:

	/**
	 * Initialize memory for ELSet.
	 */
	void lazyInit() {
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	ESet<T> *alloc(size_t sz) {
		assert_gt(sz, 0);
		ESet<T> *tmp = new ESet<T>[sz];
		gMemTally.add(cat_, sz);
		if(cat_ != 0) {
			for(size_t i = 0; i < sz; i++) {
				assert(tmp[i].ptr() == NULL);
				tmp[i].setCat(cat_);
			}
		}
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] list_;
			gMemTally.del(cat_, sz_);
			list_ = NULL;
		}
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Copy old contents into new buffer
	 * using operator=.
	 */
	void expandCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ESet<T>* tmp = alloc(newsz);
		if(list_ != NULL) {
			for(size_t i = 0; i < cur_; i++) {
				assert_eq(cat_, tmp[i].cat());
				tmp[i].xfer(list_[i]);
				assert_eq(cat_, tmp[i].cat());
			}
			free();
		}
		list_ = tmp;
		sz_ = newsz;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.  Don't copy old contents over.
	 */
	void expandNoCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		free();
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ESet<T>* tmp = alloc(newsz);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	int cat_;    // memory category, for accounting purposes
	ESet<T> *list_; // list pointer, returned from new[]
	size_t sz_;  // capacity
	size_t cur_; // occupancy (AKA size)

};

/**
 * Expandable map using a heap-allocated sorted array.
 *
 * Note that the copy constructor and operator= routines perform
 * shallow copies (w/ memcpy).
 */
template <typename K, typename V>
class EMap {

public:

	/**
	 * Allocate initial default of 128 elements.
	 */
	EMap(int cat = 0) :
		cat_(cat),
		list_(NULL),
		sz_(128),
		cur_(0)
	{
		list_ = alloc(sz_);
	}

	/**
	 * Initially allocate given number of elements; should be > 0.
	 */
	EMap(size_t isz, int cat = 0) :
		cat_(cat),
		list_(NULL),
		sz_(isz),
		cur_(0)
	{
		assert_gt(isz, 0);
		list_ = alloc(sz_);
	}

	/**
	 * Copy from another ESet.
	 */
	EMap(const EMap<K, V>& o) : list_(NULL) {
		*this = o;
	}

	/**
	 * Destructor.
	 */
	~EMap() { free(); }

	/**
	 * Copy contents of given ESet into this ESet.
	 */
	EMap& operator=(const EMap<K, V>& o) {
		sz_ = o.sz_;
		cur_ = o.cur_;
		free();
		list_ = alloc(sz_);
		memcpy(list_, o.list_, cur_ * sizeof(std::pair<K, V>));
		return *this;
	}

	/**
	 * Return number of elements.
	 */
	size_t size() const { return cur_; }

	/**
	 * Return true iff there are no elements.
	 */
	bool empty() const { return cur_ == 0; }

	/**
	 * Insert a new element into the set in sorted order.
	 */
	bool insert(const std::pair<K, V>& el) {
		size_t i = 0;
		if(cur_ == 0) {
			insert(el, 0);
			return true;
		}
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el.first);
		} else {
			// Binary search
			i = bsearchLoBound(el.first);
		}
		if(list_[i] == el) return false;
		insert(el, i);
		return true;
	}

	/**
	 * Return true iff this set contains 'el'.
	 */
	bool contains(const K& el) const {
		if(cur_ == 0) return false;
		else if(cur_ == 1) return el == list_[0].first;
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i].first == el;
	}

	/**
	 * Remove element from set.
	 */
	void remove(const K& el) {
		size_t i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		assert(i != cur_ && list_[i].first == el);
		erase(i);
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz <= cur_) return;
		if(sz_ < sz) expandCopy(sz);
	}
	
	/**
	 * Get the ith key, value pair in the map.
	 */
	const std::pair<K, V>& get(size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}
	
	/**
	 * Get the ith key, value pair in the map.
	 */
	const std::pair<K, V>& operator[](size_t i) const {
		return get(i);
	}

	/**
	 * Clear set without deallocating (or setting) anything.
	 */
	void clear() { cur_ = 0; }

private:

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	std::pair<K, V> *alloc(size_t sz) {
		assert_gt(sz, 0);
		std::pair<K, V> *tmp = new std::pair<K, V>[sz];
		gMemTally.add(cat_, sz);
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] list_;
			gMemTally.del(cat_, sz_);
			list_ = NULL;
		}
	}

	/**
	 * Simple linear scan that returns the index of the first element
	 * of list_ that is not less than el, or cur_ if all elements are
	 * less than el.
	 */
	size_t scanLoBound(const K& el) const {
		for(size_t i = 0; i < cur_; i++) {
			if(!(list_[i].first < el)) {
				// Shouldn't be equal
				return i;
			}
		}
		return cur_;
	}

	/**
	 * Perform a binary search for the first element that is not less
	 * than 'el'.  Return cur_ if all elements are less than el.
	 */
	size_t bsearchLoBound(const K& el) const {
		size_t hi = cur_;
		size_t lo = 0;
		while(true) {
			if(lo == hi) {
#ifndef NDEBUG
				if((rand() % 10) == 0) {
					assert_eq(lo, scanLoBound(el));
				}
#endif
				return lo;
			}
			size_t mid = lo + ((hi-lo)>>1);
			assert_neq(mid, hi);
			if(list_[mid].first < el) {
				if(lo == mid) {
#ifndef NDEBUG
					if((rand() % 10) == 0) {
						assert_eq(hi, scanLoBound(el));
					}
#endif
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
	}

	/**
	 * Return true if sorted, assert otherwise.
	 */
	bool sorted() const {
		if(cur_ <= 1) return true;
		for(size_t i = 0; i < cur_-1; i++) {
			assert(list_[i] < list_[i+1]);
		}
		return true;
	}

	/**
	 * Insert value 'el' at offset 'idx'.  It's OK to insert at cur_,
	 * which is equivalent to appending.
	 */
	void insert(const std::pair<K, V>& el, size_t idx) {
		assert_leq(idx, cur_);
		if(cur_ == sz_) expandCopy(sz_+1);
		for(size_t i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
		assert(sorted());
	}

	/**
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert_lt(idx, cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
		assert(sorted());
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Expansions are quadratic.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = sz_ * 2;
		while(newsz < thresh) newsz *= 2;
		std::pair<K, V>* tmp = alloc(newsz);
		free();
		list_ = tmp;
		sz_ = newsz;
	}

	int cat_;    // memory category, for accounting purposes
	std::pair<K, V> *list_; // list pointer, returned from new[]
	size_t sz_;  // capacity
	size_t cur_; // occupancy (AKA size)
};

template <typename T, int S = 128>
class EFactory {

public:

	explicit EFactory(size_t isz, int cat = 0) : l_(isz, cat) { }
	
	explicit EFactory(int cat = 0) : l_(cat) { }
	
	void clear() {
		l_.clear();
	}
	
	size_t alloc() {
		l_.expand();
		return l_.size()-1;
	}
	
	size_t size() const {
		return l_.size();
	}

	bool empty() const {
		return size() == 0;
	}
	
	T& pop() {
		T& ret = l_[l_.size()-1];
		l_.resize(l_.size()-1);
		return ret;
	}
	
	T& operator[](size_t off) {
		return l_[off];
	}

	const T& operator[](size_t off) const {
		return l_[off];
	}

protected:

	EList<T, S> l_;
};

/**
 * Dispenses pages of memory for all the lists in the cache, including
 * the sequence-to-range map, the range list, the edits list, and the
 * offsets list.  All lists contend for the same pool of memory.
 */
class Pool {
public:
	Pool(
		uint64_t bytes,
		uint32_t pagesz,
		int cat = 0) :
		cat_(cat),
		cur_(0),
		bytes_(bytes),
		pagesz_(pagesz),
		pages_(cat)
	{
		for(size_t i = 0; i < ((bytes+pagesz-1)/pagesz); i++) {
			pages_.push_back(new uint8_t[pagesz]);
			gMemTally.add(cat, pagesz);
			assert(pages_.back() != NULL);
		}
		assert(repOk());
	}
	
	/**
	 * Free each page.
	 */
	~Pool() {
		for(size_t i = 0; i < pages_.size(); i++) {
			assert(pages_[i] != NULL);
			delete[] pages_[i];
			gMemTally.del(cat_, pagesz_);
		}
	}

	/**
	 * Allocate one page, or return NULL if no pages are left.
	 */
	uint8_t * alloc() {
		assert(repOk());
		if(cur_ == pages_.size()) return NULL;
		return pages_[cur_++];
	}

	/**
	 * Clear the pool so that no pages are considered allocated.
	 */
	void clear() {
		cur_ = 0;
		assert(repOk());
	}

	/**
	 * Reset the Pool to be as though
	 */
	void free() {
		// Currently a no-op because the only freeing method supported
		// now is to clear the entire pool
	}

	/**
	 * Check that pool is internally consistent.
	 */
	bool repOk() const {
		assert_leq(cur_, pages_.size());
		assert(!pages_.empty());
		assert_gt(bytes_, 0);
		assert_gt(pagesz_, 0);
		return true;
	}

private:
	int             cat_;    // memory category, for accounting purposes
	uint32_t        cur_;    // next page to hand out
	const uint64_t  bytes_;  // total bytes in the pool
	const uint32_t  pagesz_; // size of a single page
	EList<uint8_t*> pages_;  // the pages themselves
};

/**
 * An expandable list backed by a pool.
 */
template<typename T, int S>
class PList {

#define PLIST_PER_PAGE (S / sizeof(T))

public:
	/**
	 * Initialize the current-edit pointer to 0 and set the number of
	 * edits per memory page.
	 */
	PList(int cat = 0) :
		cur_(0),
		curPage_(0),
		pages_(cat) { }

	/**
	 * Add 1 object to the list.
	 */
	bool add(Pool& p, const T& o) {
		assert(repOk());
		if(!ensure(p, 1)) return false;
		if(cur_ == PLIST_PER_PAGE) {
			cur_ = 0;
			curPage_++;
		}
		assert_lt(curPage_, pages_.size());
		assert(repOk());
		assert_lt(cur_, PLIST_PER_PAGE);
		pages_[curPage_][cur_++] = o;
		return true;
	}

	/**
	 * Add a list of objects to the list.
	 */
	bool add(Pool& p, const EList<T>& os) {
		if(!ensure(p, os.size())) return false;
		for(size_t i = 0; i < os.size(); i++) {
			if(cur_ == PLIST_PER_PAGE) {
				cur_ = 0;
				curPage_++;
			}
			assert_lt(curPage_, pages_.size());
			assert(repOk());
			assert_lt(cur_, PLIST_PER_PAGE);
			pages_[curPage_][cur_++] = os[i];
		}
		return true;
	}

	/**
	 * Add a list of objects to the list.
	 */
	bool copy(
		Pool& p,
		const PList<T, S>& src,
		size_t i,
		size_t len)
	{
		if(!ensure(p, src.size())) return false;
		for(size_t i = 0; i < src.size(); i++) {
			if(cur_ == PLIST_PER_PAGE) {
				cur_ = 0;
				curPage_++;
			}
			assert_lt(curPage_, pages_.size());
			assert(repOk());
			assert_lt(cur_, PLIST_PER_PAGE);
			pages_[curPage_][cur_++] = src[i];
		}
		return true;
	}

	/**
	 * Add 'num' objects, all equal to 'o' to the list.
	 */
	bool addFill(Pool& p, size_t num, const T& o) {
		if(!ensure(p, num)) return false;
		for(size_t i = 0; i < num; i++) {
			if(cur_ == PLIST_PER_PAGE) {
				cur_ = 0;
				curPage_++;
			}
			assert_lt(curPage_, pages_.size());
			assert(repOk());
			assert_lt(cur_, PLIST_PER_PAGE);
			pages_[curPage_][cur_++] = o;
		}
		return true;
	}

	/**
	 * Free all pages associated with the list.
	 */
	void clear() {
		pages_.clear();
		cur_ = curPage_ = 0;
	}

	/**
	 * Check that list is internally consistent.
	 */
	bool repOk() const {
		assert(pages_.size() == 0 || curPage_ < pages_.size());
		assert_leq(cur_, PLIST_PER_PAGE);
		return true;
	}

	/**
	 * Return the number of elements in the list.
	 */
	size_t size() const {
		return curPage_ * PLIST_PER_PAGE + cur_;
	}
	
	/**
	 * Return true iff the PList has no elements.
	 */
	bool empty() const {
		return size() == 0;
	}

	/**
	 * Get the ith element added to the list.
	 */
	inline const T& getConst(size_t i) const {
		assert_lt(i, size());
		size_t page = i / PLIST_PER_PAGE;
		size_t elt = i % PLIST_PER_PAGE;
		return pages_[page][elt];
	}

	/**
	 * Get the ith element added to the list.
	 */
	inline T& get(size_t i) {
		assert_lt(i, size());
		size_t page = i / PLIST_PER_PAGE;
		size_t elt = i % PLIST_PER_PAGE;
		assert_lt(page, pages_.size());
		assert(page < pages_.size()-1 || elt < cur_);
		return pages_[page][elt];
	}
	
	/**
	 * Get the most recently added element.
	 */
	inline T& back() {
		size_t page = (size()-1) / PLIST_PER_PAGE;
		size_t elt = (size()-1) % PLIST_PER_PAGE;
		assert_lt(page, pages_.size());
		assert(page < pages_.size()-1 || elt < cur_);
		return pages_[page][elt];
	}
	
	/**
	 * Get const version of the most recently added element.
	 */
	inline const T& back() const {
		size_t page = (size()-1) / PLIST_PER_PAGE;
		size_t elt = (size()-1) % PLIST_PER_PAGE;
		assert_lt(page, pages_.size());
		assert(page < pages_.size()-1 || elt < cur_);
		return pages_[page][elt];
	}

	/**
	 * Get the element most recently added to the list.
	 */
	T& last() {
		assert(!pages_.empty());
		assert_gt(PLIST_PER_PAGE, 0);
		if(cur_ == 0) {
			assert_gt(pages_.size(), 1);
			return pages_[pages_.size()-2][PLIST_PER_PAGE-1];
		} else {
			return pages_.back()[cur_-1];
		}
	}

	/**
	 * Return true iff 'num' additional objects will fit in the pages
	 * allocated to the list.  If more pages are needed, they are
	 * added if possible.
	 */
	bool ensure(Pool& p, size_t num) {
		assert(repOk());
		if(num == 0) return true;
		// Allocation of the first page
		if(pages_.size() == 0) {
			if(expand(p) == NULL) {
				return false;
			}
			assert_eq(1, pages_.size());
		}
		size_t cur = cur_;
		size_t curPage = curPage_;
		while(cur + num > PLIST_PER_PAGE) {
			assert_lt(curPage, pages_.size());
			if(curPage == pages_.size()-1 && expand(p) == NULL) {
				return false;
			}
			num -= (PLIST_PER_PAGE - cur);
			cur = 0;
			curPage++;
		}
		return true;
	}

protected:

	/**
	 * Expand our page supply by 1
	 */
	T* expand(Pool& p) {
		T* newpage = (T*)p.alloc();
		if(newpage == NULL) {
			return NULL;
		}
		pages_.push_back(newpage);
		return pages_.back();
	}

	size_t       cur_;     // current elt within page
	size_t       curPage_; // current page
	EList<T*>    pages_;   // the pages
};

/**
 * An expandable list backed by a pool.
 */
template<typename T, int S>
class PListSlice {

public:
	PListSlice() :
		i_(0),
		len_(0),
		list_()
	{ }

	PListSlice(
		PList<T, S>& list,
		uint32_t i,
		uint32_t len) :
		i_(i),
		len_(len),
		list_(&list)
	{ }
	
	/**
	 * Initialize from a piece of another PListSlice.
	 */
	void init(const PListSlice<T, S>& sl, size_t first, size_t last) {
		assert_gt(last, first);
		assert_leq(last - first, sl.len_);
		i_ = (uint32_t)(sl.i_ + first);
		len_ = (uint32_t)(last - first);
		list_ = sl.list_;
	}
	
	/**
	 * Reset state to be empty.
	 */
	void reset() {
		i_ = len_ = 0;
		list_ = NULL;
	}
	
	/**
	 * Get the ith element of the slice.
	 */
	inline const T& get(size_t i) const {
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i+i_);
	}

	/**
	 * Get the ith element of the slice.
	 */
	inline T& get(size_t i) {
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i+i_);
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline T& operator[](size_t i) {
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i+i_);
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const T& operator[](size_t i) const {
		assert(valid());
		assert_lt(i, len_);
		return list_->get(i+i_);
	}

	/**
	 * Return true iff this slice is initialized.
	 */
	bool valid() const {
		return len_ != 0;
	}
	
	/**
	 * Return number of elements in the slice.
	 */
	size_t size() const {
		return len_;
	}
	
	/**
	 * Ensure that the PListSlice is internally consistent and
	 * consistent with the backing PList.
	 */
	bool repOk() const {
		assert_leq(i_ + len_, list_->size());
		return true;
	}
	
	/**
	 * Return true iff this slice refers to the same slice of the same
	 * list as the given slice.
	 */
	bool operator==(const PListSlice& sl) const {
		return i_ == sl.i_ && len_ == sl.len_ && list_ == sl.list_;
	}

	/**
	 * Return false iff this slice refers to the same slice of the same
	 * list as the given slice.
	 */
	bool operator!=(const PListSlice& sl) const {
		return !(*this == sl);
	}
	
	/**
	 * Set the length.  This could leave things inconsistent (e.g. could
	 * include elements that fall off the end of list_).
	 */
	void setLength(size_t nlen) {
		len_ = (uint32_t)nlen;
	}
	
protected:
	uint32_t i_;
	uint32_t len_;
	PList<T, S>* list_;
};

/**
 * A Red-Black tree node.  Links to parent & left and right children.
 * Key and Payload are of types K and P.  Node total ordering is based
 * on K's total ordering.  K must implement <, == and > operators.
 */
template<typename K, typename P> // K=key, P=payload
class RedBlackNode {

	typedef RedBlackNode<K,P> TNode;

public:
	TNode *parent;  // parent
	TNode *left;    // left child
	TNode *right;   // right child
	bool   red;     // true -> red, false -> black
	K      key;     // key, for ordering
	P      payload; // payload (i.e. value)

	/**
	 * Return the parent of this node's parent, or NULL if none exists.
	 */
	RedBlackNode *grandparent() {
		return parent != NULL ? parent->parent : NULL;
	}

	/**
	 * Return the sibling of this node's parent, or NULL if none exists.
	 */
	RedBlackNode *uncle() {
		if(parent == NULL) return NULL; // no parent
		if(parent->parent == NULL) return NULL; // parent has no siblings
		return (parent->parent->left == parent) ? parent->parent->right : parent->parent->left;
	}
	
	/**
	 * Return true iff this node is its parent's left child.
	 */
	bool isLeftChild() const { assert(parent != NULL); return parent->left == this; }

	/**
	 * Return true iff this node is its parent's right child.
	 */
	bool isRightChild() const { assert(parent != NULL); return parent->right == this; }

	/**
	 * Return true iff this node is its parent's right child.
	 */
	void replaceChild(RedBlackNode* ol, RedBlackNode* nw) {
		if(left == ol) {
			left = nw;
		} else {
			assert(right == ol);
			right = nw;
		}
	}

	/**
	 * Return the number of non-null children this node has.
	 */
	int numChildren() const {
		return ((left != NULL) ? 1 : 0) + ((right != NULL) ? 1 : 0);
	}
	
	/**
	 * Check that node is internally consistent.
	 */ 
	bool repOk() const {
		if(parent != NULL) {
			assert(parent->left == this || parent->right == this);
		}
		return true;
	}

	/**
	 * True -> my key is less than than the given node's key.
	 */
	bool operator<(const TNode& o) const { return key < o.key; }

	/**
	 * True -> my key is greater than the given node's key.
	 */
	bool operator>(const TNode& o) const { return key > o.key; }

	/**
	 * True -> my key equals the given node's key.
	 */
	bool operator==(const TNode& o) const { return key == o.key; }

	/**
	 * True -> my key is less than the given key.
	 */
	bool operator<(const K& okey) const { return key < okey; }

	/**
	 * True -> my key is greater than the given key.
	 */
	bool operator>(const K& okey) const { return key > okey; }

	/**
	 * True -> my key is equal to the given key.
	 */
	bool operator==(const K& okey) const { return key == okey; }
};

/**
 * A Red-Black tree that associates keys (of type K) with payloads (of
 * type P).  Red-Black trees are self-balancing and guarantee that the
 * tree as always "balanced" to a factor of 2, i.e., the longest
 * root-to-leaf path is never more than twice as long as the shortest
 * root-to-leaf path.
 */
template<typename K, typename P> // K=key, P=payload
class RedBlack {

	typedef RedBlackNode<K,P> TNode;

public:
	/**
	 * Initialize the current-edit pointer to 0 and set the number of
	 * edits per memory page.
	 */
	RedBlack(uint32_t pageSz, int cat = 0) :
		perPage_(pageSz/sizeof(TNode)), pages_(cat) { clear(); }

	/**
	 * Given a DNA string, find the red-black node corresponding to it,
	 * if one exists.
	 */
	inline TNode* lookup(const K& key) const {
		TNode* cur = root_;
		while(cur != NULL) {
			if((*cur) == key) return cur;
			if((*cur) < key) {
				cur = cur->right;
			} else {
				cur = cur->left;
			}
		}
		return NULL;
	}

	/**
	 * Add a new key as a node in the red-black tree.
	 */
	TNode* add(
		Pool& p,      // in: pool for memory pages
		const K& key, // in: key to insert
		bool* added)  // if true, assert is thrown if key exists
	{
		// Look for key; if it's not there, get its parent
		TNode* cur = root_;
		assert(root_ == NULL || !root_->red);
		TNode* parent = NULL;
		bool leftChild = true;
		while(cur != NULL) {
			if((*cur) == key) {
				// Found it; break out of loop with cur != NULL
				break;
			}
			parent = cur;
			if((*cur) < key) {
				if((cur = cur->right) == NULL) {
					// Fell off the bottom of the tree as the right
					// child of parent 'lastCur'
					leftChild = false;
				}
			} else {
				if((cur = cur->left) == NULL) {
					// Fell off the bottom of the tree as the left
					// child of parent 'lastCur'
					leftChild = true;
				}
			}
		}
		if(cur != NULL) {
			// Found an entry; assert if we weren't supposed to
			if(added != NULL) *added = false;
		} else {
			assert(root_ == NULL || !root_->red);
			if(!addNode(p, cur)) {
				// Exhausted memory
				return NULL;
			}
			assert(cur != NULL);
			assert(cur != root_);
			assert(cur != parent);
			// Initialize new node
			cur->key = key;
			cur->left = cur->right = NULL;
			cur->red = true; // red until proven black
			keys_++;
			if(added != NULL) *added = true;
			// Put it where we know it should go
			addNode(cur, parent, leftChild);
		}
		return cur; // return the added or found node
	}

	/**
	 * Check that list is internally consistent.
	 */
	bool repOk() const {
		assert(curPage_ == 0 || curPage_ < pages_.size());
		assert_leq(cur_, perPage_);
		assert(root_ == NULL || !root_->red);
		return true;
	}
	
	/**
	 * Clear all state.
	 */
	void clear() {
		cur_ = curPage_ = 0;
		root_ = NULL;
		keys_ = 0;
		intenseRepOkCnt_ = 0;
		pages_.clear();
	}
	
	/**
	 * Return number of keys added.
	 */
	size_t size() const {
		return keys_;
	}
	
	/**
	 * Return true iff there are no keys in the map.
	 */
	bool empty() const {
		return keys_ == 0;
	}

	/**
	 * Add another node and return a pointer to it in 'node'.  A new
	 * page is allocated if necessary.  If the allocation fails, false
	 * is returned.
	 */
	bool addNode(Pool& p, TNode*& node) {
		assert_leq(cur_, perPage_);
		assert(repOk());
		assert(this != NULL);
		// Allocation of the first page
		if(pages_.size() == 0) {
			if(addPage(p) == NULL) {
				node = NULL;
				return false;
			}
			assert_eq(1, pages_.size());
		}
		if(cur_ == perPage_) {
			assert_lt(curPage_, pages_.size());
			if(curPage_ == pages_.size()-1 && addPage(p) == NULL) {
				return false;
			}
			cur_ = 0;
			curPage_++;
		}
		assert_lt(cur_, perPage_);
		assert_lt(curPage_, pages_.size());
		node = &pages_[curPage_][cur_];
		assert(node != NULL);
		cur_++;
		return true;
	}

protected:

	/**
	 * Check specifically that the red-black invariants are satistfied.
	 */
	bool redBlackRepOk(TNode* n) {
		if(n == NULL) return true;
		if(++intenseRepOkCnt_ < 500) return true;
		intenseRepOkCnt_ = 0;
		int minNodes = -1; // min # nodes along any n->leaf path
		int maxNodes = -1; // max # nodes along any n->leaf path
		// The number of black nodes along paths from n to leaf
		// (must be same for all paths)
		int blackConst = -1;
		size_t nodesTot = 0;
		redBlackRepOk(
			n,
			1, /* 1 node so far */
			n->red ? 0 : 1, /* black nodes so far */
			blackConst,
			minNodes,
			maxNodes,
			nodesTot);
		if(n == root_) {
			assert_eq(nodesTot, keys_);
		}
		assert_gt(minNodes, 0);
		assert_gt(maxNodes, 0);
		assert_leq(maxNodes, 2*minNodes);
		return true;
	}

	/**
	 * Check specifically that the red-black invariants are satistfied.
	 */
	bool redBlackRepOk(
		TNode* n,
		int nodes,
		int black,
		int& blackConst,
		int& minNodes,
		int& maxNodes,
		size_t& nodesTot) const
	{
		assert_gt(black, 0);
		nodesTot++; // account for leaf node
		if(n->left == NULL) {
			if(blackConst == -1) blackConst = black;
			assert_eq(black, blackConst);
			if(nodes+1 > maxNodes) maxNodes = nodes+1;
			if(nodes+1 < minNodes || minNodes == -1) minNodes = nodes+1;
		} else {
			if(n->red) assert(!n->left->red); // Red can't be child of a red
			redBlackRepOk(
				n->left,                         // next node
				nodes + 1,                       // # nodes so far on path
				black + (n->left->red ? 0 : 1),  // # black so far on path
				blackConst,                      // invariant # black nodes on root->leaf path
				minNodes,                        // min root->leaf len so far         
				maxNodes,                        // max root->leaf len so far
				nodesTot);                       // tot nodes so far
		}
		if(n->right == NULL) {
			if(blackConst == -1) blackConst = black;
			assert_eq(black, blackConst);
			if(nodes+1 > maxNodes) maxNodes = nodes+1;
			if(nodes+1 < minNodes || minNodes == -1) minNodes = nodes+1;
		} else {
			if(n->red) assert(!n->right->red); // Red can't be child of a red
			redBlackRepOk(
				n->right,                        // next node
				nodes + 1,                       // # nodes so far on path
				black + (n->right->red ? 0 : 1), // # black so far on path
				blackConst,                      // invariant # black nodes on root->leaf path
				minNodes,                        // min root->leaf len so far         
				maxNodes,                        // max root->leaf len so far
				nodesTot);                       // tot nodes so far
		}
		return true;
	}

	/**
	 * Rotate to the left such that n is replaced by its right child
	 * w/r/t n's current parent.
	 */
	void leftRotate(TNode* n) {
		TNode* r = n->right;
		assert(n->repOk());
		assert(r->repOk());
		n->right = r->left;
		if(n->right != NULL) {
			n->right->parent = n;
			assert(n->right->repOk());
		}
		r->parent = n->parent;
		n->parent = r;
		r->left = n;
		if(r->parent != NULL) {
			r->parent->replaceChild(n, r);
		}
		if(root_ == n) root_ = r;
		assert(!root_->red);
		assert(n->repOk());
		assert(r->repOk());
	}

	/**
	 * Rotate to the right such that n is replaced by its left child
	 * w/r/t n's current parent.  n moves down to the right and loses
	 * its left child, while its former left child moves up and gains a
	 * right child.
	 */
	void rightRotate(TNode* n) {
		TNode* r = n->left;
		assert(n->repOk());
		assert(r->repOk());
		n->left = r->right;
		if(n->left != NULL) {
			n->left->parent = n;
			assert(n->left->repOk());
		}
		r->parent = n->parent;
		n->parent = r;
		r->right = n;
		if(r->parent != NULL) {
			r->parent->replaceChild(n, r);
		}
		if(root_ == n) root_ = r;
		assert(!root_->red);
		assert(n->repOk());
		assert(r->repOk());
	}

	/**
	 * Add a node to the red-black tree, maintaining the red-black
	 * invariants.
	 */
	void addNode(TNode* n, TNode* parent, bool leftChild) {
		assert(n != NULL);
		if(parent == NULL) {
			// Case 1: inserted at root
			root_ = n;
			root_->red = false; // root must be black
			n->parent = NULL;
			assert(redBlackRepOk(root_));
			assert(n->repOk());
		} else {
			assert(!root_->red);
			// Add new node to tree
			if(leftChild) {
				assert(parent->left == NULL);
				parent->left = n;
			} else {
				assert(parent->right == NULL);
				parent->right = n;
			}
			n->parent = parent;
			int thru = 0;
			while(true) {
				thru++;
				parent = n->parent;
				if(parent != NULL) assert(parent->repOk());
				if(parent == NULL && n->red) {
					n->red = false;
				}
				if(parent == NULL || !parent->red) {
					assert(redBlackRepOk(root_));
					break;
				}
				TNode* uncle = n->uncle();
				TNode* gparent = n->grandparent();
				assert(gparent != NULL); // if parent is red, grandparent must exist
				bool uncleRed = (uncle != NULL ? uncle->red : false);
				if(uncleRed) {
					// Parent is red, uncle is red; recursive case
					assert(uncle != NULL);
					parent->red = uncle->red = false;
					gparent->red = true;
					n = gparent;
					continue;
				} else {
					if(parent->isLeftChild()) {
						// Parent is red, uncle is black, parent is
						// left child
						if(!n->isLeftChild()) {
							n = parent;
							leftRotate(n);
						}
						n = n->parent;
						n->red = false;
						n->parent->red = true;
						rightRotate(n->parent);
						assert(redBlackRepOk(n));
						assert(redBlackRepOk(root_));
					} else {
						// Parent is red, uncle is black, parent is
						// right child.
						if(!n->isRightChild()) {
							n = parent;
							rightRotate(n);
						}
						n = n->parent;
						n->red = false;
						n->parent->red = true;
						leftRotate(n->parent);
						assert(redBlackRepOk(n));
						assert(redBlackRepOk(root_));
					}
				}
				break;
			}
		}
		assert(redBlackRepOk(root_));
	}

	/**
	 * Expand our page supply by 1
	 */
	TNode* addPage(Pool& p) {
		TNode *n = (TNode *)p.alloc();
		if(n != NULL) {
			pages_.push_back(n);
		}
		return n;
	}

	size_t        keys_;    // number of keys so far
	size_t        cur_;     // current elt within page
	size_t        curPage_; // current page
	const size_t  perPage_; // # edits fitting in a page
	TNode*        root_;    // root node
	EList<TNode*> pages_;   // the pages
	int intenseRepOkCnt_;   // counter for the computationally intensive repOk function
};

/**
 * For assembling doubly-linked lists of Edits.
 */
template <typename T>
struct DoublyLinkedList {
	
	DoublyLinkedList() : payload(), prev(NULL), next(NULL) { }
	
	/**
	 * Add all elements in the doubly-linked list to the provided EList.
	 */
	void toList(EList<T>& l) {
		// Add this and all subsequent elements
		DoublyLinkedList<T> *cur = this;
		while(cur != NULL) {
			l.push_back(cur->payload);
			cur = cur->next;
		}
		// Add all previous elements
		cur = prev;
		while(cur != NULL) {
			l.push_back(cur->payload);
			cur = cur->prev;
		}
	}
	
	T                    payload;
	DoublyLinkedList<T> *prev;
	DoublyLinkedList<T> *next;
};

template <typename T1, typename T2, typename T3>
struct Triple {
	T1 first;
	T2 second;
	T3 third;
};

#endif /* DS_H_ */
