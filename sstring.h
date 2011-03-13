/*
 * sstring.h
 *
 *  Created on: Jan 16, 2010
 *      Author: Ben Langmead
 */

#ifndef SSTRING_H_
#define SSTRING_H_

#include <string.h>
#include <iostream>
#include <algorithm>
#include "assert_helpers.h"
#include "alphabet.h"
#include "random_source.h"

const static int BTString_len = 1024;

/**
 * Simple string class.
 */
template<typename T, int S = 1024, int M = 2>
class SStringExpandable {

public:

	SStringExpandable() :
		cs_(NULL),
		printcs_(NULL),
		len_(0),
		sz_(0)
	{ }

	SStringExpandable(size_t sz) :
		cs_(NULL),
		printcs_(NULL),
		len_(0),
		sz_(0)
	{
		expandNoCopy(sz);
	}

	/**
	 * Create an SStringExpandable from another SStringExpandable.
	 */
	SStringExpandable(const SStringExpandable<T, S>& o) :
		cs_(NULL),
		printcs_(NULL),
		len_(0),
		sz_(0)
	{
		*this = o;
	}

	/**
	 * Create an SStringExpandable from a std::basic_string of the
	 * appropriate type.
	 */
	SStringExpandable(const std::basic_string<T>& str) :
		cs_(NULL),
		printcs_(NULL),
		len_(0),
		sz_(0)
	{
		install(str.c_str(), str.length());
	}

	/**
	 * Create an SStringExpandable from an array and size.
	 */
	SStringExpandable(const T* b, size_t sz) :
		cs_(NULL),
		printcs_(NULL),
		len_(0),
		sz_(0)
	{
		install(b, sz);
	}

	/**
	 * Create an SStringExpandable from a zero-terminated array.
	 */
	SStringExpandable(const T* b) :
		cs_(NULL),
		printcs_(NULL),
		len_(0),
		sz_(0)
	{
		install(b, strlen(b));
	}

	/**
	 * Destroy the expandable string object.
	 */
	virtual ~SStringExpandable() {
		if(cs_ != NULL) {
			delete[] cs_;
			cs_ = NULL;
		}
		if(printcs_ != NULL) {
			delete[] printcs_;
			printcs_ = NULL;
		}
		sz_ = len_ = 0;
	}

	/**
	 * Return true iff this string is lexicographically less than the
	 * given string.
	 */
	bool operator< (const SStringExpandable<T, S, M>& b) const {
		for(size_t i = 0; i < std::min(len_, b.len_); i++) {
			if((int)cs_[i] < (int)b.cs_[i]) {
				return true;
			} else if((int)cs_[i] > (int)b.cs_[i]) {
				return false;
			}
		}
		// All are equal up to min(len, b.len)
		return len_ < b.len_;
	}

	/**
	 * Return true iff this string is lexicographically greater than
	 * the given string.
	 */
	bool operator> (const SStringExpandable<T, S, M>& b) const {
		for(size_t i = 0; i < std::min(len_, b.len_); i++) {
			if((int)cs_[i] > (int)b.cs_[i]) {
				return true;
			} else if((int)cs_[i] < (int)b.cs_[i]) {
				return false;
			}
		}
		// All are equal up to min(len, b.len)
		return len_ > b.len_;
	}

	/**
	 * Return ith character from the left of either the forward or the
	 * reverse-complement version of the read.
	 */
	T windowGet(
		size_t i,
		bool   fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		return fw ? cs_[depth+i] : cs_[depth+len-i-1];
	}

	/**
	 * Return ith character from the left of either the forward or the
	 * reverse-complement version of the read.
	 */
	void windowGet(
		T& ret,
		bool   fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = len_;
		assert_leq(len, len_ - depth);
		for(size_t i = 0; i < len; i++) {
			ret.append(fw ? cs_[depth+i] : cs_[depth+len-i-1]);
		}
	}

	/**
	 * Assignment to other SStringFixed.
	 */
	SStringExpandable<T,S>& operator=(const SStringExpandable<T,S>& o) {
		install(o.cs_, o.len_);
		return *this;
	}

	/**
	 * Return true iff all corresponding elements between this string
	 * and the given string are equal (i.e. operator== returns true).
	 * Lengths must also be equal.
	 */
	bool operator==(const SStringExpandable<T, S, M>& o) const {
		if(len_ != o.len_) return false;
		for(size_t i = 0; i < len_; i++) {
			if(cs_[i] != o.cs_[i]) return false;
		}
		return true;
	}

	/**
	 * Return true iff all corresponding elements between this string
	 * and the given string are not equal (i.e. operator== returns
	 * false).  Lengths must also be equal.
	 */
	bool operator!=(const SStringExpandable<T, S, M>& o) const {
		return !(*this == o);
	}

	/**
	 * Insert char c before position 'idx'; slide subsequent chars down.
	 */
	void insert(const T& c, size_t idx) {
		assert_lt(idx, len_);
		if(sz_ < len_ + 1) expandCopy((len_ + 1 + S) * M);
		// Move everyone down by 1
		for(int i = len_; i > idx; i--) {
			cs_[i] = cs_[i-1];
		}
		cs_[idx] = c;
		len_++;
	}

	/**
	 * Set character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, len_);
		cs_[idx] = c;
	}

	/**
	 * Append char c.
	 */
	void append(const T& c) {
		if(sz_ < len_ + 1) expandCopy((len_ + 1 + S) * M);
		cs_[len_++] = c;
	}

	/**
	 * Delete char at position 'idx'; slide subsequent chars up.
	 */
	void remove(size_t idx) {
		assert_lt(idx, len_);
		assert_gt(len_, 0);
		for(size_t i = idx; i < len_-1; i++) {
			cs_[i] = cs_[i+1];
		}
		len_--;
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const T& operator[](size_t i) const {
		assert_lt(i, len_);
		return cs_[i];
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const T& get(size_t i) const {
		assert_lt(i, len_);
		return cs_[i];
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const T* b, size_t sz) {
		if(sz_ < sz) expandNoCopy((sz + S) * M);
		memcpy(cs_, b, sz * sizeof(T));
		len_ = sz;
	}


	/**
	 * Copy all bytes from zero-terminated buffer 'b' into this string.
	 */
	void install(const T* b) { install(b, strlen(b)); }

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const char* b, size_t sz) {
		if(sz_ < sz) expandNoCopy((sz + S) * M);
		for(size_t i = 0; i < sz; i++) {
			cs_[i] = b[sz-i-1];
		}
		len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const SStringExpandable<T, S>& b) {
		if(sz_ < b.len_) expandNoCopy((b.len_ + S) * M);
		for(size_t i = 0; i < b.len_; i++) {
			cs_[i] = b.cs_[b.len_ - i - 1];
		}
		len_ = b.len_;
	}

	/**
	 * Reverse the buffer in place.
	 */
	void reverse() {
		for(size_t i = 0; i < (len_ >> 1); i++) {
			T tmp = cs_[i];
			cs_[i] = cs_[len_-i-1];
			cs_[len_-i-1] = tmp;
		}
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, the newly-added elements will contain garbage and should
	 * be initialized immediately.
	 */
	void resize(size_t len) {
		if(sz_ < len) expandCopy((len + S) * M);
		len_ = len;
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, new elements will be initialized with 'el'.
	 */
	void resize(size_t len, const T& el) {
		if(sz_ < len) expandCopy((len + S) * M);
		if(len > len_) {
			for(size_t i = len_; i < len; i++) {
				cs_[i] = el;
			}
		}
		len_ = len;
	}

	/**
	 * Set the first len elements of the buffer to el.
	 */
	void fill(size_t len, const T& el) {
		assert_leq(len, len_);
		for(size_t i = 0; i < len; i++) {
			cs_[i] = el;
		}
	}

	/**
	 * Set all elements of the buffer to el.
	 */
	void fill(const T& el) {
		fill(len_, el);
	}

	/**
	 * Trim len characters from the beginning of the string.
	 */
	void trimBegin(size_t len) {
		assert_leq(len, len_);
		if(len == len_) {
			len_ = 0; return;
		}
		for(size_t i = 0; i < len_-len; i++) {
			cs_[i] = cs_[i+len];
		}
		len_ -= len;
	}

	/**
	 * Trim len characters from the end of the string.
	 */
	void trimEnd(size_t len) {
		if(len >= len_) len_ = 0;
		else len_ -= len;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	void append(const T* b, size_t sz) {
		if(sz_ < len_ + sz) expandCopy((len_ + sz + S) * M);
		memcpy(cs_ + len_, b, sz * sizeof(T));
		len_ += sz;
	}

	/**
	 * Copy bytes from zero-terminated buffer 'b' into this string.
	 */
	void append(const T* b) {
		append(b, strlen(b));
	}

	/**
	 * Return the length of the string.
	 */
	size_t length() const { return len_; }

	/**
	 * Clear the buffer.
	 */
	void clear() { len_ = 0; }

	/**
	 * Return true iff the buffer is empty.
	 */
	bool empty() const { return len_ == 0; }

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	const T* toZBufXForm(const char *xform) const {
		ASSERT_ONLY(size_t xformElts = strlen(xform));
		T* printcs = const_cast<T*>(printcs_);
		for(size_t i = 0; i < len_; i++) {
			assert_lt(cs_[i], (int)xformElts);
			printcs[i] = xform[cs_[i]];
		}
		printcs[len_] = 0;
		return printcs_;
	}

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const T* toZBuf() const {
		assert_lt(len_, sz_);
		const_cast<T*>(cs_)[len_] = 0;
		return cs_;
	}

	/**
	 * Return true iff this DNA string matches the given nucleotide
	 * character string.
	 */
	bool eq(const char *str) const {
		const char *self = toZBuf();
		return strcmp(str, self) == 0;
	}

	/**
	 * Return a const version of the raw buffer.
	 */
	const T* buf() const { return cs_; }

	/**
	 * Return a writeable version of the raw buffer.
	 */
	T* wbuf() { return cs_; }

protected:
	/**
	 * Allocate new, bigger buffer and copy old contents into it.  If
	 * requested size can be accommodated by current buffer, do nothing.
	 */
	void expandCopy(size_t sz) {
		if(sz_ >= sz) return; // done!
		T *tmp  = new T[sz + 1];
		T *ptmp = new T[sz + 1];
		if(cs_ != NULL) {
			memcpy(tmp, cs_, sizeof(T)*len_);
			delete[] cs_;
		}
		if(printcs_ != NULL) {
			memcpy(ptmp, printcs_, sizeof(T)*len_);
			delete[] printcs_;
		}
		cs_ = tmp;
		printcs_ = ptmp;
		sz_ = sz;
	}

	/**
	 * Allocate new, bigger buffer.  If requested size can be
	 * accommodated by current buffer, do nothing.
	 */
	void expandNoCopy(size_t sz) {
		if(sz_ >= sz) return; // done!
		if(cs_      != NULL) delete[] cs_;
		if(printcs_ != NULL) delete[] printcs_;
		cs_ = new T[sz + 1];
		printcs_ = new T[sz + 1];
		sz_ = sz;
	}

	T *cs_;      // +1 so that we have the option of dropping in a terminating "\0"
	T *printcs_; // +1 so that we have the option of dropping in a terminating "\0"
	size_t len_; // # filled-in elements
	size_t sz_;  // size capacity of cs_
};

/**
 * Simple string class with in-object storage.
 *
 * All copies induced by, e.g., operator=, the copy constructor,
 * install() and append(), are shallow (using memcpy/sizeof).  If deep
 * copies are needed, use a different class.
 *
 * Reading from an uninitialized element results in an assert as long
 * as NDEBUG is not defined.  If NDEBUG is defined, the result is
 * undefined.
 */
template<typename T, int S = BTString_len>
class SStringFixed {
public:
	SStringFixed() : len_(0) { }

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SStringFixed(const SStringFixed<T, S>& o) {
		*this = o;
	}

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SStringFixed(const std::basic_string<T>& str) {
		install(str.c_str(), str.length());
	}

	/**
	 * Create an SStringFixed from an array and size.
	 */
	SStringFixed(const T* b, size_t sz) {
		install(b, sz);
	}

	/**
	 * Create an SStringFixed from a zero-terminated string.
	 */
	SStringFixed(const T* b) {
		install(b, strlen(b));
	}

	virtual ~SStringFixed() { } // C++ needs this

	/**
	 * Retrieve constant version of element i.
	 */
	const T& operator[](size_t i) const {
		return get(i);
	}

	/**
	 * Return true iff this string is lexicographically less than the
	 * given string.
	 */
	bool operator< (const SStringFixed<T, S>& b) const {
		for(size_t i = 0; i < std::min(len_, b.len_); i++) {
			if((int)cs_[i] < (int)b.cs_[i]) {
				return true;
			} else if((int)cs_[i] > (int)b.cs_[i]) {
				return false;
			}
		}
		// All are equal up to min(len, b.len)
		return len_ < b.len_;
	}

	/**
	 * Return true iff this string is lexicographically greater than
	 * the given string.
	 */
	bool operator> (const SStringFixed<T, S>& b) const {
		for(size_t i = 0; i < std::min(len_, b.len_); i++) {
			if((int)cs_[i] > (int)b.cs_[i]) {
				return true;
			} else if((int)cs_[i] < (int)b.cs_[i]) {
				return false;
			}
		}
		// All are equal up to min(len, b.len)
		return len_ > b.len_;
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const T& get(size_t i) const {
		assert_lt(i, len_);
		return cs_[i];
	}

	/**
	 * Return ith character from the left of either the forward or the
	 * reverse-complement version of the read.
	 */
	T windowGet(
		size_t i,
		bool   fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = len_;
		assert_lt(i, len);
		assert_leq(len, len_ - depth);
		return fw ? cs_[depth+i] : cs_[depth+len-i-1];
	}

	/**
	 * Return ith character from the left of either the forward or the
	 * reverse-complement version of the read.
	 */
	void windowGet(
		T& ret,
		bool   fw,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = len_;
		assert_leq(len, len_ - depth);
		for(size_t i = 0; i < len; i++) {
			ret.append(fw ? cs_[depth+i] : cs_[depth+len-i-1]);
		}
	}

	/**
	 * Assignment to other SStringFixed.
	 */
	SStringFixed<T,S>& operator=(const SStringFixed<T,S>& o) {
		install(o.cs_, o.len_);
		return *this;
	}

	/**
	 * Return true iff all corresponding elements between this string
	 * and the given string are equal (i.e. operator== returns true).
	 * Lengths must also be equal.
	 */
	bool operator==(const SStringFixed<T,S>& o) const {
		if(len_ != o.len_) return false;
		for(size_t i = 0; i < len_; i++) {
			if(cs_[i] != o.cs_[i]) return false;
		}
		return true;
	}

	/**
	 * Return the inverse of operator==.
	 */
	bool operator!=(const SStringFixed<T,S>& o) const {
		return !operator==(o);
	}

	/**
	 * Insert char c before position 'idx'; slide subsequent chars down.
	 */
	void insert(const T& c, size_t idx) {
		assert_lt(len_, S);
		assert_lt(idx, len_);
		// Move everyone down by 1
		for(int i = len_; i > idx; i--) {
			cs_[i] = cs_[i-1];
		}
		cs_[idx] = c;
		len_++;
	}

	/**
	 * Set character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, len_);
		cs_[idx] = c;
	}

	/**
	 * Append char c.
	 */
	void append(const T& c) {
		assert_lt(len_, S);
		cs_[len_++] = c;
	}

	/**
	 * Delete char at position 'idx'; slide subsequent chars up.
	 */
	void remove(size_t idx) {
		assert_lt(idx, len_);
		assert_gt(len_, 0);
		for(size_t i = idx; i < len_-1; i++) {
			cs_[i] = cs_[i+1];
		}
		len_--;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const T* b, size_t sz) {
		assert_leq(sz, S);
		memcpy(cs_, b, sz * sizeof(T));
		len_ = sz;
	}

	/**
	 * Copy all bytes from zero-terminated buffer 'b' into this string.
	 */
	void install(const T* b) { install(b, strlen(b)); }

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			cs_[i] = b[sz-i-1];
		}
		len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reversing them
	 * in the process.
	 */
	void installReverse(const SStringFixed<T, S>& b) {
		assert_leq(b.len_, S);
		for(size_t i = 0; i < b.len_; i++) {
			cs_[i] = b.cs_[b.len_ - i - 1];
		}
		len_ = b.len_;
	}

	/**
	 * Reverse the buffer in place.
	 */
	void reverse() {
		for(size_t i = 0; i < (len_ >> 1); i++) {
			T tmp = cs_[i];
			cs_[i] = cs_[len_-i-1];
			cs_[len_-i-1] = tmp;
		}
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, the newly-added elements will contain garbage and should
	 * be initialized immediately.
	 */
	void resize(size_t len) {
		assert_lt(len, S);
		len_ = len;
	}

	/**
	 * Simply resize the buffer.  If the buffer is resized to be
	 * longer, new elements will be initialized with 'el'.
	 */
	void resize(size_t len, const T& el) {
		assert_lt(len, S);
		if(len > len_) {
			for(size_t i = len_; i < len; i++) {
				cs_[i] = el;
			}
		}
		len_ = len;
	}

	/**
	 * Set the first len elements of the buffer to el.
	 */
	void fill(size_t len, const T& el) {
		assert_leq(len, len_);
		for(size_t i = 0; i < len; i++) {
			cs_[i] = el;
		}
	}

	/**
	 * Set all elements of the buffer to el.
	 */
	void fill(const T& el) {
		fill(len_, el);
	}

	/**
	 * Trim len characters from the beginning of the string.
	 */
	void trimBegin(size_t len) {
		assert_leq(len, len_);
		if(len == len_) {
			len_ = 0; return;
		}
		for(size_t i = 0; i < len_-len; i++) {
			cs_[i] = cs_[i+len];
		}
		len_ -= len;
	}

	/**
	 * Trim len characters from the end of the string.
	 */
	void trimEnd(size_t len) {
		if(len >= len_) len_ = 0;
		else len_ -= len;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	void append(const T* b, size_t sz) {
		assert_leq(sz + len_, S);
		memcpy(cs_ + len_, b, sz * sizeof(T));
		len_ += sz;
	}

	/**
	 * Copy bytes from zero-terminated buffer 'b' into this string.
	 */
	void append(const T* b) {
		append(b, strlen(b));
	}

	/**
	 * Return the length of the string.
	 */
	size_t length() const { return len_; }

	/**
	 * Clear the buffer.
	 */
	void clear() { len_ = 0; }

	/**
	 * Return true iff the buffer is empty.
	 */
	bool empty() const { return len_ == 0; }

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const T* toZBuf() const {
		const_cast<T*>(cs_)[len_] = 0;
		return cs_;
	}

	/**
	 * Return true iff this DNA string matches the given nucleotide
	 * character string.
	 */
	bool eq(const char *str) const {
		const char *self = toZBuf();
		return strcmp(str, self) == 0;
	}
	
	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	const T* toZBufXForm(const char *xform) const {
		ASSERT_ONLY(size_t xformElts = strlen(xform));
		T* printcs = const_cast<T*>(printcs_);
		for(size_t i = 0; i < len_; i++) {
			assert_lt(cs_[i], (int)xformElts);
			printcs[i] = xform[cs_[i]];
		}
		printcs[len_] = 0;
		return printcs_;
	}

	/**
	 * Return a const version of the raw buffer.
	 */
	const T* buf() const { return cs_; }

	/**
	 * Return a writeable version of the raw buffer.
	 */
	T* wbuf() { return cs_; }

protected:
	T cs_[S+1]; // +1 so that we have the option of dropping in a terminating "\0"
	T printcs_[S+1]; // +1 so that we have the option of dropping in a terminating "\0"
	size_t len_;
};

//
// Stream put operators
//

template <typename T, int S, int M>
std::ostream& operator<< (std::ostream& os, const SStringExpandable<T, S, M>& str) {
	os << str.toZBuf();
	return os;
}

template <typename T, int S>
std::ostream& operator<< (std::ostream& os, const SStringFixed<T, S>& str) {
	os << str.toZBuf();
	return os;
}

extern uint8_t asc2dna[];
extern uint8_t asc2col[];

/**
 * Encapsulates a fixed-length DNA string with characters encoded as
 * chars.  Only capable of encoding A, C, G, T and N.  The length is
 * specified via the template parameter S.
 */
template<int S = BTString_len>
class SDnaStringFixed : public SStringFixed<char, S> {
public:

	SDnaStringFixed() : SStringFixed<char, S>() { }

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SDnaStringFixed(const SDnaStringFixed<S>& o) :
		SStringFixed<char, S>(o) { }

	/**
	 * Create an SStringFixed from a C++ basic_string.
	 */
	SDnaStringFixed(const std::basic_string<char>& str) :
		SStringFixed<char, S>(str) { }

	/**
	 * Create an SStringFixed from an array and size.
	 */
	SDnaStringFixed(const char* b, size_t sz) :
		SStringFixed<char, S>(b, sz) { }

	/**
	 * Create an SStringFixed from a zero-terminated string.
	 */
	SDnaStringFixed(const char* b, bool chars = false, bool colors = false) :
		SStringFixed<char, S>()
	{
		if(chars) {
			if(colors) {
				installColors(b, strlen(b));
			} else {
				installChars(b, strlen(b));
			}
		} else {
			install(b, strlen(b));
		}
	}

	virtual ~SDnaStringFixed() { } // C++ needs this

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			this->cs_[i] = (b[sz-i-1] == 4 ? 4 : b[sz-i-1] ^ 3);
		}
		this->len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const SDnaStringFixed<S>& b) {
		assert_leq(b.len_, S);
		for(size_t i = 0; i < b.len_; i++) {
			this->cs_[i] = (b.cs_[b.len_-i-1] == 4 ? 4 : b.cs_[b.len_-i-1] ^ 3);
		}
		this->len_ = b.len_;
	}

	/**
	 * Either reverse or reverse-complement (depending on "color") this
	 * DNA buffer in-place.
	 */
	void reverseComp(bool color = false) {
		if(color) {
			this->reverse();
		} else {
			for(size_t i = 0; i < (this->len_ >> 1); i++) {
				char tmp1 = (this->cs_[i] == 4 ? 4 : this->cs_[i] ^ 3);
				char tmp2 = (this->cs_[this->len_-i-1] == 4 ? 4 : this->cs_[this->len_-i-1] ^ 3);
				this->cs_[i] = tmp2;
				this->cs_[this->len_-i-1] = tmp1;
			}
			// Do middle element iff there are an odd number
			if((this->len_ & 1) != 0) {
				char tmp = this->cs_[this->len_ >> 1];
				tmp = (tmp == 4 ? 4 : tmp ^ 3);
				this->cs_[this->len_ >> 1] = tmp;
			}
		}
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const char* b, size_t sz) {
		assert_leq(sz, S);
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for(size_t i = 0; i < sz; i++) {
			assert_leq(this->cs_[i], 4);
			assert_geq(this->cs_[i], 0);
		}
#endif
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			assert_in(toupper(b[i]), "ACGTN-");
			this->cs_[i] = asc2dna[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII color characters into normal DNA
	 * characters.
	 */
	virtual void installColors(const char* b, size_t sz) {
		assert_leq(sz, S);
		for(size_t i = 0; i < sz; i++) {
			assert_in(b[i], "0123.");
			this->cs_[i] = asc2col[(int)b[i]];
			assert_geq(this->cs_[i], 0);
			assert_leq(this->cs_[i], 4);
		}
		this->len_ = sz;
	}

	/**
	 * Copy C++ string of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const std::basic_string<char>& str) {
		installChars(str.c_str(), str.length());
	}

	/**
	 * Copy C++ string of ASCII color characters into normal DNA
	 * characters.
	 */
	virtual void installColors(const std::basic_string<char>& str) {
		installColors(str.c_str(), str.length());
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_leq(c, 4);
		assert_geq(c, 0);
		this->cs_[idx] = c;
	}

	/**
	 * Append DNA char c.
	 */
	void append(const char& c) {
		assert_lt(this->len_, S);
		assert_leq(c, 4);
		assert_geq(c, 0);
		this->cs_[this->len_++] = c;
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void setChar(char c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_in(toupper(c), "ACGTN");
		this->cs_[idx] = asc2dna[(int)c];
	}

	/**
	 * Append DNA character.
	 */
	void appendChar(char c) {
		assert_lt(this->len_, S);
		assert_in(toupper(c), "ACGTN");
		this->cs_[this->len_++] = asc2dna[(int)c];
	}

	/**
	 * Return DNA character corresponding to element 'idx'.
	 */
	char toChar(size_t idx) const {
		assert_geq((int)this->cs_[idx], 0);
		assert_leq((int)this->cs_[idx], 4);
		return "ACGTN"[(int)this->cs_[idx]];
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const char& operator[](size_t i) const {
		return this->get(i);
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const char& get(size_t i) const {
		assert_lt(i, this->len_);
		assert_leq(this->cs_[i], 4);
		assert_geq(this->cs_[i], 0);
		return this->cs_[i];
	}

	/**
	 * Return the ith character in the window defined by fw, color,
	 * depth and len.
	 */
	char windowGetDna(
		size_t i,
		bool   fw,
		bool   color,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = this->len_;
		assert_lt(i, len);
		assert_leq(len, this->len_ - depth);
		if(fw) return this->cs_[depth+i];
		else   return color ? this->cs_[depth+len-i-1] :
		                      compDna(this->cs_[depth+len-i-1]);
	}

	/**
	 * Fill the given DNA buffer with the substring specified by fw,
	 * color, depth and len.
	 */
	void windowGetDna(
		SDnaStringFixed<S>& buf,
		bool   fw,
		bool   color,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = this->len_;
		assert_leq(len, this->len_ - depth);
		for(size_t i = 0; i < len; i++) {
			buf.append(fw ? this->cs_[depth+i] :
			                (color ? this->cs_[depth+len-i-1] :
			                         compDna(this->cs_[depth+len-i-1])));
		}
	}

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const char* toZBuf() const { return this->toZBufXForm("ACGTN"); }
};

/**
 * Encapsulates a fixed-length DNA string with characters encoded as
 * chars.  Only capable of encoding A, C, G, T and N.  The length is
 * specified via the template parameter S.
 */

template<int S = 1024, int M = 2>
class SDnaStringExpandable : public SStringExpandable<char, S, M> {
public:

	SDnaStringExpandable() : SStringExpandable<char, S, M>() { }

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SDnaStringExpandable(const SDnaStringExpandable<S, M>& o) :
		SStringExpandable<char, S, M>(o) { }

	/**
	 * Create an SStringFixed from a C++ basic_string.
	 */
	SDnaStringExpandable(const std::basic_string<char>& str) :
		SStringExpandable<char, S, M>(str) { }

	/**
	 * Create an SStringFixed from an array and size.
	 */
	SDnaStringExpandable(const char* b, size_t sz) :
		SStringExpandable<char, S, M>(b, sz) { }

	/**
	 * Create an SStringFixed from a zero-terminated string.
	 */
	SDnaStringExpandable(
		const char* b,
		bool chars = false,
		bool colors = false) :
		SStringExpandable<char, S, M>()
	{
		if(chars) {
			if(colors) {
				installColors(b, strlen(b));
			} else {
				installChars(b, strlen(b));
			}
		} else {
			install(b, strlen(b));
		}
	}

	virtual ~SDnaStringExpandable() { } // C++ needs this

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const char* b, size_t sz) {
		if(this->sz_ < sz) this->expandCopy((sz + S) * M);
		for(size_t i = 0; i < sz; i++) {
			this->cs_[i] = (b[sz-i-1] == 4 ? 4 : b[sz-i-1] ^ 3);
		}
		this->len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const SDnaStringExpandable<S, M>& b) {
		if(this->sz_ < b.len_) this->expandCopy((b.len_ + S) * M);
		for(size_t i = 0; i < b.len_; i++) {
			this->cs_[i] = (b.cs_[b.len_-i-1] == 4 ? 4 : b.cs_[b.len_-i-1] ^ 3);
		}
		this->len_ = b.len_;
	}

	/**
	 * Either reverse or reverse-complement (depending on "color") this
	 * DNA buffer in-place.
	 */
	void reverseComp(bool color = false) {
		if(color) {
			this->reverse();
		} else {
			for(size_t i = 0; i < (this->len_ >> 1); i++) {
				char tmp1 = (this->cs_[i] == 4 ? 4 : this->cs_[i] ^ 3);
				char tmp2 = (this->cs_[this->len_-i-1] == 4 ? 4 : this->cs_[this->len_-i-1] ^ 3);
				this->cs_[i] = tmp2;
				this->cs_[this->len_-i-1] = tmp1;
			}
			// Do middle element iff there are an odd number
			if((this->len_ & 1) != 0) {
				char tmp = this->cs_[this->len_ >> 1];
				tmp = (tmp == 4 ? 4 : tmp ^ 3);
				this->cs_[this->len_ >> 1] = tmp;
			}
		}
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const char* b, size_t sz) {
		if(this->sz_ < sz) this->expandCopy((sz + S) * M);
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for(size_t i = 0; i < sz; i++) {
			assert_range(0, 4, (int)this->cs_[i]);
		}
#endif
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const char* b, size_t sz) {
		if(this->sz_ < sz) this->expandCopy((sz + S) * M);
		for(size_t i = 0; i < sz; i++) {
			assert_in(toupper(b[i]), "ACGTN-");
			this->cs_[i] = asc2dna[(int)b[i]];
			assert_range(0, 4, (int)this->cs_[i]);
		}
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII color characters into normal DNA
	 * characters.
	 */
	virtual void installColors(const char* b, size_t sz) {
		if(this->sz_ < sz) this->expandCopy((sz + S) * M);
		for(size_t i = 0; i < sz; i++) {
			assert_in(b[i], "0123.");
			this->cs_[i] = asc2col[(int)b[i]];
			assert_range(0, 4, (int)this->cs_[i]);
		}
		this->len_ = sz;
	}

	/**
	 * Copy C++ string of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const std::basic_string<char>& str) {
		installChars(str.c_str(), str.length());
	}

	/**
	 * Copy C++ string of ASCII color characters into normal DNA
	 * characters.
	 */
	virtual void installColors(const std::basic_string<char>& str) {
		installColors(str.c_str(), str.length());
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_range(0, 4, c);
		this->cs_[idx] = c;
	}

	/**
	 * Append DNA char c.
	 */
	void append(const char& c) {
		if(this->sz_ < this->len_ + 1) {
			this->expandCopy((this->len_ + 1 + S) * M);
		}
		assert_range(0, 4, (int)c);
		this->cs_[this->len_++] = c;
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void setChar(char c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_in(toupper(c), "ACGTN");
		this->cs_[idx] = asc2dna[(int)c];
	}

	/**
	 * Append DNA character.
	 */
	void appendChar(char c) {
		if(this->sz_ < this->len_ + 1) {
			this->expandCopy((this->len_ + 1 + S) * M);
		}
		assert_in(toupper(c), "ACGTN");
		this->cs_[this->len_++] = asc2dna[(int)c];
	}

	/**
	 * Return DNA character corresponding to element 'idx'.
	 */
	char toChar(size_t idx) const {
		assert_range(0, 4, (int)this->cs_[idx]);
		return "ACGTN"[(int)this->cs_[idx]];
	}

	/**
	 * Retrieve constant version of element i.
	 */
	inline const char& operator[](size_t i) const {
		return this->get(i);
	}

	/**
	 * Retrieve constant version of element i.
	 */
	inline const char& get(size_t i) const {
		assert_lt(i, this->len_);
		assert_range(0, 4, (int)this->cs_[i]);
		return this->cs_[i];
	}

	/**
	 * Return the ith character in the window defined by fw, color,
	 * depth and len.
	 */
	char windowGetDna(
		size_t i,
		bool   fw,
		bool   color,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = this->len_;
		assert_lt(i, len);
		assert_leq(len, this->len_ - depth);
		if(fw) return this->cs_[depth+i];
		else   return color ? this->cs_[depth+len-i-1] :
		                      compDna(this->cs_[depth+len-i-1]);
	}

	/**
	 * Fill the given DNA buffer with the substring specified by fw,
	 * color, depth and len.
	 */
	void windowGetDna(
		SDnaStringExpandable<S, M>& buf,
		bool   fw,
		bool   color,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = this->len_;
		assert_leq(len, this->len_ - depth);
		for(size_t i = 0; i < len; i++) {
			buf.append(fw ? this->cs_[depth+i] :
			                (color ? this->cs_[depth+len-i-1] :
			                         compDna(this->cs_[depth+len-i-1])));
		}
	}

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const char* toZBuf() const { return this->toZBufXForm("ACGTN"); }
};

/**
 * Encapsulates an expandable DNA string with characters encoded as
 * char-sized masks.  Encodes A, C, G, T, and all IUPAC, as well as the
 * empty mask indicating "matches nothing."
 */
template<int S = 16, int M = 2>
class SDnaMaskString : public SStringExpandable<char, S, M> {
public:

	SDnaMaskString() : SStringExpandable<char, S, M>() { }

	/**
	 * Create an SStringFixed from another SStringFixed.
	 */
	SDnaMaskString(const SDnaMaskString<S, M>& o) :
		SStringExpandable<char, S, M>(o) { }

	/**
	 * Create an SStringFixed from a C++ basic_string.
	 */
	SDnaMaskString(const std::basic_string<char>& str) :
		SStringExpandable<char, S, M>(str) { }

	/**
	 * Create an SStringFixed from an array and size.
	 */
	SDnaMaskString(const char* b, size_t sz) :
		SStringExpandable<char, S, M>(b, sz) { }

	/**
	 * Create an SStringFixed from a zero-terminated string.
	 */
	SDnaMaskString(const char* b, bool chars = false) :
		SStringExpandable<char, S, M>()
	{
		if(chars) {
			installChars(b, strlen(b));
		} else {
			install(b, strlen(b));
		}
	}

	virtual ~SDnaMaskString() { } // C++ needs this

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const char* b, size_t sz) {
		while(this->sz_ < sz) {
			this->expandNoCopy((sz + S) * M);
		}
		for(size_t i = 0; i < sz; i++) {
			this->cs_[i] = maskcomp[(int)b[sz-i-1]];
		}
		this->len_ = sz;
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string, reverse-
	 * complementing them in the process, assuming an encoding where
	 * 0=A, 1=C, 2=G, 3=T, 4=N.
	 */
	void installReverseComp(const SDnaMaskString<S, M>& b) {
		while(this->sz_ < b.len_) {
			this->expandNoCopy((b.len_ + S) * M);
		}
		for(size_t i = 0; i < b.len_; i++) {
			this->cs_[i] = maskcomp[(int)b.cs_[b.len_-i-1]];
		}
		this->len_ = b.len_;
	}

	/**
	 * Either reverse or reverse-complement (depending on "color") this
	 * DNA buffer in-place.
	 */
	void reverseComp(bool color = false) {
		if(color) {
			this->reverse();
		} else {
			for(size_t i = 0; i < (this->len_ >> 1); i++) {
				char tmp1 = maskcomp[(int)this->cs_[i]];
				char tmp2 = maskcomp[(int)this->cs_[this->len_-i-1]];
				this->cs_[i] = tmp2;
				this->cs_[this->len_-i-1] = tmp1;
			}
			// Do middle element iff there are an odd number
			if((this->len_ & 1) != 0) {
				char tmp = this->cs_[this->len_ >> 1];
				tmp = maskcomp[(int)tmp];
				this->cs_[this->len_ >> 1] = tmp;
			}
		}
	}

	/**
	 * Copy 'sz' bytes from buffer 'b' into this string.
	 */
	virtual void install(const char* b, size_t sz) {
		while(this->sz_ < sz) {
			this->expandNoCopy((sz + S) * M);
		}
		memcpy(this->cs_, b, sz);
#ifndef NDEBUG
		for(size_t i = 0; i < sz; i++) {
			assert_range((int)this->cs_[i], 0, 15);
		}
#endif
		this->len_ = sz;
	}

	/**
	 * Copy buffer 'b' of ASCII DNA characters into DNA masks.
	 */
	virtual void installChars(const char* b, size_t sz) {
		while(this->sz_ < sz) {
			this->expandNoCopy((sz + S) * M);
		}
		for(size_t i = 0; i < sz; i++) {
			assert_in(b[i], iupacs);
			this->cs_[i] = asc2dnamask[(int)b[i]];
			assert_range((int)this->cs_[i], 0, 15);
		}
		this->len_ = sz;
	}

	/**
	 * Copy C++ string of ASCII DNA characters into normal DNA
	 * characters.
	 */
	virtual void installChars(const std::basic_string<char>& str) {
		installChars(str.c_str(), str.length());
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void set(int c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_range(c, 0, 15);
		this->cs_[idx] = c;
	}

	/**
	 * Append DNA char c.
	 */
	void append(const char& c) {
		while(this->sz_ < this->len_+1) {
			this->expandNoCopy((this->len_ + 1 + S) * M);
		}
		assert_range((int)c, 0, 15);
		this->cs_[this->len_++] = c;
	}

	/**
	 * Set DNA character at index 'idx' to 'c'.
	 */
	void setChar(char c, size_t idx) {
		assert_lt(idx, this->len_);
		assert_in(toupper(c), iupacs);
		this->cs_[idx] = asc2dnamask[(int)c];
	}

	/**
	 * Append DNA character.
	 */
	void appendChar(char c) {
		while(this->sz_ < this->len_+1) {
			expandNoCopy((this->len_ + 1 + S) * M);
		}
		assert_in(toupper(c), iupacs);
		this->cs_[this->len_++] = asc2dnamask[(int)c];
	}

	/**
	 * Return DNA character corresponding to element 'idx'.
	 */
	char toChar(size_t idx) const {
		assert_range((int)this->cs_[idx], 0, 15);
		return mask2iupac[(int)this->cs_[idx]];
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const char& operator[](size_t i) const {
		return this->get(i);
	}

	/**
	 * Retrieve mutable version of element i.
	 */
	char& operator[](size_t i) {
		return this->get(i);
	}

	/**
	 * Retrieve constant version of element i.
	 */
	const char& get(size_t i) const {
		assert_lt(i, this->len_);
		assert_range((int)this->cs_[i], 0, 15);
		return this->cs_[i];
	}

	/**
	 * Retrieve mutable version of element i.
	 */
	char& get(size_t i) {
		assert_lt(i, this->len_);
		assert_range((int)this->cs_[i], 0, 15);
		return this->cs_[i];
	}

	/**
	 * Return the ith character in the window defined by fw, color,
	 * depth and len.
	 */
	char windowGetDna(
		size_t i,
		bool   fw,
		bool   color,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = this->len_;
		assert_lt(i, len);
		assert_leq(len, this->len_ - depth);
		if(fw) return this->cs_[depth+i];
		else   return color ? this->cs_[depth+len-i-1] :
		                      maskcomp[this->cs_[depth+len-i-1]];
	}

	/**
	 * Fill the given DNA buffer with the substring specified by fw,
	 * color, depth and len.
	 */
	void windowGetDna(
		SDnaStringFixed<S>& buf,
		bool   fw,
		bool   color,
		size_t depth = 0,
		size_t len = 0) const
	{
		if(len == 0) len = this->len_;
		assert_leq(len, this->len_ - depth);
		for(size_t i = 0; i < len; i++) {
			buf.append(fw ? this->cs_[depth+i] :
			                (color ? this->cs_[depth+len-i-1] :
			                         maskcomp[this->cs_[depth+len-i-1]]));
		}
	}

	/**
	 * Sample a random substring of the given length from this DNA
	 * string and install the result in 'dst'.
	 */
	template<typename T>
	void randSubstr(
		RandomSource& rnd,  // pseudo-random generator
		T& dst,             // put sampled substring here
		size_t len,         // length of substring to extract
		bool watson = true, // true -> possibly extract from Watson strand
		bool crick = true)  // true -> possibly extract from Crick strand
	{
		assert(watson || crick);
		assert_geq(this->len_, len);
		size_t poss = this->len_ - len + 1;
		assert_gt(poss, 0);
		uint32_t rndoff = rnd.nextU32() % poss;
		bool fw;
		if     (watson && !crick) fw = true;
		else if(!watson && crick) fw = false;
		else {
			fw = (rnd.nextU2() == 0) ? true : false;
		}
		if(fw) {
			// Install Watson substring
			for(size_t i = 0; i < len; i++) {
				dst[i] = this->cs_[i + rndoff];
			}
		} else {
			// Install Crick substring
			for(size_t i = 0; i < len; i++) {
				dst[i] = maskcomp[(int)this->cs_[i + rndoff + (len - i - 1)]];
			}
		}
	}

	/**
	 * Put a terminator in the 'len_'th element and then return a
	 * pointer to the buffer.  Useful for printing.
	 */
	virtual const char* toZBuf() const { return this->toZBufXForm(iupacs); }
};

typedef SStringExpandable<char, 1024, 2> BTString;
typedef SDnaStringExpandable<1024, 2>    BTDnaString;
typedef SDnaMaskString<32, 2>            BTDnaMask;

#endif /* SSTRING_H_ */
