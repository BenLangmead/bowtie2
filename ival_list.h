/**
 * ival_list.h
 */

#ifndef IVAL_LIST_H_
#define IVAL_LIST_H_

#include "ds.h"
#include <algorithm>
#include <utility>

/**
 * Simple interval.
 */
struct EIvalMergeListPair {
	EIvalMergeListPair() { one = two = 0; }
	EIvalMergeListPair(const int64_t& o, const int64_t& t) { one = o; two = t; }
	bool operator<(const EIvalMergeListPair& o) const {
		if(one < o.one) {
			return true;
		}
		if(one > o.one) {
			return false;
		}
		return two < o.two;
	}
	int64_t one, two;
}; 

/**
 * Encapsulates the "union" of a collection of intervals.  Intervals are stored
 * in a sorted list.  Intervals can be added but not removed.  Supports just
 * one type of query for now: locusPresent().
 */
class EIvalMergeList {
public:

	typedef EIvalMergeListPair TPair;
	static const size_t DEFAULT_UNSORTED_SZ = 16;

	explicit EIvalMergeList(int cat = 0) :
		sorted_(cat),
		sortedLhs_(cat),
		unsorted_(cat),
		unsortedSz_(DEFAULT_UNSORTED_SZ)
	{ }

	explicit EIvalMergeList(size_t unsortedSz, int cat = 0) :
		sorted_(cat),
		sortedLhs_(cat),
		unsorted_(cat),
		unsortedSz_(unsortedSz)
	{ }
	
	/**
	 * Add a new interval to the list.
	 */
	void add(const int64_t& i, const int64_t& f) {
		assert_leq(unsorted_.size(), unsortedSz_);
		if(unsorted_.size() < unsortedSz_) {
			unsorted_.push_back(TPair(i, f));
		}
		if(unsorted_.size() == unsortedSz_) {
			flush();
		}
	}

	/**
	 * Move all unsorted interval information into the sorted list and re-sort.
	 * Merge overlapping intervals.
	 */
	void flush() {
		for(size_t i = 0; i < unsorted_.size(); i++) {
			sorted_.push_back(unsorted_[i]);
		}
		sorted_.sort();
		merge();
		sortedLhs_.clear();
		for(size_t i = 0; i < sorted_.size(); i++) {
			sortedLhs_.push_back(sorted_[i].one);
		}
		assert(sortedLhs_.sorted());
		unsorted_.clear();
	}
	
	/**
	 * Check that this interval list is internally consistent.
	 */
	bool repOk() const {
		assert_eq(sorted_.size(), sortedLhs_.size());
		return true;
	}
	
	/**
	 * Remove all ranges from the list.
	 */
	void reset() { clear(); }
	
	/**
	 * Remove all ranges from the list.
	 */
	void clear() {
		sorted_.clear();
		sortedLhs_.clear();
		unsorted_.clear();
	}
	
	/**
	 * Return true iff this locus is present in one of the intervals in the
	 * list.
	 */
	bool locusPresent(const int64_t& loc) const {
		return locusPresentUnsorted(loc) || locusPresentSorted(loc);
	}
	
	/**
	 * Return the number of intervals added since the last call to reset() or
	 * clear().
	 */
	size_t size() const {
		return sorted_.size() + unsorted_.size();
	}
	
	/**
	 * Return true iff list is empty.
	 */
	bool empty() const {
		return sorted_.empty() && unsorted_.empty();
	}
	
protected:
	
	/**
	 * Go through the sorted interval list and merge adjacent entries that
	 * overlap.
	 */
	void merge() {
		size_t nmerged = 0;
		for(size_t i = 1; i < sorted_.size(); i++) {
			if(sorted_[i-1].two >= sorted_[i].one) {
				nmerged++;
				assert_leq(sorted_[i-1].one, sorted_[i].one);
				sorted_[i].one = sorted_[i-1].one;
				sorted_[i].two = std::max(sorted_[i-1].two, sorted_[i].two);
				sorted_[i-1].one = std::numeric_limits<int64_t>::max();
				sorted_[i-1].two = std::numeric_limits<int64_t>::max();
			}
		}
		sorted_.sort();
		assert_lt(nmerged, sorted_.size());
		sorted_.resize(sorted_.size()-nmerged);
	}

	/**
	 * Return true iff the given locus is present in one of the intervals in
	 * the sorted list.
	 */
	bool locusPresentSorted(int64_t loc) const {
		assert(repOk());
		if(sorted_.empty()) {
			return false;
		}
		size_t beg = sortedLhs_.bsearchLoBound(loc);
		if(beg == sortedLhs_.size() || sortedLhs_[beg] > loc) {
			// Check element before
			if(beg == 0) {
				return false;
			}
			return (sorted_[beg-1].one <= loc && sorted_[beg-1].two > loc);
		} else {
			assert_eq(loc, sortedLhs_[beg]);
			return true;
		}
	}

	/**
	 * Return true iff the given locus is present in one of the intervals in
	 * the unsorted list.
	 */
	bool locusPresentUnsorted(int64_t loc) const {
		for(size_t i = 0; i < unsorted_.size(); i++) {
			if(loc >= unsorted_[i].one && loc < unsorted_[i].two) { 
				return true;
			}
		}
		return false;
	}

	EList<TPair>   sorted_;     // LHS, RHS sorted
	EList<int64_t> sortedLhs_;  // LHS, index into sorted_, sorted
	EList<TPair>   unsorted_;   // unsorted
	size_t         unsortedSz_; // max allowed size of unsorted_
};

/**
 * Added strandedness to EIvalMergeList by maintaining separate lists for
 * Watson and Crick strands.
 */
class EIvalMergeStrandedList {
public:

	explicit EIvalMergeStrandedList(int cat = 0) :
		listfw_(cat),
		listrc_(cat)
	{ }

	explicit EIvalMergeStrandedList(size_t unsortedSz, int cat = 0) :
		listfw_(unsortedSz, cat),
		listrc_(unsortedSz, cat)
	{ }

	/**
	 * Add a new interval to the list.
	 */
	void add(const int64_t& i, const int64_t& f, bool fw) {
		if(fw) {
			listfw_.add(i, f);
		} else {
			listrc_.add(i, f);
		}
	}
	
	/**
	 * Check that this interval list is internally consistent.
	 */
	bool repOk() const {
		assert(listfw_.repOk());
		assert(listrc_.repOk());
		return true;
	}
	
	/**
	 * Remove all ranges from the list.
	 */
	void reset() { clear(); }
	
	/**
	 * Remove all ranges from the list.
	 */
	void clear() {
		listfw_.clear();
		listrc_.clear();
	}
	
	/**
	 * Return true iff this locus is present in one of the intervals in the
	 * list.
	 */
	bool locusPresent(const int64_t& loc, bool fw) const {
		if(fw) {
			return listfw_.locusPresent(loc);
		} else {
			return listrc_.locusPresent(loc);
		}
	}
	
	/**
	 * Return the number of intervals added since the last call to reset() or
	 * clear().
	 */
	size_t size() const {
		return listfw_.size() + listrc_.size();
	}
	
	/**
	 * Return true iff list is empty.
	 */
	bool empty() const {
		return listfw_.empty() && listrc_.empty();
	}
	
protected:
	
	EIvalMergeList listfw_, listrc_;
};

#endif /*ndef IVAL_LIST_H_*/
