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

#ifndef IVAL_LIST_H_
#define IVAL_LIST_H_

#include "ds.h"
#include "ref_coord.h"
#include <utility>

/**
 * Encapsulates the "union" of a collection of intervals.  Intervals are stored
 * in a sorted list.  Intervals can be added but not removed.  Supports just
 * one type of query for now: locusPresent().
 */
class EIvalMergeList {
public:

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
	 * Set the maximum size of the unsorted list.
	 */
	void setUnsortedSize(size_t usz) {
		unsortedSz_ = usz;
	} 
	
	/**
	 * Add a new interval to the list.
	 */
	void add(const Interval& i) {
		assert_leq(unsorted_.size(), unsortedSz_);
		if(unsorted_.size() < unsortedSz_) {
			unsorted_.push_back(i);
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
			sortedLhs_.push_back(sorted_[i].upstream());
		}
		assert(sortedLhs_.sorted());
		unsorted_.clear();
	}
	
#ifndef NDEBUG
	/**
	 * Check that this interval list is internally consistent.
	 */
	bool repOk() const {
		assert_eq(sorted_.size(), sortedLhs_.size());
		return true;
	}
#endif
	
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
	bool locusPresent(const Coord& loc) const {
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
			if(sorted_[i-1].downstream() >= sorted_[i].upstream()) {
				nmerged++;
				assert_leq(sorted_[i-1].upstream(), sorted_[i].upstream());
				Coord up = std::min(sorted_[i-1].upstream(), sorted_[i].upstream());
				Coord dn = std::max(sorted_[i-1].downstream(), sorted_[i].downstream());
				sorted_[i].setUpstream(up);
				sorted_[i].setLength(dn.off() - up.off());
				sorted_[i-1].invalidate();
			}
		}
		sorted_.sort();
		assert_lt(nmerged, sorted_.size());
		sorted_.resize(sorted_.size()-nmerged);
#ifndef NDEBUG
		for(size_t i = 0; i < sorted_.size(); i++) {
			assert(sorted_[i].valid());
		}
#endif
	}

	/**
	 * Return true iff the given locus is present in one of the intervals in
	 * the sorted list.
	 */
	bool locusPresentSorted(const Coord& loc) const {
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
			return sorted_[beg-1].contains(loc);
		} else {
			assert_eq(loc, sortedLhs_[beg]);
			return true;
		}
	}

	/**
	 * Return true iff the given locus is present in one of the intervals in
	 * the unsorted list.
	 */
	bool locusPresentUnsorted(const Coord& loc) const {
		for(size_t i = 0; i < unsorted_.size(); i++) {
			if(unsorted_[i].contains(loc)) {
				return true;
			}
		}
		return false;
	}

	EList<Interval> sorted_;     // LHS, RHS sorted
	EList<Coord>    sortedLhs_;  // LHS, index into sorted_, sorted
	EList<Interval> unsorted_;   // unsorted
	size_t          unsortedSz_; // max allowed size of unsorted_
};

/**
 * Binned version of the above.  We bin using the low bits of the reference
 * sequence.
 */
class EIvalMergeListBinned {
public:

	static const size_t NBIN = 7;

	explicit EIvalMergeListBinned(int cat = 0) : bins_(1 << NBIN, cat) {
		bins_.resize(1 << NBIN);
	}

	explicit EIvalMergeListBinned(
		size_t unsortedSz,
		int cat = 0) : bins_(1 << NBIN, cat)
	{
		bins_.resize(1 << NBIN);
		for(size_t i = 0; i < (1 << NBIN); i++) {
			bins_[i].setUnsortedSize(unsortedSz);
		}
	}

	/**
	 * Add a new interval to the list.
	 */
	void add(const Interval& i) {
		size_t bin = i.ref() & ~(0xffffffff << NBIN);
		assert_lt(bin, bins_.size());
		bins_[bin].add(i);
	}

#ifndef NDEBUG
	/**
	 * Check that this interval list is internally consistent.
	 */
	bool repOk() const {
		for(size_t i = 0; i < bins_.size(); i++) {
			assert(bins_[i].repOk());
		}
		return true;
	}
#endif
	
	/**
	 * Remove all ranges from the list.
	 */
	void reset() { clear(); }
	
	/**
	 * Remove all ranges from the list.
	 */
	void clear() {
		for(size_t i = 0; i < bins_.size(); i++) {
			bins_[i].clear();
		}
	}
	
	/**
	 * Return true iff this locus is present in one of the intervals in the
	 * list.
	 */
	bool locusPresent(const Coord& loc) const {
		size_t bin = loc.ref() & ~(0xffffffff << NBIN);
		assert_lt(bin, bins_.size());
		return bins_[bin].locusPresent(loc);
	}
	
	/**
	 * Return the number of intervals added since the last call to reset() or
	 * clear().
	 */
	size_t size() const {
		// TODO: Keep track of size
		size_t sz = 0;
		for(size_t i = 0; i < bins_.size(); i++) {
			sz += bins_[i].size();
		}
		return sz;
	}
	
	/**
	 * Return true iff list is empty.
	 */
	bool empty() const {
		return size() == 0;
	}
	
protected:
	
	EList<EIvalMergeList> bins_;
};

#endif /*ndef IVAL_LIST_H_*/
