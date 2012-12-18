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

#ifndef REF_COORD_H_
#define REF_COORD_H_

#include <stdint.h>
#include <iostream>
#include <limits>
#include "assert_helpers.h"

typedef int64_t TRefId;
typedef int64_t TRefOff;

/**
 * Encapsulates a reference coordinate; i.e. identifiers for (a) a
 * reference sequence, and (b) a 0-based offset into that sequence.
 */
class Coord {

public:

	Coord() { invalidate(); }

	Coord(const Coord& c) { init(c); }
	
	Coord(TRefId rf, TRefOff of, bool fw) { init(rf, of, fw); }

	/**
	 * Copy given fields into this Coord.
	 */
	void init(TRefId rf, TRefOff of, bool fw) {
		ref_ = rf;
		off_ = of;
		orient_ = (fw ? 1 : 0);
	}

	/**
	 * Copy contents of given Coord into this one.
	 */
	void init(const Coord& c) {
		ref_ = c.ref_;
		off_ = c.off_;
		orient_ = c.orient_;
	}
	
	/**
	 * Return true iff this Coord is identical to the given Coord.
	 */
	bool operator==(const Coord& o) const {
		assert(valid());
		assert(o.valid());
		return ref_ == o.ref_ && off_ == o.off_ && fw() == o.fw();
	}

	/**
	 * Return true iff this Coord is less than the given Coord.  One Coord is
	 * less than another if (a) its reference id is less, (b) its orientation is
	 * less, or (c) its offset is less.
	 */
	bool operator<(const Coord& o) const {
		//assert(valid());
		//assert(o.valid());
		if(ref_ < o.ref_) return true;
		if(ref_ > o.ref_) return false;
		if(orient_ < o.orient_) return true;
		if(orient_ > o.orient_) return false;
		if(off_ < o.off_) return true;
		if(off_ > o.off_) return false;
		return false;
	}
	
	/**
	 * Return the opposite result from operator<.
	 */
	bool operator>=(const Coord& o) const {
		return !((*this) < o);
	}
	
	/**
	 * Return true iff this Coord is greater than the given Coord.  One Coord
	 * is greater than another if (a) its reference id is greater, (b) its
	 * orientation is greater, or (c) its offset is greater.
	 */
	bool operator>(const Coord& o) const {
		//assert(valid());
		//assert(o.valid());
		if(ref_ > o.ref_) return true;
		if(ref_ < o.ref_) return false;
		if(orient_ > o.orient_) return true;
		if(orient_ < o.orient_) return false;
		if(off_ > o.off_) return true;
		if(off_ < o.off_) return false;
		return false;
	}
	
	/**
	 * Return the opposite result from operator>.
	 */
	bool operator<=(const Coord& o) const {
		return !((*this) > o);
	}
	
	/**
	 * Make this coord invalid.
	 */
	void invalidate() {
		ref_ = std::numeric_limits<TRefId>::max();
		off_ = std::numeric_limits<TRefOff>::max();
		orient_ = -1;
	}
	
	/**
	 * Return true iff this Coord is valid (i.e. ref and off have both
	 * been set since the last call to invalidate()).
	 */
	bool valid() const {
		if(ref_ != std::numeric_limits<TRefId>::max() &&
		   off_ != std::numeric_limits<TRefOff>::max())
		{
			assert(orient_ == 0 || orient_ == 1);
			return true;
		}
		return false;
	}
	
	/**
	 * Get orientation of the Coord.
	 */
	bool fw() const {
		assert(valid());
		assert(orient_ == 0 || orient_ == 1);
		return orient_ == 1;
	}
	
#ifndef NDEBUG
	/**
	 * Check that coord is internally consistent.
	 */
	bool repOk() const {
		if(ref_ != std::numeric_limits<TRefId>::max() &&
		   off_ != std::numeric_limits<TRefOff>::max())
		{
			assert(orient_ == 0 || orient_ == 1);
		}
		return true;
	}
#endif
	
	/**
	 * Check whether an interval defined by this coord and having
	 * length 'len' is contained within an interval defined by
	 * 'inbegin' and 'inend'.
	 */
	bool within(int64_t len, int64_t inbegin, int64_t inend) const {
		return off_ >= inbegin && off_ + len <= inend;
	}
	
	inline TRefId  ref()    const { return ref_; }
	inline TRefOff off()    const { return off_; }
	inline int     orient() const { return orient_; }
	
	inline void setRef(TRefId  id)  { ref_ = id;  }
	inline void setOff(TRefOff off) { off_ = off; }

	inline void adjustOff(TRefOff off) { off_ += off; }

protected:

	TRefId  ref_;    // which reference?
	TRefOff off_;    // 0-based offset into reference
	int     orient_; // true -> Watson strand
};

std::ostream& operator<<(std::ostream& out, const Coord& c);

/**
 * Encapsulates a reference interval, which consists of a Coord and a length.
 */
class Interval {

public:
	
	Interval() { invalidate(); }
	
	explicit Interval(const Coord& upstream, TRefOff len) {
		init(upstream, len);
	}

	explicit Interval(TRefId rf, TRefOff of, bool fw, TRefOff len) {
		init(rf, of, fw, len);
	}

	void init(const Coord& upstream, TRefOff len) {
		upstream_ = upstream;
		len_ = len;
	}
	
	void init(TRefId rf, TRefOff of, bool fw, TRefOff len) {
		upstream_.init(rf, of, fw);
		len_ = len;
	}
	
	/**
	 * Set offset.
	 */
	void setOff(TRefOff of) {
		upstream_.setOff(of);
	}

	/**
	 * Set length.
	 */
	void setLen(TRefOff len) {
		len_ = len;
	}

	/**
	 * Make this coord invalid.
	 */
	void invalidate() {
		upstream_.invalidate();
		len_ = 0;
	}
	
	/**
	 * Return true iff this Interval is valid.
	 */
	bool valid() const {
		if(upstream_.valid()) {
			assert_gt(len_, 0);
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Return true iff this Interval is equal to the given Interval,
	 * i.e. if they cover the same set of positions.
	 */
	bool operator==(const Interval& o) const {
		return upstream_ == o.upstream_ &&
		       len_ == o.len_;
	}

	/**
	 * Return true iff this Interval is less than the given Interval.
	 * One interval is less than another if its upstream location is
	 * prior to the other's or, if their upstream locations are equal,
	 * if its length is less than the other's.
	 */
	bool operator<(const Interval& o) const {
		if(upstream_ < o.upstream_) return true;
		if(upstream_ > o.upstream_) return false;
		if(len_ < o.len_) return true;
		return false;
	}
	
	/**
	 * Return opposite result from operator<.
	 */
	bool operator>=(const Interval& o) const {
		return !((*this) < o);
	}

	/**
	 * Return true iff this Interval is greater than than the given
	 * Interval.  One interval is greater than another if its upstream
	 * location is after the other's or, if their upstream locations
	 * are equal, if its length is greater than the other's.
	 */
	bool operator>(const Interval& o) const {
		if(upstream_ > o.upstream_) return true;
		if(upstream_ < o.upstream_) return false;
		if(len_ > o.len_) return true;
		return false;
	}

	/**
	 * Return opposite result from operator>.
	 */
	bool operator<=(const Interval& o) const {
		return !((*this) > o);
	}
	
	/**
	 * Set upstream Coord.
	 */
	void setUpstream(const Coord& c) {
		upstream_ = c;
	}

	/**
	 * Set length.
	 */
	void setLength(TRefOff l) {
		len_ = l;
	}
	
	inline TRefId  ref()    const { return upstream_.ref(); }
	inline TRefOff off()    const { return upstream_.off(); }
	inline TRefOff dnoff()  const { return upstream_.off() + len_; }
	inline int     orient() const { return upstream_.orient(); }

	/**
	 * Return a Coord encoding the coordinate just past the downstream edge of
	 * the interval.
	 */
	inline Coord downstream() const {
		return Coord(
			upstream_.ref(),
			upstream_.off() + len_,
			upstream_.orient());
	}
	
	/**
	 * Return true iff the given Coord is inside this Interval.
	 */
	inline bool contains(const Coord& c) const {
		return
			c.ref()    == ref() &&
			c.orient() == orient() &&
			c.off()    >= off() &&
			c.off()    <  dnoff();
	}

	/**
	 * Return true iff the given Coord is inside this Interval, without
	 * requiring orientations to match.
	 */
	inline bool containsIgnoreOrient(const Coord& c) const {
		return
			c.ref()    == ref() &&
			c.off()    >= off() &&
			c.off()    <  dnoff();
	}

	/**
	 * Return true iff the given Interval is inside this Interval.
	 */
	inline bool contains(const Interval& c) const {
		return
			c.ref()    == ref() &&
			c.orient() == orient() &&
			c.off()    >= off() &&
			c.dnoff()  <= dnoff();
	}

	/**
	 * Return true iff the given Interval is inside this Interval, without
	 * requiring orientations to match.
	 */
	inline bool containsIgnoreOrient(const Interval& c) const {
		return
			c.ref()    == ref() &&
			c.off()    >= off() &&
			c.dnoff()  <= dnoff();
	}

	/**
	 * Return true iff the given Interval overlaps this Interval.
	 */
	inline bool overlaps(const Interval& c) const {
		return
			c.ref()    == upstream_.ref() &&
			c.orient() == upstream_.orient() &&
			((off() <= c.off()   && dnoff() > c.off())   ||
			 (off() <= c.dnoff() && dnoff() > c.dnoff()) ||
			 (c.off() <= off()   && c.dnoff() > off())   ||
			 (c.off() <= dnoff() && c.dnoff() > dnoff()));
	}

	/**
	 * Return true iff the given Interval overlaps this Interval, without
	 * requiring orientations to match.
	 */
	inline bool overlapsIgnoreOrient(const Interval& c) const {
		return
			c.ref()    == upstream_.ref() &&
			((off() <= c.off()   && dnoff() > c.off())   ||
			 (off() <= c.dnoff() && dnoff() > c.dnoff()) ||
			 (c.off() <= off()   && c.dnoff() > off())   ||
			 (c.off() <= dnoff() && c.dnoff() > dnoff()));
	}
	
	inline const Coord&  upstream()   const { return upstream_; }
	inline TRefOff       len()      const { return len_;      }

#ifndef NDEBUG
	/**
	 * Check that the Interval is internally consistent.
	 */
	bool repOk() const {
		assert(upstream_.repOk());
		assert_geq(len_, 0);
		return true;
	}
#endif

	inline void adjustOff(TRefOff off) { upstream_.adjustOff(off); }

protected:

	Coord   upstream_;
	TRefOff len_;
};

std::ostream& operator<<(std::ostream& out, const Interval& c);

#endif /*ndef REF_COORD_H_*/
