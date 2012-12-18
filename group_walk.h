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

/*
 * group_walk.h
 *
 * Classes and routines for walking a set of BW ranges backwards from the edge
 * of a seed hit with the goal of resolving the offset of each row in each
 * range.  Here "offset" means offset into the concatenated string of all
 * references.  The main class is 'GroupWalk' and an important helper is
 * 'GWState'.
 *
 * For each combination of seed offset and orientation, there is an associated
 * QVal.  Each QVal describes a (possibly empty) set of suffix array ranges.
 * Call these "seed range sets."  Each range in the set is "backed" by a range
 * of the salist, represented as a PListSlice. Such a range is the origin of a
 * walk.
 *
 * When an offset is resolved, it is entered into the salist via the
 * PListSlice.  Note that other routines in this same thread might also be
 * setting elements of the salist, so routines here should expect that elements
 * can go from unresolved to resolved at any time.
 *
 * What bookkeeping do we have to do as we walk?  Before the first step, we
 * convert the initial QVal into a list of SATuples; the SATuples are our link
 * to the correpsonding ranges in the suffix array.  The list of SATuples is
 * then converted to a list of GWState objects; these keep track of where we
 * are in our walk (e.g. what 'top' and 'bot' are, how many steps have we gone,
 * etc) as well as how the elements in the current range correspond to elements
 * from the original range.
 *
 * The user asks the GroupWalk to resolve another offset by calling advance().
 * advance() can be called in various ways:
 *
 * (a) The user can request that the GroupWalk proceed until a
 *     *particular* element is resolved, then return that resolved
 *     element.  Other elements may be resolved along the way, but
 *     those results are buffered and may be dispensed in future calls
 *     to advance().
 *
 * (b) The user can request that the GroupWalk select an as-yet-
 *     unreported element at random and and proceed until that element
 *     is resolved and report it.  Again, other elements may be
 *     resolved along the way but they are buffered.
 *
 * (c) The user can request that the GroupWalk resolve elements in a
 *     particular BW range (with a particular offset and orientation)
 *     in an order of its choosing.  The GroupWalk in this case
 *     attempts to resolve as many offsets as possible as quickly as
 *     possible, and returns them as soon as they're found.  The res_
 *     buffer is used in this case.
 *
 * (d) Like (c) but resolving elements at a paritcular offset and
 *     orientation instead of at a specific BW range.  The res_ buffer
 *     is used in this case, since there's a chance that the 
 *
 * There are simple ways to heuristically reduce the problem size while
 * maintaining randomness.  For instance, the user put a ceiling on the
 * number of elements that we walk from any given seed offset or range.
 * We can then trim away random subranges to reduce the size of the
 * problem.  There is no need for the caller to do this for us.
 */

#ifndef GROUP_WALK_H_
#define GROUP_WALK_H_

#include <stdint.h>
#include <limits>
#include "ds.h"
#include "bt2_idx.h"
#include "read.h"
#include "reference.h"
#include "mem_ids.h"

typedef uint32_t TIndexOff;

/**
 * Encapsulate an SA range and an associated list of slots where the resolved
 * offsets can be placed.
 */
template<typename T>
class SARangeWithOffs {

public:

	SARangeWithOffs() { reset(); };

	SARangeWithOffs(TIndexOff tf, size_t len, const T& o) {
		init(tf, len, o);
	}
	
	void init(TIndexOff tf, size_t len_, const T& o) {
		topf = tf; len = len_, offs = o;
	}

	/**
	 * Reset to uninitialized state.
	 */
	void reset() { topf = std::numeric_limits<TIndexOff>::max(); }
	
	/**
	 * Return true if this is initialized.
	 */
	bool inited() const {
		return topf != std::numeric_limits<TIndexOff>::max();
	}
	
	/**
	 * Return the number of times this reference substring occurs in the
	 * reference, which is also the size of the 'offs' TSlice.
	 */
	size_t size() const { return offs.size(); }

	TIndexOff topf; // top in BWT index
	size_t    len;  // length of the reference sequence involved
	T         offs; // offsets
};

/**
 * A group of per-thread state that can be shared between all the GroupWalks
 * used in that thread.
 */
struct GroupWalkState {

	GroupWalkState(int cat) : map(cat) {
		masks[0].setCat(cat);
		masks[1].setCat(cat);
		masks[2].setCat(cat);
		masks[3].setCat(cat);
	}

	EList<bool> masks[4];      // temporary list for masks; used in GWState
	EList<uint32_t, 16> map;   // temporary list of GWState maps
};

/**
 * Encapsulates counters that encode how much work the walk-left logic
 * has done.
 */
struct WalkMetrics {

	WalkMetrics() { reset(); MUTEX_INIT(lock); }

	/**
	 * Sum each across this object and 'm'.  This is the only safe way
	 * to update a WalkMetrics shared by many threads.
	 */
	void merge(const WalkMetrics& m, bool getLock = false) {
		ThreadSafe ts(&lock, getLock);
		bwops += m.bwops;
		branches += m.branches;
		resolves += m.resolves;
		refresolves += m.refresolves;
		reports += m.reports;
	}
	
	/**
	 * Set all to 0.
	 */
	void reset() {
		bwops = branches = resolves = refresolves = reports = 0;
	}

	uint64_t bwops;       // Burrows-Wheeler operations
	uint64_t branches;    // BW range branch-offs
	uint64_t resolves;    // # offs resolved with BW walk-left
	uint64_t refresolves; // # resolutions caused by reference scanning
	uint64_t reports;     // # offs reported (1 can be reported many times)
	MUTEX_T lock;
};

/**
 * Coordinates for a BW element that the GroupWalk might resolve.
 */
struct GWElt {

	GWElt() { reset(); }
	
	/**
	 * Reset GWElt to uninitialized state.
	 */
	void reset() {
		offidx = range = elt = len = 0xffffffff;
		fw = false;
	}

	/**
	 * Initialize this WalkResult.
	 */
	void init(
		uint32_t oi,
		bool f,
		uint32_t r,
		uint32_t e,
		uint32_t l)
	{
		offidx = oi;
		fw = f;
		range = r;
		elt = e;
		len = l;
	}

	/**
	 * Return true iff this GWElt and the given GWElt refer to the same
	 * element.
	 */
	bool operator==(const GWElt& o) const {
		return offidx == o.offidx &&
		       fw == o.fw &&
		       range == o.range &&
		       elt == o.elt &&
		       len == o.len;
	}
	
	/**
	 * Return true iff this GWElt and the given GWElt refer to
	 * different elements.
	 */
	bool operator!=(const GWElt& o) const {
		return !(*this == o);
	}

	uint32_t offidx; // seed offset index
	bool     fw;     // strand
	uint32_t range;  // range
	uint32_t elt;    // element
	uint32_t len;    // length
};

/**
 * A record encapsulating the result of looking up one BW element in
 * the Bowtie index.
 */
struct WalkResult {

	WalkResult() { reset(); }
	
	/**
	 * Reset GWElt to uninitialized state.
	 */
	void reset() {
		elt.reset();
		bwrow = toff = 0xffffffff;
	}

	/**
	 * Initialize this WalkResult.
	 */
	void init(
		uint32_t oi,  // seed offset index
		bool f,       // strand
		uint32_t r,   // range
		uint32_t e,   // element
		uint32_t bwr, // BW row
		uint32_t len, // length
		uint32_t to)  // text offset
	{
		elt.init(oi, f, r, e, len);
		bwrow = bwr;
		toff = to;
	}

	GWElt    elt;   // element resolved
	uint32_t bwrow; // SA row resolved
	uint32_t toff;  // resolved offset from SA sample
};

/**
 * A GW hit encapsulates an SATuple describing a reference substring
 * in the cache, along with a bool indicating whether each element of
 * the hit has been reported yet.
 */
template<typename T>
class GWHit {

public:

	GWHit() :
		fmap(0, GW_CAT),
		offidx(0xffffffff),
		fw(false),
		range(0xffffffff),
		len(0xffffffff),
		reported_(0, GW_CAT),
		nrep_(0)
	{
		assert(repOkBasic());
	}

	/**
	 * Initialize with a new SA range.  Resolve the done vector so that
	 * there's one bool per suffix array element.
	 */
	void init(
		SARangeWithOffs<T>& sa,
		uint32_t oi,
		bool f,
		uint32_t r)
	{
		nrep_ = 0;
		offidx = oi;
		fw = f;
		range = r;
		len = (uint32_t)sa.len;
		reported_.resize(sa.offs.size());
		reported_.fill(false);
		fmap.resize(sa.offs.size());
		fmap.fill(make_pair(0xffffffff, 0xffffffff));
	}
	
	/**
	 * Clear contents of sat and done.
	 */
	void reset() {
		reported_.clear();
		fmap.clear();
		nrep_ = 0;
		offidx = 0xffffffff;
		fw = false;
		range = 0xffffffff;
		len = 0xffffffff;
	}
	
#ifndef NDEBUG
	/**
	 * Check that GWHit is internally consistent.  If a pointer to an
	 * EList of GWStates is given, we assume that it is the EList
	 * corresponding to this GWHit and check whether the forward and
	 * reverse mappings match up for the as-yet-unresolved elements.
	 */
	bool repOk(const SARangeWithOffs<T>& sa) const {
		assert_eq(reported_.size(), sa.offs.size());
		assert_eq(fmap.size(), sa.offs.size());
		// Shouldn't be any repeats among as-yet-unresolveds
		size_t nrep = 0;
		for(size_t i = 0; i < fmap.size(); i++) {
			if(reported_[i]) nrep++;
			if(sa.offs[i] != 0xffffffff) {
				continue;
			}
			for(size_t j = i+1; j < fmap.size(); j++) {
				if(sa.offs[j] != 0xffffffff) {
					continue;
				}
				assert(fmap[i] != fmap[j]);
			}
		}
		assert_eq(nrep_, nrep);
		return true;
	}

	/**
	 * Return true iff this GWHit is not obviously corrupt.
	 */
	bool repOkBasic() {
		return true;
	}
#endif
	
	/**
	 * Set the ith element to be reported.
	 */
	void setReported(size_t i) {
		assert(!reported_[i]);
		assert_lt(i, reported_.size());
		reported_[i] = true;
		nrep_++;
	}
	
	/**
	 * Return true iff element i has been reported.
	 */
	bool reported(size_t i) const {
		assert_lt(i, reported_.size());
		return reported_[i];
	}
	
	/**
	 * Return true iff all elements have been reported.
	 */
	bool done() const {
		assert_leq(nrep_, reported_.size());
		return nrep_ == reported_.size();
	}

	EList<std::pair<uint32_t, uint32_t>, 16> fmap; // forward map; to GWState & elt
	uint32_t offidx; // offset idx
	bool fw;         // orientation
	uint32_t range;  // original range index
	uint32_t len;    // length of hit

protected:

	EList<bool, 16> reported_; // per-elt bool indicating whether it's been reported
	size_t nrep_;
};

/**
 * Encapsulates the progress made along a particular path from the original
 * range.
 */
template<typename T>
class GWState {
	
public:

	GWState() : map_(0, GW_CAT) {
		reset(); assert(repOkBasic());
	}
	
	/**
	 * Initialize this GWState with new ebwt, top, bot, step, and sat.
	 *
	 * We assume map is already set up.
	 *
	 * Returns true iff at least one elt was resolved.
	 */
	template<int S>
	pair<int, int> init(
		const Ebwt& ebwt,             // index to walk left in
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T>& sa,       // SA range with offsets
		EList<GWState, S>& sts,       // EList of GWStates for range being advanced
		GWHit<T>& hit,                // Corresponding hit structure
		uint32_t range,               // which range is this?
		bool reportList,              // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult, 16>* res,   // EList where resolved offsets should be appended
		uint32_t tp,                  // top of range at this step
		uint32_t bt,                  // bot of range at this step
		uint32_t st,                  // # steps taken to get to this step
		WalkMetrics& met)
	{
		assert_gt(bt, tp);
		assert_lt(range, sts.size());
		top = tp;
		bot = bt;
		step = st;
		assert(!inited_);
		ASSERT_ONLY(inited_ = true);
		ASSERT_ONLY(lastStep_ = step-1);
		return init(ebwt, ref, sa, sts, hit, range, reportList, res, met);
	}

	/**
	 * Initialize this GWState.
	 *
	 * We assume map is already set up, and that 'step' is equal to the
	 * number of steps taken to get to the new top/bot pair *currently*
	 * in the top and bot fields.
	 *
	 * Returns a pair of numbers, the first being the number of
	 * resolved but unreported offsets found during this advance, the
	 * second being the number of as-yet-unresolved offsets.
	 */
	template<int S>
	pair<int, int> init(
		const Ebwt& ebwt,             // forward Bowtie index
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T>& sa,       // SA range with offsets
		EList<GWState, S>& st,        // EList of GWStates for advancing range
		GWHit<T>& hit,                // Corresponding hit structure
		uint32_t range,               // range being inited
		bool reportList,              // report resolutions, adding to 'res' list?
		EList<WalkResult, 16>* res,   // EList to append resolutions
		WalkMetrics& met)             // update these metrics
	{
		assert(inited_);
		assert_eq(step, lastStep_+1);
		ASSERT_ONLY(lastStep_++);
		assert_leq((uint32_t)step, ebwt.eh().len());
		assert_lt(range, st.size());
		pair<int, int> ret = make_pair(0, 0);
		uint32_t trimBegin = 0, trimEnd = 0;
		bool empty = true; // assume all resolved until proven otherwise
		// Commit new information, if any, to the PListSlide.  Also,
		// trim and check if we're done.
		for(size_t i = mapi_; i < map_.size(); i++) {
			bool resolved = (off(i, sa) != 0xffffffff);
			if(!resolved) {
				// Elt not resolved yet; try to resolve it now
				uint32_t bwrow = (uint32_t)(top - mapi_ + i);
				uint32_t toff = ebwt.tryOffset(bwrow);
				ASSERT_ONLY(uint32_t origBwRow = sa.topf + map(i));
				assert_eq(bwrow, ebwt.walkLeft(origBwRow, step));
				if(toff != 0xffffffff) {
					// Yes, toff was resolvable
					assert_eq(toff, ebwt.getOffset(bwrow));
					met.resolves++;
					toff += step;
					assert_eq(toff, ebwt.getOffset(origBwRow));
					setOff(i, toff, sa, met);
					if(!reportList) ret.first++;
#if 0
// used to be #ifndef NDEBUG, but since we no longer require that the reference
// string info be included, this is no longer relevant.

					// Sanity check that the reference characters under this
					// hit match the seed characters in hit.satup->key.seq.
					// This is NOT a check that we associated the exact right
					// text offset with the BW row.  This is an important
					// distinction because when resolved offsets are filled in
					// via refernce scanning, they are not necessarily the
					// exact right text offsets to associate with the
					// respective BW rows but they WILL all be correct w/r/t
					// the reference sequence underneath, which is what really
					// matters here.
					uint32_t tidx = 0xffffffff, tof, tlen;
					bool straddled = false;
					ebwt.joinedToTextOff(
						hit.len, // length of seed
						toff,    // offset in joined reference string
						tidx,    // reference sequence id
						tof,     // offset in reference coordinates
						tlen,    // length of reference sequence
						true,    // don't reject straddlers
						straddled);
					if(tidx != 0xffffffff &&
					   hit.satup->key.seq != std::numeric_limits<uint64_t>::max())
					{
						// key: 2-bit characters packed into a 64-bit word with
						// the least significant bitpair corresponding to the
						// rightmost character on the Watson reference strand.
						uint64_t key = hit.satup->key.seq;
						for(int64_t j = tof + hit.len-1; j >= tof; j--) {
							// Get next reference base to the left
							int c = ref.getBase(tidx, j);
							assert_range(0, 3, c);
							// Must equal least significant bitpair of key
							if(c != (int)(key & 3)) {
								// Oops; when we jump to the piece of the
								// reference where the seed hit is, it doesn't
								// match the seed hit.  Before dying, check
								// whether we have the right spot in the joined
								// reference string
								SString<char> jref;
								ebwt.restore(jref);
								uint64_t key2 = hit.satup->key.seq;
								for(int64_t k = toff + hit.len-1; k >= toff; k--) {
									int c = jref[k];
									assert_range(0, 3, c);
									assert_eq(c, (int)(key2 & 3));
									key2 >>= 2;
								}
								assert(false);
							}
							key >>= 2;
						}
					}
#endif
				}
			}
			// Is the element resolved?  We ask this regardless of how it was
			// resolved (whether this function did it just now, whether it did
			// it a while ago, or whether some other function outside GroupWalk
			// did it).
			if(off(i, sa) != 0xffffffff) {
				if(reportList && !hit.reported(map(i))) {
					// Report it
					uint32_t toff = off(i, sa);
					assert(res != NULL);
					res->expand();
					uint32_t origBwRow = sa.topf + map(i);
					res->back().init(
						hit.offidx, // offset idx
						hit.fw,     // orientation
						hit.range,  // original range index
						map(i),     // original element offset
						origBwRow,  // BW row resolved
						hit.len,    // hit length
						toff);      // text offset
					hit.setReported(map(i));
					met.reports++;
				}
				// Offset resolved
				if(empty) {
					// Haven't seen a non-empty entry yet, so we
					// can trim this from the beginning.
					trimBegin++;
				} else {
					trimEnd++;
				}
			} else {
				// Offset not yet resolved
				ret.second++;
				trimEnd = 0;
				empty = false;
				// Set the forward map in the corresponding GWHit
				// object to point to the appropriate element of our
				// range
				assert_geq(i, mapi_);
				uint32_t bmap = map(i);
				hit.fmap[bmap].first = range;
				hit.fmap[bmap].second = (uint32_t)i;
#ifndef NDEBUG
				for(size_t j = 0; j < bmap; j++) {
					if(sa.offs[j] == 0xffffffff &&
					   hit.fmap[j].first == range)
					{
						assert_neq(i, hit.fmap[j].second);
					}
				}
#endif
			}
		}
		// Trim from beginning
		assert_geq(trimBegin, 0);
		mapi_ += trimBegin;
		top += trimBegin;
		if(trimEnd > 0) {
			// Trim from end
			map_.resize(map_.size() - trimEnd);
			bot -= trimEnd;
		}
		if(empty) {
			assert(done());
#ifndef NDEBUG
			// If range is done, all elements from map should be
			// resolved
			for(size_t i = mapi_; i < map_.size(); i++) {
				assert_neq(0xffffffff, off(i, sa));
			}
			// If this range is done, then it should be the case that
			// all elements in the corresponding GWHit that point to
			// this range are resolved.
			for(size_t i = 0; i < hit.fmap.size(); i++) {
				if(sa.offs[i] == 0xffffffff) {
					assert_neq(range, hit.fmap[i].first);
				}
			}
#endif
			return ret;
		} else {
			assert(!done());
		}
		// Is there a dollar sign in the middle of the range?
		assert_neq(top, ebwt._zOff);
		assert_neq(bot-1, ebwt._zOff);
		if(ebwt._zOff > top && ebwt._zOff < bot-1) {
			// Yes, the dollar sign is in the middle of this range.  We
			// must split it into the two ranges on either side of the
			// dollar.  Let 'bot' and 'top' delimit the portion of the
			// range prior to the dollar.
			uint32_t oldbot = bot;
			bot = ebwt._zOff;
			// Note: might be able to do additional trimming off the
			// end.
			// Create a new range for the portion after the dollar.
			st.expand();
			st.back().reset();
			uint32_t ztop = ebwt._zOff+1;
			st.back().initMap(oldbot - ztop);
			assert_eq(map_.size(), oldbot-top+mapi_);
			for(size_t i = ztop; i < oldbot; i++) {
				st.back().map_[i - ztop] = map(i-top+mapi_);
			}
			map_.resize(bot - top + mapi_);
			st.back().init(
				ebwt,
				ref,
				sa,
				st,
				hit,
				(uint32_t)st.size()-1,
				reportList,
				res,
				ztop,
				oldbot,
				step,
				met);
		}
		assert_gt(bot, top);
		// Prepare SideLocus's for next step
		if(bot-top > 1) {
			SideLocus::initFromTopBot(top, bot, ebwt.eh(), ebwt.ebwt(), tloc, bloc);
			assert(tloc.valid()); assert(tloc.repOk(ebwt.eh()));
			assert(bloc.valid()); assert(bloc.repOk(ebwt.eh()));
		} else {
			tloc.initFromRow(top, ebwt.eh(), ebwt.ebwt());
			assert(tloc.valid()); assert(tloc.repOk(ebwt.eh()));
			bloc.invalidate();
		}
		return ret;
	}
	
#ifndef NDEBUG
	/**
	 * Check if this GWP is internally consistent.
	 */
	bool repOk(
		const Ebwt& ebwt,
		GWHit<T>& hit,
		uint32_t range) const
	{
		assert(done() || bot > top);
		assert(doneResolving(hit) || (tloc.valid() && tloc.repOk(ebwt.eh())));
		assert(doneResolving(hit) || bot == top+1 || (bloc.valid() && bloc.repOk(ebwt.eh())));
		assert_eq(map_.size()-mapi_, bot-top);
		// Make sure that 'done' is compatible with whether we have >=
		// 1 elements left to resolve.
		int left = 0;
		for(size_t i = mapi_; i < map_.size(); i++) {
			ASSERT_ONLY(uint32_t row = (uint32_t)(top + i - mapi_));
			ASSERT_ONLY(uint32_t origRow = hit.satup->topf + map(i));
			assert(step == 0 || row != origRow);
			assert_eq(row, ebwt.walkLeft(origRow, step));
			assert_lt(map_[i], hit.satup->offs.size());
			if(off(i, hit) == 0xffffffff) left++;
		}
		assert(repOkMapRepeats());
		assert(repOkMapInclusive(hit, range));
		return true;
	}
	
	/**
	 * Return true iff this GWState is not obviously corrupt.
	 */
	bool repOkBasic() {
		assert_geq(bot, top);
		return true;
	}

	/**
	 * Check that the fmap elements pointed to by our map_ include all
	 * of the fmap elements that point to this range.
	 */
	bool repOkMapInclusive(GWHit<T>& hit, uint32_t range) const {
		for(size_t i = 0; i < hit.fmap.size(); i++) {
			if(hit.satup->offs[i] == 0xffffffff) {
				if(range == hit.fmap[i].first) {
					ASSERT_ONLY(bool found = false);
					for(size_t j = mapi_; j < map_.size(); j++) {
						if(map(j) == i) {
							ASSERT_ONLY(found = true);
							break;
						}
					}
					assert(found);
				}
			}
		}
		return true;
	}
	
	/**
	 * Check that no two elements in map_ are the same.
	 */
	bool repOkMapRepeats() const {
		for(size_t i = mapi_; i < map_.size(); i++) {
			for(size_t j = i+1; j < map_.size(); j++) {
				assert_neq(map_[i], map_[j]);
			}
		}
		return true;
	}
#endif
	
	/**
	 * Return the offset currently assigned to the ith element.  If it
	 * has not yet been resolved, return 0xffffffff.
	 */
	uint32_t off(
		size_t i,
		const SARangeWithOffs<T>& sa)
	{
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		assert_lt(map_[i], sa.offs.size());
		return sa.offs.get(map_[i]);
	}

	/**
	 * Return the offset of the element within the original range's
	 * PListSlice that the ith element of this range corresponds to.
	 */
	uint32_t map(size_t i) const {
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		return map_[i];
	}

	/**
	 * Return the offset of the first untrimmed offset in the map.
	 */
	uint32_t mapi() const {
		return mapi_;
	}

	/**
	 * Return number of active elements in the range being tracked by
	 * this GWState.
	 */
	size_t size() const {
		return map_.size() - mapi_;
	}
	
	/**
	 * Return true iff all elements in this leaf range have been
	 * resolved.
	 */
	bool done() const {
		return size() == 0;
	}

	/**
	 * Set the PListSlice element that corresponds to the ith element
	 * of 'map' to the specified offset.
	 */
	void setOff(
		size_t i,
		uint32_t off,
		SARangeWithOffs<T>& sa,
		WalkMetrics& met)
	{
		assert_lt(i + mapi_, map_.size());
		assert_lt(map_[i + mapi_], sa.offs.size());
		size_t saoff = map_[i + mapi_];
		sa.offs[saoff] = off;
		assert_eq(off, sa.offs[saoff]);
	}

	/**
	 * Advance this GWState by one step (i.e. one BW operation).  In
	 * the event of a "split", more elements are added to the EList
	 * 'st', which must have room for at least 3 more elements without
	 * needing another expansion.  If an expansion of 'st' is
	 * triggered, this GWState object becomes invalid.
	 *
	 * Returns a pair of numbers, the first being the number of
	 * resolved but unreported offsets found during this advance, the
	 * second being the number of as-yet-unresolved offsets.
	 */
	template <int S>
	pair<int, int> advance(
		const Ebwt& ebwt,            // the forward Bowtie index, for stepping left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T>& sa,      // SA range with offsets
		GWHit<T>& hit,               // the associated GWHit object
		uint32_t range,              // which range is this?
		bool reportList,             // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult, 16>* res,  // EList where resolved offsets should be appended
		EList<GWState, S>& st,       // EList of GWStates for range being advanced
		GroupWalkState& gws,         // temporary storage for masks
		WalkMetrics& met,
		PerReadMetrics& prm)
	{
		ASSERT_ONLY(uint32_t origTop = top);
		ASSERT_ONLY(uint32_t origBot = bot);
		assert_geq(step, 0);
		assert_eq(step, lastStep_);
		assert_geq(st.capacity(), st.size() + 4);
		assert(tloc.valid()); assert(tloc.repOk(ebwt.eh()));
		assert_eq(bot-top, map_.size()-mapi_);
		pair<int, int> ret = make_pair(0, 0);
		assert_eq(top, tloc.toBWRow());
		if(bloc.valid()) {
			// Still multiple elements being tracked
			assert_lt(top+1, bot);
			uint32_t upto[4], in[4];
			upto[0] = in[0] = upto[1] = in[1] =
			upto[2] = in[2] = upto[3] = in[3] = 0;
			assert_eq(bot, bloc.toBWRow());
			met.bwops++;
			prm.nExFmops++;
			// Assert that there's not a dollar sign in the middle of
			// this range
			assert(bot <= ebwt._zOff || top > ebwt._zOff);
			ebwt.mapLFRange(tloc, bloc, bot-top, upto, in, gws.masks);
#ifndef NDEBUG
			for(int i = 0; i < 4; i++) {
				assert_eq(bot-top, gws.masks[i].size());
			}
#endif
			bool first = true;
			ASSERT_ONLY(uint32_t sum = 0);
			uint32_t newtop = 0, newbot = 0;
			gws.map.clear();
			for(int i = 0; i < 4; i++) {
				if(in[i] > 0) {
					// Non-empty range resulted
					if(first) {
						// For the first one, 
						first = false;
						newtop = upto[i];
						newbot = newtop + in[i];
						assert_leq(newbot-newtop, bot-top);
						// Range narrowed so we have to look at the masks
						for(size_t j = 0; j < gws.masks[i].size(); j++) {
							assert_lt(j+mapi_, map_.size());
							if(gws.masks[i][j]) {
								gws.map.push_back(map_[j+mapi_]);
								assert(gws.map.size() <= 1 || gws.map.back() != gws.map[gws.map.size()-2]);
#ifndef NDEBUG
								// If this element is not yet resolved,
								// then check that it really is the
								// expected number of steps to the left
								// of the corresponding element in the
								// root range
								assert_lt(gws.map.back(), sa.size());
								if(sa.offs[gws.map.back()] == 0xffffffff) {
									assert_eq(newtop + gws.map.size() - 1,
											  ebwt.walkLeft(sa.topf + gws.map.back(), step+1));
								}
#endif
							}
						}
 						assert_eq(newbot-newtop, gws.map.size());
					} else {
						// For each beyond the first, create a new
						// GWState and add it to the GWState list. 
						// NOTE: this can cause the underlying list to
						// be expanded which in turn might leave 'st'
						// pointing to bad memory.
						st.expand();
						st.back().reset();
						uint32_t ntop = upto[i];
						uint32_t nbot = ntop + in[i];
						assert_lt(nbot-ntop, bot-top);
						st.back().mapi_ = 0;
						st.back().map_.clear();
						met.branches++;
						// Range narrowed so we have to look at the masks
						for(size_t j = 0; j < gws.masks[i].size(); j++) {
							if(gws.masks[i][j]) st.back().map_.push_back(map_[j+mapi_]);
						}
						pair<int, int> rret =
						st.back().init(
							ebwt,        // forward Bowtie index
							ref,         // bitpair-encodede reference
							sa,          // SA range with offsets
							st,          // EList of all GWStates associated with original range
							hit,         // associated GWHit object
							(uint32_t)st.size()-1, // range offset
							reportList,  // if true, report hits to 'res' list
							res,         // report hits here if reportList is true
							ntop,        // BW top of new range
							nbot,        // BW bot of new range
							step+1,      // # steps taken to get to this new range
							met);        // update these metrics
						ret.first += rret.first;
						ret.second += rret.second;
					}
					ASSERT_ONLY(sum += in[i]);
				}
			}
			mapi_ = 0;
			assert_eq(bot-top, sum);
			assert_gt(newbot, newtop);
			assert_leq(newbot-newtop, bot-top);
			assert(top != newtop || bot != newbot);
			//assert(!(newtop < top && newbot > top));
			top = newtop;
			bot = newbot;
			if(!gws.map.empty()) {
				map_ = gws.map;
			}
			//assert(repOkMapRepeats());
			//assert(repOkMapInclusive(hit, range));
			assert_eq(bot-top, map_.size());
		} else {
			// Down to one element
			assert_eq(bot, top+1);
			assert_eq(1, map_.size()-mapi_);
			// Sets top, returns char walked through (which we ignore)
			ASSERT_ONLY(uint32_t oldtop = top);
			met.bwops++;
			prm.nExFmops++;
			ebwt.mapLF1(top, tloc);
			assert_neq(top, oldtop);
			bot = top+1;
			if(mapi_ > 0) {
				map_[0] = map_[mapi_];
				mapi_ = 0;
			}
			map_.resize(1);
		}
		assert(top != origTop || bot != origBot);
		step++;
		assert_gt(step, 0);
		assert_leq((uint32_t)step, ebwt.eh().len());
		pair<int, int> rret =
		init<S>(
			ebwt,       // forward Bowtie index
			ref,        // bitpair-encodede reference
			sa,         // SA range with offsets
			st,         // EList of all GWStates associated with original range
			hit,        // associated GWHit object
			range,      // range offset
			reportList, // if true, report hits to 'res' list
			res,        // report hits here if reportList is true
			met);       // update these metrics
		ret.first += rret.first;
		ret.second += rret.second;
		return ret;
	}

	/**
	 * Clear all state in preparation for the next walk.
	 */
	void reset() {
		top = bot = step = mapi_ = 0;
		ASSERT_ONLY(lastStep_ = -1);
		ASSERT_ONLY(inited_ = false);
		tloc.invalidate();
		bloc.invalidate();
		map_.clear();
	}
	
	/**
	 * Resize the map_ field to the given size.
	 */
	void initMap(size_t newsz) {
		mapi_ = 0;
		map_.resize(newsz);
		for(size_t i = 0; i < newsz; i++) {
			map_[i] = (uint32_t)i;
		}
	}

	/**
	 * Return true iff all rows corresponding to this GWState have been
	 * resolved and reported.
	 */
	bool doneReporting(const GWHit<T>& hit) const {
		for(size_t i = mapi_; i < map_.size(); i++) {
			if(!hit.reported(map(i))) return false;
		}
		return true;
	}

	/**
	 * Return true iff all rows corresponding to this GWState have been
	 * resolved (but not necessarily reported).
	 */
	bool doneResolving(const SARangeWithOffs<T>& sa) const {
		for(size_t i = mapi_; i < map_.size(); i++) {
			if(sa.offs[map(i)] == 0xffffffff) return false;
		}
		return true;
	}

	SideLocus tloc;      // SideLocus for top
	SideLocus bloc;      // SideLocus for bottom
	uint32_t  top;       // top elt of range in BWT
	uint32_t  bot;       // bot elt of range in BWT
	int       step;      // how many steps have we walked to the left so far

protected:
	
	ASSERT_ONLY(bool inited_);
	ASSERT_ONLY(int lastStep_);
	EList<uint32_t, 16> map_; // which elts in range 'range' we're tracking
	uint32_t mapi_;           // first untrimmed element of map
};

template<typename T, int S>
class GroupWalk2S {
public:
	typedef EList<GWState<T>, S> TStateV;

	GroupWalk2S() : st_(8, GW_CAT) {
		reset();
	}
	
	/**
	 * Reset the GroupWalk in preparation for the next SeedResults.
	 */
	void reset() {
		elt_ = rep_ = 0;
		ASSERT_ONLY(inited_ = false);
	}

	/**
	 * Initialize a new group walk w/r/t a QVal object.
	 */
	void init(
		const Ebwt& ebwtFw,         // forward Bowtie index for walking left
		const BitPairReference& ref,// bitpair-encoded reference
		SARangeWithOffs<T>& sa,     // SA range with offsets
		RandomSource& rnd,          // pseudo-random generator for sampling rows
		WalkMetrics& met)           // update metrics here
	{
		reset();
#ifndef NDEBUG
		inited_ = true;
#endif
		// Init GWHit
		hit_.init(sa, 0, false, 0);
		// Init corresponding GWState
		st_.resize(1);
		st_.back().reset();
		assert(st_.back().repOkBasic());
		uint32_t top = sa.topf;
		uint32_t bot = (uint32_t)(top + sa.size());
		st_.back().initMap(bot-top);
		st_.ensure(4);
		st_.back().init(
			ebwtFw,             // Bowtie index
			ref,                // bitpair-encoded reference
			sa,                 // SA range with offsets
			st_,                // EList<GWState>
			hit_,               // GWHit
			0,                  // range 0
			false,              // put resolved elements into res_?
			NULL,               // put resolved elements here
			top,                // BW row at top
			bot,                // BW row at bot
			0,                  // # steps taken
			met);               // update metrics here
		elt_ += sa.size();
		assert(hit_.repOk(sa));
	}

	//
	// ELEMENT-BASED
	//

	/**
	 * Advance the GroupWalk until all elements have been resolved.
	 */
	void resolveAll(WalkMetrics& met, PerReadMetrics& prm) {
		WalkResult res; // ignore results for now
		for(size_t i = 0; i < elt_; i++) {
			advanceElement((uint32_t)i, res, met, prm);
		}
	}

	/**
	 * Advance the GroupWalk until the specified element has been
	 * resolved.
	 */
	bool advanceElement(
		uint32_t elt,                // element within the range
		const Ebwt& ebwtFw,          // forward Bowtie index for walking left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T>& sa,      // SA range with offsets
		GroupWalkState& gws,         // GroupWalk state; scratch space
		WalkResult& res,             // put the result here
		WalkMetrics& met,            // metrics
		PerReadMetrics& prm)         // per-read metrics
	{
		assert(inited_);
		assert(!done());
		assert(hit_.repOk(sa));
		assert_lt(elt, sa.size()); // elt must fall within range
		// Until we've resolved our element of interest...
		while(sa.offs[elt] == 0xffffffff) {
			// Get the GWState that contains our element of interest
			size_t range = hit_.fmap[elt].first;
			st_.ensure(4);
			GWState<T>& st = st_[range];
			assert(!st.doneResolving(sa));
			// Returns a pair of numbers, the first being the number of
			// resolved but unreported offsets found during this advance, the
			// second being the number of as-yet-unresolved offsets.
			st.advance(
				ebwtFw,
				ref,
				sa,
				hit_,
				(uint32_t)range,
				false,
				NULL,
				st_,
				gws,
				met,
				prm);
			assert(sa.offs[elt] != 0xffffffff ||
			       !st_[hit_.fmap[elt].first].doneResolving(sa));
		}
		assert_neq(0xffffffff, sa.offs[elt]);
		// Report it!
		if(!hit_.reported(elt)) {
			hit_.setReported(elt);
		}
		met.reports++;
		res.init(
			0,              // seed offset
			false,          // orientation
			0,              // range
			elt,            // element
			sa.topf + elt,  // bw row
			(uint32_t)sa.len, // length of hit
			sa.offs[elt]);  // resolved text offset
		rep_++;
		return true;
	}

	/**
	 * Return true iff all elements have been resolved and reported.
	 */
	bool done() const { return rep_ == elt_; }
	
#ifndef NDEBUG
	/**
	 * Check that GroupWalk is internally consistent.
	 */
	bool repOk(const SARangeWithOffs<T>& sa) const {
		assert(hit_.repOk(sa));
		assert_leq(rep_, elt_);
		// This is a lot of work
		size_t resolved = 0, reported = 0;
		// For each element
		const size_t sz = sa.size();
		for(size_t m = 0; m < sz; m++) {
			// Is it resolved?
			if(sa.offs[m] != 0xffffffff) {
				resolved++;
			} else {
				assert(!hit_.reported(m));
			}
			// Is it reported?
			if(hit_.reported(m)) {
				reported++;
			}
			assert_geq(resolved, reported);
		}
		assert_geq(resolved, reported);
		assert_eq(rep_, reported);
		assert_eq(elt_, sz);
		return true;
	}
#endif

	/**
	 * Return the number of BW elements that we can resolve.
	 */
	size_t numElts() const { return elt_; }
	
	/**
	 * Return the size occupied by this GroupWalk and all its constituent
	 * objects.
	 */
	size_t totalSizeBytes() const {
		return 2 * sizeof(size_t) + st_.totalSizeBytes() + sizeof(GWHit<T>);
	}
	/**
	 * Return the capacity of this GroupWalk and all its constituent objects.
	 */
	size_t totalCapacityBytes() const {
		return 2 * sizeof(size_t) + st_.totalCapacityBytes() + sizeof(GWHit<T>);
	}
	
#ifndef NDEBUG
	bool initialized() const { return inited_; }
#endif
	
protected:

	ASSERT_ONLY(bool inited_);    // initialized?
	
	size_t elt_;    // # BW elements under the control of the GropuWalk
	size_t rep_;    // # BW elements reported

	// For each orientation and seed offset, keep a GWState object that
	// holds the state of the walk so far.
	TStateV st_;

	// For each orientation and seed offset, keep an EList of GWHit.
	GWHit<T> hit_;
};

#endif /*GROUP_WALK_H_*/
