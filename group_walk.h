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

/**
 * Encapsulate an SA range and an associated list of slots where the resolved
 * offsets can be placed.
 */
template<typename T>
class SARangeWithOffs {

public:

	SARangeWithOffs() { reset(); };

	SARangeWithOffs(TIndexOffU tf, size_t len, const T& o) {
		init(tf, len, o);
	}

	void init(TIndexOffU tf, const T& o) {
		topf = tf, offs = o;
	}

	/**
	 * Reset to uninitialized state.
	 */
	void reset() { topf = std::numeric_limits<TIndexOffU>::max(); }

	/**
	 * Return true if this is initialized.
	 */
	bool inited() const {
		return topf != std::numeric_limits<TIndexOffU>::max();
	}

	/**
	 * Return the number of times this reference substring occurs in the
	 * reference, which is also the size of the 'offs' TSlice.
	 */
	size_t size() const { return offs.size(); }

	TIndexOffU topf; // top in BWT index
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
	EList<TIndexOffU, 16> map;   // temporary list of GWState maps
};

/**
 * Encapsulates counters that encode how much work the walk-left logic
 * has done.
 */
struct WalkMetrics {

	WalkMetrics() {
		reset();
	}

	/**
	 * Sum each across this object and 'm'.  This is the only safe way
	 * to update a WalkMetrics shared by many threads.
	 */
	void merge(const WalkMetrics& m, bool getLock = false) {
		if(getLock) {
			ThreadSafe ts(mutex_m);
			mergeImpl(m);
		} else {
			mergeImpl(m);
		}
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
	MUTEX_T mutex_m;

private:
	/**
	 * Sum each across this object and 'm'.  This is the only safe way
	 * to update a WalkMetrics shared by many threads.
	 */
	void mergeImpl(const WalkMetrics& m) {
		bwops += m.bwops;
		branches += m.branches;
		resolves += m.resolves;
		refresolves += m.refresolves;
		reports += m.reports;
	}
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
		offidx = range = elt = len = OFF_MASK;
		fw = false;
	}

	/**
	 * Initialize this WalkResult.
	 */
	void init(
		TIndexOffU oi,
		bool f,
		TIndexOffU r,
		TIndexOffU e,
		TIndexOffU l)
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

	TIndexOffU offidx; // seed offset index
	bool     fw;     // strand
	TIndexOffU range;  // range
	TIndexOffU elt;    // element
	TIndexOffU len;    // length
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
		bwrow = toff = OFF_MASK;
	}

	/**
	 * Initialize this WalkResult.
	 */
	void init(
		TIndexOffU oi,  // seed offset index
		bool f,       // strand
		TIndexOffU r,   // range
		TIndexOffU e,   // element
		TIndexOffU bwr, // BW row
		TIndexOffU len, // length
		TIndexOffU to)  // text offset
	{
		elt.init(oi, f, r, e, len);
		bwrow = bwr;
		toff = to;
	}

	GWElt    elt;   // element resolved
	TIndexOffU bwrow; // SA row resolved
	TIndexOffU toff;  // resolved offset from SA sample
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
		offidx(OFF_MASK),
		fw(false),
		range(OFF_MASK),
		len(OFF_MASK),
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
		TIndexOffU oi,
		bool f,
		TIndexOffU r)
	{
		nrep_ = 0;
		offidx = oi;
		fw = f;
		range = r;
		len = (TIndexOffU)sa.len;
		reported_.resize(sa.offs.size());
		reported_.fill(false);
		fmap.resize(sa.offs.size());
		fmap.fill(make_pair(OFF_MASK, OFF_MASK));
	}

	/**
	 * Clear contents of sat and done.
	 */
	void reset() {
		reported_.clear();
		fmap.clear();
		nrep_ = 0;
		offidx = OFF_MASK;
		fw = false;
		range = OFF_MASK;
		len = OFF_MASK;
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
			if(sa.offs[i] != OFF_MASK) {
				continue;
			}
			for(size_t j = i+1; j < fmap.size(); j++) {
				if(sa.offs[j] != OFF_MASK) {
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

	EList<std::pair<TIndexOffU, TIndexOffU>, 16> fmap; // forward map; to GWState & elt
	TIndexOffU offidx; // offset idx
	bool fw;         // orientation
	TIndexOffU range;  // original range index
	TIndexOffU len;    // length of hit

protected:

	EList<bool, 16> reported_; // per-elt bool indicating whether it's been reported
	size_t nrep_;
};

// kSpecMask[r]: the low r bit-pairs of a 64-bit word
static const uint64_t kSpecMask[33] = {
	0x0ull, 0x3ull, 0xfull, 0x3full, 0xffull, 0x3ffull, 0xfffull, 0x3fffull,
	0xffffull, 0x3ffffull, 0xfffffull, 0x3fffffull, 0xffffffull, 0x3ffffffull,
	0xfffffffull, 0x3fffffffull, 0xffffffffull, 0x3ffffffffull, 0xfffffffffull,
	0x3fffffffffull, 0xffffffffffull, 0x3ffffffffffull, 0xfffffffffffull,
	0x3fffffffffffull, 0xffffffffffffull, 0x3ffffffffffffull, 0xfffffffffffffull,
	0x3fffffffffffffull, 0xffffffffffffffull, 0x3ffffffffffffffull,
	0xfffffffffffffffull, 0x3fffffffffffffffull, 0xffffffffffffffffull };

#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
// GWState::advance() steps a range of several rows a 64-bit word of the BWT
// (32 characters) at a time, reading the side's bytes as little-endian words
#define GW_RANGE_WORDS 1
#else
#define GW_RANGE_WORDS 0
#endif

/**
 * Bit 2k of the result is set iff the kth character (bit-pair) of 'dw' is
 * 'c' (as countInU64() counts them).
 */
static inline uint64_t gwCharBits(uint64_t dw, int c) {
	const uint64_t x0 = dw ^ c_table[c];
	return x0 & (x0 >> 1) & 0x5555555555555555ull;
}

/**
 * Call f(dw, rm, j, lo) for each 64-bit word of the BWT holding the last
 * characters of 'n' consecutive rows, the first of which is character 'k' of
 * the side at 'side': dw is the word, rm the bits (as gwCharBits() sets them)
 * of the characters that belong to the rows, j the index (from 0) of the row
 * of the first of them, and lo its character position in the word.  Sides
 * hold 48*OFF_SIZE characters, as SideLocus::initFromRow() assumes.
 */
template<typename F>
static inline void gwRangeWords(
	const uint8_t* side,
	int k,
	TIndexOffU n,
	size_t sideSz,
	F f)
{
	const int CPS = 48 * OFF_SIZE;
	TIndexOffU j = 0;
	while(j < n) {
		const int lo = k & 31;
		const TIndexOffU left = n - j;
		const int take = (TIndexOffU)(32 - lo) < left ? 32 - lo : (int)left;
		uint64_t dw;
		memcpy(&dw, side + 8 * (k >> 5), 8);
		f(dw, (kSpecMask[lo + take] ^ kSpecMask[lo]) & 0x5555555555555555ull, j, lo);
		j += take;
		k += take;
		if(k == CPS) {
			k = 0;
			side += sideSz;
		}
	}
}

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
	pair<TIndexOff, TIndexOff> init(
		const Ebwt& ebwt,             // index to walk left in
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T>& sa,       // SA range with offsets
		EList<GWState, S>& sts,       // EList of GWStates for range being advanced
		GWHit<T>& hit,                // Corresponding hit structure
		TIndexOffU range,               // which range is this?
		bool reportList,              // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult, 16>* res,   // EList where resolved offsets should be appended
		TIndexOffU tp,                  // top of range at this step
		TIndexOffU bt,                  // bot of range at this step
		TIndexOffU st,                  // # steps taken to get to this step
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
	pair<TIndexOff, TIndexOff> init(
		const Ebwt& ebwt,             // forward Bowtie index
		const BitPairReference& ref,  // bitpair-encoded reference
		SARangeWithOffs<T>& sa,       // SA range with offsets
		EList<GWState, S>& st,        // EList of GWStates for advancing range
		GWHit<T>& hit,                // Corresponding hit structure
		TIndexOffU range,               // range being inited
		bool reportList,              // report resolutions, adding to 'res' list?
		EList<WalkResult, 16>* res,   // EList to append resolutions
		WalkMetrics& met)             // update these metrics
	{
		assert(inited_);
		assert_eq(step, lastStep_+1);
		ASSERT_ONLY(lastStep_++);
		assert_leq((TIndexOffU)step, ebwt.eh().len());
		assert_lt(range, st.size());
		pair<TIndexOff, TIndexOff> ret = make_pair(0, 0);
		TIndexOffU trimBegin = 0, trimEnd = 0;
		bool empty = true; // assume all resolved until proven otherwise
		// Commit new information, if any, to the PListSlide.  Also,
		// trim and check if we're done.  (initScan() does this loop's work
		// when it can, faster.)
		if(!initScan(ebwt, sa, hit, range, reportList, ret, trimBegin, trimEnd, empty, met))
		for(size_t i = mapi_; i < map_.size(); i++) {
			bool resolved = (off(i, sa) != OFF_MASK);
			if(!resolved) {
				// Elt not resolved yet; try to resolve it now
				TIndexOffU bwrow = (TIndexOff)(top - mapi_ + i);
				TIndexOffU toff = ebwt.tryOffset(bwrow);
				ASSERT_ONLY(TIndexOffU origBwRow = sa.topf + map(i));
				assert_eq(bwrow, ebwt.walkLeft(origBwRow, step));
				if(toff != OFF_MASK) {
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
					TIndexOffU tidx = OFF_MASK, tof, tlen;
					bool straddled = false;
					ebwt.joinedToTextOff(
						hit.len, // length of seed
						toff,    // offset in joined reference string
						tidx,    // reference sequence id
						tof,     // offset in reference coordinates
						tlen,    // length of reference sequence
						true,    // don't reject straddlers
						straddled);
					if(tidx != OFF_MASK &&
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
			if(off(i, sa) != OFF_MASK) {
				if(reportList && !hit.reported(map(i))) {
					// Report it
					TIndexOffU toff = off(i, sa);
					assert(res != NULL);
					res->expand();
					TIndexOffU origBwRow = sa.topf + map(i);
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
				TIndexOffU bmap = map(i);
				hit.fmap[bmap].first = range;
				hit.fmap[bmap].second = (TIndexOffU)i;
#ifndef NDEBUG
				for(size_t j = 0; j < bmap; j++) {
					if(sa.offs[j] == OFF_MASK &&
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
				assert_neq(OFF_MASK, off(i, sa));
			}
			// If this range is done, then it should be the case that
			// all elements in the corresponding GWHit that point to
			// this range are resolved.
			for(size_t i = 0; i < hit.fmap.size(); i++) {
				if(sa.offs[i] == OFF_MASK) {
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
			TIndexOffU oldbot = bot;
			bot = ebwt._zOff;
			// Note: might be able to do additional trimming off the
			// end.
			// Create a new range for the portion after the dollar.
			st.expand();
			st.back().reset();
			TIndexOffU ztop = ebwt._zOff+1;
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
				(TIndexOffU)st.size()-1,
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

	/**
	 * What the loop of init() over the elements does, for when nothing is
	 * reported as it is resolved (reportList false), mapi_ is 0 and the
	 * range's offset slots are contiguous in memory: the same offsets
	 * resolved, forward-map entries written, counts and trims, with the
	 * lists' pointers and the index's fields held in locals (the build's
	 * -fno-strict-aliasing otherwise reloads them after every store).  Only
	 * the rows whose offsets tryOffset() knows are visited to resolve; the
	 * elements are independent, so the order of the writes does not matter.
	 * Returns false, having done nothing, when those conditions do not hold.
	 */
	bool initScan(
		const Ebwt& ebwt,             // forward Bowtie index
		SARangeWithOffs<T>& sa,       // SA range with offsets
		GWHit<T>& hit,                // Corresponding hit structure
		TIndexOffU range,             // range being inited
		bool reportList,              // report resolutions, adding to 'res' list?
		pair<TIndexOff, TIndexOff>& ret,
		TIndexOffU& trimBegin,
		TIndexOffU& trimEnd,
		bool& empty,
		WalkMetrics& met)             // update these metrics
	{
		if(reportList || mapi_ != 0 || map_.empty()) {
			return false;
		}
		TIndexOffU *offs = sa.offs.contiguous();
		if(offs == NULL) {
			return false;
		}
		const TIndexOffU *mp = map_.ptr();
		const size_t sz = map_.size();
		std::pair<TIndexOffU, TIndexOffU> *fmap = hit.fmap.ptr();
		const TIndexOffU zOff = ebwt._zOff;
		const TIndexOffU offMask = ebwt.eh()._offMask;
		const int offRate = ebwt.eh()._offRate;
		const TIndexOffU *saOffs = ebwt.offs();
		const TIndexOffU tp = top;
		const TIndexOffU stp = (TIndexOffU)step;
		TIndexOff nres = 0, nleft = 0;
		TIndexOffU tb = 0, te = 0;
		bool emp = true;
		if(sz < 32) {
			// Few rows: one pass, as init()'s loop
			for(size_t i = 0; i < sz; i++) {
				const TIndexOffU e = mp[i];
				assert_lt(e, sa.offs.size());
				TIndexOffU o = offs[e];
				if(o == OFF_MASK) {
					// Elt not resolved yet; try to resolve it now, as
					// Ebwt::tryOffset() does
					const TIndexOffU bwrow = (TIndexOffU)(TIndexOff)(tp + i);
					assert_eq(bwrow, ebwt.walkLeft(sa.topf + e, step));
					if(bwrow == zOff || (bwrow & offMask) == bwrow) {
						const TIndexOffU toff = (bwrow == zOff) ? 0 : saOffs[bwrow >> offRate];
						assert_eq(toff, ebwt.tryOffset(bwrow));
						if(toff != OFF_MASK) {
							o = toff + stp;
							assert_eq(o, ebwt.getOffset(sa.topf + e));
							offs[e] = o;
							nres++;
						}
					}
				}
				if(o != OFF_MASK) {
					// Offset resolved
					if(emp) {
						tb++;
					} else {
						te++;
					}
				} else {
					// Offset not yet resolved
					nleft++;
					te = 0;
					emp = false;
					fmap[e].first = range;
					fmap[e].second = (TIndexOffU)i;
				}
			}
		} else {
			const uint64_t tend = (uint64_t)tp + sz;
			// Resolve the unresolved elements whose rows' offsets are
			// known, as Ebwt::tryOffset() finds them: the '$' row's, and
			// those of the rows that are multiples of 2^offRate (sampled)
			const uint64_t per = (uint64_t)1 << offRate;
			assert_eq((TIndexOffU)~(per - 1), offMask);
			for(uint64_t r = ((uint64_t)tp + per - 1) & ~(per - 1); r < tend; r += per) {
				if((TIndexOffU)r == zOff) {
					continue;
				}
				const TIndexOffU e = mp[r - tp];
				if(offs[e] == OFF_MASK) {
					const TIndexOffU toff = saOffs[r >> offRate];
					assert_eq(toff, ebwt.tryOffset((TIndexOffU)r));
					if(toff != OFF_MASK) {
						offs[e] = toff + stp;
						assert_eq(offs[e], ebwt.getOffset(sa.topf + e));
						nres++;
					}
				}
			}
			if((uint64_t)zOff >= tp && (uint64_t)zOff < tend) {
				const TIndexOffU e = mp[zOff - tp];
				if(offs[e] == OFF_MASK) {
					offs[e] = 0 + stp;
					nres++;
				}
			}
			// Point the forward map at each element still unresolved; trim
			// the resolved ones before the first and after the last of
			// those.  The unresolved are listed a block of rows at a time
			// without branching on them (which rows are resolved follows no
			// pattern, so such branches are mispredicted).
			size_t first = 0, last = 0;
			const size_t BLK = 256;
			TIndexOffU blk[BLK];
			for(size_t i0 = 0; i0 < sz; i0 += BLK) {
				const size_t i1 = (sz - i0 < BLK) ? sz : i0 + BLK;
				size_t k = 0;
				for(size_t i = i0; i < i1; i++) {
					blk[k] = (TIndexOffU)i;
					k += (offs[mp[i]] == OFF_MASK) ? 1 : 0;
				}
				if(k > 0) {
					if(nleft == 0) {
						first = blk[0];
					}
					last = blk[k-1];
					for(size_t j = 0; j < k; j++) {
						const TIndexOffU i = blk[j];
						assert_eq(tp + i, ebwt.walkLeft(sa.topf + mp[i], step));
						std::pair<TIndexOffU, TIndexOffU>& f = fmap[mp[i]];
						f.first = range;
						f.second = i;
					}
					nleft += (TIndexOff)k;
				}
			}
			if(nleft > 0) {
				tb = (TIndexOffU)first;
				te = (TIndexOffU)(sz - 1 - last);
				emp = false;
			} else {
				tb = (TIndexOffU)sz;
			}
		}
		met.resolves += (uint64_t)nres;
		ret.first += nres;
		ret.second += nleft;
		trimBegin = tb;
		trimEnd = te;
		empty = emp;
		return true;
	}

#ifndef NDEBUG
	/**
	 * Check if this GWP is internally consistent.
	 */
	bool repOk(
		const Ebwt& ebwt,
		GWHit<T>& hit,
		TIndexOffU range) const
	{
		assert(done() || bot > top);
		assert(doneResolving(hit) || (tloc.valid() && tloc.repOk(ebwt.eh())));
		assert(doneResolving(hit) || bot == top+1 || (bloc.valid() && bloc.repOk(ebwt.eh())));
		assert_eq(map_.size()-mapi_, bot-top);
		// Make sure that 'done' is compatible with whether we have >=
		// 1 elements left to resolve.
		TIndexOff left = 0;
		for(size_t i = mapi_; i < map_.size(); i++) {
			ASSERT_ONLY(TIndexOffU row = (TIndexOffU)(top + i - mapi_));
			ASSERT_ONLY(TIndexOffU origRow = hit.satup->topf + map(i));
			assert(step == 0 || row != origRow);
			assert_eq(row, ebwt.walkLeft(origRow, step));
			assert_lt(map_[i], hit.satup->offs.size());
			if(off(i, hit) == OFF_MASK) left++;
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
	bool repOkMapInclusive(GWHit<T>& hit, TIndexOffU range) const {
		for(size_t i = 0; i < hit.fmap.size(); i++) {
			if(hit.satup->offs[i] == OFF_MASK) {
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
	TIndexOffU off(
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
	TIndexOffU map(size_t i) const {
		assert_geq(i, mapi_);
		assert_lt(i, map_.size());
		return map_[i];
	}

	/**
	 * Return the offset of the first untrimmed offset in the map.
	 */
	TIndexOffU mapi() const {
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
		TIndexOffU off,
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
	pair<TIndexOff, TIndexOff> advance(
		const Ebwt& ebwt,            // the forward Bowtie index, for stepping left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T>& sa,      // SA range with offsets
		GWHit<T>& hit,               // the associated GWHit object
		TIndexOffU range,              // which range is this?
		bool reportList,             // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult, 16>* res,  // EList where resolved offsets should be appended
		EList<GWState, S>& st,       // EList of GWStates for range being advanced
		GroupWalkState& gws,         // temporary storage for masks
		WalkMetrics& met,
		PerReadMetrics& prm)
	{
		ASSERT_ONLY(TIndexOffU origTop = top);
		ASSERT_ONLY(TIndexOffU origBot = bot);
		assert_geq(step, 0);
		assert_eq(step, lastStep_);
		assert_geq(st.capacity(), st.size() + 4);
		assert(tloc.valid()); assert(tloc.repOk(ebwt.eh()));
		assert_eq(bot-top, map_.size()-mapi_);
		pair<TIndexOff, TIndexOff> ret = make_pair(0, 0);
		assert_eq(top, tloc.toBWRow());
		if(bloc.valid()) {
			// Still multiple elements being tracked
			assert_lt(top+1, bot);
			assert_eq(bot, bloc.toBWRow());
			met.bwops++;
			prm.nExFmops++;
			// Assert that there's not a dollar sign in the middle of
			// this range
			assert(bot <= ebwt._zOff || top > ebwt._zOff);
#if GW_RANGE_WORDS
			(void)gws;
			stepRange<S>(ebwt, ref, sa, hit, reportList, res, st, met, ret);
#else
			TIndexOffU upto[4], in[4];
			upto[0] = in[0] = upto[1] = in[1] =
			upto[2] = in[2] = upto[3] = in[3] = 0;
			ebwt.mapLFRange(tloc, bloc, bot-top, upto, in, gws.masks);
#ifndef NDEBUG
			for(int i = 0; i < 4; i++) {
				assert_eq(bot-top, gws.masks[i].size());
			}
#endif
			bool first = true;
			ASSERT_ONLY(TIndexOffU sum = 0);
			TIndexOffU newtop = 0, newbot = 0;
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
								if(sa.offs[gws.map.back()] == OFF_MASK) {
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
						TIndexOffU ntop = upto[i];
						TIndexOffU nbot = ntop + in[i];
						assert_lt(nbot-ntop, bot-top);
						st.back().mapi_ = 0;
						st.back().map_.clear();
						met.branches++;
						// Range narrowed so we have to look at the masks
						for(size_t j = 0; j < gws.masks[i].size(); j++) {
							if(gws.masks[i][j]) st.back().map_.push_back(map_[j+mapi_]);
						}
						pair<TIndexOff, TIndexOff> rret =
						st.back().init(
							ebwt,        // forward Bowtie index
							ref,         // bitpair-encodede reference
							sa,          // SA range with offsets
							st,          // EList of all GWStates associated with original range
							hit,         // associated GWHit object
							(TIndexOffU)st.size()-1, // range offset
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
#endif
		} else {
			// Down to one element
			assert_eq(bot, top+1);
			assert_eq(1, map_.size()-mapi_);
			// Sets top, returns char walked through (which we ignore)
			ASSERT_ONLY(TIndexOffU oldtop = top);
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
		assert_leq((TIndexOffU)step, ebwt.eh().len());
		pair<TIndexOff, TIndexOff> rret =
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

#if GW_RANGE_WORDS
	/**
	 * For the 'n' rows starting at the row at locus 'l': 'upto' gets, for
	 * each character, the row the first of them that ends in the character
	 * moves to with an LF step, and 'in' how many end in the character;
	 * what Ebwt::mapLFRange() puts in its cntsUpto and cntsIn.  Counted a
	 * 64-bit word of the side at a time.
	 */
	static void countRange(
		const Ebwt& ebwt,
		const SideLocus& l,
		TIndexOffU n,
		TIndexOffU* upto,
		TIndexOffU* in)
	{
		assert_gt(n, 0);
		const EbwtParams& eh = ebwt.eh();
		const uint8_t *side = l.side(ebwt.ebwt());
		const int co = (int)l._charOff;
		const size_t sideSz = (size_t)eh._sideSz;
		// Occurrences of each character before the top row in its side, then
		// as Ebwt::countBt2SideRange() completes them
		upto[0] = upto[1] = upto[2] = upto[3] = 0;
		const int NW = (48 * OFF_SIZE) / 32;
		for(int w = 0; w < NW; w++) {
			uint64_t dw;
			memcpy(&dw, side + 8 * w, 8);
			int r = co - 32 * w;
			r = r < 0 ? 0 : r;
			r = r > 32 ? 32 : r;
			const uint64_t m = kSpecMask[r];
			for(int c = 0; c < 3; c++) {
				upto[c] += (TIndexOffU)__builtin_popcountll(gwCharBits(dw, c) & m);
			}
		}
		upto[3] = (TIndexOffU)co - upto[0] - upto[1] - upto[2];
		if(l._sideByteOff <= ebwt._zEbwtByteOff &&
		   l._sideByteOff + l._by >= ebwt._zEbwtByteOff)
		{
			// '$' is represented as an 'A' (as in Ebwt::countBt2SideRange())
			if((l._sideByteOff + l._by > ebwt._zEbwtByteOff) ||
			   (l._sideByteOff + l._by == ebwt._zEbwtByteOff && l._bp > ebwt._zEbwtBpOff))
			{
				upto[0]--;
			}
		}
		const TIndexOffU *acgt = reinterpret_cast<const TIndexOffU*>(side + eh._sideBwtSz);
		for(int c = 0; c < 4; c++) {
			upto[c] += endianizeU<TIndexOffU>(acgt[c], ebwt._switchEndian) + ebwt.fchr()[c];
		}
		// Occurrences of each character in the range
		in[0] = in[1] = in[2] = in[3] = 0;
		gwRangeWords(side, co, n, sideSz,
			[&](uint64_t dw, uint64_t rm, TIndexOffU, int) {
				for(int c = 0; c < 3; c++) {
					in[c] += (TIndexOffU)__builtin_popcountll(gwCharBits(dw, c) & rm);
				}
			});
		in[3] = n - in[0] - in[1] - in[2];
	}

	/**
	 * The step of advance() for a range of several rows (bloc valid): what
	 * mapLFRange() and the loops over its masks do, a word of the BWT at a
	 * time.  The rows whose last character is the first (lowest) character
	 * present stay in this state; those of each other character present, in
	 * increasing order, go to a new state appended to 'st' and initialized
	 * as soon as it is made, as advance() makes them.  Each list keeps the
	 * order of the rows.  top, bot, mapi_ and map_ are left as advance()
	 * leaves them before its final init().  The range never holds the '$'
	 * row (init() trims or splits it off).
	 */
	template <int S>
	void stepRange(
		const Ebwt& ebwt,            // the forward Bowtie index, for stepping left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T>& sa,      // SA range with offsets
		GWHit<T>& hit,               // the associated GWHit object
		bool reportList,             // if true, "report" resolved offsets immediately by adding them to 'res' list
		EList<WalkResult, 16>* res,  // EList where resolved offsets should be appended
		EList<GWState, S>& st,       // EList of GWStates for range being advanced
		WalkMetrics& met,
		pair<TIndexOff, TIndexOff>& ret)
	{
		const TIndexOffU n = bot - top;
		assert_gt(n, 1);
		const uint8_t *side = tloc.side(ebwt.ebwt());
		const int co = (int)tloc._charOff;
		const size_t sideSz = (size_t)ebwt.eh()._sideSz;
		TIndexOffU upto[4], in[4];
		countRange(ebwt, tloc, n, upto, in);
#ifndef NDEBUG
		{
			TIndexOffU upto2[4] = { 0, 0, 0, 0 }, in2[4] = { 0, 0, 0, 0 };
			GroupWalkState gws2(GW_CAT);
			SideLocus tl = tloc, bl = bloc;
			ebwt.mapLFRange(tl, bl, n, upto2, in2, gws2.masks);
			for(int c = 0; c < 4; c++) {
				assert_eq(upto2[c], upto[c]);
				assert_eq(in2[c], in[c]);
			}
		}
#endif
		// cf: the first character present.  The rows of each later one go
		// to a new state, made and initialized in character order.
		int cf = 0;
		while(in[cf] == 0) cf++;
		const TIndexOffU *src = map_.ptr() + mapi_;
		for(int c = cf + 1; c < 4; c++) {
			if(in[c] == 0) continue;
			// NOTE: this can cause the underlying list to be expanded
			// which in turn might leave 'st' pointing to bad memory.
			st.expand();
			GWState& nst = st.back();
			nst.reset();
			nst.map_.resize(in[c]);
			gather(side, co, n, sideSz, c, src, nst.map_.ptr());
			met.branches++;
			pair<TIndexOff, TIndexOff> rret =
			nst.init(
				ebwt,        // forward Bowtie index
				ref,         // bitpair-encodede reference
				sa,          // SA range with offsets
				st,          // EList of all GWStates associated with original range
				hit,         // associated GWHit object
				(TIndexOffU)st.size()-1, // range offset
				reportList,  // if true, report hits to 'res' list
				res,         // report hits here if reportList is true
				upto[c],     // BW top of new range
				upto[c] + in[c], // BW bot of new range
				step+1,      // # steps taken to get to this new range
				met);        // update these metrics
			ret.first += rret.first;
			ret.second += rret.second;
		}
		// This state keeps the rows of character cf, moved to the front of
		// map_ (each moves to the same or an earlier slot)
		if(in[cf] < n) {
			gather(side, co, n, sideSz, cf, src, map_.ptr());
		} else if(mapi_ > 0) {
			memmove(map_.ptr(), src, n * sizeof(TIndexOffU));
		}
		map_.resize(in[cf]);
		mapi_ = 0;
		assert(top != upto[cf] || bot != upto[cf] + in[cf]);
		top = upto[cf];
		bot = top + in[cf];
	}

	/**
	 * Write the entries of 'src' of the rows (of the n starting at character
	 * 'k' of 'side') whose last character is 'c', in order, to 'dst'.  'dst'
	 * may be 'src' or before it in the same list.
	 */
	static void gather(
		const uint8_t* side,
		int k,
		TIndexOffU n,
		size_t sideSz,
		int c,
		const TIndexOffU* src,
		TIndexOffU* dst)
	{
		gwRangeWords(side, k, n, sideSz,
			[&](uint64_t dw, uint64_t rm, TIndexOffU j, int lo) {
				// One row per set bit, in order (each read before any
				// write that could reach it, as dst never passes src + j)
				uint64_t b = gwCharBits(dw, c) & rm;
				while(b != 0) {
					const int p = __builtin_ctzll(b) >> 1;
					*dst++ = src[j + (TIndexOffU)(p - lo)];
					b &= b - 1;
				}
			});
	}
#endif

	/**
	 * Leave this state as init() leaves it at step 0 for a one-row range
	 * (map_ holding element 0) whose element is unresolved and whose row is
	 * neither SA-sampled nor the '$' row: init() resolves nothing, points
	 * the forward map at it and sets tloc at the row.  The caller vouches
	 * that the offset slot was unresolved when the walk was set up; it is
	 * not read here, so a slot resolved since leaves the state as init()
	 * left it then.
	 */
	void initFreshSingle(
		const Ebwt& ebwt,
		GWHit<T>& hit,
		TIndexOffU row)
	{
		assert_eq(1, map_.size());
		assert_eq(0, mapi_);
		assert(!inited_);
		assert_neq(row, ebwt._zOff);
		assert_neq(row, row & ebwt.eh()._offMask);
		ASSERT_ONLY(inited_ = true);
		ASSERT_ONLY(lastStep_ = 0);
		top = row;
		bot = row + 1;
		step = 0;
		hit.fmap[0].first = 0;
		hit.fmap[0].second = 0;
		tloc.initFromRow(top, ebwt.eh(), ebwt.ebwt());
		bloc.invalidate();
	}

	/**
	 * Leave this single-element state as advance() and init() leave it when
	 * its element is resolved after 'k' (>= 1) steps at row 'rowk', the last
	 * step having started from row 'rowkm1': top and bot trimmed past the
	 * resolved row, tloc at the last unresolved row.  For committing a walk
	 * that SpecWalk took instead of advance().
	 */
	void finishSingle(
		const Ebwt& ebwt,
		TIndexOffU k,
		TIndexOffU rowk,
		TIndexOffU rowkm1)
	{
		assert_eq(1, map_.size());
		assert_eq(0, mapi_);
		assert_eq(0, step);
		assert_gt(k, 0);
		top = bot = rowk + 1;
		step = (TIndexOff)k;
		mapi_ = 1;
		tloc.initFromRow(rowkm1, ebwt.eh(), ebwt.ebwt(), false);
		bloc.invalidate();
		ASSERT_ONLY(lastStep_ = step);
	}

	/**
	 * Do what 'k' (>= 1) calls of advance() do to this state, which tracks
	 * one row (bloc invalid) whose element is unresolved, when the walk from
	 * its row takes k LF steps to 'rowk', the first row whose offset
	 * tryOffset() knows, the last step starting from 'rowkm1': the same
	 * counts (bwops, nExFmops, resolves), the offset, the forward-map entry
	 * (each step but the last points it at slot 0), and the state as the
	 * last step's init() leaves it.  The walk is SpecWalk's; nothing else
	 * writes the element's offset slot between those k steps.
	 */
	void commitWalk(
		const Ebwt& ebwt,
		SARangeWithOffs<T>& sa,
		GWHit<T>& hit,
		TIndexOffU range,
		TIndexOffU k,
		TIndexOffU rowk,
		TIndexOffU rowkm1,
		WalkMetrics& met,
		PerReadMetrics& prm)
	{
		assert(!bloc.valid());
		assert_eq(1, size());
		assert_gt(k, 0);
		assert_eq(top, tloc.toBWRow());
		assert_eq(rowk, ebwt.walkLeft(top, k));
		const TIndexOffU e = map_[mapi_];
		assert_eq(OFF_MASK, sa.offs[e]);
		// As advance() moves the element to the front
		map_[0] = e;
		map_.resize(1);
		met.bwops += k;
		prm.nExFmops += k;
		if(k > 1) {
			// init() of each step that leaves the element unresolved
			hit.fmap[e].first = range;
			hit.fmap[e].second = 0;
			tloc.initFromRow(rowkm1, ebwt.eh(), ebwt.ebwt(), false);
		}
		step += (TIndexOff)k;
		const TIndexOffU toff = ebwt.tryOffset(rowk);
		assert_neq(OFF_MASK, toff);
		met.resolves++;
		sa.offs[e] = toff + (TIndexOffU)step;
		assert_eq(sa.offs[e], ebwt.getOffset(sa.topf + e));
		// The last step's init() trims the resolved row
		mapi_ = 1;
		top = bot = rowk + 1;
		ASSERT_ONLY(lastStep_ = step);
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
			map_[i] = (TIndexOffU)i;
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
			if(sa.offs[map(i)] == OFF_MASK) return false;
		}
		return true;
	}

	SideLocus tloc;      // SideLocus for top
	SideLocus bloc;      // SideLocus for bottom
	TIndexOffU  top;       // top elt of range in BWT
	TIndexOffU  bot;       // bot elt of range in BWT
	TIndexOff   step;      // how many steps have we walked to the left so far

protected:

	ASSERT_ONLY(bool inited_);
	ASSERT_ONLY(TIndexOff lastStep_);
	EList<TIndexOffU, 16> map_; // which elts in range 'range' we're tracking
	TIndexOffU mapi_;           // first untrimmed element of map
};

class SpecWalk;

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
		TIndexOffU top = sa.topf;
		TIndexOffU bot = (TIndexOffU)(top + sa.size());
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

	/**
	 * Do what init() does for a one-element range whose row is neither
	 * SA-sampled nor the '$' row and whose offset slot was unresolved at
	 * the time an unbatched walk would call init(): such an init() reads the
	 * slot, resolves nothing and counts nothing.  The slot is not read, so
	 * the walk is left as init() left it then even if another walk sharing
	 * the slot has resolved it since (advanceElement() then only reports
	 * it, as it would have).
	 */
	void initFreshSingle(
		const Ebwt& ebwtFw,         // forward Bowtie index for walking left
		SARangeWithOffs<T>& sa)     // SA range with offsets, one element
	{
		assert_eq(1, sa.size());
		reset();
#ifndef NDEBUG
		inited_ = true;
#endif
		hit_.init(sa, 0, false, 0);
		st_.resize(1);
		st_.back().reset();
		st_.back().initMap(1);
		st_.ensure(4);
		st_.back().initFreshSingle(ebwtFw, hit_, sa.topf);
		elt_ += 1;
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
			advanceElement((TIndexOffU)i, res, met, prm);
		}
	}

	/**
	 * Advance the GroupWalk until the specified element has been
	 * resolved.
	 */
	bool advanceElement(
		TIndexOffU elt,                // element within the range
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
		while(sa.offs[elt] == OFF_MASK) {
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
				(TIndexOffU)range,
				false,
				NULL,
				st_,
				gws,
				met,
				prm);
			assert(sa.offs[elt] != OFF_MASK ||
			       !st_[hit_.fmap[elt].first].doneResolving(sa));
		}
		assert_neq(OFF_MASK, sa.offs[elt]);
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
			(TIndexOffU)sa.len, // length of hit
			sa.offs[elt]);  // resolved text offset
		rep_++;
		return true;
	}

	/**
	 * Do what advanceElement() does, but when the state holding the element
	 * tracks one row, take its remaining LF steps from 'spec' (starting the
	 * walk there if it is not already finished) and commit them at once.
	 * After each step of a state of several rows, the walks of the
	 * one-row states the step leaves are started in 'spec' while it has
	 * free lanes, so that they proceed alongside the walks needed first.
	 * Defined after SpecWalk.
	 */
	bool advanceElementSpec(
		TIndexOffU elt,              // element within the range
		const Ebwt& ebwtFw,          // forward Bowtie index for walking left
		const BitPairReference& ref, // bitpair-encoded reference
		SARangeWithOffs<T>& sa,      // SA range with offsets
		GroupWalkState& gws,         // GroupWalk state; scratch space
		WalkResult& res,             // put the result here
		WalkMetrics& met,            // metrics
		PerReadMetrics& prm,         // per-read metrics
		SpecWalk& spec);             // single-row walks

	/**
	 * Return true iff this walk has a single element and advanceElement()
	 * would step to resolve it: nothing reported yet, no step taken, and
	 * its offset not resolved, neither at init() nor since (by a walk that
	 * shares the offset's slot).  Such an element's walk depends only on
	 * its BW row, sa.topf.
	 */
	bool freshSingle(const SARangeWithOffs<T>& sa) const {
		return elt_ == 1 && rep_ == 0 && st_.size() == 1 &&
		       st_[0].step == 0 && st_[0].size() == 1 &&
		       sa.offs.get(0) == OFF_MASK;
	}

	/**
	 * Do what advanceElement(0, ...) does for a freshSingle() walk, given the
	 * walk's outcome as SpecWalk computed it: 'k' LF steps from sa.topf to
	 * 'rowk', the first row with a known offset, the last step starting from
	 * 'rowkm1'.  Every count advanceElement() keeps (bwops, resolves,
	 * reports, nExFmops), the offset, the reported flag and the state are
	 * left as advanceElement() leaves them.
	 */
	void commitSingle(
		const Ebwt& ebwtFw,          // forward Bowtie index
		SARangeWithOffs<T>& sa,      // SA range with offsets
		TIndexOffU k,                // # LF steps
		TIndexOffU rowk,             // row reached, offset known
		TIndexOffU rowkm1,           // row the last step left
		WalkResult& res,             // put the result here
		WalkMetrics& met,            // metrics
		PerReadMetrics& prm)         // per-read metrics
	{
		assert(inited_);
		assert(freshSingle(sa));
		assert_eq(st_[0].top, sa.topf);
		assert_gt(k, 0);
		assert_eq(rowk, ebwtFw.walkLeft(sa.topf, k));
		TIndexOffU toff = ebwtFw.tryOffset(rowk);
		assert_neq(OFF_MASK, toff);
		met.bwops += k;
		prm.nExFmops += k;
		met.resolves++;
		toff += k;
		assert_eq(toff, ebwtFw.getOffset(sa.topf));
		sa.offs[0] = toff;
		st_[0].finishSingle(ebwtFw, k, rowk, rowkm1);
		// As advanceElement() does once the element is resolved
		if(!hit_.reported(0)) {
			hit_.setReported(0);
		}
		met.reports++;
		res.init(
			0,              // seed offset
			false,          // orientation
			0,              // range
			0,              // element
			sa.topf,        // bw row
			(TIndexOffU)sa.len, // length of hit
			sa.offs[0]);    // resolved text offset
		rep_++;
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
			if(sa.offs[m] != OFF_MASK) {
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
	 * Return the number of BW elements reported so far.
	 */
	size_t numReported() const { return rep_; }

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

/**
 * Walks to resolve the offsets of single BW rows, several stepped in turn so
 * that their memory waits overlap.  A walk starts at a row and takes LF steps
 * until it reaches a row whose offset is known (sampled, or the '$' row),
 * exactly as GroupWalk2S::advanceElement() walks a single-element range.  It
 * reads only the index and writes only this object: nothing is counted or
 * resolved here.  A finished walk is kept, keyed by its starting row, until
 * GroupWalk2S::commitSingle() applies it where an unbatched walk would have
 * walked; walks that are never needed are dropped.  The outcome of a walk
 * depends only on its starting row and the index.
 */
class SpecWalk {

public:

	static const int LANES = 16;   // walks stepped in turn
	static const int TBITS = 10;   // finished walks kept: 2^TBITS, direct-mapped

	SpecWalk() : ebwt_(NULL), nlive_(0) {
		for(int i = 0; i < LANES; i++) live_[i] = false;
		for(int i = 0; i < (1 << TBITS); i++) key_[i] = OFF_MASK;
	}

	/**
	 * Drop the walks in progress; finished ones are kept.
	 */
	void dropLive() {
		for(int i = 0; i < LANES; i++) live_[i] = false;
		nlive_ = 0;
	}

	/**
	 * Find the finished walk from row 'row0'; return false if there is none.
	 */
	bool lookup(
		const Ebwt& ebwt,
		TIndexOffU row0,
		TIndexOffU& k,
		TIndexOffU& rowk,
		TIndexOffU& rowkm1) const
	{
		if(&ebwt != ebwt_) return false;
		const size_t h = slot(row0);
		if(key_[h] != row0) return false;
		k = k_[h]; rowk = rowk_[h]; rowkm1 = rowkm1_[h];
		return true;
	}

	/**
	 * Return the number of walks in progress.
	 */
	int live() const { return nlive_; }

	/**
	 * Return true iff a walk from 'row0' is finished or in progress.
	 */
	bool have(const Ebwt& ebwt, TIndexOffU row0) const {
		if(&ebwt != ebwt_) return false;
		if(key_[slot(row0)] == row0) return true;
		return lane(row0) >= 0;
	}

	/**
	 * Start a walk from 'row0', a row whose offset is not known.  If every
	 * lane is busy, step the walks in progress until one finishes.
	 */
	void start(const Ebwt& ebwt, TIndexOffU row0) {
		if(&ebwt != ebwt_) {
			// Results belong to one index
			for(int i = 0; i < (1 << TBITS); i++) key_[i] = OFF_MASK;
			dropLive();
			ebwt_ = &ebwt;
		}
		assert(row0 != ebwt._zOff);
		assert_neq((row0 & ebwt.eh()._offMask), row0);
		while(nlive_ == LANES) {
			round(ebwt);
		}
		for(int i = 0; i < LANES; i++) {
			if(!live_[i]) {
				Lane& x = lanes_[i];
				x.row0 = x.row = x.prev = row0;
				x.k = 0;
				x.loc.initFromRow(row0, ebwt.eh(), ebwt.ebwt());
				live_[i] = true;
				nlive_++;
				return;
			}
		}
		assert(false);
	}

	/**
	 * Step the walks in progress until the walk from 'row0' finishes.
	 */
	void finish(const Ebwt& ebwt, TIndexOffU row0) {
		const int l = lane(row0);
		if(l < 0) return;
		while(live_[l]) {
			round(ebwt);
		}
	}

	/**
	 * Step the walks in progress until all have finished.
	 */
	void finishAll(const Ebwt& ebwt) {
		while(nlive_ > 0) {
			round(ebwt);
		}
	}

protected:

	struct Lane {
		TIndexOffU row0;   // starting row
		TIndexOffU row;    // current row
		TIndexOffU prev;   // row the last step started from
		TIndexOffU k;      // # steps taken
		SideLocus loc;     // locus of 'row'; set (and prefetched) as the row is reached
	};

	size_t slot(TIndexOffU row0) const {
		return (size_t)(((uint64_t)row0 * 0x9E3779B97F4A7C15ull) >> (64 - TBITS));
	}

	/**
	 * LF(row) for the row at locus 'l' (not the '$' row): what
	 * Ebwt::mapLF1(row, l) sets the row to.  The occurrences of the row's
	 * character before it in its side are counted over every 64-bit word of
	 * the side with masks, rather than by Ebwt::countUpTo()'s loops, whose
	 * trip counts depend on the row and are mispredicted.
	 */
	static TIndexOffU lf(const Ebwt& ebwt, const SideLocus& l) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
		const uint8_t *side = l.side(ebwt.ebwt());
		const int c = unpack_2b_from_8b(side[l._by], l._bp);
		const uint64_t pat = c_table[c];
		const int n = (int)l._charOff;    // characters before the row in its side
		const int NW = (48 * OFF_SIZE) / 32; // 32 characters per word
		TIndexOffU cnt = 0;
		for(int w = 0; w < NW; w++) {
			uint64_t dw;
			memcpy(&dw, side + 8 * w, 8);
			const uint64_t x0 = dw ^ pat;
			const uint64_t x3 = x0 & (x0 >> 1) & 0x5555555555555555ull;
			int r = n - 32 * w;
			r = r < 0 ? 0 : r;
			r = r > 32 ? 32 : r;
			cnt += (TIndexOffU)__builtin_popcountll(x3 & kSpecMask[r]);
		}
		assert_eq(cnt, ebwt.countUpTo(l, c));
		if(c == 0 && l._sideByteOff <= ebwt._zEbwtByteOff &&
		   l._sideByteOff + l._by >= ebwt._zEbwtByteOff)
		{
			// '$' is represented as an 'A' (as in Ebwt::countBt2Side())
			if((l._sideByteOff + l._by > ebwt._zEbwtByteOff) ||
			   (l._sideByteOff + l._by == ebwt._zEbwtByteOff && l._bp > ebwt._zEbwtBpOff))
			{
				cnt--;
			}
		}
		const TIndexOffU *acgt = reinterpret_cast<const TIndexOffU*>(side + ebwt.eh()._sideBwtSz);
		const TIndexOffU ret = endianizeU<TIndexOffU>(acgt[c], ebwt._switchEndian) + cnt + ebwt.fchr()[c];
		assert_eq(ret, ebwt.countBt2Side(l, c));
		return ret;
#else
		TIndexOffU row = l.toBWRow();
		ebwt.mapLF1(row, l);
		return row;
#endif
	}

	int lane(TIndexOffU row0) const {
		for(int i = 0; i < LANES; i++) {
			if(live_[i] && lanes_[i].row0 == row0) return i;
		}
		return -1;
	}

	/**
	 * One LF step of every walk in progress.  A walk ends at the first row
	 * whose offset Ebwt::tryOffset() knows, as GWState::init() checks.
	 */
	void round(const Ebwt& ebwt) {
		const EbwtParams& eh = ebwt.eh();
		for(int i = 0; i < LANES; i++) {
			if(!live_[i]) continue;
			Lane& x = lanes_[i];
			x.prev = x.row;
			x.row = lf(ebwt, x.loc);
			x.k++;
			if(x.row == ebwt._zOff || (x.row & eh._offMask) == x.row) {
				if(x.row != ebwt._zOff) {
					__builtin_prefetch(ebwt.offs() + (x.row >> eh._offRate));
				}
				const size_t h = slot(x.row0);
				key_[h] = x.row0;
				k_[h] = x.k;
				rowk_[h] = x.row;
				rowkm1_[h] = x.prev;
				live_[i] = false;
				nlive_--;
			} else {
				x.loc.initFromRow(x.row, eh, ebwt.ebwt());
			}
		}
	}

	const Ebwt* ebwt_;          // index the walks belong to
	Lane lanes_[LANES];
	bool live_[LANES];
	int nlive_;
	TIndexOffU key_[1 << TBITS];    // starting row of a finished walk, OFF_MASK if none
	TIndexOffU k_[1 << TBITS];      // # steps
	TIndexOffU rowk_[1 << TBITS];   // row reached
	TIndexOffU rowkm1_[1 << TBITS]; // row the last step started from
};

template<typename T, int S>
bool GroupWalk2S<T, S>::advanceElementSpec(
	TIndexOffU elt,              // element within the range
	const Ebwt& ebwtFw,          // forward Bowtie index for walking left
	const BitPairReference& ref, // bitpair-encoded reference
	SARangeWithOffs<T>& sa,      // SA range with offsets
	GroupWalkState& gws,         // GroupWalk state; scratch space
	WalkResult& res,             // put the result here
	WalkMetrics& met,            // metrics
	PerReadMetrics& prm,         // per-read metrics
	SpecWalk& spec)              // single-row walks
{
	assert(inited_);
	assert(!done());
	assert(hit_.repOk(sa));
	assert_lt(elt, sa.size()); // elt must fall within range
	// Until we've resolved our element of interest...
	while(sa.offs[elt] == OFF_MASK) {
		// Get the GWState that contains our element of interest
		size_t range = hit_.fmap[elt].first;
		st_.ensure(4);
		if(!st_[range].bloc.valid()) {
			// One row left, the element's: its remaining steps are the walk
			// from the state's row, or the rest of the walk from the
			// element's own row (resolveElement() walks the rows of small
			// ranges)
			GWState<T>& st = st_[range];
			assert_eq(elt, st.map(st.mapi()));
			TIndexOffU k = 0, rowk = 0, rowkm1 = 0;
			bool have = spec.lookup(ebwtFw, st.top, k, rowk, rowkm1);
			if(!have &&
			   spec.lookup(ebwtFw, sa.topf + elt, k, rowk, rowkm1) &&
			   k > (TIndexOffU)st.step)
			{
				k -= (TIndexOffU)st.step;
				have = true;
			}
			if(!have) {
				if(!spec.have(ebwtFw, st.top)) {
					spec.start(ebwtFw, st.top);
				}
				spec.finish(ebwtFw, st.top);
				have = spec.lookup(ebwtFw, st.top, k, rowk, rowkm1);
			}
			if(have) {
				st.commitWalk(ebwtFw, sa, hit_, (TIndexOffU)range, k, rowk, rowkm1, met, prm);
				continue;
			}
			// The finished walk was displaced by another's: step as before
		}
		const size_t n0 = st_.size();
		st_[range].advance(
			ebwtFw,
			ref,
			sa,
			hit_,
			(TIndexOffU)range,
			false,
			NULL,
			st_,
			gws,
			met,
			prm);
		// Start the walks of the one-row states this step left (this one
		// and those it made) while lanes are free.  Most are never needed
		// (the extension usually stops after a few elements of a large
		// range), so none is started by stepping the others to free a lane:
		// they only progress alongside walks that are needed.
		for(size_t j = n0; j <= st_.size() && spec.live() < SpecWalk::LANES; j++) {
			const GWState<T>& s = st_[j == st_.size() ? range : j];
			if(s.done() || s.bloc.valid()) {
				continue;
			}
			const TIndexOffU e = s.map(s.mapi());
			TIndexOffU k = 0, rowk = 0, rowkm1 = 0;
			if(sa.offs[e] != OFF_MASK ||
			   spec.have(ebwtFw, s.top) ||
			   spec.lookup(ebwtFw, sa.topf + e, k, rowk, rowkm1))
			{
				continue;
			}
			spec.start(ebwtFw, s.top);
		}
		assert(sa.offs[elt] != OFF_MASK ||
		       !st_[hit_.fmap[elt].first].doneResolving(sa));
	}
	assert_neq(OFF_MASK, sa.offs[elt]);
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
		(TIndexOffU)sa.len, // length of hit
		sa.offs[elt]);  // resolved text offset
	rep_++;
	return true;
}

#endif /*GROUP_WALK_H_*/
