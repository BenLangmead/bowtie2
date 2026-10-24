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

#include <array>
#include <atomic>
#include <cstring>
#include "aligner_cache.h"
#include "aligner_seed.h"
#include "search_globals.h"
#include "bt2_idx.h"

using namespace std;

/**
 * Construct a constraint with no edits of any kind allowed.
 */
Constraint Constraint::exact() {
	Constraint c;
	c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::penaltyBased(int pen) {
	Constraint c;
	c.penalty = pen;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint related to the length of the read.
 */
Constraint Constraint::penaltyFuncBased(const SimpleFunc& f) {
	Constraint c;
	c.penFunc = f;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::mmBased(int mms) {
	Constraint c;
	c.mms = mms;
	c.edits = c.dels = c.ins = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::editBased(int edits) {
	Constraint c;
	c.edits = edits;
	c.dels = c.ins = c.mms = 0;
	return c;
}

// Input to seachSeedBi
class SeedAligner::SeedAlignerSearchParams {
public:
	class CacheAndSeed {
	public:
		CacheAndSeed(
			SeedSearchCache &_cache,         // local seed alignment cache
			const InstantiatedSeed& _seed,   // current instantiated seed
		        const Ebwt* ebwtFw, 	         // forward index (BWT)
			const Ebwt* ebwtBw               // backward/mirror index (BWT')

		) : cache(_cache)
		  , seed(_seed)
		  , hasi0(false), fwi0(0), bwi0(0) // just set a default
		{
	                int off = seed.steps[0];
	                bool ltr = off > 0;
        	        off = abs(off)-1;
			int ftabLen = ebwtFw->eh().ftabChars();
			hasi0 = (ftabLen > 1 && ftabLen <= seed.maxjump);
			if(hasi0) {
				if(!ltr) {
					assert_geq(off+1, ftabLen-1);
					off = off - ftabLen + 1;
				}
				// startSearchSeedBi will need them, start prefetching now
				fwi0 = ebwtFw->ftabSeqToInt( cache.getSeq(), off, false);
				ebwtFw->ftabLoHiPrefetch(fwi0);
				if(ebwtBw!=NULL) {
					bwi0 = ebwtBw->ftabSeqToInt( cache.getSeq(), off, false);
					ebwtBw->ftabLoHiPrefetch(bwi0);
				}
			}
		}

		CacheAndSeed(CacheAndSeed &other) = default;
		CacheAndSeed(CacheAndSeed &&other) = default;

		SeedSearchCache &cache;        // local seed alignment cache
		const InstantiatedSeed& seed;  // current instantiated seed
		bool hasi0;
		TIndexOffU fwi0;                // Idx of fw ftab
		TIndexOffU bwi0;                // Idx of bw ftab
	};

	SeedAlignerSearchParams(
		CacheAndSeed &_cs,
		const int _step,              // depth into steps[] array
		const int _depth,             // recursion depth
		const BwtTopBot &_bwt,        // The 4 BWT idxs
		const SideLocus &_tloc,       // locus for top (perhaps unititialized)
		const SideLocus &_bloc,       // locus for bot (perhaps unititialized)
		const std::array<Constraint,3> _cv,        // constraints to enforce in seed zones
		const Constraint &_overall,   // overall constraints to enforce
		DoublyLinkedList<Edit> *_prevEdit)  // previous edit
	: cs(_cs)
	, step(_step)
	, depth(_depth)
	, bwt(_bwt)
	, tloc(_tloc)
	, bloc(_bloc)
	, cv(_cv)
	, overall(_overall)
	, prevEdit(_prevEdit)
	{}

	SeedAlignerSearchParams(
		CacheAndSeed &_cs,
		const int _step,              // depth into steps[] array
		const int _depth,             // recursion depth
		const BwtTopBot &_bwt,        // The 4 BWT idxs
		const SideLocus &_tloc,       // locus for top (perhaps unititialized)
		const SideLocus &_bloc,       // locus for bot (perhaps unititialized)
		const Constraint &_c0,        // constraints to enforce in seed zone 0
		const Constraint &_c1,        // constraints to enforce in seed zone 1
		const Constraint &_c2,        // constraints to enforce in seed zone 2
		const Constraint &_overall,   // overall constraints to enforce
		DoublyLinkedList<Edit> *_prevEdit)  // previous edit
	: cs(_cs)
	, step(_step)
	, depth(_depth)
	, bwt(_bwt)
	, tloc(_tloc)
	, bloc(_bloc)
	, cv{ _c0, _c1, _c2 }
	, overall(_overall)
	, prevEdit(_prevEdit)
	{}

	// create an empty bwt, tloc and bloc, with step=0
	// and constratins from seed, for initial searchSeedBi invocation
	SeedAlignerSearchParams(
		SeedSearchCache &cache,         // local seed alignment cache
		const InstantiatedSeed& seed,   // current instantiated seed
	        const Ebwt* ebwtFw, 	        // forward index (BWT)
		const Ebwt* ebwtBw)             // backward/mirror index (BWT')
	: cs(cache, seed, ebwtFw, ebwtBw)
	, step(0)
	, depth(0)
	, bwt()
	, tloc()
	, bloc()
	, cv{ seed.cons[0], seed.cons[1], seed.cons[2]  }
	, overall(seed.overall)
	, prevEdit(NULL)
	{}

	void checkCV() const {
			assert(cv[0].acceptable());
			assert(cv[1].acceptable());
			assert(cv[2].acceptable());
	}

	CacheAndSeed cs;      // local seed alignment cache and associated instatiated seed
	int step;             // depth into steps[] array
	int depth;            // recursion depth
	BwtTopBot bwt;        // The 4 BWT idxs
	SideLocus tloc;       // locus for top (perhaps unititialized)
	SideLocus bloc;       // locus for bot (perhaps unititialized)
	std::array<Constraint,3> cv;        // constraints to enforce in seed zones
	Constraint overall;   // overall constraints to enforce
	DoublyLinkedList<Edit> *prevEdit;  // previous edit
};


//
// Some static methods for constructing some standard SeedPolicies
//

/**
 * Given a read, depth and orientation, extract a seed data structure
 * from the read and fill in the steps & zones arrays.  The Seed
 * contains the sequence and quality values.
 */
bool
Seed::instantiate(
	const Read& read,
	const BTDnaString& seq, // seed read sequence
	const BTString& qual,   // seed quality sequence
	const Scoring& pens,
	int depth,
	int seedoffidx,
	int seedtypeidx,
	bool fw,
	InstantiatedSeed& is) const
{
	assert(overall != NULL);
	int seedlen = len;
	if((int)read.length() < seedlen) {
		// Shrink seed length to fit read if necessary
		seedlen = (int)read.length();
	}
	assert_gt(seedlen, 0);
	is.steps.resize(seedlen);
	is.zones.resize(seedlen);
	// Fill in 'steps' and 'zones'
	//
	// The 'steps' list indicates which read character should be
	// incorporated at each step of the search process.  Often we will
	// simply proceed from one end to the other, in which case the
	// 'steps' list is ascending or descending.  In some cases (e.g.
	// the 2mm case), we might want to switch directions at least once
	// during the search, in which case 'steps' will jump in the
	// middle.  When an element of the 'steps' list is negative, this
	// indicates that the next
	//
	// The 'zones' list indicates which zone constraint is active at
	// each step.  Each element of the 'zones' list is a pair; the
	// first pair element indicates the applicable zone when
	// considering either mismatch or delete (ref gap) events, while
	// the second pair element indicates the applicable zone when
	// considering insertion (read gap) events.  When either pair
	// element is a negative number, that indicates that we are about
	// to leave the zone for good, at which point we may need to
	// evaluate whether we have reached the zone's budget.
	//
	switch(type) {
		case SEED_TYPE_EXACT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = -(seedlen - k);
				// Zone 0 all the way
				is.zones[k].first = is.zones[k].second = 0;
			}
			break;
		}
		case SEED_TYPE_LEFT_TO_RIGHT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = k+1;
				// Zone 0 from 0 up to ceil(len/2), then 1
				is.zones[k].first = is.zones[k].second = ((k < (seedlen+1)/2) ? 0 : 1);
			}
			// Zone 1 ends at the RHS
			is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
			break;
		}
		case SEED_TYPE_RIGHT_TO_LEFT: {
			for(int k = 0; k < seedlen; k++) {
				is.steps[k] = -(seedlen - k);
				// Zone 0 from 0 up to floor(len/2), then 1
				is.zones[k].first  = ((k < seedlen/2) ? 0 : 1);
				// Inserts: Zone 0 from 0 up to ceil(len/2)-1, then 1
				is.zones[k].second = ((k < (seedlen+1)/2+1) ? 0 : 1);
			}
			is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
			break;
		}
		case SEED_TYPE_INSIDE_OUT: {
			// Zone 0 from ceil(N/4) up to N-floor(N/4)
			int step = 0;
			for(int k = (seedlen+3)/4; k < seedlen - (seedlen/4); k++) {
				is.zones[step].first = is.zones[step].second = 0;
				is.steps[step++] = k+1;
			}
			// Zone 1 from N-floor(N/4) up
			for(int k = seedlen - (seedlen/4); k < seedlen; k++) {
				is.zones[step].first = is.zones[step].second = 1;
				is.steps[step++] = k+1;
			}
			// No Zone 1 if seedlen is short (like 2)
			//assert_eq(1, is.zones[step-1].first);
			is.zones[step-1].first = is.zones[step-1].second = -1;
			// Zone 2 from ((seedlen+3)/4)-1 down to 0
			for(int k = ((seedlen+3)/4)-1; k >= 0; k--) {
				is.zones[step].first = is.zones[step].second = 2;
				is.steps[step++] = -(k+1);
			}
			assert_eq(2, is.zones[step-1].first);
			is.zones[step-1].first = is.zones[step-1].second = -2;
			assert_eq(seedlen, step);
			break;
		}
		default:
			throw 1;
	}
	// Instantiate constraints
	for(int i = 0; i < 3; i++) {
		is.cons[i] = zones[i];
		is.cons[i].instantiate(read.length());
	}
	is.overall = *overall;
	is.overall.instantiate(read.length());
	// Take a sweep through the seed sequence.  Consider where the Ns
	// occur and how zones are laid out.  Calculate the maximum number
	// of positions we can jump over initially (e.g. with the ftab) and
	// perhaps set this function's return value to false, indicating
	// that the arrangements of Ns prevents the seed from aligning.
	bool streak = true;
	is.maxjump = 0;
	bool ret = true;
	bool ltr = (is.steps[0] > 0); // true -> left-to-right
	for(size_t i = 0; i < is.steps.size(); i++) {
		assert_neq(0, is.steps[i]);
		int off = is.steps[i];
		off = abs(off)-1;
		Constraint& cons = is.cons[abs(is.zones[i].first)];
		int c = seq[off];  assert_range(0, 4, c);
		int q = qual[off];
		if(ltr != (is.steps[i] > 0) || // changed direction
		   is.zones[i].first != 0 ||   // changed zone
		   is.zones[i].second != 0)    // changed zone
		{
			streak = false;
		}
		if(c == 4) {
			// Induced mismatch
			if(cons.canN(q, pens)) {
				cons.chargeN(q, pens);
			} else {
				// Seed disqualified due to arrangement of Ns
				return false;
			}
		}
		if(streak) is.maxjump++;
	}
	is.seedoff = depth;
	is.seedoffidx = seedoffidx;
	is.fw = fw;
	is.s = *this;
	return ret;
}

/**
 * Bitwise equality of two doubles.
 */
static inline bool sameDouble(double a, double b) {
	uint64_t x, y;
	memcpy(&x, &a, sizeof(x));
	memcpy(&y, &b, sizeof(y));
	return x == y;
}

/**
 * Return true iff two penalty functions have the same type and parameters.
 */
static bool sameFunc(const SimpleFunc& a, const SimpleFunc& b) {
	return a.getType() == b.getType() &&
	       sameDouble(a.getMin(),   b.getMin()) &&
	       sameDouble(a.getMax(),   b.getMax()) &&
	       sameDouble(a.getConst(), b.getConst()) &&
	       sameDouble(a.getCoeff(), b.getCoeff());
}

/**
 * Return true iff every field of two constraints is the same.
 */
static bool sameConstraint(const Constraint& a, const Constraint& b) {
	return a.edits       == b.edits &&
	       a.mms         == b.mms &&
	       a.ins         == b.ins &&
	       a.dels        == b.dels &&
	       a.penalty     == b.penalty &&
	       a.editsCeil   == b.editsCeil &&
	       a.mmsCeil     == b.mmsCeil &&
	       a.insCeil     == b.insCeil &&
	       a.delsCeil    == b.delsCeil &&
	       a.penaltyCeil == b.penaltyCeil &&
	       a.instantiated == b.instantiated &&
	       sameFunc(a.penFunc, b.penFunc);
}

/**
 * Return true iff every field of two seed policies is the same (the
 * overall constraint is compared by pointer).
 */
static bool sameSeed(const Seed& a, const Seed& b) {
	return a.len == b.len &&
	       a.type == b.type &&
	       a.overall == b.overall &&
	       sameConstraint(a.zones[0], b.zones[0]) &&
	       sameConstraint(a.zones[1], b.zones[1]) &&
	       sameConstraint(a.zones[2], b.zones[2]);
}

#ifndef NDEBUG
/**
 * Return true iff two InstantiatedSeeds hold the same fields as
 * Seed::instantiate writes them.
 */
static bool sameInstantiation(const InstantiatedSeed& a, const InstantiatedSeed& b) {
	if(a.steps.size() != b.steps.size() || a.zones.size() != b.zones.size()) {
		return false;
	}
	for(size_t k = 0; k < a.steps.size(); k++) {
		if(a.steps[k] != b.steps[k]) return false;
	}
	for(size_t k = 0; k < a.zones.size(); k++) {
		if(a.zones[k] != b.zones[k]) return false;
	}
	return sameConstraint(a.cons[0], b.cons[0]) &&
	       sameConstraint(a.cons[1], b.cons[1]) &&
	       sameConstraint(a.cons[2], b.cons[2]) &&
	       sameConstraint(a.overall, b.overall) &&
	       a.maxjump == b.maxjump &&
	       a.seedoff == b.seedoff &&
	       a.seedoffidx == b.seedoffidx &&
	       a.fw == b.fw &&
	       sameSeed(a.s, b.s);
}
#endif

/**
 * Source of SeedInstKey generation numbers; 0 is never handed out.
 */
static std::atomic<uint64_t> seedInstGens(0);

bool SeedInstKey::matches(const Seed& sd, size_t rdl) const {
	assert(sd.overall != NULL);
	return gen != 0 &&
	       rdlen == rdl &&
	       sameSeed(seed, sd) &&
	       sameConstraint(overall, *sd.overall);
}

void SeedInstKey::set(const Seed& sd, size_t rdl) {
	assert(sd.overall != NULL);
	seed = sd;
	overall = *sd.overall;
	rdlen = rdl;
	gen = seedInstGens.fetch_add(1, std::memory_order_relaxed) + 1;
}

/**
 * Return a set consisting of 1 seed encapsulating an exact matching
 * strategy.
 */
void
Seed::zeroMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_EXACT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::exact();
	pols.back().zones[2] = Constraint::exact();
	pols.back().overall = &oall;
}

/**
 * Return a set of 2 seeds encapsulating a half-and-half 1mm strategy.
 */
void
Seed::oneMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0;
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
}

/**
 * Return a set of 3 seeds encapsulating search roots for:
 *
 * 1. Starting from the left-hand side and searching toward the
 *    right-hand side allowing 2 mismatches in the right half.
 * 2. Starting from the right-hand side and searching toward the
 *    left-hand side allowing 2 mismatches in the left half.
 * 3. Starting (effectively) from the center and searching out toward
 *    both the left and right-hand sides, allowing one mismatch on
 *    either side.
 *
 * This is not exhaustive.  There are 2 mismatch cases mised; if you
 * imagine the seed as divided into four successive quarters A, B, C
 * and D, the cases we miss are when mismatches occur in A and C or B
 * and D.
 */
void
Seed::twoMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[1].mmsCeil = 1; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 3: inside-out search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_INSIDE_OUT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::mmBased(1);
	pols.back().zones[2].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().overall = &oall;
}

/**
 * Types of actions that can be taken by the SeedAligner.
 */
enum {
	SA_ACTION_TYPE_RESET = 1,
	SA_ACTION_TYPE_SEARCH_SEED, // 2
	SA_ACTION_TYPE_FTAB,        // 3
	SA_ACTION_TYPE_FCHR,        // 4
	SA_ACTION_TYPE_MATCH,       // 5
	SA_ACTION_TYPE_EDIT         // 6
};

/**
 * Given a read and a few coordinates that describe a substring of the read (or
 * its reverse complement), fill in 'seq' and 'qual' objects with the seed
 * sequence and qualities.
 *
 * The seq field is filled with the sequence as it would align to the Watson
 * reference strand.  I.e. if fw is false, then the sequence that appears in
 * 'seq' is the reverse complement of the raw read substring.
 */
bool
SeedAligner::instantiateSeq(
	const Read& read, // input read
	BTDnaString& seq, // output sequence
	BTString& qual,   // output qualities
	int len,          // seed length
	int depth,        // seed's 0-based offset from 5' end
	bool fw) const    // seed's orientation
{
	// Fill in 'seq' and 'qual'
	int seedlen = len;
	if((int)read.length() < seedlen) seedlen = (int)read.length();
	seq.resize(len);
	qual.resize(len);
	// If fw is false, we take characters starting at the 3' end of the
	// reverse complement of the read.
	int cor = 0; // OR of the characters (0-3, N = 4)
	for(int i = 0; i < len; i++) {
		const int c = read.patFw.windowGetDna(i, fw, depth, len);
		cor |= c;
		seq.set(c, i);
		qual.set(read.qual.windowGet(i, fw, depth, len), i);
	}
	return (cor & 4) != 0;
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed, instantiate the seed, retracting if necessary.
 */
pair<int, int> SeedAligner::instantiateSeeds(
	const EList<Seed>& seeds,  // search seeds
	size_t off,                // offset into read to start extracting
	int per,                   // interval between seeds
	const Read& read,          // read to align
	const Scoring& pens,       // scoring scheme
	bool nofw,                 // don't align forward read
	bool norc,                 // don't align revcomp read
	AlignmentCacheIface& cache,// holds some seed hits from previous reads
	SeedResults& sr,           // holds all the seed hits
	SeedSearchMetrics& met,    // metrics
	pair<int, int>& instFw,
	pair<int, int>& instRc)
{
	assert(!seeds.empty());
	assert_gt(read.length(), 0);
	// Check whether read has too many Ns
	offIdx2off_.clear();
	int len = seeds[0].len; // assume they're all the same length
#ifndef NDEBUG
	for(size_t i = 1; i < seeds.size(); i++) {
		assert_eq(len, seeds[i].len);
	}
#endif
	// Calc # seeds within read interval
	int nseeds = 1;
	if((int)read.length() - (int)off > len) {
		nseeds += ((int)read.length() - (int)off - len) / per;
	}
	for(int i = 0; i < nseeds; i++) {
		offIdx2off_.push_back(per * i + (int)off);
	}
	pair<int, int> ret;
	ret.first = 0;  // # seeds that require alignment
	ret.second = 0; // # seeds that hit in cache with non-empty results
	sr.reset(read, offIdx2off_, nseeds);
	assert(sr.repOk(&cache.current(), true)); // require that SeedResult be initialized
	// The instantiation keys of sr's InstantiatedSeeds, one per seed policy
	std::vector<SeedInstKey>& keys = sr.instKeys();
	if(keys.size() < seeds.size()) keys.resize(seeds.size());
	for(size_t j = 0; j < seeds.size(); j++) {
		if(!keys[j].matches(seeds[j], read.length())) {
			keys[j].set(seeds[j], read.length());
		}
	}
	// For each seed position
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			// Skip this orientation b/c user specified --nofw or --norc
			continue;
		}
		// For each seed position
		for(int i = 0; i < nseeds; i++) {
			int depth = i * per + (int)off;
			int seedlen = seeds[0].len;
			// Extract the seed sequence at this offset
			// If fw == true, we extract the characters from i*per to
			// i*(per-1) (exclusive).  If fw == false, 
			const bool hasN = instantiateSeq(
				read,
				sr.seqs(fw)[i],
				sr.quals(fw)[i],
				std::min<int>((int)seedlen, (int)read.length()),
				depth,
				fw);
			QKey qk(sr.seqs(fw)[i] ASSERT_ONLY(, tmpdnastr_));
			// For each search strategy
			EList<InstantiatedSeed>& iss = sr.instantiatedSeeds(fw, i);
			for(int j = 0; j < (int)seeds.size(); j++) {
				iss.expand();
				assert_eq(seedlen, seeds[j].len);
				InstantiatedSeed* is = &iss.back();
				const uint64_t gen = keys[j].gen;
				if(!hasN && is->instGen == gen) {
					// 'is' already holds this policy's instantiation for
					// this read length, and there are no Ns to charge to
					// its constraints: only the offset fields differ
					is->seedoff = depth;
					is->seedoffidx = i;
					is->fw = fw;
#ifndef NDEBUG
					InstantiatedSeed chk;
					const bool chkok = seeds[j].instantiate(
						read, sr.seqs(fw)[i], sr.quals(fw)[i], pens,
						depth, i, j, fw, chk);
					assert(chkok);
					assert(sameInstantiation(chk, *is));
#endif
					ret.first++;
					if(fwi == 0) { instFw.first++; } else { instRc.first++; }
				} else if(seeds[j].instantiate(
					read,
					sr.seqs(fw)[i],
					sr.quals(fw)[i],
					pens,
					depth,
					i,
					j,
					fw,
					*is))
				{
					// Without Ns the constraints are uncharged and 'is'
					// holds the key's instantiation
					is->instGen = hasN ? 0 : gen;
					// Can we fill this seed hit in from the cache?
					ret.first++;
					if(fwi == 0) { instFw.first++; } else { instRc.first++; }
				} else {
					// Seed may fail to instantiate if there are Ns
					// that prevent it from matching
					is->instGen = 0;
					met.filteredseed++;
					iss.pop_back();
				}
			}
		}
	}
	return ret;
}

class SeedAlignerSearchState {
public:
	TIndexOffU tp[4], bp[4]; // dest BW ranges for "prime" index
	TIndexOffU t[4], b[4];   // dest BW ranges
	TIndexOffU *tf, *tb, *bf, *bb; // depend on ltr
	const Ebwt* ebwt;

	TIndexOffU ntop;
	int off;
	bool ltr;
	bool done;
public:
	SeedAlignerSearchState()
	: tp{0,0,0,0}, bp{0,0,0,0}
	, t{0,0,0,0}, b{0,0,0,0}
	, tf(NULL), tb(NULL), bf(NULL), bb(NULL)
	, ebwt(NULL)
	, ntop(0)
	, off(0)
	, ltr(false)
	, done(false)
	{}

	void setOff(
		size_t _off,
		const BwtTopBot &bwt,      // The 4 BWT idxs
		const Ebwt* ebwtFw_,       // forward index (BWT)
		const Ebwt* ebwtBw_)       // backward/mirror index (BWT')
	{
		off = _off;
		ltr = off > 0;
		t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
		off = abs(off)-1;
		if(ltr) {
			ebwt = ebwtBw_;
			tp[0] = tp[1] = tp[2] = tp[3] = bwt.topf;
			bp[0] = bp[1] = bp[2] = bp[3] = bwt.botf;
			tf = tp; tb = t;
			bf = bp; bb = b;
			ntop = bwt.topb;
		} else {
			ebwt = ebwtFw_;
			tp[0] = tp[1] = tp[2] = tp[3] = bwt.topb;
			bp[0] = bp[1] = bp[2] = bp[3] = bwt.botb;
			tf = t; tb = tp;
			bf = b; bb = bp;
			ntop = bwt.topf;
		}
		assert(ebwt != NULL);
	}

public:
#ifndef NDEBUG
	TIndexOffU lasttot;

	void initLastTot(TIndexOffU tot) { lasttot = tot;}
	void assertLeqAndSetLastTot(TIndexOffU tot) {
		assert_leq(tot, lasttot);
		lasttot = tot;
	}
#else
	// noop in production code
	void initLastTot(TIndexOffU tot) {}
	void assertLeqAndSetLastTot(TIndexOffU tot) {};
#endif

};

// searchAllSeeds scratch, kept across calls so that it allocates only
// when a read needs more seed searches than any read before it
struct SeedAligner::SeedSearchScratch {
	std::vector<SeedAlignerSearchParams> params; // one per seed search
	std::vector<SeedAlignerSearchState> states;  // one per seed search
	std::vector<size_t> first;                   // per cache: its first search
};

SeedAligner::~SeedAligner() {
	delete sss_;
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed:
 *
 * 1. Instantiate all seeds, retracting them if necessary.
 * 2. Calculate zone boundaries for each seed
 *
 * The seed offsets of both strands are searched together, up to
 * SEED_BATCH_CACHES offsets per searchSeedBi batch, so that more searches
 * are in flight than with one batch per strand.  The searches write only
 * their SeedSearchCaches and the bwops_/bwedits_ counts, and the searches
 * of one offset stay in one batch, in their order; the results are then
 * handed to the alignment cache and to sr in the order the per-strand
 * batches used: forward seed offsets in order, then reverse-complement ones.
 */
void SeedAligner::searchAllSeeds(
	const EList<Seed>& seeds,    // search seeds
	const Ebwt* ebwtFw,          // BWT index
	const Ebwt* ebwtBw,          // BWT' index
	const Read& read,            // read to align
	const Scoring& pens,         // scoring scheme
	size_t ncut_f,               // min seed quality for BWT 
	AlignmentCacheIface& cache,  // local cache for seed alignments
	SeedResults& sr,             // holds all the seed hits
	SeedSearchMetrics& met,      // metrics
	PerReadMetrics& prm)         // per-read metrics
{
	assert(!seeds.empty());
	assert(ebwtFw != NULL);
	assert(ebwtFw->isInMemory());
	assert(sr.repOk(&cache.current()));
	ebwtFw_ = ebwtFw;
	ebwtBw_ = ebwtBw;
	sc_ = &pens;
	read_ = &read;
	ncut_f_ = ncut_f;
	bwops_ = bwedits_ = 0;
	uint64_t possearches = 0, seedsearches = 0, intrahits = 0, interhits = 0, ooms = 0;

	if(sss_ == NULL) sss_ = new SeedSearchScratch();
	SeedSearchMultiCache& mcache = mcache_;
	std::vector<SeedAlignerSearchParams>& paramVec = sss_->params;
	std::vector<SeedAlignerSearchState>& sstateVec = sss_->states;
	std::vector<size_t>& first = sss_->first;
	mcache.clear();
	paramVec.clear();
	first.clear();

	// One cache per seed offset with something to search, forward offsets
	// first.  All caches are taken before the searches refer to them.
	for(int fwi = 0; fwi < 2; fwi++) {
		const bool fw = (fwi == 0);
		for(int i = 0; i < (int)sr.numOffs(); i++) {
			assert(sr.repOk(&cache.current()));
			if(sr.instantiatedSeeds(fw, i).empty()) {
				// Cache hit in an across-read cache
				continue;
			}
			mcache.add(sr.seqs(fw)[i], sr.quals(fw)[i], i, fw);
			possearches++;
		}
	}
	// The searches, in the same order
	for(size_t mnr = 0; mnr < mcache.size(); mnr++) {
		const bool fw = mcache.getFw(mnr);
		const int i = mcache.getSeedOffIdx(mnr);
		EList<InstantiatedSeed>& iss = sr.instantiatedSeeds(fw, i);
		SeedSearchCache &srcache = mcache[mnr];
		first.push_back(paramVec.size());
		for(size_t j = 0; j < iss.size(); j++) {
			// Set seq and qual appropriately, using the seed sequences
			// and qualities already installed in SeedResults
			assert_eq(fw, iss[j].fw);
			assert_eq(i, (int)iss[j].seedoffidx);
			paramVec.emplace_back(srcache, iss[j], ebwtFw_, ebwtBw_);
			seedsearches++;
		}
	}
	first.push_back(paramVec.size());

	// do the searches
	if(sstateVec.size() < paramVec.size()) sstateVec.resize(paramVec.size());
	for(size_t c0 = 0; c0 < mcache.size(); c0 += SEED_BATCH_CACHES) {
		const size_t c1 = std::min(c0 + SEED_BATCH_CACHES, mcache.size());
		const size_t p0 = first[c0], p1 = first[c1];
		if(p1 > p0) {
			searchSeedBi(p1 - p0, &(paramVec[p0]), &(sstateVec[p0]));
		}
	}

	// finish aligning and add to SeedResult
	for (size_t mnr=0; mnr<mcache.size(); mnr++) {
		SeedSearchCache &srcache = mcache[mnr];
		// Tell the cache that we've started aligning, so the cache can
		// expect a series of on-the-fly updates
		int ret = srcache.beginAlign(cache);
		if(ret == -1) {
			// Out of memory when we tried to add key to map
			ooms++;
			continue;
		}
		assert(srcache.aligning());
		if(!srcache.addAllCached()){
			// Memory exhausted during copy
			ooms++;
			continue;
		}
		srcache.finishAlign();
		assert(!srcache.aligning());
		if(srcache.qvValid()) {
			sr.add(
				srcache.getQv(),   // range of ranges in cache
				cache.current(), // cache
				mcache.getSeedOffIdx(mnr),     // seed index (from 5' end)
				mcache.getFw(mnr));   // whether seed is from forward read
		}
	} // mnr loop
	prm.nSeedRanges = sr.numRanges();
	prm.nSeedElts = sr.numElts();
	prm.nSeedRangesFw = sr.numRangesFw();
	prm.nSeedRangesRc = sr.numRangesRc();
	prm.nSeedEltsFw = sr.numEltsFw();
	prm.nSeedEltsRc = sr.numEltsRc();
	prm.seedMedian = (uint64_t)(sr.medianHitsPerSeed() + 0.5);
	prm.seedMean = (uint64_t)sr.averageHitsPerSeed();

	prm.nSdFmops += bwops_;
	met.seedsearch += seedsearches;
	met.nrange += sr.numRanges();
	met.nelt += sr.numElts();
	met.possearch += possearches;
	met.intrahit += intrahits;
	met.interhit += interhits;
	met.ooms += ooms;
	met.bwops += bwops_;
	met.bweds += bwedits_;
}

bool SeedAligner::sanityPartial(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	const BTDnaString& seq,
	size_t dep,
	size_t len,
	bool do1mm,
	TIndexOffU topfw,
	TIndexOffU botfw,
	TIndexOffU topbw,
	TIndexOffU botbw)
{
	tmpdnastr_.clear();
	for(size_t i = dep; i < len; i++) {
		tmpdnastr_.append(seq[i]);
	}
	TIndexOffU top_fw = 0, bot_fw = 0;
	ebwtFw->contains(tmpdnastr_, &top_fw, &bot_fw);
	assert_eq(top_fw, topfw);
	assert_eq(bot_fw, botfw);
	if(do1mm && ebwtBw != NULL) {
		tmpdnastr_.reverse();
		TIndexOffU top_bw = 0, bot_bw = 0;
		ebwtBw->contains(tmpdnastr_, &top_bw, &bot_bw);
		assert_eq(top_bw, topbw);
		assert_eq(bot_bw, botbw);
	}
	return true;
}

inline void exactSweepInit(
	const Ebwt&        ebwt,
	const BTDnaString& seq,
	const int          ftabLen,
	const size_t       len,
	size_t            &dep,
	TIndexOffU        &top, 
	TIndexOffU        &bot
	)
{
	top = bot = 0;

	const size_t left = len - dep;
	assert_gt(left, 0);
	bool doFtab = ftabLen > 1 && left >= (size_t)ftabLen;
	if(doFtab) {
		const size_t endi = len-dep-1;
		// Does N interfere with use of Ftab?
		for(size_t i = 0; i < (size_t)ftabLen; i++) {
			int c = seq[endi-i];
			if(c > 3) {
				doFtab = false;
				break;
			}
		}
	}
	if(doFtab) {
		// Use ftab
		ebwt.ftabLoHi(seq, left - ftabLen, false, top, bot);
		dep += (size_t)ftabLen;
	} else {
		// Use fchr
		int c = seq[len-dep-1];
		if(c < 4) {
			top = ebwt.fchr()[c];
			bot = ebwt.fchr()[c+1];
		}
		dep++;
	}
}

/**
 * The first part of exactSweepInit for a sweep at depth dep: true if it
 * uses the ftab, and then the ftab index in fi (computed as ftabLoHi does).
 */
inline bool exactSweepFtabIdx(
	const Ebwt&        ebwt,
	const BTDnaString& seq,
	const int          ftabLen,
	const size_t       len,
	const size_t       dep,
	TIndexOffU        &fi)
{
	const size_t left = len - dep;
	assert_gt(left, 0);
	if(ftabLen <= 1 || left < (size_t)ftabLen) {
		return false;
	}
	const size_t endi = len-dep-1;
	// Does N interfere with use of Ftab?
	for(size_t i = 0; i < (size_t)ftabLen; i++) {
		int c = seq[endi-i];
		if(c > 3) {
			return false;
		}
	}
	fi = ebwt.ftabSeqToInt(seq, left - ftabLen, false);
	return true;
}

/**
 * The rest of exactSweepInit, given exactSweepFtabIdx's result for the
 * same depth: the same top, bot and dep.
 */
inline void exactSweepInitIdx(
	const Ebwt&        ebwt,
	const BTDnaString& seq,
	const int          ftabLen,
	const size_t       len,
	const bool         doFtab,
	const TIndexOffU   fi,
	size_t            &dep,
	TIndexOffU        &top,
	TIndexOffU        &bot
	)
{
	top = bot = 0;
	if(doFtab) {
		// Use ftab
		if(fi != std::numeric_limits<TIndexOffU>::max()) {
			ebwt.ftabLoHi(fi, top, bot);
		}
		dep += (size_t)ftabLen;
	} else {
		// Use fchr
		int c = seq[len-dep-1];
		if(c < 4) {
			top = ebwt.fchr()[c];
			bot = ebwt.fchr()[c+1];
		}
		dep++;
	}
}

inline void exactSweepMapLF(
	const Ebwt&        ebwt,
	const BTDnaString& seq,
	const size_t       len,
	const size_t       dep,
	const SideLocus   &tloc, 
	const SideLocus  &bloc,
	TIndexOffU        &top, 
	TIndexOffU        &bot,
	uint64_t          &bwops           // Burrows-Wheeler operations
)
{
	int c = seq[len-dep-1];
	if(c > 3) {
		top = bot = 0;
	} else {
		if(bloc.valid()) {
			bwops += 2;
			top = ebwt.mapLF(tloc, c);
			bot = ebwt.mapLF(bloc, c);
		} else {
			bwops++;
			top = ebwt.mapLF1(top, tloc, c);
			if(top == OFF_MASK) {
				top = bot = 0;
			} else {
				bot = top+1;
			}
		}
	}
}


inline bool exactSweepStep(
	const Ebwt&        ebwt,    // BWT index
	const TIndexOffU   top, 
	const TIndexOffU   bot,
	const size_t       mineMax, // don't care about edit bounds > this
	SideLocus         &tloc, 
	SideLocus         &bloc,
	size_t            &mineCnt, // minimum # edits
	size_t            &nedit,
	bool              &done
	)
{
	if(bot <= top) {
		nedit++;
		if(nedit >= mineMax) {
			mineCnt = nedit;
			done = true;
		}
		return true;
	}
	INIT_LOCS(top, bot, tloc, bloc, ebwt);
	return false;
}

/**
 * Sweep right-to-left and left-to-right using exact matching.  Remember all
 * the SA ranges encountered along the way.  Report exact matches if there are
 * any.  Calculate a lower bound on the number of edits in an end-to-end
 * alignment.
 */
size_t SeedAligner::exactSweep(
	const Ebwt&        ebwt,    // BWT index
	const Read&        read,    // read to align
	const Scoring&     sc,      // scoring scheme
	bool               nofw,    // don't align forward read
	bool               norc,    // don't align revcomp read
	size_t             mineMax, // don't care about edit bounds > this
	size_t&            mineFw,  // minimum # edits for forward read
	size_t&            mineRc,  // minimum # edits for revcomp read
	bool               repex,   // report 0mm hits?
	SeedResults&       hits,    // holds all the seed hits (and exact hit)
	SeedSearchMetrics& met)     // metrics
{
	assert_gt(mineMax, 0);
	const size_t len = read.length();
	const int ftabLen = ebwt.eh().ftabChars();

	size_t nelt = 0;

	std::array<SideLocus,2> tloc;
	std::array<SideLocus,2> bloc;
	TIndexOffU top[2] = {0, 0};
	TIndexOffU bot[2] = {0, 0};

	size_t dep[2] = {0, 0};
	size_t nedit[2] = {0, 0};
	bool doInit[2] = {true, true};

	size_t prefetch_count = 0;
	bool done[2] = {nofw, norc};

	for(int fwi = 0; fwi < 2; fwi++) {
		if (!done[fwi]) {
			bool fw = (fwi == 0);
			const BTDnaString& seq = fw ? read.patFw : read.patRc;
			assert(!seq.empty());
			__builtin_prefetch(&(seq[len-1]));
			if (len>48) __builtin_prefetch(&(seq[len-49])); // HW prefetch prediction assumes forward, help it
		}
	}

	while( (dep[0] < len && !done[0]) || (dep[1] < len && !done[1]) ) {
		prefetch_count++;
		if (prefetch_count>=48) { // cache line is 64 bytes, but we may skip some deps
			for(int fwi = 0; fwi < 2; fwi++) {
				if (dep[fwi] < len && !done[fwi]) {
					bool fw = (fwi == 0);
					const BTDnaString& seq = fw ? read.patFw : read.patRc;
					const size_t left = len-dep[fwi];
					if (left>48) {
						__builtin_prefetch(&(seq[left-49])); // HW prefetch prediction assumes forward, help it
					}
				}
			}
			prefetch_count=0;
		}
		// by doing both fw in the internal loop, I give the prefetch in exactSweepStep to be effective
		for(int fwi = 0; fwi < 2; fwi++) {
			if (dep[fwi] < len && !done[fwi]) {
				bool fw = (fwi == 0);
				const BTDnaString& seq = fw ? read.patFw : read.patRc;

				if (doInit[fwi]) {
					exactSweepInit(ebwt, seq, ftabLen, len,            // in
							dep[fwi], top[fwi], bot[fwi]);          // out
					if ( exactSweepStep(ebwt, top[fwi], bot[fwi], mineMax,
							tloc[fwi], bloc[fwi],
							fw ? mineFw : mineRc,
							nedit[fwi], done[fwi]) ) {
						continue;
					}
					doInit[fwi]=false;
				}

				if (dep[fwi]< len) {
					exactSweepMapLF(ebwt, seq, len, dep[fwi], tloc[fwi], bloc[fwi],
							top[fwi], bot[fwi], bwops_);

					if ( exactSweepStep(ebwt, top[fwi], bot[fwi], mineMax,
								tloc[fwi], bloc[fwi],
								fw ? mineFw : mineRc,
								nedit[fwi], done[fwi]) ) {
						doInit[fwi]=true;
					}
					dep[fwi]++;
				}
			}
		}
	}

	for(int fwi = 0; fwi < 2; fwi++) {
		if( (!done[fwi]) && (dep[fwi] >= len) ) {
			const bool fw = (fwi == 0);

			// Set the minimum # edits
			if(fw) { mineFw = nedit[fwi]; } else { mineRc = nedit[fwi]; }
			// Done
			if(nedit[fwi] == 0 && bot[fwi] > top[fwi]) {
				if(repex) {
					// This is an exact hit
					int64_t score = len * sc.match();
					if(fw) {
						hits.addExactEeFw(top[fwi], bot[fwi], NULL, NULL, fw, score);
						assert(ebwt.contains(fw ? read.patFw : read.patRc, NULL, NULL));
					} else {
						hits.addExactEeRc(top[fwi], bot[fwi], NULL, NULL, fw, score);
						assert(ebwt.contains(fw ? read.patFw : read.patRc, NULL, NULL));
					}
				}
				nelt += (bot[fwi] - top[fwi]);
			}
		}
	}
	return nelt;
}

/**
 * exactSweep for many reads at once.  exactSweep steps the sweeps of the read
 * and of its reverse complement in turn, so that each step's prefetch of the
 * sides its sweep needs next is hidden by the other sweep's step.  Usually
 * only one of the two survives past the first few steps, so for most of a
 * sweep one load is in flight.  Here the sweeps of all the jobs (each read,
 * then its reverse complement, job after job) are stepped in turn in
 * EXACT_SWEEP_LANES lanes; when a sweep ends, the next one takes its lane.
 * The sweeps share no state, and each takes exactly the steps it takes in
 * exactSweep, so each job's edit bounds, exact hits, element count (nelt)
 * and LF operations are the ones exactSweep gives for its read.
 */
void SeedAligner::exactSweepBatch(
	const Ebwt&        ebwt,    // BWT index
	const Scoring&     sc,      // scoring scheme
	size_t             mineMax, // don't care about edit bounds > this
	bool               repex,   // report 0mm hits?
	ExactSweepJob* const* jobs, // the reads to sweep
	size_t             njobs)   // # jobs
{
	assert_gt(mineMax, 0);
	const int ftabLen = ebwt.eh().ftabChars();

	// One sweep in progress
	struct Lane {
		ExactSweepJob*     job;    // NULL: lane is idle
		int                fwi;    // 0: the read, 1: its reverse complement
		const BTDnaString* seq;
		size_t             len;
		size_t             dep;
		size_t             nedit;
		TIndexOffU         top;
		TIndexOffU         bot;
		SideLocus          tloc;
		SideLocus          bloc;
		bool               doInit;
		bool               done;
		size_t             prefetch_count;
		bool               doFtab; // the next (re)start uses the ftab
		TIndexOffU         fi;     // at this ftab index
#ifndef NDEBUG
		size_t             fidep;  // the depth they were computed for
#endif
	};
	Lane lanes[EXACT_SWEEP_LANES];

	// Before a lane's turn that (re)starts its sweep: compute the ftab
	// index the start reads and prefetch its entries, so that they arrive
	// while the other lanes take their turns
	const TIndexOffU* const ftab = ebwt.ftab();
	auto prepInit = [&](Lane& ln) {
		ln.doFtab = exactSweepFtabIdx(ebwt, *ln.seq, ftabLen, ln.len, ln.dep, ln.fi);
		if(ln.doFtab && ln.fi != std::numeric_limits<TIndexOffU>::max()) {
			__builtin_prefetch(ftab + ln.fi);     // ftabHi(fi)
			__builtin_prefetch(ftab + ln.fi + 1); // ftabLo(fi+1)
		}
#ifndef NDEBUG
		ln.fidep = ln.dep;
#endif
	};

	// Sweep s is job s/2's read (s even) or its reverse complement (s odd)
	const size_t nsweeps = 2 * njobs;
	size_t next = 0;

	// The end of a sweep: what exactSweep does after its loop
	auto finish = [&](Lane& ln) {
		if( (!ln.done) && (ln.dep >= ln.len) ) {
			ExactSweepJob& job = *ln.job;
			// Set the minimum # edits
			job.mine[ln.fwi] = ln.nedit;
			// Done
			if(ln.nedit == 0 && ln.bot > ln.top) {
				if(repex) {
					// This is an exact hit
					int64_t score = ln.len * sc.match();
					job.hit[ln.fwi] = true;
					job.top[ln.fwi] = ln.top;
					job.bot[ln.fwi] = ln.bot;
					job.score[ln.fwi] = score;
					assert(ebwt.contains(*ln.seq, NULL, NULL));
				}
				job.nelt += (ln.bot - ln.top);
			}
		}
	};

	// Start the next sweep in a lane; false if none is left
	auto start = [&](Lane& ln) {
		while(next < nsweeps) {
			ExactSweepJob& job = *jobs[next >> 1];
			const int fwi = (int)(next & 1);
			next++;
			if(fwi == 0 ? job.nofw : job.norc) {
				continue;
			}
			ln.job = &job;
			ln.fwi = fwi;
			ln.seq = fwi == 0 ? &job.read->patFw : &job.read->patRc;
			ln.len = job.read->length();
			ln.dep = 0;
			ln.nedit = 0;
			ln.top = ln.bot = 0;
			ln.doInit = true;
			ln.done = false;
			ln.prefetch_count = 0;
			if(ln.len == 0) {
				finish(ln);
				continue;
			}
			const BTDnaString& seq = *ln.seq;
			__builtin_prefetch(&(seq[ln.len-1]));
			if (ln.len>48) __builtin_prefetch(&(seq[ln.len-49])); // HW prefetch prediction assumes forward, help it
			prepInit(ln);
			return true;
		}
		ln.job = NULL;
		return false;
	};

	size_t nlanes = 0;
	for(size_t l = 0; l < EXACT_SWEEP_LANES; l++) {
		if(start(lanes[l])) {
			nlanes++;
		}
	}
	while(nlanes > 0) {
		for(size_t l = 0; l < EXACT_SWEEP_LANES; l++) {
			Lane& ln = lanes[l];
			if(ln.job == NULL) {
				continue;
			}
			const BTDnaString& seq = *ln.seq;
			if (++ln.prefetch_count >= 48) { // cache line is 64 bytes, but we may skip some deps
				const size_t left = ln.len-ln.dep;
				if (left>48) {
					__builtin_prefetch(&(seq[left-49])); // HW prefetch prediction assumes forward, help it
				}
				ln.prefetch_count = 0;
			}
			// One turn of exactSweep's loop for this sweep
			bool step = true;
			if (ln.doInit) {
				assert_eq(ln.fidep, ln.dep);
				exactSweepInitIdx(ebwt, seq, ftabLen, ln.len,    // in
						ln.doFtab, ln.fi,                       // in
						ln.dep, ln.top, ln.bot);                // out
				if ( exactSweepStep(ebwt, ln.top, ln.bot, mineMax,
						ln.tloc, ln.bloc,
						ln.job->mine[ln.fwi],
						ln.nedit, ln.done) ) {
					step = false;
				} else {
					ln.doInit=false;
					// Its first LF step reads the side INIT_LOCS has just
					// prefetched: leave it to the lane's next turn
					step = false;
				}
			}
			if (step && ln.dep < ln.len) {
				exactSweepMapLF(ebwt, seq, ln.len, ln.dep, ln.tloc, ln.bloc,
						ln.top, ln.bot, ln.job->bwops);

				if ( exactSweepStep(ebwt, ln.top, ln.bot, mineMax,
							ln.tloc, ln.bloc,
							ln.job->mine[ln.fwi],
							ln.nedit, ln.done) ) {
					ln.doInit=true;
				}
				ln.dep++;
			}
			if (ln.done || ln.dep >= ln.len) {
				// This sweep has ended; the next one takes the lane
				finish(ln);
				if(!start(ln)) {
					nlanes--;
				}
			} else if (ln.doInit) {
				// The sweep restarts on the lane's next turn
				prepInit(ln);
			}
		}
	}
}

/**
 * Compute the ftab indices oneMmStart reads for the job's directions (each
 * as ftabLoHi computes it from the read), store them in the job, and
 * prefetch their ftab entries.  Directions oneMmStart skips get an index
 * too; it is not read.
 */
void SeedAligner::oneMmStartPrefetch(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	bool               rep1mm, // report 1mm hits?
	OneMmJob&          job)
{
	const Read& read = *job.read;
	const size_t len = read.length();
	const int ftabLen = ebwtFw->eh().ftabChars();
	job.ftabPre = true;
	job.ftabPreRep1mm = rep1mm;
	if(ftabLen <= 1) {
		return; // oneMmStart uses fchr
	}
	const TIndexOffU nomatch = std::numeric_limits<TIndexOffU>::max();
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if( fw && job.nofw) continue;
		if(!fw && job.norc) continue;
		int lim = rep1mm ? 2 : 1;
		for(int ebwtfwi = 0; ebwtfwi < lim; ebwtfwi++) {
			bool ebwtfw = (ebwtfwi == 0);
			// The near half's length, as in oneMmStart
			const size_t nea = (len >> 1) + (ebwtfw ? 0 : (len & 1));
			if((size_t)ftabLen > nea) {
				continue; // oneMmStart uses fchr
			}
			const Ebwt* ebwt  = (ebwtfw ? ebwtFw : ebwtBw);
			const Ebwt* ebwtp = (ebwtfw ? ebwtBw : ebwtFw);
			const BTDnaString& seq =
				(fw ? (ebwtfw ? read.patFw : read.patFwRev) :
				      (ebwtfw ? read.patRc : read.patRcRev));
			const bool rev = !ebwtfw;
			const int k = 2 * fwi + ebwtfwi;
			const TIndexOffU fi = ebwt->ftabSeqToInt(seq, len - ftabLen, rev);
			job.fi[k] = fi;
			if(fi != nomatch) {
				__builtin_prefetch(ebwt->ftab() + fi);
				__builtin_prefetch(ebwt->ftab() + fi + 1);
			}
			if(rep1mm) {
				const TIndexOffU fip = ebwtp->ftabSeqToInt(seq, len - ftabLen, rev);
				job.fip[k] = fip;
				if(fip != nomatch) {
					__builtin_prefetch(ebwtp->ftab() + fip);
					__builtin_prefetch(ebwtp->ftab() + fip + 1);
				}
			}
		}
	}
}

/**
 * Set up the directions of a 1-mismatch search and their main walks, which
 * are the first job.ndir walks.  Returns false if the read cannot have a hit
 * (too many Ns); the job then has no walks.
 */
bool SeedAligner::oneMmStart(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	const Scoring&     sc,     // scoring
	bool               repex,  // report 0mm hits?
	bool               rep1mm, // report 1mm hits?
	OneMmJob&          job)
{
	assert(!rep1mm || ebwtBw != NULL);
	const Read& read = *job.read;
	const size_t len = read.length();
	job.repex = repex;
	job.rep1mm = rep1mm;
	job.len = len;
	job.nceil = sc.nCeil.f<int>((double)len);
	size_t ns = read.ns();
	job.ns = ns;
	if(ns > 1) {
		// Can't align this with <= 1 mismatches
		return false;
	} else if(ns == 1 && !rep1mm) {
		// Can't align this with 0 mismatches
		return false;
	}
	assert_geq(len, 2);
	assert(!rep1mm || ebwtBw->eh().ftabChars() == ebwtFw->eh().ftabChars());
#ifndef NDEBUG
	if(ebwtBw != NULL) {
		for(int i = 0; i < 4; i++) {
			assert_eq(ebwtBw->fchr()[i], ebwtFw->fchr()[i]);
		}
	}
#endif
	size_t halfFw = len >> 1;
	size_t halfBw = len >> 1;
	if((len & 1) != 0) {
		halfBw++;
	}
	assert_geq(halfFw, 1);
	assert_geq(halfBw, 1);
	int ndir = 0;
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if( fw && job.nofw) continue;
		if(!fw && job.norc) continue;
		// Align going right-to-left, left-to-right
		int lim = rep1mm ? 2 : 1;
		for(int ebwtfwi = 0; ebwtfwi < lim; ebwtfwi++) {
			bool ebwtfw = (ebwtfwi == 0);
			const Ebwt* ebwt  = (ebwtfw ? ebwtFw : ebwtBw);
			const Ebwt* ebwtp = (ebwtfw ? ebwtBw : ebwtFw);
			assert(rep1mm || ebwt->fw());
			const BTDnaString& seq =
				(fw ? (ebwtfw ? read.patFw : read.patFwRev) :
				      (ebwtfw ? read.patRc : read.patRcRev));
			assert(!seq.empty());
			const BTString& qual =
				(fw ? (ebwtfw ? read.qual    : read.qualRev) :
				      (ebwtfw ? read.qualRev : read.qual));
			int ftabLen = ebwt->eh().ftabChars();
			size_t nea = ebwtfw ? halfFw : halfBw;
			// Check if there's an N in the near portion
			bool skip = false;
			for(size_t dep = 0; dep < nea; dep++) {
				if(seq[len-dep-1] > 3) {
					skip = true;
					break;
				}
			}
			if(skip) {
				continue;
			}
			size_t dep = 0;
			TIndexOffU top = 0, bot = 0, topp = 0, botp = 0;
			// Align near half
			if(ftabLen > 1 && (size_t)ftabLen <= nea) {
				// Use ftab to jump partway into near half
				bool rev = !ebwtfw;
				if(job.ftabPre && job.ftabPreRep1mm == rep1mm) {
					// The indices ftabLoHi would compute, computed earlier
					const TIndexOffU nomatch = std::numeric_limits<TIndexOffU>::max();
					const int k = 2 * fwi + ebwtfwi;
					assert_eq(job.fi[k], ebwt->ftabSeqToInt(seq, len - ftabLen, rev));
					if(job.fi[k] != nomatch) {
						ebwt->ftabLoHi(job.fi[k], top, bot);
					}
					if(rep1mm) {
						assert_eq(job.fip[k], ebwtp->ftabSeqToInt(seq, len - ftabLen, rev));
						if(job.fip[k] != nomatch) {
							ebwtp->ftabLoHi(job.fip[k], topp, botp);
						}
						assert_eq(bot - top, botp - topp);
					}
				} else {
					ebwt->ftabLoHi(seq, len - ftabLen, rev, top, bot);
					if(rep1mm) {
						ebwtp->ftabLoHi(seq, len - ftabLen, rev, topp, botp);
						assert_eq(bot - top, botp - topp);
					}
				}
				if(bot - top == 0) {
					continue;
				}
				dep = ftabLen;
			} else {
				// Use fchr to jump in by 1 pos
				int c = seq[len-1];
				assert_range(0, 3, c);
				top = topp = ebwt->fchr()[c];
				bot = botp = ebwt->fchr()[c+1];
				if(bot - top == 0) {
					continue;
				}
				dep = 1;
			}
			assert(sanityPartial(ebwt, ebwtp, seq, len-dep, len, rep1mm, top, bot, topp, botp));
			OneMmDir& d = job.dirs[ndir];
			d.ebwt = ebwt;
			d.ebwtp = ebwtp;
			d.seq = &seq;
			d.qual = &qual;
			d.nea = nea;
			d.fw = fw;
			d.ebwtfw = ebwtfw;
			d.exact = false;
			d.etop = d.ebot = 0;
			job.walks.expand();
			job.walks.back().init(ndir, false, top, bot, topp, botp, dep);
			ndir++;
		}
	}
	job.ndir = ndir;
	job.mainLeft = (size_t)ndir;
	return true;
}

/**
 * One LF step of main walk i of a job from loci tloc and bloc, which it then
 * sets for the next step; returns false when the walk has ended.  It may
 * append mismatch walks to the job's walks, so it works on a copy.
 */
inline bool SeedAligner::oneMmStepMain(
	OneMmJob& job,
	size_t i,
	SideLocus& tloc,
	SideLocus& bloc,
	bool prefetch)
{
	OneMmWalk m = job.walks[i];
	OneMmDir& d = job.dirs[m.dir];
	const Ebwt* ebwt = d.ebwt;
	const BTDnaString& seq = *d.seq;
	const size_t len = job.len;
	const bool rep1mm = job.rep1mm;
	TIndexOffU t[4], b[4];   // dest BW ranges for BWT
	t[0] = t[1] = t[2] = t[3] = 0;
	b[0] = b[1] = b[2] = b[3] = 0;
	TIndexOffU tp[4], bp[4]; // dest BW ranges for BWT'
	tp[0] = tp[1] = tp[2] = tp[3] = m.topp;
	bp[0] = bp[1] = bp[2] = bp[3] = m.botp;
	if(m.dep < d.nea) {
		// Near half: exact matching
		assert_lt(m.dep, len);
		int rdc = seq[len - m.dep - 1];
		if(bloc.valid()) {
			job.bwops++;
			ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
			SANITY_CHECK_4TUP(t, b, tp, bp);
			m.top = t[rdc]; m.bot = b[rdc];
			if(m.bot <= m.top) {
				return false;
			}
			m.topp = tp[rdc]; m.botp = bp[rdc];
			assert(!rep1mm || m.bot - m.top == m.botp - m.topp);
		} else {
			assert_eq(m.bot, m.top+1);
			assert(!rep1mm || m.botp == m.topp+1);
			job.bwops++;
			m.top = ebwt->mapLF1(m.top, tloc, rdc);
			if(m.top == OFF_MASK) {
				return false;
			}
			m.bot = m.top + 1;
			// topp/botp stay the same
		}
		if(prefetch) {
			INIT_LOCS(m.top, m.bot, tloc, bloc, *ebwt);
		} else {
			INIT_LOCS_NOPF(m.top, m.bot, tloc, bloc, *ebwt);
		}
		assert(sanityPartial(ebwt, d.ebwtp, seq, len - m.dep - 1, len, rep1mm, m.top, m.bot, m.topp, m.botp));
		m.dep++;
		job.walks[i] = m;
		return true;
	}
	// Far half: exact matching, and at each position a mismatch walk
	// for each other character that occurs
	if(m.dep >= len) {
		return false;
	}
	int rdc = seq[len - m.dep - 1];
	int quc = (*d.qual)[len - m.dep - 1];
	if(rdc > 3 && job.nceil == 0) {
		return false;
	}
	int clo = 0, chi = 3;
	bool match = true;
	if(bloc.valid()) {
		job.bwops++;
		ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
		SANITY_CHECK_4TUP(t, b, tp, bp);
		match = rdc < 4;
		if(match) {
			m.top = t[rdc]; m.bot = b[rdc];
			m.topp = tp[rdc]; m.botp = bp[rdc];
		}
	} else {
		assert_eq(m.bot, m.top+1);
		assert(!rep1mm || m.botp == m.topp+1);
		job.bwops++;
		clo = ebwt->mapLF1(m.top, tloc);
		match = (clo == rdc);
		assert_range(-1, 3, clo);
		if(clo < 0) {
			return false; // Hit the $
		}
		t[clo] = m.top;
		b[clo] = m.bot = m.top + 1;
		bp[clo] = m.botp;
		tp[clo] = m.topp;
		assert(!rep1mm || m.bot - m.top == m.botp - m.topp);
		assert(!rep1mm || b[clo] - t[clo] == bp[clo] - tp[clo]);
		chi = clo;
	}
	if(rep1mm && (job.ns == 0 || rdc > 3)) {
		for(int j = clo; j <= chi; j++) {
			if(j == rdc || b[j] == t[j]) {
				// Either matches read or isn't a possibility
				continue;
			}
			// Potential mismatch - next, try
			assert_eq(b[j] - t[j], bp[j] - tp[j]);
			job.walks.expand();
			OneMmWalk& w = job.walks.back();
			w.init(m.dir, true, t[j], b[j], tp[j], bp[j], m.dep + 1);
			w.mmdep = (uint32_t)m.dep;
			w.mmc = (int8_t)j;
			w.rdc = (int8_t)rdc;
			w.quc = quc;
			if(w.dep == len) {
				// Mismatch at the last position: nothing left to match
				w.done = w.hit = true;
			}
		}
	}
	if(m.bot > m.top && match) {
		assert_lt(rdc, 4);
		if(m.dep == len-1) {
			// Success; this is an exact hit
			if(d.ebwtfw && job.repex) {
				d.exact = true;
				d.etop = m.top;
				d.ebot = m.bot;
			}
			return false; // End of far loop
		}
		if(prefetch) {
			INIT_LOCS(m.top, m.bot, tloc, bloc, *ebwt);
		} else {
			INIT_LOCS_NOPF(m.top, m.bot, tloc, bloc, *ebwt);
		}
		assert(sanityPartial(ebwt, d.ebwtp, seq, len - m.dep - 1, len, rep1mm, m.top, m.bot, m.topp, m.botp));
		m.dep++;
		job.walks[i] = m;
		return true;
	}
	return false; // End of far loop
}

/**
 * One LF step of mismatch walk i of a job, as oneMmStepMain.
 */
inline bool SeedAligner::oneMmStepMm(
	OneMmJob& job,
	size_t i,
	SideLocus& tloc,
	SideLocus& bloc,
	bool prefetch)
{
	OneMmWalk& w = job.walks[i];
	const OneMmDir& d = job.dirs[w.dir];
	const Ebwt* ebwt = d.ebwt;
	const size_t len = job.len;
	int rdcm = (*d.seq)[len - w.dep - 1];
	assert_lt(rdcm, 4);
	if(bloc.valid()) {
		job.bwops++;
		TIndexOffU tm[4], bm[4];   // dest BW ranges for BWT
		tm[0] = tm[1] = tm[2] = tm[3] = 0;
		bm[0] = bm[1] = bm[2] = bm[3] = 0;
		TIndexOffU tmp[4], bmp[4]; // dest BW ranges for BWT'
		tmp[0] = tmp[1] = tmp[2] = tmp[3] = w.topp;
		bmp[0] = bmp[1] = bmp[2] = bmp[3] = w.botp;
		ebwt->mapBiLFEx(tloc, bloc, tm, bm, tmp, bmp);
		SANITY_CHECK_4TUP(tm, bm, tmp, bmp);
		w.top = tm[rdcm]; w.bot = bm[rdcm];
		w.topp = tmp[rdcm]; w.botp = bmp[rdcm];
		if(w.bot <= w.top) {
			return false;
		}
	} else {
		assert_eq(w.bot, w.top+1);
		assert_eq(w.botp, w.topp+1);
		job.bwops++;
		w.top = ebwt->mapLF1(w.top, tloc, rdcm);
		if(w.top == OFF_MASK) {
			return false;
		}
		w.bot = w.top + 1;
		// topp/botp stay the same
	}
	if(++w.dep == len) {
		w.hit = true;
		return false;
	}
	if(prefetch) {
		INIT_LOCS(w.top, w.bot, tloc, bloc, *ebwt);
	} else {
		INIT_LOCS_NOPF(w.top, w.bot, tloc, bloc, *ebwt);
	}
	return true;
}

/**
 * Step the walks of the jobs in turn, at most LANES at a time.  The walks of
 * each job are started in the order they were made, job after job; starting
 * a walk sets its loci, which prefetches the sides of its first step.  A
 * job's main walks may make more walks for it until they end, so every job
 * from the first one that may still have a walk to start is looked at
 * whenever a lane is free.  The walks share no state: every walk takes the
 * same steps as it would alone.
 *
 * With GROUP, setting a lane's loci does not prefetch their sides; the
 * prefetches of every lane are issued together once per round, after the
 * walks that take free lanes are started and before the lanes are stepped,
 * so that their page walks overlap.  The same lines are prefetched.
 * oneMmSearchBatch's 32-lane runs group; oneMmSearch's one-read runs do
 * not (with 4 KB pages, grouping made the per-read runs slower).
 */
template<size_t LANES, bool GROUP>
void SeedAligner::oneMmRun(OneMmJob* const* jobs, size_t njobs)
{
	OneMmJob* actj[LANES];
	size_t act[LANES];
	SideLocus tlocs[LANES], blocs[LANES];
	const uint8_t* ebs[LANES]; // with GROUP: the BWT each lane's walk steps in
	size_t nact = 0, jfirst = 0;
	while(true) {
		while(jfirst < njobs && jobs[jfirst]->mainLeft == 0 &&
		      jobs[jfirst]->nextw == jobs[jfirst]->walks.size())
		{
			jfirst++;
		}
		for(size_t j = jfirst; j < njobs && nact < LANES; j++) {
			OneMmJob& job = *jobs[j];
			while(nact < LANES && job.nextw < job.walks.size()) {
				const OneMmWalk& w = job.walks[job.nextw];
				if(!w.done) {
					const Ebwt& e = *job.dirs[w.dir].ebwt;
					if(GROUP) {
						INIT_LOCS_NOPF(w.top, w.bot, tlocs[nact], blocs[nact], e);
						ebs[nact] = e.ebwt();
					} else {
						INIT_LOCS(w.top, w.bot, tlocs[nact], blocs[nact], e);
					}
					actj[nact] = &job;
					act[nact++] = job.nextw;
				}
				job.nextw++;
			}
		}
		if(nact == 0) {
			break;
		}
		if(GROUP) {
			// Every lane's loci were set, without prefetches, by its step
			// in the last round or by the start of its walk above
			for(size_t k = 0; k < nact; k++) {
				SIDE_LOCI_PREFETCH(tlocs[k], blocs[k], ebs[k]);
			}
		}
		for(size_t k = 0; k < nact;) {
			OneMmJob& job = *actj[k];
			const size_t i = act[k];
			const bool mm = job.walks[i].mm;
			if(mm ? oneMmStepMm(job, i, tlocs[k], blocs[k], !GROUP) :
			        oneMmStepMain(job, i, tlocs[k], blocs[k], !GROUP))
			{
				k++;
			} else {
				job.walks[i].done = true;
				if(!mm) {
					job.mainLeft--;
				}
				nact--;
				actj[k] = actj[nact];
				act[k] = act[nact];
				tlocs[k] = tlocs[nact];
				blocs[k] = blocs[nact];
				if(GROUP) {
					ebs[k] = ebs[nact];
				}
			}
		}
	}
}

/**
 * Report the hits of a job whose walks have all ended, in the order of a
 * one-walk-at-a-time search: by direction, then mismatch position, then
 * character.  Mismatch walks are made in that order within a direction.
 */
bool SeedAligner::oneMmReport(
	const OneMmJob&    job,
	const Ebwt*        ebwtFw, // BWT index
	const Scoring&     sc,     // scoring
	int64_t            minsc,  // minimum score
	bool               local,  // 1mm hits must be legal local alignments
	SeedResults&       hits)   // holds all the seed hits (and exact hit)
{
	const Read& read = *job.read;
	(void)read;
	(void)ebwtFw;  // read only by the debug checks
	const size_t len = job.len;
	const int ndir = job.ndir;
	const bool repex = job.repex;
	(void)repex;
	bool results = false;
	for(int di = 0; di < ndir; di++) {
		const OneMmDir& d = job.dirs[di];
		const bool fw = d.fw;
		const bool ebwtfw = d.ebwtfw;
		for(size_t i = (size_t)ndir; i < job.walks.size(); i++) {
			const OneMmWalk& w = job.walks[i];
			if(w.dir != di || !w.hit) {
				continue;
			}
			assert(w.mm);
			const size_t dep = w.mmdep;
			const int j = w.mmc;
			const int rdc = w.rdc;
			const int quc = w.quc;
			const TIndexOffU topm = w.top, botm = w.bot;
			const TIndexOffU topmp = w.topp, botmp = w.botp;
			// Success; this is a 1MM hit
			size_t off5p = dep;  // offset from 5' end of read
			if(fw == ebwtfw) {
				off5p = len - off5p - 1;
			}
			Edit e((uint32_t)off5p, j, rdc, EDIT_TYPE_MM, false);
			results = true;
			int64_t score = (len - 1) * sc.match();
			// In --local mode, need to double-check that
			// end-to-end alignment doesn't violate  local
			// alignment principles.  Specifically, it
			// shouldn't to or below 0 anywhere in the middle.
			int pen = sc.score(rdc, (int)(1 << j), quc - 33);
			score += pen;
			bool valid = true;
			if(local) {
				int64_t locscore_fw = 0, locscore_bw = 0;
				for(size_t i = 0; i < len; i++) {
					if(i == dep) {
						if(locscore_fw + pen <= 0) {
							valid = false;
							break;
						}
						locscore_fw += pen;
					} else {
						locscore_fw += sc.match();
					}
					if(len-i-1 == dep) {
						if(locscore_bw + pen <= 0) {
							valid = false;
							break;
						}
						locscore_bw += pen;
					} else {
						locscore_bw += sc.match();
					}
				}
			}
			if(valid) {
				valid = score >= minsc;
			}
			if(valid) {
#ifndef NDEBUG
				BTDnaString& rf = tmprfdnastr_;
				rf.clear();
				edits_.clear();
				edits_.push_back(e);
				if(!fw) Edit::invertPoss(edits_, len, false);
				Edit::toRef(fw ? read.patFw : read.patRc, edits_, rf);
				if(!fw) Edit::invertPoss(edits_, len, false);
				assert_eq(len, rf.length());
				for(size_t i = 0; i < len; i++) {
					assert_lt((int)rf[i], 4);
				}
				ASSERT_ONLY(TIndexOffU toptmp = 0);
				ASSERT_ONLY(TIndexOffU bottmp = 0);
				assert(ebwtFw->contains(rf, &toptmp, &bottmp));
#endif
				TIndexOffU toprep = ebwtfw ? topm : topmp;
				TIndexOffU botrep = ebwtfw ? botm : botmp;
				assert_eq(toprep, toptmp);
				assert_eq(botrep, bottmp);
				hits.add1mmEe(toprep, botrep, &e, NULL, fw, score);
			}
		}
		if(d.exact) {
			assert(ebwtfw && repex);
			const BTDnaString& seq = *d.seq;
			(void)seq;
			results = true;
			int64_t score = len * sc.match();
			if(fw) {
				hits.addExactEeFw(d.etop, d.ebot, NULL, NULL, fw, score);
				assert(ebwtFw->contains(seq, NULL, NULL));
			} else {
				hits.addExactEeRc(d.etop, d.ebot, NULL, NULL, fw, score);
				assert(ebwtFw->contains(seq, NULL, NULL));
			}
		}
	}
	return results;
}

/**
 * Search for end-to-end alignments with up to 1 mismatch.  Return true iff
 * one is found.
 *
 * There are up to four directions: each strand of the read, matched
 * right-to-left in the forward index and, when looking for 1-mismatch hits,
 * in the mirror index.  Each direction has a main walk that matches the read
 * exactly, and every mismatch the main walk tries in the far half starts
 * another walk that matches the rest of the read exactly.  These walks are
 * independent of each other, and each LF step waits on a load from a random
 * place in the index.  So the walks are stepped in turn, one LF step each:
 * each step prefetches the sides its walk needs next, and the steps of the
 * other walks hide that load.  Every walk takes the same steps as it would
 * alone, and the hits are reported afterwards in the order in which a
 * one-walk-at-a-time search finds them.
 */
bool SeedAligner::oneMmSearch(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	const Read&        read,   // read to align
	const Scoring&     sc,     // scoring
	int64_t            minsc,  // minimum score
	bool               nofw,   // don't align forward read
	bool               norc,   // don't align revcomp read
	bool               local,  // 1mm hits must be legal local alignments
	bool               repex,  // report 0mm hits?
	bool               rep1mm, // report 1mm hits?
	SeedResults&       hits,   // holds all the seed hits (and exact hit)
	SeedSearchMetrics& met)    // metrics
{
	OneMmJob& job = mmJob_;
	job.init(&read, nofw, norc);
	if(!oneMmStart(ebwtFw, ebwtBw, sc, repex, rep1mm, job)) {
		return false;
	}
	OneMmJob* jobs[1] = { &job };
	oneMmRun<ONEMM_WALKS, false>(jobs, 1);
	bwops_ += job.bwops;
	return oneMmReport(job, ebwtFw, sc, minsc, local, hits);
}

/**
 * oneMmSearch's walks for many reads at once.  oneMmSearch steps the walks
 * of one read in turn, but most reads have only two or three walks going at
 * a time.  Here the walks of all the jobs are stepped in turn in
 * ONEMM_BATCH_LANES lanes.  Each walk takes the same steps as in
 * oneMmSearch, so each job ends with the walks, hits and LF operations
 * oneMmSearch has for its read; oneMmSearchApply reports them.
 */
void SeedAligner::oneMmSearchBatch(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	const Scoring&     sc,     // scoring
	bool               repex,  // report 0mm hits?
	bool               rep1mm, // report 1mm hits?
	OneMmJob* const*   jobs,   // the reads to search
	size_t             njobs)  // # jobs
{
	for(size_t j = 0; j < njobs; j++) {
		oneMmStartPrefetch(ebwtFw, ebwtBw, rep1mm, *jobs[j]);
	}
	for(size_t j = 0; j < njobs; j++) {
		oneMmStart(ebwtFw, ebwtBw, sc, repex, rep1mm, *jobs[j]);
	}
	oneMmRun<ONEMM_BATCH_LANES, true>(jobs, njobs);
}

inline void
SeedAligner::prefetchNextLocsBi(
	const InstantiatedSeed& seed, // current instantiated seed
        TIndexOffU topf,              // top in BWT
        TIndexOffU botf,              // bot in BWT
        TIndexOffU topb,              // top in BWT'
        TIndexOffU botb,              // bot in BWT'
        int step                    // step to get ready for
        )
{
	if(step == (int)seed.steps.size()) return; // no more steps!
	// Which direction are we going in next?
	if(seed.steps[step] > 0) {
		// Left to right; use BWT'
		if(botb - topb == 1) {
			// Already down to 1 row; just init top locus
			SideLocus::prefetchFromRow(
				topb, ebwtBw_->eh(), ebwtBw_->ebwt());
		} else {
			SideLocus::prefetchFromTopBot(
				topb, botb, ebwtBw_->eh(), ebwtBw_->ebwt());
		}
	} else {
		// Right to left; use BWT
		if(botf - topf == 1) {
			// Already down to 1 row; just init top locus
			SideLocus::prefetchFromRow(
				topf, ebwtFw_->eh(), ebwtFw_->ebwt());
		} else {
			SideLocus::prefetchFromTopBot(
				topf, botf, ebwtFw_->eh(), ebwtFw_->ebwt());
		}
	}
}

/**
 * Get tloc, bloc ready for the next step.  If the new range is under
 * the ceiling.
 */
inline void
SeedAligner::nextLocsBi(
	const InstantiatedSeed& seed, // current instantiated seed
	SideLocus& tloc,              // top locus
	SideLocus& bloc,              // bot locus
	TIndexOffU topf,              // top in BWT
	TIndexOffU botf,              // bot in BWT
	TIndexOffU topb,              // top in BWT'
	TIndexOffU botb,              // bot in BWT'
	int step                    // step to get ready for
#if 0
	, const SABWOffTrack* prevOt, // previous tracker
	SABWOffTrack& ot            // current tracker
#endif
	)
{
	assert_gt(botf, 0);
	assert(ebwtBw_ == NULL || botb > 0);
	assert_geq(step, 0); // next step can't be first one
	assert(ebwtBw_ == NULL || botf-topf == botb-topb);
	if(step == (int)seed.steps.size()) return; // no more steps!
	// Which direction are we going in next?
	if(seed.steps[step] > 0) {
		// Left to right; use BWT'
		if(botb - topb == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topb, ebwtBw_->eh(), ebwtBw_->ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topb, botb, ebwtBw_->eh(), ebwtBw_->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	} else {
		// Right to left; use BWT
		if(botf - topf == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topf, ebwtFw_->eh(), ebwtFw_->ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topf, botf, ebwtFw_->eh(), ebwtFw_->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	// Check if we should update the tracker with this refinement
#if 0
	if(botf-topf <= BW_OFF_TRACK_CEIL) {
		if(ot.size() == 0 && prevOt != NULL && prevOt->size() > 0) {
			// Inherit state from the predecessor
			ot = *prevOt;
		}
		bool ltr = seed.steps[step-1] > 0;
		int adj = abs(seed.steps[step-1])-1;
		const Ebwt* ebwt = ltr ? ebwtBw_ : ebwtFw_;
		ot.update(
			ltr ? topb : topf,    // top
			ltr ? botb : botf,    // bot
			adj,                  // adj (to be subtracted from offset)
			ebwt->offs(),         // offs array
			ebwt->eh().offRate(), // offrate (sample = every 1 << offrate elts)
			NULL                  // dead
		);
		assert_gt(ot.size(), 0);
	}
#endif
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
}

/**
 * Report a seed hit found by searchSeedBi(), but first try to extend it out in
 * either direction as far as possible without hitting any edits.  This will
 * allow us to prioritize the seed hits better later on.  Call reportHit() when
 * we're done, which actually adds the hit to the cache.  Returns result from
 * calling reportHit().
 */
void
SeedAligner::extendAndReportHit(
	SeedSearchCache &cache,              // local seed alignment cache
	size_t off,                          // offset of seed currently being searched
	bool fw,                             // orientation of seed currently being searched
	TIndexOffU topf,                     // top in BWT
	TIndexOffU botf,                     // bot in BWT
	TIndexOffU topb,                     // top in BWT'
	TIndexOffU botb,                     // bot in BWT'
	uint16_t len,                      // length of hit
	DoublyLinkedList<Edit> *prevEdit)  // previous edit
{
	size_t nlex = 0, nrex = 0;
	TIndexOffU t[4], b[4];
	TIndexOffU tp[4], bp[4];
	SideLocus tloc, bloc;
	if(off > 0) {
		const Ebwt *ebwt = ebwtFw_;
		assert(ebwt != NULL);
		// Extend left using forward index
		const BTDnaString& seq = fw ? read_->patFw : read_->patRc;
		// See what we get by extending 
		TIndexOffU top = topf, bot = botf;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topb;
		bp[0] = bp[1] = bp[2] = bp[3] = botb;
		SideLocus tloc, bloc;
		INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		for(size_t ii = off; ii > 0; ii--) {
			size_t i = ii-1;
			// Get char from read
			int rdc = seq.get(i);
			// See what we get by extending 
			if(bloc.valid()) {
				bwops_++;
				t[0] = t[1] = t[2] = t[3] =
				b[0] = b[1] = b[2] = b[3] = 0;
				ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				for(int j = 0; j < 4; j++) {
					if(b[i] > t[i]) {
						if(nonz >= 0) {
							abort = true;
							break;
						}
						nonz = j;
						top = t[i]; bot = b[i];
					}
				}
				if(abort || nonz != rdc) {
					break;
				}
			} else {
				assert_eq(bot, top+1);
				bwops_++;
				int c = ebwt->mapLF1(top, tloc);
				if(c != rdc) {
					break;
				}
				bot = top + 1;
			}
			if(++nlex == 255) {
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		}
	}
	size_t rdlen = read_->length();
	size_t nright = rdlen - off - len;
	if(nright > 0 && ebwtBw_ != NULL) {
		const Ebwt *ebwt = ebwtBw_;
		assert(ebwt != NULL);
		// Extend right using backward index
		const BTDnaString& seq = fw ? read_->patFw : read_->patRc;
		// See what we get by extending 
		TIndexOffU top = topb, bot = botb;
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		tp[0] = tp[1] = tp[2] = tp[3] = topb;
		bp[0] = bp[1] = bp[2] = bp[3] = botb;
		INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		for(size_t i = off + len; i < rdlen; i++) {
			// Get char from read
			int rdc = seq.get(i);
			// See what we get by extending 
			if(bloc.valid()) {
				bwops_++;
				t[0] = t[1] = t[2] = t[3] =
				b[0] = b[1] = b[2] = b[3] = 0;
				ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
				SANITY_CHECK_4TUP(t, b, tp, bp);
				int nonz = -1;
				bool abort = false;
				for(int j = 0; j < 4; j++) {
					if(b[i] > t[i]) {
						if(nonz >= 0) {
							abort = true;
							break;
						}
						nonz = j;
						top = t[i]; bot = b[i];
					}
				}
				if(abort || nonz != rdc) {
					break;
				}
			} else {
				assert_eq(bot, top+1);
				bwops_++;
				int c = ebwt->mapLF1(top, tloc);
				if(c != rdc) {
					break;
				}
				bot = top + 1;
			}
			if(++nrex == 255) {
				break;
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
		}
	}
	assert_lt(nlex, rdlen);
	assert_leq(nlex, off);
	assert_lt(nrex, rdlen);
	reportHit(cache, topf, botf, topb, botb, len, prevEdit);
}

/**
 * Report a seed hit found by searchSeedBi() by adding it to the cache.  Return
 * false if the hit could not be reported because of, e.g., cache exhaustion.
 */
void
SeedAligner::reportHit(
	SeedSearchCache &cache,              // local seed alignment cache
	TIndexOffU topf,                     // top in BWT
	TIndexOffU botf,                     // bot in BWT
	TIndexOffU topb,                     // top in BWT'
	TIndexOffU botb,                     // bot in BWT'
	uint16_t len,                      // length of hit
	DoublyLinkedList<Edit> *prevEdit)  // previous edit
{
	if ( (botf - topf) > ncut_f_) return; // low quality seed, just ignore

	const BTDnaString& seq = cache.getSeq();

	// Add information about the seed hit to AlignmentCache.  This
	// information eventually makes its way back to the SeedResults
	// object when we call finishAlign(...).
	BTDnaString& rf = tmprfdnastr_;
	rf.clear();
	edits_.clear();
	if(prevEdit != NULL) {
		prevEdit->toList(edits_);
		Edit::sort(edits_);
		assert(Edit::repOk(edits_, seq));
		Edit::toRef(seq, edits_, rf);
	} else {
		rf = seq;
	}
	// Note: Disabled, as we now use memory cache
	// Sanity check: shouldn't add the same hit twice.  If this
	// happens, it may be because our zone Constraints are not set up
	// properly and erroneously return true from acceptable() when they
	// should return false in some cases.
	//assert_eq(hits_.size(), cache.curNumRanges());
	//assert(hits_.insert(rf));
	cache.addOnTheFly(rf, topf, botf, topb, botb);
	//assert_eq(hits_.size(), cache.curNumRanges());
#ifndef NDEBUG
	// Sanity check that the topf/botf and topb/botb ranges really
	// correspond to the reference sequence aligned to
	{
		BTDnaString rfr;
		TIndexOffU tpf, btf, tpb, btb;
		tpf = btf = tpb = btb = 0;
		assert(ebwtFw_->contains(rf, &tpf, &btf));
		if(ebwtBw_ != NULL) {
			rfr = rf;
			rfr.reverse();
			assert(ebwtBw_->contains(rfr, &tpb, &btb));
			assert_eq(tpf, topf);
			assert_eq(btf, botf);
			assert_eq(tpb, topb);
			assert_eq(btb, botb);
		}
	}
#endif
	return;
}

// return true, if we are already done
bool
SeedAligner::startSearchSeedBi(SeedAligner::SeedAlignerSearchParams &p)
{
	SeedSearchCache &cache = p.cs.cache;
	const InstantiatedSeed& seed = p.cs.seed;
	const BTDnaString& seq = cache.getSeq();

	assert_gt(seed.steps.size(), 0);
	assert(ebwtBw_ == NULL || ebwtBw_->eh().ftabChars() == ebwtFw_->eh().ftabChars());
#ifndef NDEBUG
	for(int i = 0; i < 4; i++) {
		assert(ebwtBw_ == NULL || ebwtBw_->fchr()[i] == ebwtFw_->fchr()[i]);
	}
#endif
	if(p.step == (int)seed.steps.size()) {
		// Finished aligning seed
		p.checkCV();
		reportHit(cache, p.bwt, seq.length(), p.prevEdit);
		return true;
	}
#ifndef NDEBUG
	if(p.depth > 0) {
		assert(p.bwt.botf - p.bwt.topf == 1 ||  p.bloc.valid());
		assert(p.bwt.botf - p.bwt.topf > 1  || !p.bloc.valid());
	}
#endif
	if(p.step == 0) {
		// Just starting
		assert(p.prevEdit == NULL);
		assert(!p.tloc.valid());
		assert(!p.bloc.valid());
		int off = seed.steps[0];
		off = abs(off)-1;
		// Check whether/how far we can jump using ftab or fchr
		int ftabLen = ebwtFw_->eh().ftabChars();
		if(p.cs.hasi0) { //if(ftabLen > 1 && ftabLen <= seed.maxjump)
			ebwtFw_->ftabLoHi(p.cs.fwi0, p.bwt.topf, p.bwt.botf);
			#ifdef NDEBUG
			if(p.bwt.botf - p.bwt.topf == 0) return true;
			#endif
			#ifdef NDEBUG
			if(ebwtBw_ != NULL) {
				p.bwt.topb = ebwtBw_->ftabHi(p.cs.bwi0);
				p.bwt.botb = p.bwt.topb + (p.bwt.botf-p.bwt.topf);
			}
			#else
			if(ebwtBw_ != NULL) {
				ebwtBw_->ftabLoHi(p.cs.bwi0, p.bwt.topb, p.bwt.botb);
				assert_eq(p.bwt.botf-p.bwt.topf, p.bwt.botb-p.bwt.topb);
			}
			if(p.bwt.botf - p.bwt.topf == 0) return true;
			#endif
			p.step += ftabLen;
		} else if(seed.maxjump > 0) {
			// Use fchr
			int c = seq[off];
			assert_range(0, 3, c);
			p.bwt.topf = p.bwt.topb = ebwtFw_->fchr()[c];
			p.bwt.botf = p.bwt.botb = ebwtFw_->fchr()[c+1];
			if(p.bwt.botf - p.bwt.topf == 0) return true;
			p.step++;
		} else {
			assert_eq(0, seed.maxjump);
			p.bwt.topf = p.bwt.topb = 0;
			p.bwt.botf = p.bwt.botb = ebwtFw_->fchr()[4];
		}
		if(p.step == (int)seed.steps.size()) {
			// Finished aligning seed
			p.checkCV();
			reportHit(cache, p.bwt, seq.length(), p.prevEdit);
			return true;
		}
		nextLocsBi(seed, p.tloc, p.bloc, p.bwt, p.step);
		assert(p.tloc.valid());
	} else assert(p.prevEdit != NULL);
	assert(p.tloc.valid());
	assert(p.bwt.botf - p.bwt.topf == 1 ||  p.bloc.valid());
	assert(p.bwt.botf - p.bwt.topf > 1  || !p.bloc.valid());
	assert_geq(p.step, 0);

	return false;
}

class SeedAlignerSearchSave {
public:
	SeedAlignerSearchSave(
		Constraint &cons, Constraint &ovCons,
		SideLocus &tloc, SideLocus &bloc)
	: orgCons(cons), orgOvCons(ovCons)
	, orgTloc(tloc), orgBloc(bloc)
	, oldCons(cons), oldOvCons(ovCons)
	, oldTloc(tloc), oldBloc(bloc)
	{}

	~SeedAlignerSearchSave() {
		orgCons = oldCons; orgOvCons = oldOvCons;
		orgTloc = oldTloc; orgBloc = oldBloc;
	}

private:
	// reference to the original variables
	Constraint &orgCons;
	Constraint &orgOvCons;
	SideLocus &orgTloc;
	SideLocus &orgBloc;

	// copy of the originals
	const Constraint oldCons;
	const Constraint oldOvCons;
	const SideLocus oldTloc;
	const SideLocus oldBloc;
};

// State used for recursion
class SeedAlignerSearchRecState {
public:
	SeedAlignerSearchRecState(
			int j,
			int c,
			const SeedAlignerSearchState &sstate,
			DoublyLinkedList<Edit> *prevEdit)
	: bwt(sstate.tf[j], sstate.bf[j], sstate.tb[j], sstate.bb[j])
	, edit(sstate.off, j, c, EDIT_TYPE_MM, false)
	, editl()
	, _prevEdit(prevEdit)
	{
		assert(_prevEdit == NULL || _prevEdit->next == NULL);
		editl.payload = edit;
		if(_prevEdit != NULL) {
			_prevEdit->next = &editl;
			editl.prev = _prevEdit;
		}
		assert(editl.next == NULL);
	}

	~SeedAlignerSearchRecState() {
		if(_prevEdit != NULL) _prevEdit->next = NULL;
	}

	BwtTopBot bwt;
	Edit edit;
	DoublyLinkedList<Edit> editl;
private:
	DoublyLinkedList<Edit> *_prevEdit;
};

/**
 * Given a seed, search.  Assumes zone 0 = no backtracking.
 *
 * Return a list of Seed hits.
 * 1. Edits
 * 2. Bidirectional BWT range(s) on either end
 */
void
SeedAligner::searchSeedBi(
	const size_t nparams,
	SeedAligner::SeedAlignerSearchParams paramVec[],
	SeedAlignerSearchState sstateVec[])
{
	size_t nleft = nparams; // will keep track of how many are not done yet

	for (size_t n=0; n<nparams; n++) {
		SeedAlignerSearchParams& p= paramVec[n];
		const bool done = startSearchSeedBi(p);

		SeedAlignerSearchState& sstate = sstateVec[n];
		sstate.done = done;
		if(done) {
			nleft--;
		} else {
			sstate.initLastTot(p.bwt.botf - p.bwt.topf);
		}
	}

	while (nleft>0) {
	   // Note: We can do the params in any order we want
	   // but we must do the steps inside the same param in order
	   // Will loop over all of them, and just check which ones are invalid
           for (size_t n=0; n<nparams; n++) {
                SeedAlignerSearchState& sstate = sstateVec[n];
		if (sstate.done) continue;

		SeedAlignerSearchParams& p= paramVec[n];
		const InstantiatedSeed& seed = p.cs.seed;
		if (p.step >= (int) seed.steps.size()) {
			sstate.done = true;
			nleft--;
			continue;
		}
		size_t i = p.step; // call the stepIdx i for historical reasons
		p.step++; // get ready for the next iteration

		SeedSearchCache &cache = p.cs.cache;
		const BTDnaString& seq = cache.getSeq();
		const BTString& qual = cache.getQual();

		assert_gt(p.bwt.botf, p.bwt.topf);
		assert(p.bwt.botf - p.bwt.topf == 1 ||  p.bloc.valid());
		assert(p.bwt.botf - p.bwt.topf > 1  || !p.bloc.valid());
		assert(ebwtBw_ == NULL || p.bwt.botf-p.bwt.topf == p.bwt.botb-p.bwt.topb);
		assert(p.tloc.valid());
		sstate.setOff(seed.steps[i], p.bwt, ebwtFw_, ebwtBw_);
		__builtin_prefetch(&(seq[sstate.off]));
		__builtin_prefetch(&(qual[sstate.off]));
		if(p.bloc.valid()) {
			// Range delimited by tloc/bloc has size >1.  If size == 1,
			// we use a simpler query (see if(!bloc.valid()) blocks below)
			bwops_++;
			sstate.ebwt->mapBiLFEx(p.tloc, p.bloc, sstate.t, sstate.b, sstate.tp, sstate.bp);
			ASSERT_ONLY(TIndexOffU tot = (sstate.b[0]-sstate.t[0])+(sstate.b[1]-sstate.t[1])+(sstate.b[2]-sstate.t[2])+(sstate.b[3]-sstate.t[3]));
			ASSERT_ONLY(TIndexOffU totp = (sstate.bp[0]-sstate.tp[0])+(sstate.bp[1]-sstate.tp[1])+(sstate.bp[2]-sstate.tp[2])+(sstate.bp[3]-sstate.tp[3]));
			assert_eq(tot, totp);
#ifndef NDEBUG
			sstate.assertLeqAndSetLastTot(tot);
#endif
		}
		int c = seq[sstate.off];  assert_range(0, 4, c);
		//
		bool leaveZone = seed.zones[i].first < 0;
		//bool leaveZoneIns = zones_[i].second < 0;
		Constraint& cons    = p.cv[abs(seed.zones[i].first)];
		//Constraint& insCons = p.cv[abs(seed.zones[i].second)];
		// Is it legal for us to advance on characters other than 'c'?
		if(!(cons.mustMatch() && !p.overall.mustMatch()) || c == 4) {
			// There may be legal edits
			bool bail = false;
			if(!p.bloc.valid()) {
				// Range delimited by tloc/bloc has size 1
				bwops_++;
				int cc = sstate.ebwt->mapLF1(sstate.ntop, p.tloc);
				assert_range(-1, 3, cc);
				if(cc < 0) bail = true;
				else { sstate.t[cc] = sstate.ntop; sstate.b[cc] = sstate.ntop+1; }
			}
			if(!bail) {
				int q = qual[sstate.off];
				if((cons.canMismatch(q, *sc_) && p.overall.canMismatch(q, *sc_)) || c == 4) {
					SeedAlignerSearchSave save(cons, p.overall, p.tloc, p.bloc);
					if(c != 4) {
						cons.chargeMismatch(q, *sc_);
						p.overall.chargeMismatch(q, *sc_);
					}
					// Can leave the zone as-is
					if(!leaveZone || (cons.acceptable() && p.overall.acceptable())) {
						for(int j = 0; j < 4; j++) {
							if(j == c || sstate.b[j] == sstate.t[j]) continue;
							// Potential mismatch
							SeedAlignerSearchRecState rstate(j, c, sstate, p.prevEdit);
							nextLocsBi(seed, p.tloc, p.bloc, rstate.bwt, i+1);
							bwedits_++;
							SeedAlignerSearchParams p2(
								p.cs,
								i+1,             // depth into steps[] array
								p.depth+1,       // recursion depth
								rstate.bwt,      // The 4 BWT idxs
								p.tloc,          // locus for top (perhaps unititialized)
								p.bloc,          // locus for bot (perhaps unititialized)
								p.cv,            // constraints to enforce in seed zones
								p.overall,       // overall constraints to enforce
								&rstate.editl);  // latest edit
							// recursion is rare, so just do one at a time
							SeedAlignerSearchState sstate2;
							searchSeedBi(1, &p2, &sstate2);
							// as rstate gets out of scope, p.prevEdit->next is updated
						}
					} else {
						// Not enough edits to make this path
						// non-redundant with other seeds
					}

					// as olds gets out of scope,
					// restores cons, p.overall, p.tloc, p.bloc
				}
				if(cons.canGap() && p.overall.canGap()) {
					throw 1; // TODO
//					int delEx = 0;
//					if(cons.canDelete(delEx, *sc_) && overall.canDelete(delEx, *sc_)) {
//						// Try delete
//					}
//					int insEx = 0;
//					if(insCons.canInsert(insEx, *sc_) && overall.canInsert(insEx, *sc_)) {
//						// Try insert
//					}
				}
			} // if(!bail)
		}
		if(c == 4) { // couldn't handle the N
			sstate.done = true;
			nleft--;
			continue;
		}
		if(leaveZone && (!cons.acceptable() || !p.overall.acceptable())) {
			// Not enough edits to make this path non-redundant with
			// other seeds
			sstate.done = true;
			nleft--;
			continue;
		}
		if(!p.bloc.valid()) {
			assert(ebwtBw_ == NULL || sstate.bp[c] == sstate.tp[c]+1);
			// Range delimited by tloc/bloc has size 1
			bwops_++;
			sstate.t[c] = sstate.ebwt->mapLF1(sstate.ntop, p.tloc, c);
			if(sstate.t[c] == OFF_MASK) {
				sstate.done = true;
				nleft--;
				continue;
			}
			assert_geq(sstate.t[c], sstate.ebwt->fchr()[c]);
			assert_lt(sstate.t[c],  sstate.ebwt->fchr()[c+1]);
			sstate.b[c] = sstate.t[c]+1;
			assert_gt(sstate.b[c], 0);
		}
		assert(ebwtBw_ == NULL || sstate.bf[c]-sstate.tf[c] == sstate.bb[c]-sstate.tb[c]);
#ifndef NDEBUG
		sstate.assertLeqAndSetLastTot(sstate.bf[c]-sstate.tf[c]);
#endif
		if(sstate.b[c] == sstate.t[c]) {
			sstate.done = true;
			nleft--;
			continue;
		}
		p.bwt.set(sstate.tf[c], sstate.bf[c], sstate.tb[c], sstate.bb[c]);
		if(i+1 == seed.steps.size()) {
			// Finished aligning seed
			p.checkCV();
			reportHit(cache, p.bwt, seq.length(), p.prevEdit);
			sstate.done = true;
			nleft--;
			continue;
		}
		nextLocsBi(seed, p.tloc, p.bloc, p.bwt, i+1);
	   } // for n
	} // while
	return;
}

#ifdef ALIGNER_SEED_MAIN

#include <getopt.h>
#include <string>

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(const char *errmsg, const char *arg) {
	long l;
	char *endPtr = NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}

enum {
	ARG_NOFW = 256,
	ARG_NORC,
	ARG_MM,
	ARG_SHMEM,
	ARG_TESTS,
	ARG_RANDOM_TESTS,
	ARG_SEED
};

static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char*)"verbose",  no_argument,       0, 'v'},
	{(char*)"timing",   no_argument,       0, 't'},
	{(char*)"nofw",     no_argument,       0, ARG_NOFW},
	{(char*)"norc",     no_argument,       0, ARG_NORC},
	{(char*)"mm",       no_argument,       0, ARG_MM},
	{(char*)"shmem",    no_argument,       0, ARG_SHMEM},
	{(char*)"tests",    no_argument,       0, ARG_TESTS},
	{(char*)"random",   required_argument, 0, ARG_RANDOM_TESTS},
	{(char*)"seed",     required_argument, 0, ARG_SEED},
};

static void printUsage(ostream& os) {
	os << "Usage: ac [options]* <index> <patterns>" << endl;
	os << "Options:" << endl;
	os << "  --mm                memory-mapped mode" << endl;
	os << "  --shmem             shared memory mode" << endl;
	os << "  --nofw              don't align forward-oriented read" << endl;
	os << "  --norc              don't align reverse-complemented read" << endl;
	os << "  -t/--timing         show timing information" << endl;
	os << "  -v/--verbose        talkative mode" << endl;
}

bool gNorc = false;
bool gNofw = false;
int gVerbose = 0;
int gGapBarrier = 1;
int gSnpPhred = 30;
bool gReportOverhangs = true;

extern void aligner_seed_tests();
extern void aligner_random_seed_tests(
	int num_tests,
	TIndexOffU qslo,
	TIndexOffU qshi,
	uint32_t seed);

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {
	bool useMm = false;
	bool useShmem = false;
	bool mmSweep = false;
	bool noRefNames = false;
	bool sanity = false;
	bool timing = false;
	int option_index = 0;
	int seed = 777;
	int next_option;
	do {
		next_option = getopt_long(
			argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'v':       gVerbose = true; break;
			case 't':       timing   = true; break;
			case ARG_NOFW:  gNofw    = true; break;
			case ARG_NORC:  gNorc    = true; break;
			case ARG_MM:    useMm    = true; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_SEED:  seed = parseInt("", optarg); break;
			case ARG_TESTS: {
				aligner_seed_tests();
				aligner_random_seed_tests(
					100,     // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					18);   // pseudo-random seed
				return 0;
			}
			case ARG_RANDOM_TESTS: {
				seed = parseInt("", optarg);
				aligner_random_seed_tests(
					100,   // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					seed); // pseudo-random seed
				return 0;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	char *reffn;
	if(optind >= argc) {
		cerr << "No reference; quitting..." << endl;
		return 1;
	}
	reffn = argv[optind++];
	if(optind >= argc) {
		cerr << "No reads; quitting..." << endl;
		return 1;
	}
	string ebwtBase(reffn);
	BitPairReference ref(
		ebwtBase,    // base path
		false,       // whether we expect it to be colorspace
		sanity,      // whether to sanity-check reference as it's loaded
		NULL,        // fasta files to sanity check reference against
		NULL,        // another way of specifying original sequences
		false,       // true -> infiles (2 args ago) contains raw seqs
		useMm,       // use memory mapping to load index?
		useShmem,    // use shared memory (not memory mapping)
		mmSweep,     // touch all the pages after memory-mapping the index
		gVerbose,    // verbose
		gVerbose);   // verbose but just for startup messages
	Timer *t = new Timer(cerr, "Time loading fw index: ", timing);
	Ebwt ebwtFw(
		ebwtBase,
		false,       // index is colorspace
		0,           // don't need entireReverse for fw index
		true,        // index is for the forward direction
		-1,          // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	t = new Timer(cerr, "Time loading bw index: ", timing);
	Ebwt ebwtBw(
		ebwtBase + ".rev",
		false,       // index is colorspace
		1,           // need entireReverse
		false,       // index is for the backward direction
		-1,          // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		false,       // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	for(int i = optind; i < argc; i++) {
	}
}
#endif
