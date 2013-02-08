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

#ifndef ALIGNER_SEED_H_
#define ALIGNER_SEED_H_

#include <iostream>
#include <utility>
#include <limits>
#include "qual.h"
#include "ds.h"
#include "sstring.h"
#include "alphabet.h"
#include "edit.h"
#include "read.h"
// Threading is necessary to synchronize the classes that dump
// intermediate alignment results to files.  Otherwise, all data herein
// is constant and shared, or per-thread.
#include "threading.h"
#include "aligner_result.h"
#include "aligner_cache.h"
#include "scoring.h"
#include "mem_ids.h"
#include "simple_func.h"

/**
 * A constraint to apply to an alignment zone, or to an overall
 * alignment.
 *
 * The constraint can put both caps and ceilings on the number and
 * types of edits allowed.
 */
struct Constraint {
	
	Constraint() { init(); }
	
	/**
	 * Initialize Constraint to be fully permissive.
	 */
	void init() {
		edits = mms = ins = dels = penalty = editsCeil = mmsCeil =
		insCeil = delsCeil = penaltyCeil = MAX_I;
		penFunc.reset();
		instantiated = false;
	}
	
	/**
	 * Return true iff penalities and constraints prevent us from
	 * adding any edits.
	 */
	bool mustMatch() {
		assert(instantiated);
		return (mms == 0 && edits == 0) ||
		        penalty == 0 ||
		       (mms == 0 && dels == 0 && ins == 0);
	}
	
	/**
	 * Return true iff a mismatch of the given quality is permitted.
	 */
	bool canMismatch(int q, const Scoring& cm) {
		assert(instantiated);
		return (mms > 0 || edits > 0) &&
		       penalty >= cm.mm(q);
	}

	/**
	 * Return true iff a mismatch of the given quality is permitted.
	 */
	bool canN(int q, const Scoring& cm) {
		assert(instantiated);
		return (mms > 0 || edits > 0) &&
		       penalty >= cm.n(q);
	}
	
	/**
	 * Return true iff a mismatch of *any* quality (even qual=1) is
	 * permitted.
	 */
	bool canMismatch() {
		assert(instantiated);
		return (mms > 0 || edits > 0) && penalty > 0;
	}

	/**
	 * Return true iff a mismatch of *any* quality (even qual=1) is
	 * permitted.
	 */
	bool canN() {
		assert(instantiated);
		return (mms > 0 || edits > 0);
	}
	
	/**
	 * Return true iff a deletion of the given extension (0=open, 1=1st
	 * extension, etc) is permitted.
	 */
	bool canDelete(int ex, const Scoring& cm) {
		assert(instantiated);
		return (dels > 0 && edits > 0) &&
		       penalty >= cm.del(ex);
	}

	/**
	 * Return true iff a deletion of any extension is permitted.
	 */
	bool canDelete() {
		assert(instantiated);
		return (dels > 0 || edits > 0) &&
		       penalty > 0;
	}
	
	/**
	 * Return true iff an insertion of the given extension (0=open,
	 * 1=1st extension, etc) is permitted.
	 */
	bool canInsert(int ex, const Scoring& cm) {
		assert(instantiated);
		return (ins > 0 || edits > 0) &&
		       penalty >= cm.ins(ex);
	}

	/**
	 * Return true iff an insertion of any extension is permitted.
	 */
	bool canInsert() {
		assert(instantiated);
		return (ins > 0 || edits > 0) &&
		       penalty > 0;
	}
	
	/**
	 * Return true iff a gap of any extension is permitted
	 */
	bool canGap() {
		assert(instantiated);
		return ((ins > 0 || dels > 0) || edits > 0) && penalty > 0;
	}
	
	/**
	 * Charge a mismatch of the given quality.
	 */
	void chargeMismatch(int q, const Scoring& cm) {
		assert(instantiated);
		if(mms == 0) { assert_gt(edits, 0); edits--; }
		else mms--;
		penalty -= cm.mm(q);
		assert_geq(mms, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	
	/**
	 * Charge an N mismatch of the given quality.
	 */
	void chargeN(int q, const Scoring& cm) {
		assert(instantiated);
		if(mms == 0) { assert_gt(edits, 0); edits--; }
		else mms--;
		penalty -= cm.n(q);
		assert_geq(mms, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	
	/**
	 * Charge a deletion of the given extension.
	 */
	void chargeDelete(int ex, const Scoring& cm) {
		assert(instantiated);
		dels--;
		edits--;
		penalty -= cm.del(ex);
		assert_geq(dels, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}

	/**
	 * Charge an insertion of the given extension.
	 */
	void chargeInsert(int ex, const Scoring& cm) {
		assert(instantiated);
		ins--;
		edits--;
		penalty -= cm.ins(ex);
		assert_geq(ins, 0);
		assert_geq(edits, 0);
		assert_geq(penalty, 0);
	}
	
	/**
	 * Once the constrained area is completely explored, call this
	 * function to check whether there were *at least* as many
	 * dissimilarities as required by the constraint.  Bounds like this
	 * are helpful to resolve instances where two search roots would
	 * otherwise overlap in what alignments they can find.
	 */
	bool acceptable() {
		assert(instantiated);
		return edits   <= editsCeil &&
		       mms     <= mmsCeil   &&
		       ins     <= insCeil   &&
		       dels    <= delsCeil  &&
		       penalty <= penaltyCeil;
	}
	
	/**
	 * Instantiate a constraint w/r/t the read length and the constant
	 * and linear coefficients for the penalty function.
	 */
	static int instantiate(size_t rdlen, const SimpleFunc& func) {
		return func.f<int>((double)rdlen);
	}
	
	/**
	 * Instantiate this constraint w/r/t the read length.
	 */
	void instantiate(size_t rdlen) {
		assert(!instantiated);
		if(penFunc.initialized()) {
			penalty = Constraint::instantiate(rdlen, penFunc);
		}
		instantiated = true;
	}
	
	int edits;      // # edits permitted
	int mms;        // # mismatches permitted
	int ins;        // # insertions permitted
	int dels;       // # deletions permitted
	int penalty;    // penalty total permitted
	int editsCeil;  // <= this many edits can be left at the end
	int mmsCeil;    // <= this many mismatches can be left at the end
	int insCeil;    // <= this many inserts can be left at the end
	int delsCeil;   // <= this many deletions can be left at the end
	int penaltyCeil;// <= this much leftover penalty can be left at the end
	SimpleFunc penFunc;// penalty function; function of read len
	bool instantiated; // whether constraint is instantiated w/r/t read len
	
	//
	// Some static methods for constructing some standard Constraints
	//

	/**
	 * Construct a constraint with no edits of any kind allowed.
	 */
	static Constraint exact();
	
	/**
	 * Construct a constraint where the only constraint is a total
	 * penalty constraint.
	 */
	static Constraint penaltyBased(int pen);

	/**
	 * Construct a constraint where the only constraint is a total
	 * penalty constraint related to the length of the read.
	 */
	static Constraint penaltyFuncBased(const SimpleFunc& func);

	/**
	 * Construct a constraint where the only constraint is a total
	 * penalty constraint.
	 */
	static Constraint mmBased(int mms);

	/**
	 * Construct a constraint where the only constraint is a total
	 * penalty constraint.
	 */
	static Constraint editBased(int edits);
};

/**
 * We divide seed search strategies into three categories:
 *
 * 1. A left-to-right search where the left half of the read is
 *    constrained to match exactly and the right half is subject to
 *    some looser constraint (e.g. 1mm or 2mm)
 * 2. Same as 1, but going right to left with the exact matching half
 *    on the right.
 * 3. Inside-out search where the center half of the read is
 *    constrained to match exactly, and the extreme quarters of the
 *    read are subject to a looser constraint.
 */
enum {
	SEED_TYPE_EXACT = 1,
	SEED_TYPE_LEFT_TO_RIGHT,
	SEED_TYPE_RIGHT_TO_LEFT,
	SEED_TYPE_INSIDE_OUT
};

struct InstantiatedSeed;

/**
 * Policy dictating how to size and arrange seeds along the length of
 * the read, and what constraints to force on the zones of the seed.
 * We assume that seeds are plopped down at regular intervals from the
 * 5' to 3' ends, with the first seed flush to the 5' end.
 *
 * If the read is shorter than a single seed, one seed is used and it
 * is shrunk to accommodate the read.
 */
struct Seed {

	int len;             // length of a seed
	int type;            // dictates anchor portion, direction of search
	Constraint *overall; // for the overall alignment

	Seed() { init(0, 0, NULL); }

	/**
	 * Construct and initialize this seed with given length and type.
	 */
	Seed(int ln, int ty, Constraint* oc) {
		init(ln, ty, oc);
	}

	/**
	 * Initialize this seed with given length and type.
	 */
	void init(int ln, int ty, Constraint* oc) {
		len = ln;
		type = ty;
		overall = oc;
	}
	
	// If the seed is split into halves, we just use zones[0] and
	// zones[1]; 0 is the near half and 1 is the far half.  If the seed
	// is split into thirds (i.e. inside-out) then 0 is the center, 1
	// is the far portion on the left, and 2 is the far portion on the
	// right.
	Constraint zones[3];

	/**
	 * Once the constrained seed is completely explored, call this
	 * function to check whether there were *at least* as many
	 * dissimilarities as required by all constraints.  Bounds like this
	 * are helpful to resolve instances where two search roots would
	 * otherwise overlap in what alignments they can find.
	 */
	bool acceptable() {
		assert(overall != NULL);
		return zones[0].acceptable() &&
		       zones[1].acceptable() &&
		       zones[2].acceptable() &&
		       overall->acceptable();
	}

	/**
	 * Given a read, depth and orientation, extract a seed data structure
	 * from the read and fill in the steps & zones arrays.  The Seed
	 * contains the sequence and quality values.
	 */
	bool instantiate(
		const Read& read,
		const BTDnaString& seq, // already-extracted seed sequence
		const BTString& qual,   // already-extracted seed quality sequence
		const Scoring& pens,
		int depth,
		int seedoffidx,
		int seedtypeidx,
		bool fw,
		InstantiatedSeed& si) const;

	/**
	 * Return a list of Seed objects encapsulating
	 */
	static void mmSeeds(
		int mms,
		int ln,
		EList<Seed>& pols,
		Constraint& oall)
	{
		if(mms == 0) {
			zeroMmSeeds(ln, pols, oall);
		} else if(mms == 1) {
			oneMmSeeds(ln, pols, oall);
		} else if(mms == 2) {
			twoMmSeeds(ln, pols, oall);
		} else throw 1;
	}
	
	static void zeroMmSeeds(int ln, EList<Seed>&, Constraint&);
	static void oneMmSeeds (int ln, EList<Seed>&, Constraint&);
	static void twoMmSeeds (int ln, EList<Seed>&, Constraint&);
};

/**
 * An instantiated seed is a seed (perhaps modified to fit the read)
 * plus all data needed to conduct a search of the seed.
 */
struct InstantiatedSeed {

	InstantiatedSeed() : steps(AL_CAT), zones(AL_CAT) { }

	// Steps map.  There are as many steps as there are positions in
	// the seed.  The map is a helpful abstraction because we sometimes
	// visit seed positions in an irregular order (e.g. inside-out
	// search).
	EList<int> steps;

	// Zones map.  For each step, records what constraint to charge an
	// edit to.  The first entry in each pair gives the constraint for
	// non-insert edits and the second entry in each pair gives the
	// constraint for insert edits.  If the value stored is negative,
	// this indicates that the zone is "closed out" after this
	// position, so zone acceptility should be checked.
	EList<pair<int, int> > zones;

	// Nucleotide sequence covering the seed, extracted from read
	BTDnaString *seq;
	
	// Quality sequence covering the seed, extracted from read
	BTString *qual;
	
	// Initial constraints governing zones 0, 1, 2.  We precalculate
	// the effect of Ns on these.
	Constraint cons[3];
	
	// Overall constraint, tailored to the read length.
	Constraint overall;
	
	// Maximum number of positions that the aligner may advance before
	// its first step.  This lets the aligner know whether it can use
	// the ftab or not.
	int maxjump;
	
	// Offset of seed from 5' end of read
	int seedoff;

	// Id for seed offset; ids are such that the smallest index is the
	// closest to the 5' end and consecutive ids are adjacent (i.e.
	// there are no intervening offsets with seeds)
	int seedoffidx;
	
	// Type of seed (left-to-right, etc)
	int seedtypeidx;
	
	// Seed comes from forward-oriented read?
	bool fw;
	
	// Filtered out due to the pattern of Ns present.  If true, this
	// seed should be ignored by searchAllSeeds().
	bool nfiltered;
	
	// Seed this was instantiated from
	Seed s;
	
#ifndef NDEBUG
	/**
	 * Check that InstantiatedSeed is internally consistent.
	 */
	bool repOk() const {
		return true;
	}
#endif
};

/**
 * Simple struct for holding a end-to-end alignments for the read with at most
 * 2 edits.
 */
struct EEHit {
	
	EEHit() { reset(); }
	
	void reset() {
		top = bot = 0;
		fw = false;
		e1.reset();
		e2.reset();
		score = MIN_I64;
	}
	
	void init(
		uint32_t top_,
		uint32_t bot_,
		const Edit* e1_,
		const Edit* e2_,
		bool fw_,
		int64_t score_)
	{
		top = top_; bot = bot_;
		if(e1_ != NULL) {
			e1 = *e1_;
		} else {
			e1.reset();
		}
		if(e2_ != NULL) {
			e2 = *e2_;
		} else {
			e2.reset();
		}
		fw = fw_;
		score = score_;
	}
	
	/**
	 * Return number of mismatches in the alignment.
	 */
	int mms() const {
		if     (e2.inited()) return 2;
		else if(e1.inited()) return 1;
		else                 return 0;
	}
	
	/**
	 * Return the number of Ns involved in the alignment.
	 */
	int ns() const {
		int ns = 0;
		if(e1.inited() && e1.hasN()) {
			ns++;
			if(e2.inited() && e2.hasN()) {
				ns++;
			}
		}
		return ns;
	}

	/**
	 * Return the number of Ns involved in the alignment.
	 */
	int refns() const {
		int ns = 0;
		if(e1.inited() && e1.chr == 'N') {
			ns++;
			if(e2.inited() && e2.chr == 'N') {
				ns++;
			}
		}
		return ns;
	}
	
	/**
	 * Return true iff there is no hit.
	 */
	bool empty() const {
		return bot <= top;
	}
	
	/**
	 * Higher score = higher priority.
	 */
	bool operator<(const EEHit& o) const {
		return score > o.score;
	}
	
	/**
	 * Return the size of the alignments SA range.s
	 */
	uint32_t size() const { return bot - top; }
	
#ifndef NDEBUG
	/**
	 * Check that hit is sane w/r/t read.
	 */
	bool repOk(const Read& rd) const {
		assert_gt(bot, top);
		if(e1.inited()) {
			assert_lt(e1.pos, rd.length());
			if(e2.inited()) {
				assert_lt(e2.pos, rd.length());
			}
		}
		return true;
	}
#endif
	
	uint32_t top;
	uint32_t bot;
	Edit     e1;
	Edit     e2;
	bool     fw;
	int64_t  score;
};

/**
 * Data structure for holding all of the seed hits associated with a read.  All
 * the seed hits for a given read are encapsulated in a single QVal object.  A
 * QVal refers to a range of values in the qlist, where each qlist value is a 
 * BW range and a slot to hold the hit's suffix array offset.  QVals are kept
 * in two lists (hitsFw_ and hitsRc_), one for seeds on the forward read strand,
 * one for seeds on the reverse read strand.  The list is indexed by read
 * offset index (e.g. 0=closest-to-5', 1=second-closest, etc).
 *
 * An assumption behind this data structure is that all the seeds are found
 * first, then downstream analyses try to extend them.  In between finding the
 * seed hits and extending them, the sort() member function is called, which
 * ranks QVals according to the order they should be extended.  Right now the
 * policy is that QVals with fewer elements (hits) should be tried first.
 */
class SeedResults {

public:
	SeedResults() :
		seqFw_(AL_CAT),
		seqRc_(AL_CAT),
		qualFw_(AL_CAT),
		qualRc_(AL_CAT),
		hitsFw_(AL_CAT),
		hitsRc_(AL_CAT),
		isFw_(AL_CAT),
		isRc_(AL_CAT),
		sortedFw_(AL_CAT),
		sortedRc_(AL_CAT),
		offIdx2off_(AL_CAT),
		rankOffs_(AL_CAT),
		rankFws_(AL_CAT),
		mm1Hit_(AL_CAT)
	{
		clear();
	}
	
	/**
	 * Set the current read.
	 */
	void nextRead(const Read& read) {
		read_ = &read;
	}

	/**
	 * Set the appropriate element of either hitsFw_ or hitsRc_ to the given
	 * QVal.  A QVal encapsulates all the BW ranges for reference substrings 
	 * that are within some distance of the seed string.
	 */
	void add(
		const QVal& qv,           // range of ranges in cache
		const AlignmentCache& ac, // cache
		uint32_t seedIdx,         // seed index (from 5' end)
		bool     seedFw)          // whether seed is from forward read
	{
		assert(qv.repOk(ac));
		assert(repOk(&ac));
		assert_lt(seedIdx, hitsFw_.size());
		assert_gt(numOffs_, 0); // if this fails, probably failed to call reset
		if(qv.empty()) return;
		if(seedFw) {
			assert(!hitsFw_[seedIdx].valid());
			hitsFw_[seedIdx] = qv;
			numEltsFw_ += qv.numElts();
			numRangesFw_ += qv.numRanges();
			if(qv.numRanges() > 0) nonzFw_++;
		} else {
			assert(!hitsRc_[seedIdx].valid());
			hitsRc_[seedIdx] = qv;
			numEltsRc_ += qv.numElts();
			numRangesRc_ += qv.numRanges();
			if(qv.numRanges() > 0) nonzRc_++;
		}
		numElts_ += qv.numElts();
		numRanges_ += qv.numRanges();
		if(qv.numRanges() > 0) {
			nonzTot_++;
		}
		assert(repOk(&ac));
	}

	/**
	 * Clear buffered seed hits and state.  Set the number of seed
	 * offsets and the read.
	 */
	void reset(
		const Read& read,
		const EList<uint32_t>& offIdx2off,
		size_t numOffs)
	{
		assert_gt(numOffs, 0);
		clearSeeds();
		numOffs_ = numOffs;
		seqFw_.resize(numOffs_);
		seqRc_.resize(numOffs_);
		qualFw_.resize(numOffs_);
		qualRc_.resize(numOffs_);
		hitsFw_.resize(numOffs_);
		hitsRc_.resize(numOffs_);
		isFw_.resize(numOffs_);
		isRc_.resize(numOffs_);
		sortedFw_.resize(numOffs_);
		sortedRc_.resize(numOffs_);
		offIdx2off_ = offIdx2off;
		for(size_t i = 0; i < numOffs_; i++) {
			sortedFw_[i] = sortedRc_[i] = false;
			hitsFw_[i].reset();
			hitsRc_[i].reset();
			isFw_[i].clear();
			isRc_[i].clear();
		}
		read_ = &read;
		sorted_ = false;
	}
	
	/**
	 * Clear seed-hit state.
	 */
	void clearSeeds() {
		sortedFw_.clear();
		sortedRc_.clear();
		rankOffs_.clear();
		rankFws_.clear();
		offIdx2off_.clear();
		hitsFw_.clear();
		hitsRc_.clear();
		isFw_.clear();
		isRc_.clear();
		seqFw_.clear();
		seqRc_.clear();
		nonzTot_ = 0;
		nonzFw_ = 0;
		nonzRc_ = 0;
		numOffs_ = 0;
		numRanges_ = 0;
		numElts_ = 0;
		numRangesFw_ = 0;
		numEltsFw_ = 0;
		numRangesRc_ = 0;
		numEltsRc_ = 0;
	}
	
	/**
	 * Clear seed-hit state and end-to-end alignment state.
	 */
	void clear() {
		clearSeeds();
		read_ = NULL;
		exactFwHit_.reset();
		exactRcHit_.reset();
		mm1Hit_.clear();
		mm1Sorted_ = false;
		mm1Elt_ = 0;
		assert(empty());
	}
	
	/**
	 * Extract key summaries from this SeedResults and put into 'ssum'.
	 */
	void toSeedAlSumm(SeedAlSumm& ssum) const {
		// Number of positions with at least 1 range
		ssum.nonzTot   = nonzTot_;
		ssum.nonzFw    = nonzFw_;
		ssum.nonzRc    = nonzRc_;

		// Number of ranges
		ssum.nrangeTot = numRanges_;
		ssum.nrangeFw  = numRangesFw_;
		ssum.nrangeRc  = numRangesRc_;

		// Number of elements
		ssum.neltTot   = numElts_;
		ssum.neltFw    = numEltsFw_;
		ssum.neltRc    = numEltsRc_;
		
		// Other summaries
		ssum.maxNonzRangeFw = ssum.minNonzRangeFw = 0;
		ssum.maxNonzRangeRc = ssum.minNonzRangeRc = 0;
		ssum.maxNonzEltFw = ssum.minNonzEltFw = 0;
		ssum.maxNonzEltRc = ssum.minNonzEltRc = 0;
		for(size_t i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid()) {
				if(ssum.minNonzEltFw == 0 || hitsFw_[i].numElts() < ssum.minNonzEltFw) {
					ssum.minNonzEltFw = hitsFw_[i].numElts();
				}
				if(ssum.maxNonzEltFw == 0 || hitsFw_[i].numElts() > ssum.maxNonzEltFw) {
					ssum.maxNonzEltFw = hitsFw_[i].numElts();
				}
				if(ssum.minNonzRangeFw == 0 || hitsFw_[i].numRanges() < ssum.minNonzRangeFw) {
					ssum.minNonzRangeFw = hitsFw_[i].numRanges();
				}
				if(ssum.maxNonzRangeFw == 0 || hitsFw_[i].numRanges() > ssum.maxNonzRangeFw) {
					ssum.maxNonzRangeFw = hitsFw_[i].numRanges();
				}
			}
			if(hitsRc_[i].valid()) {
				if(ssum.minNonzEltRc == 0 || hitsRc_[i].numElts() < ssum.minNonzEltRc) {
					ssum.minNonzEltRc = hitsRc_[i].numElts();
				}
				if(ssum.maxNonzEltRc == 0 || hitsRc_[i].numElts() > ssum.maxNonzEltRc) {
					ssum.maxNonzEltRc = hitsRc_[i].numElts();
				}
				if(ssum.minNonzRangeRc == 0 || hitsRc_[i].numRanges() < ssum.minNonzRangeRc) {
					ssum.minNonzRangeRc = hitsRc_[i].numRanges();
				}
				if(ssum.maxNonzRangeRc == 0 || hitsRc_[i].numRanges() > ssum.maxNonzRangeRc) {
					ssum.maxNonzRangeRc = hitsRc_[i].numRanges();
				}
			}
		}
	}
	
	/**
	 * Return average number of hits per seed.
	 */
	float averageHitsPerSeed() const {
		return (float)numElts_ / (float)nonzTot_;
	}
	
	/**
	 * Return median of all the non-zero per-seed # hits
	 */
	float medianHitsPerSeed() const {
		EList<size_t>& median = const_cast<EList<size_t>&>(tmpMedian_);
		median.clear();
		for(size_t i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid() && hitsFw_[i].numElts() > 0) {
				median.push_back(hitsFw_[i].numElts());
			}
			if(hitsRc_[i].valid() && hitsRc_[i].numElts() > 0) {
				median.push_back(hitsRc_[i].numElts());
			}
		}
		if(tmpMedian_.empty()) {
			return 0.0f;
		}
		median.sort();
		float med1 = (float)median[tmpMedian_.size() >> 1];
		float med2 = med1;
		if((median.size() & 1) == 0) {
			med2 = (float)median[(tmpMedian_.size() >> 1) - 1];
		}
		return med1 + med2 * 0.5f;
	}
	
	/**
	 * Return a number that's meant to quantify how hopeful we are that this
	 * set of seed hits will lead to good alignments.
	 */
	double uniquenessFactor() const {
		double result = 0.0;
		for(size_t i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid()) {
				size_t nelt = hitsFw_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
			if(hitsRc_[i].valid()) {
				size_t nelt = hitsRc_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
		}
		return result;
	}

	/**
	 * Return the number of ranges being held.
	 */
	size_t numRanges() const { return numRanges_; }

	/**
	 * Return the number of elements being held.
	 */
	size_t numElts() const { return numElts_; }

	/**
	 * Return the number of ranges being held for seeds on the forward
	 * read strand.
	 */
	size_t numRangesFw() const { return numRangesFw_; }

	/**
	 * Return the number of elements being held for seeds on the
	 * forward read strand.
	 */
	size_t numEltsFw() const { return numEltsFw_; }

	/**
	 * Return the number of ranges being held for seeds on the
	 * reverse-complement read strand.
	 */
	size_t numRangesRc() const { return numRangesRc_; }

	/**
	 * Return the number of elements being held for seeds on the
	 * reverse-complement read strand.
	 */
	size_t numEltsRc() const { return numEltsRc_; }
	
	/**
	 * Given an offset index, return the offset that has that index.
	 */
	size_t idx2off(size_t off) const {
		return offIdx2off_[off];
	}
	
	/**
	 * Return true iff there are 0 hits being held.
	 */
	bool empty() const { return numRanges() == 0; }
	
	/**
	 * Get the QVal representing all the reference hits for the given
	 * orientation and seed offset index.
	 */
	const QVal& hitsAtOffIdx(bool fw, size_t seedoffidx) const {
		assert_lt(seedoffidx, numOffs_);
		assert(repOk(NULL));
		return fw ? hitsFw_[seedoffidx] : hitsRc_[seedoffidx];
	}

	/**
	 * Get the Instantiated seeds for the given orientation and offset.
	 */
	EList<InstantiatedSeed>& instantiatedSeeds(bool fw, size_t seedoffidx) {
		assert_lt(seedoffidx, numOffs_);
		assert(repOk(NULL));
		return fw ? isFw_[seedoffidx] : isRc_[seedoffidx];
	}
	
	/**
	 * Return the number of different seed offsets possible.
	 */
	size_t numOffs() const { return numOffs_; }
	
	/**
	 * Return the read from which seeds were extracted, aligned.
	 */
	const Read& read() const { return *read_; }
	
#ifndef NDEBUG
	/**
	 * Check that this SeedResults is internally consistent.
	 */
	bool repOk(
		const AlignmentCache* ac,
		bool requireInited = false) const
	{
		if(requireInited) {
			assert(read_ != NULL);
		}
		if(numOffs_ > 0) {
			assert_eq(numOffs_, hitsFw_.size());
			assert_eq(numOffs_, hitsRc_.size());
			assert_leq(numRanges_, numElts_);
			assert_leq(nonzTot_, numRanges_);
			size_t nonzs = 0;
			for(int fw = 0; fw <= 1; fw++) {
				const EList<QVal>& rrs = (fw ? hitsFw_ : hitsRc_);
				for(size_t i = 0; i < numOffs_; i++) {
					if(rrs[i].valid()) {
						if(rrs[i].numRanges() > 0) nonzs++;
						if(ac != NULL) {
							assert(rrs[i].repOk(*ac));
						}
					}
				}
			}
			assert_eq(nonzs, nonzTot_);
			assert(!sorted_ || nonzTot_ == rankFws_.size());
			assert(!sorted_ || nonzTot_ == rankOffs_.size());
		}
		return true;
	}
#endif
	
	/**
	 * Populate rankOffs_ and rankFws_ with the list of QVals that need to be
	 * examined for this SeedResults, in order.  The order is ascending by
	 * number of elements, so QVals with fewer elements (i.e. seed sequences
	 * that are more unique) will be tried first and QVals with more elements
	 * (i.e. seed sequences
	 */
	void rankSeedHits(RandomSource& rnd) {
		while(rankOffs_.size() < nonzTot_) {
			uint32_t minsz = 0xffffffff;
			uint32_t minidx = 0;
			bool minfw = true;
			// Rank seed-hit positions in ascending order by number of elements
			// in all BW ranges
			bool rb = rnd.nextBool();
			assert(rb == 0 || rb == 1);
			for(int fwi = 0; fwi <= 1; fwi++) {
				bool fw = (fwi == (rb ? 1 : 0));
				EList<QVal>& rrs = (fw ? hitsFw_ : hitsRc_);
				EList<bool>& sorted = (fw ? sortedFw_ : sortedRc_);
				uint32_t i = (rnd.nextU32() % (uint32_t)numOffs_);
				for(uint32_t ii = 0; ii < numOffs_; ii++) {
					if(rrs[i].valid() &&         // valid QVal
					   rrs[i].numElts() > 0 &&   // non-empty
					   !sorted[i] &&             // not already sorted
					   rrs[i].numElts() < minsz) // least elts so far?
					{
						minsz = rrs[i].numElts();
						minidx = i;
						minfw = (fw == 1);
					}
					if((++i) == numOffs_) {
						i = 0;
					}
				}
			}
			assert_neq(0xffffffff, minsz);
			if(minfw) {
				sortedFw_[minidx] = true;
			} else {
				sortedRc_[minidx] = true;
			}
			rankOffs_.push_back(minidx);
			rankFws_.push_back(minfw);
		}
		assert_eq(rankOffs_.size(), rankFws_.size());
		sorted_ = true;
	}

	/**
	 * Return the number of orientation/offsets into the read that have
	 * at least one seed hit.
	 */
	size_t nonzeroOffsets() const {
		assert(!sorted_ || nonzTot_ == rankFws_.size());
		assert(!sorted_ || nonzTot_ == rankOffs_.size());
		return nonzTot_;
	}
	
	/**
	 * Return true iff all seeds hit for forward read.
	 */
	bool allFwSeedsHit() const {
		return nonzFw_ == numOffs();
	}

	/**
	 * Return true iff all seeds hit for revcomp read.
	 */
	bool allRcSeedsHit() const {
		return nonzRc_ == numOffs();
	}
	
	/**
	 * Return the minimum number of edits that an end-to-end alignment of the
	 * fw read could have.  Uses knowledge of how many seeds have exact hits
	 * and how the seeds overlap.
	 */
	size_t fewestEditsEE(bool fw, int seedlen, int per) const {
		assert_gt(seedlen, 0);
		assert_gt(per, 0);
		size_t nonz = fw ? nonzFw_ : nonzRc_;
		if(nonz < numOffs()) {
			int maxdepth = (seedlen + per - 1) / per;
			int missing = (int)(numOffs() - nonz);
			return (missing + maxdepth - 1) / maxdepth;
		} else {
			// Exact hit is possible (not guaranteed)
			return 0;
		}
	}

	/**
	 * Return the number of offsets into the forward read that have at
	 * least one seed hit.
	 */
	size_t nonzeroOffsetsFw() const {
		return nonzFw_;
	}
	
	/**
	 * Return the number of offsets into the reverse-complement read
	 * that have at least one seed hit.
	 */
	size_t nonzeroOffsetsRc() const {
		return nonzRc_;
	}

	/**
	 * Return a QVal of seed hits of the given rank 'r'.  'offidx' gets the id
	 * of the offset from 5' from which it was extracted (0 for the 5-most
	 * offset, 1 for the next closes to 5', etc).  'off' gets the offset from
	 * the 5' end.  'fw' gets true iff the seed was extracted from the forward
	 * read.
	 */
	const QVal& hitsByRank(
		size_t    r,       // in
		uint32_t& offidx,  // out
		uint32_t& off,     // out
		bool&     fw,      // out
		uint32_t& seedlen) // out
	{
		assert(sorted_);
		assert_lt(r, nonzTot_);
		if(rankFws_[r]) {
			fw = true;
			offidx = rankOffs_[r];
			assert_lt(offidx, offIdx2off_.size());
			off = offIdx2off_[offidx];
			seedlen = (uint32_t)seqFw_[rankOffs_[r]].length();
			return hitsFw_[rankOffs_[r]];
		} else {
			fw = false;
			offidx = rankOffs_[r];
			assert_lt(offidx, offIdx2off_.size());
			off = offIdx2off_[offidx];
			seedlen = (uint32_t)seqRc_[rankOffs_[r]].length();
			return hitsRc_[rankOffs_[r]];
		}
	}

	/**
	 * Return an EList of seed hits of the given rank.
	 */
	const BTDnaString& seqByRank(size_t r) {
		assert(sorted_);
		assert_lt(r, nonzTot_);
		return rankFws_[r] ? seqFw_[rankOffs_[r]] : seqRc_[rankOffs_[r]];
	}

	/**
	 * Return an EList of seed hits of the given rank.
	 */
	const BTString& qualByRank(size_t r) {
		assert(sorted_);
		assert_lt(r, nonzTot_);
		return rankFws_[r] ? qualFw_[rankOffs_[r]] : qualRc_[rankOffs_[r]];
	}
	
	/**
	 * Return the list of extracted seed sequences for seeds on either
	 * the forward or reverse strand.
	 */
	EList<BTDnaString>& seqs(bool fw) { return fw ? seqFw_ : seqRc_; }

	/**
	 * Return the list of extracted quality sequences for seeds on
	 * either the forward or reverse strand.
	 */
	EList<BTString>& quals(bool fw) { return fw ? qualFw_ : qualRc_; }

	/**
	 * Return exact end-to-end alignment of fw read.
	 */
	EEHit exactFwEEHit() const { return exactFwHit_; }

	/**
	 * Return exact end-to-end alignment of rc read.
	 */
	EEHit exactRcEEHit() const { return exactRcHit_; }
	
	/**
	 * Return const ref to list of 1-mismatch end-to-end alignments.
	 */
	const EList<EEHit>& mm1EEHits() const { return mm1Hit_; }
	
	/**
	 * Sort the end-to-end 1-mismatch alignments, prioritizing by score (higher
	 * score = higher priority).
	 */
	void sort1mmEe(RandomSource& rnd) {
		assert(!mm1Sorted_);
		mm1Hit_.sort();
		size_t streak = 0;
		for(size_t i = 1; i < mm1Hit_.size(); i++) {
			if(mm1Hit_[i].score == mm1Hit_[i-1].score) {
				if(streak == 0) { streak = 1; }
				streak++;
			} else {
				if(streak > 1) {
					assert_geq(i, streak);
					mm1Hit_.shufflePortion(i-streak, streak, rnd);
				}
				streak = 0;
			}
		}
		if(streak > 1) {
			mm1Hit_.shufflePortion(mm1Hit_.size() - streak, streak, rnd);
		}
		mm1Sorted_ = true;
	}
	
	/**
	 * Add an end-to-end 1-mismatch alignment.
	 */
	void add1mmEe(
		uint32_t top,
		uint32_t bot,
		const Edit* e1,
		const Edit* e2,
		bool fw,
		int64_t score)
	{
		mm1Hit_.expand();
		mm1Hit_.back().init(top, bot, e1, e2, fw, score);
		mm1Elt_ += (bot - top);
	}

	/**
	 * Add an end-to-end exact alignment.
	 */
	void addExactEeFw(
		uint32_t top,
		uint32_t bot,
		const Edit* e1,
		const Edit* e2,
		bool fw,
		int64_t score)
	{
		exactFwHit_.init(top, bot, e1, e2, fw, score);
	}

	/**
	 * Add an end-to-end exact alignment.
	 */
	void addExactEeRc(
		uint32_t top,
		uint32_t bot,
		const Edit* e1,
		const Edit* e2,
		bool fw,
		int64_t score)
	{
		exactRcHit_.init(top, bot, e1, e2, fw, score);
	}
	
	/**
	 * Clear out the end-to-end exact alignments.
	 */
	void clearExactE2eHits() {
		exactFwHit_.reset();
		exactRcHit_.reset();
	}
	
	/**
	 * Clear out the end-to-end 1-mismatch alignments.
	 */
	void clear1mmE2eHits() {
		mm1Hit_.clear();     // 1-mismatch end-to-end hits
		mm1Elt_ = 0;         // number of 1-mismatch hit rows
		mm1Sorted_ = false;  // true iff we've sorted the mm1Hit_ list
	}

	/**
	 * Return the number of distinct exact and 1-mismatch end-to-end hits
	 * found.
	 */
	size_t numE2eHits() const {
		return exactFwHit_.size() + exactRcHit_.size() + mm1Elt_;
	}

	/**
	 * Return the number of distinct exact end-to-end hits found.
	 */
	size_t numExactE2eHits() const {
		return exactFwHit_.size() + exactRcHit_.size();
	}

	/**
	 * Return the number of distinct 1-mismatch end-to-end hits found.
	 */
	size_t num1mmE2eHits() const {
		return mm1Elt_;
	}
	
	/**
	 * Return the length of the read that yielded the seed hits.
	 */
	size_t readLength() const {
		assert(read_ != NULL);
		return read_->length();
	}

protected:

	// As seed hits and edits are added they're sorted into these
	// containers
	EList<BTDnaString>  seqFw_;       // seqs for seeds from forward read
	EList<BTDnaString>  seqRc_;       // seqs for seeds from revcomp read
	EList<BTString>     qualFw_;      // quals for seeds from forward read
	EList<BTString>     qualRc_;      // quals for seeds from revcomp read
	EList<QVal>         hitsFw_;      // hits for forward read
	EList<QVal>         hitsRc_;      // hits for revcomp read
	EList<EList<InstantiatedSeed> > isFw_; // hits for forward read
	EList<EList<InstantiatedSeed> > isRc_; // hits for revcomp read
	EList<bool>         sortedFw_;    // true iff fw QVal was sorted/ranked
	EList<bool>         sortedRc_;    // true iff rc QVal was sorted/ranked
	size_t              nonzTot_;     // # offsets with non-zero size
	size_t              nonzFw_;      // # offsets into fw read with non-0 size
	size_t              nonzRc_;      // # offsets into rc read with non-0 size
	size_t              numRanges_;   // # ranges added
	size_t              numElts_;     // # elements added
	size_t              numRangesFw_; // # ranges added for fw seeds
	size_t              numEltsFw_;   // # elements added for fw seeds
	size_t              numRangesRc_; // # ranges added for rc seeds
	size_t              numEltsRc_;   // # elements added for rc seeds

	EList<uint32_t>     offIdx2off_;// map from offset indexes to offsets from 5' end

	// When the sort routine is called, the seed hits collected so far
	// are sorted into another set of containers that allow easy access
	// to hits from the lowest-ranked offset (the one with the fewest
	// BW elements) to the greatest-ranked offset.  Offsets with 0 hits
	// are ignored.
	EList<uint32_t>     rankOffs_;  // sorted offests of seeds to try
	EList<bool>         rankFws_;   // sorted orientations assoc. with rankOffs_
	bool                sorted_;    // true if sort() called since last reset
	
	// These fields set once per read
	size_t              numOffs_;   // # different seed offsets possible
	const Read*         read_;      // read from which seeds were extracted
	
	EEHit               exactFwHit_; // end-to-end exact hit for fw read
	EEHit               exactRcHit_; // end-to-end exact hit for rc read
	EList<EEHit>        mm1Hit_;     // 1-mismatch end-to-end hits
	size_t              mm1Elt_;     // number of 1-mismatch hit rows
	bool                mm1Sorted_;  // true iff we've sorted the mm1Hit_ list

	EList<size_t> tmpMedian_; // temporary storage for calculating median
};


// Forward decl
class Ebwt;
struct SideLocus;

/**
 * Encapsulates a sumamry of what the searchAllSeeds aligner did.
 */
struct SeedSearchMetrics {

	SeedSearchMetrics() : mutex_m() {
	    reset();
	}

	/**
	 * Merge this metrics object with the given object, i.e., sum each
	 * category.  This is the only safe way to update a
	 * SeedSearchMetrics object shread by multiple threads.
	 */
	void merge(const SeedSearchMetrics& m, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);
		seedsearch   += m.seedsearch;
		possearch    += m.possearch;
		intrahit     += m.intrahit;
		interhit     += m.interhit;
		filteredseed += m.filteredseed;
		ooms         += m.ooms;
		bwops        += m.bwops;
		bweds        += m.bweds;
		bestmin0     += m.bestmin0;
		bestmin1     += m.bestmin1;
		bestmin2     += m.bestmin2;
	}
	
	/**
	 * Set all counters to 0.
	 */
	void reset() {
		seedsearch =
		possearch =
		intrahit =
		interhit =
		filteredseed =
		ooms =
		bwops =
		bweds =
		bestmin0 =
		bestmin1 =
		bestmin2 = 0;
	}

	uint64_t seedsearch;   // # times we executed strategy in InstantiatedSeed
	uint64_t possearch;    // # offsets where aligner executed >= 1 strategy
	uint64_t intrahit;     // # offsets where current-read cache gave answer
	uint64_t interhit;     // # offsets where across-read cache gave answer
	uint64_t filteredseed; // # seed instantiations skipped due to Ns
	uint64_t ooms;         // out-of-memory errors
	uint64_t bwops;        // Burrows-Wheeler operations
	uint64_t bweds;        // Burrows-Wheeler edits
	uint64_t bestmin0;     // # times the best min # edits was 0
	uint64_t bestmin1;     // # times the best min # edits was 1
	uint64_t bestmin2;     // # times the best min # edits was 2
	MUTEX_T  mutex_m;
};

/**
 * Given an index and a seeding scheme, searches for seed hits.
 */
class SeedAligner {

public:
	
	/**
	 * Initialize with index.
	 */
	SeedAligner() : edits_(AL_CAT), offIdx2off_(AL_CAT) { }

	/**
	 * Given a read and a few coordinates that describe a substring of the
	 * read (or its reverse complement), fill in 'seq' and 'qual' objects
	 * with the seed sequence and qualities.
	 */
	void instantiateSeq(
		const Read& read, // input read
		BTDnaString& seq, // output sequence
		BTString& qual,   // output qualities
		int len,          // seed length
		int depth,        // seed's 0-based offset from 5' end
		bool fw) const;   // seed's orientation

	/**
	 * Iterate through the seeds that cover the read and initiate a
	 * search for each seed.
	 */
	std::pair<int, int> instantiateSeeds(
		const EList<Seed>& seeds,   // search seeds
		size_t off,                 // offset into read to start extracting
		int per,                    // interval between seeds
		const Read& read,           // read to align
		const Scoring& pens,        // scoring scheme
		bool nofw,                  // don't align forward read
		bool norc,                  // don't align revcomp read
		AlignmentCacheIface& cache, // holds some seed hits from previous reads
		SeedResults& sr,            // holds all the seed hits
		SeedSearchMetrics& met);    // metrics

	/**
	 * Iterate through the seeds that cover the read and initiate a
	 * search for each seed.
	 */
	void searchAllSeeds(
		const EList<Seed>& seeds,   // search seeds
		const Ebwt* ebwtFw,         // BWT index
		const Ebwt* ebwtBw,         // BWT' index
		const Read& read,           // read to align
		const Scoring& pens,        // scoring scheme
		AlignmentCacheIface& cache, // local seed alignment cache
		SeedResults& hits,          // holds all the seed hits
		SeedSearchMetrics& met,     // metrics
		PerReadMetrics& prm);       // per-read metrics

	/**
	 * Sanity-check a partial alignment produced during oneMmSearch.
	 */
	bool sanityPartial(
		const Ebwt*        ebwtFw, // BWT index
		const Ebwt*        ebwtBw, // BWT' index
		const BTDnaString& seq,
		size_t             dep,
		size_t             len,
		bool               do1mm,
		uint32_t           topfw,
		uint32_t           botfw,
		uint32_t           topbw,
		uint32_t           botbw);

	/**
	 * Do an exact-matching sweet to establish a lower bound on number of edits
	 * and to find exact alignments.
	 */
	size_t exactSweep(
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
		SeedSearchMetrics& met);    // metrics

	/**
	 * Search for end-to-end alignments with up to 1 mismatch.
	 */
	bool oneMmSearch(
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
		SeedSearchMetrics& met);   // metrics

protected:

	/**
	 * Report a seed hit found by searchSeedBi(), but first try to extend it out in
	 * either direction as far as possible without hitting any edits.  This will
	 * allow us to prioritize the seed hits better later on.  Call reportHit() when
	 * we're done, which actually adds the hit to the cache.  Returns result from
	 * calling reportHit().
	 */
	bool extendAndReportHit(
		uint32_t topf,                     // top in BWT
		uint32_t botf,                     // bot in BWT
		uint32_t topb,                     // top in BWT'
		uint32_t botb,                     // bot in BWT'
		uint16_t len,                      // length of hit
		DoublyLinkedList<Edit> *prevEdit); // previous edit

	/**
	 * Report a seed hit found by searchSeedBi() by adding it to the cache.  Return
	 * false if the hit could not be reported because of, e.g., cache exhaustion.
	 */
	bool reportHit(
		uint32_t topf,         // top in BWT
		uint32_t botf,         // bot in BWT
		uint32_t topb,         // top in BWT'
		uint32_t botb,         // bot in BWT'
		uint16_t len,          // length of hit
		DoublyLinkedList<Edit> *prevEdit);  // previous edit
	
	/**
	 * Given an instantiated seed (in s_ and other fields), search
	 */
	bool searchSeedBi();
	
	/**
	 * Main, recursive implementation of the seed search.
	 */
	bool searchSeedBi(
		int step,              // depth into steps_[] array
		int depth,             // recursion depth
		uint32_t topf,         // top in BWT
		uint32_t botf,         // bot in BWT
		uint32_t topb,         // top in BWT'
		uint32_t botb,         // bot in BWT'
		SideLocus tloc,        // locus for top (perhaps unititialized)
		SideLocus bloc,        // locus for bot (perhaps unititialized)
		Constraint c0,         // constraints to enforce in seed zone 0
		Constraint c1,         // constraints to enforce in seed zone 1
		Constraint c2,         // constraints to enforce in seed zone 2
		Constraint overall,    // overall constraints
		DoublyLinkedList<Edit> *prevEdit);  // previous edit
	
	/**
	 * Get tloc and bloc ready for the next step.
	 */
	inline void nextLocsBi(
		SideLocus& tloc,            // top locus
		SideLocus& bloc,            // bot locus
		uint32_t topf,              // top in BWT
		uint32_t botf,              // bot in BWT
		uint32_t topb,              // top in BWT'
		uint32_t botb,              // bot in BWT'
		int step);                  // step to get ready for
	
	// Following are set in searchAllSeeds then used by searchSeed()
	// and other protected members.
	const Ebwt* ebwtFw_;       // forward index (BWT)
	const Ebwt* ebwtBw_;       // backward/mirror index (BWT')
	const Scoring* sc_;        // scoring scheme
	const InstantiatedSeed* s_;// current instantiated seed
	
	const Read* read_;         // read whose seeds are currently being aligned
	
	// The following are set just before a call to searchSeedBi()
	const BTDnaString* seq_;   // sequence of current seed
	const BTString* qual_;     // quality string for current seed
	size_t off_;               // offset of seed currently being searched
	bool fw_;                  // orientation of seed currently being searched
	
	EList<Edit> edits_;        // temporary place to sort edits
	AlignmentCacheIface *ca_;  // local alignment cache for seed alignments
	EList<uint32_t> offIdx2off_;// offset idx to read offset map, set up instantiateSeeds()
	uint64_t bwops_;           // Burrows-Wheeler operations
	uint64_t bwedits_;         // Burrows-Wheeler edits
	BTDnaString tmprfdnastr_;  // used in reportHit
	
	ASSERT_ONLY(ESet<BTDnaString> hits_); // Ref hits so far for seed being aligned
	BTDnaString tmpdnastr_;
};

#define INIT_LOCS(top, bot, tloc, bloc, e) { \
	if(bot - top == 1) { \
		tloc.initFromRow(top, (e).eh(), (e).ebwt()); \
		bloc.invalidate(); \
	} else { \
		SideLocus::initFromTopBot(top, bot, (e).eh(), (e).ebwt(), tloc, bloc); \
		assert(bloc.valid()); \
	} \
}

#define SANITY_CHECK_4TUP(t, b, tp, bp) { \
	ASSERT_ONLY(uint32_t tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3])); \
	ASSERT_ONLY(uint32_t totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3])); \
	assert_eq(tot, totp); \
}

#endif /*ALIGNER_SEED_H_*/
