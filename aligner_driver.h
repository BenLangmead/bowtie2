/*
 * Copyright 2012, Ben Langmead <langmea@cs.jhu.edu>
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
 * aligner_driver.h
 *
 * REDUNDANT SEED HITS
 *
 * We say that two seed hits are redundant if they trigger identical
 * seed-extend dynamic programming problems.  Put another way, they both lie on
 * the same diagonal of the overall read/reference dynamic programming matrix.
 * Detecting redundant seed hits is simple when the seed hits are ungapped.  We
 * do this after offset resolution but before the offset is converted to genome
 * coordinates (see uses of the seenDiags1_/seenDiags2_ fields for examples).
 *
 * REDUNDANT ALIGNMENTS
 *
 * In an unpaired context, we say that two alignments are redundant if they
 * share any cells in the global DP table.  Roughly speaking, this is like
 * saying that two alignments are redundant if any read character aligns to the
 * same reference character (same reference sequence, same strand, same offset)
 * in both alignments.
 *
 * In a paired-end context, we say that two paired-end alignments are redundant
 * if the mate #1s are redundant and the mate #2s are redundant.
 *
 * How do we enforce this?  In the unpaired context, this is relatively simple:
 * the cells from each alignment are checked against a set containing all cells
 * from all previous alignments.  Given a new alignment, for each cell in the
 * new alignment we check whether it is in the set.  If there is any overlap,
 * the new alignment is rejected as redundant.  Otherwise, the new alignment is
 * accepted and its cells are added to the set.
 *
 * Enforcement in a paired context is a little trickier.  Consider the
 * following approaches:
 *
 * 1. Skip anchors that are redundant with any previous anchor or opposite
 *    alignment.  This is sufficient to ensure no two concordant alignments
 *    found are redundant.
 *
 * 2. Same as scheme 1, but with a "transitive closure" scheme for finding all
 *    concordant pairs in the vicinity of an anchor.  Consider the AB/AC
 *    scenario from the previous paragraph.  If B is the anchor alignment, we
 *    will find AB but not AC.  But under this scheme, once we find AB we then
 *    let B be a new anchor and immediately look for its opposites.  Likewise,
 *    if we find any opposite, we make them anchors and continue searching.  We
 *    don't stop searching until every opposite is used as an anchor.
 *
 * 3. Skip anchors that are redundant with any previous anchor alignment (but
 *    allow anchors that are redundant with previous opposite alignments).
 *    This isn't sufficient to avoid redundant concordant alignments.  To avoid
 *    redundant concordants, we need an additional procedure that checks each
 *    new concordant alignment one-by-one against a list of previous concordant
 *    alignments to see if it is redundant.
 *
 * We take approach 1.
 */

#ifndef ALIGNER_DRIVER_H_
#define ALIGNER_DRIVER_H_

#include "aligner_seed2.h"
#include "simple_func.h"
#include "aln_sink.h"

/**
 * Concrete subclass of DescentRootSelector.  Puts a root every 'ival' chars,
 * where 'ival' is determined by user-specified parameters.  A root is filtered
 * out if the end of the read is less than 'landing' positions away, in the
 * direction of the search.
 */
class IntervalRootSelector : public DescentRootSelector {

public:

	IntervalRootSelector(
		double consExp,
		const SimpleFunc& rootIval,
		size_t landing)
	{
		consExp_ = consExp;
		rootIval_ = rootIval;
		landing_ = landing;
	}
	
	virtual ~IntervalRootSelector() { }

	virtual void select(
		const Read& q,                 // read that we're selecting roots for
		const Read* qo,                // opposite mate, if applicable
		bool nofw,                     // don't add roots for fw read
		bool norc,                     // don't add roots for rc read
		EList<DescentConfig>& confs,   // put DescentConfigs here
		EList<DescentRoot>& roots);    // put DescentRoot here

protected:

	double consExp_;
	SimpleFunc rootIval_;
	size_t landing_;
};

/**
 * Concrete subclass of DescentRootSelector.  Puts a root every 'ival' chars,
 * where 'ival' is determined by user-specified parameters.  A root is filtered
 * out if the end of the read is less than 'landing' positions away, in the
 * direction of the search.
 */
class PrioritizedRootSelector : public DescentRootSelector {

public:

	PrioritizedRootSelector(
		double consExp,
		const SimpleFunc& rootIval,
		size_t landing)
	{
		consExp_ = consExp;
		rootIval_ = rootIval;
		landing_ = landing;
	}
	
	virtual ~PrioritizedRootSelector() { }

	virtual void select(
		const Read& q,                 // read that we're selecting roots for
		const Read* qo,                // opposite mate, if applicable
		bool nofw,                     // don't add roots for fw read
		bool norc,                     // don't add roots for rc read
		EList<DescentConfig>& confs,   // put DescentConfigs here
		EList<DescentRoot>& roots);    // put DescentRoot here

protected:

	double consExp_;
	SimpleFunc rootIval_;
	size_t landing_;
	EHeap<DescentRoot> rootHeap_;
	EList<int> scoresOrig_[2];
	EList<int> scores_[2];
};

/**
 * Return values from extendSeeds and extendSeedsPaired.
 */
enum {
	// Candidates were examined exhaustively
	ALDRIVER_EXHAUSTED_CANDIDATES = 1,
	// The policy does not need us to look any further
	ALDRIVER_POLICY_FULFILLED,
	// We stopped because we ran up against a limit on how much work we should
	// do for one set of seed ranges, e.g. the limit on number of consecutive
	// unproductive DP extensions
	ALDRIVER_EXCEEDED_LIMIT
};

/**
 * This class is the glue between a DescentDriver and the dynamic programming
 * implementations in Bowtie 2.  The DescentDriver is used to find some very
 * high-scoring alignments, but is additionally used to rank partial alignments
 * so that they can be extended using dynamic programming.
 *
 * It is also the glue between the DescentDrivers and the DescentRootSelector
 * concrete subclasses that decide where to put the search roots.
 */
class AlignerDriver {

public:

	AlignerDriver(
		double consExp,
		bool prioritizeRoots,
		const SimpleFunc& rootIval,
		size_t landing,
		bool veryVerbose,
		const SimpleFunc& totsz,
		const SimpleFunc& totfmops) :
		alsel_(),
		dr1_(veryVerbose),
		dr2_(veryVerbose)
	{
		assert_gt(landing, 0);
		totsz_ = totsz;
		totfmops_ = totfmops;
		if(prioritizeRoots) {
			// Prioritize roots according the quality info & Ns
			sel_ = new PrioritizedRootSelector(consExp, rootIval, landing);
		} else {
			// Take a root every so many positions
			sel_ = new IntervalRootSelector(consExp, rootIval, landing);
		}
	}

	/**
	 * Destroy this AlignerDriver.
	 */
	virtual ~AlignerDriver() {
		delete sel_;
	}

	/**
	 * Initialize driver with respect to a new read or pair.
	 */
	void initRead(
		const Read& q1,
		bool nofw,
		bool norc,
		TAlScore minsc,
		TAlScore maxpen,
		const Read* q2)
	{
		// Initialize search for mate 1.  This includes instantiating and
		// prioritizing all the search roots.
		dr1_.initRead(q1, nofw, norc, minsc, maxpen, q2, sel_);
		red1_.init(q1.length());
		paired_ = false;
		if(q2 != NULL) {
			// Initialize search for mate 1.  This includes instantiating and
			// prioritizing all the search roots.
			dr2_.initRead(*q2, nofw, norc, minsc, maxpen, &q1, sel_);
			red2_.init(q2->length());
			paired_ = true;
		} else {
			dr2_.reset();
		}
		// Initialize stopping conditions.  We use two conditions:
		// totsz: when memory footprint exceeds this many bytes
		// totfmops: when we've exceeded this many FM Index ops
		size_t totsz = totsz_.f<size_t>(q1.length());
		size_t totfmops = totfmops_.f<size_t>(q1.length());
		stop_.init(
			totsz,
			0,
			true,
			totfmops);
	}
	
	/**
	 * Start the driver.  The driver will begin by conducting a best-first,
	 * index-assisted search through the space of possible full and partial
	 * alignments.  This search may be followed up with a dynamic programming
	 * extension step, taking a prioritized set of partial SA ranges found
	 * during the search and extending each with DP.  The process might also be
	 * iterated, with the search being occasioanally halted so that DPs can be
	 * tried, then restarted, etc.
	 */
	int go(
		const Scoring& sc,
		const Ebwt& ebwtFw,
		const Ebwt& ebwtBw,
		const BitPairReference& ref,
		DescentMetrics& met,
		WalkMetrics& wlm,
		PerReadMetrics& prm,
		RandomSource& rnd,
		AlnSinkWrap& sink);
	
	/**
	 * Reset state of all DescentDrivers.
	 */
	void reset() {
		dr1_.reset();
		dr2_.reset();
		red1_.reset();
		red2_.reset();
	}
	
	const DescentDriver& dr1() { return dr1_; }
	const DescentDriver& dr2() { return dr2_; }

protected:

	DescentRootSelector *sel_;        // selects where roots should go
	DescentAlignmentSelector alsel_;  // one selector can deal with >1 drivers
	DescentDriver dr1_;               // driver for mate 1/unpaired reads
	DescentDriver dr2_;               // driver for paired-end reads
	DescentStoppingConditions stop_;  // when to pause index-assisted BFS
	bool paired_;                     // current read is paired?

	SimpleFunc totsz_;      // memory limit on best-first search data
	SimpleFunc totfmops_;   // max # FM ops for best-first search

	// For detecting redundant alignments
	RedundantAlns  red1_;   // database of cells used for mate 1 alignments
	RedundantAlns  red2_;   // database of cells used for mate 2 alignments

	// For AlnRes::matchesRef
	ASSERT_ONLY(SStringExpandable<char> raw_refbuf_);
	ASSERT_ONLY(SStringExpandable<uint32_t> raw_destU32_);
	ASSERT_ONLY(EList<bool> raw_matches_);
	ASSERT_ONLY(BTDnaString tmp_rf_);
	ASSERT_ONLY(BTDnaString tmp_rdseq_);
	ASSERT_ONLY(BTString tmp_qseq_);
};

#endif /* defined(ALIGNER_DRIVER_H_) */
