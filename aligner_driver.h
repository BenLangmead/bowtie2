/*
 * Copyright 2012, Ben Langmead <blangmea@jhsph.edu>
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

typedef DescentQuery AlignerDriverQuery;

/**
 * Concrete subclass of DescentRootSelector.  Puts a root every 'ival' chars,
 * where 'ival' is determined by user-specified parameters.  A root is filtered
 * out if the end of the read is less than 'landing' positions away, in the
 * direction of the search.
 */
class AlignerDriverRootSelector : public DescentRootSelector {

public:

	AlignerDriverRootSelector(
		const SimpleFunc& descCons,
		const SimpleFunc& rootIval,
		size_t landing)
	{
		descCons_ = descCons;
		rootIval_ = rootIval;
		landing_ = landing;
	}

	virtual void select(
		const DescentQuery& q,         // read that we're selecting roots for
		const DescentQuery* qo,        // opposite mate, if applicable
		EList<DescentConfig>& confs,   // put DescentConfigs here
		EList<DescentRoot> roots);     // put DescentRoot here

protected:

	SimpleFunc descCons_;
	SimpleFunc rootIval_;
	size_t landing_;
};

/**
 * This class is the glue between a DescentDriver and the dynamic programming
 * implementations in Bowtie 2.  The DescentDriver is used to find some very
 * high-scoring alignments, but is additionally used to rank partial alignments
 * so that they can be extended using dynamic programming.
 */
class AlignerDriver {

public:

	AlignerDriver(
		const SimpleFunc& descCons,
		const SimpleFunc& rootIval,
		size_t landing,
		size_t totsz) :
		sel_(descCons, rootIval, landing),
		stop_(totsz, 0, 0)
	{
	}
	
	/**
	 * Initialize driver with respect to a new read or pair.
	 */
	void initRead(
		const AlignerDriverQuery& q1,
		const AlignerDriverQuery* q2)
	{
		dr1_.initRead(q1, q2, &sel_);
		paired_ = false;
		if(q2 != NULL) {
			dr2_.initRead(*q2, &q1, &sel_);
			paired_ = true;
		} else {
			dr2_.reset();
		}
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
	void go(
		const Scoring& sc,
		const Ebwt& ebwtFw,
		const Ebwt& ebwtBw,
		DescentMetrics& met,
		RandomSource& rnd,
		AlnSinkWrap& sink)
	{
		if(paired_) {
			// Paired-end - alternate between advancing dr1_ / dr2_ whenever a
			// new full alignment is discovered in the one currently being
			// advanced.  Whenever a new full alignment is found, check to see
			// if it pairs with a previously discovered alignment.
			bool first1 = (rnd.nextU2() == 0);
			bool first = true;
			while(true) {
				if(first && first1) {
					dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met);
				}
				dr2_.advance(stop_, sc, ebwtFw, ebwtBw, met);
				first = false;
			}
		} else {
			// Unpaired
		}
	}
	
	/**
	 * Reset state of all DescentDrivers.
	 */
	void reset() {
		dr1_.reset();
		dr2_.reset();
	}

protected:

	AlignerDriverRootSelector sel_;
	DescentDriver dr1_;
	DescentDriver dr2_;
	DescentStoppingConditions stop_;
	bool paired_;
};

#endif /* defined(ALIGNER_DRIVER_H_) */
