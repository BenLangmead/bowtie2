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

#include "aligner_driver.h"

void PrioritizedRootSelector::select(
	const Read& q,
	const Read* qo,
	bool nofw,
	bool norc,
	EList<DescentConfig>& confs,
	EList<DescentRoot>& roots)
{
	assert_gt(landing_, 0);
	// To specify a search root, we must specify an offset from the 5' end,
	// whether it is left-to-right, and whether it searchers over the read or
	// its reverse complement.
	
	// Note that it's not very sensible to pick a search root going
	// left-to-right but where its offset puts it very close to the
	// right-hand-side of the read.  I.e. that root will "bounce" almost
	// immediately and go in the other direction.
	
	// How to pick these roots?  One idea is to simply lay down roots every N
	// positions along the read and its reverse-complement.  We would do this
	// both for left-pointing and right-pointing roots.
	
	
	// Another way is to consider every possible root, then rank them according
	// to how optimistic we are that picking that root will be productive.
	// Things that make us more optimistic are:
	//
	// 1. First several read characters to align are high quality
	// 2. First several read characters to align are free of Ns
	// 3. First several read characters do not form a simple repeat
	// 4. First several k-mers are well represented both in other reads and in
	//    the reference genome
	// 5. Characters, k-mers just before the root are "bad"
	// 6. Root is flush with one end of the read or the other
	
	// Go left-to-right along the forward and reverse-complement reads,
	// compiling info about the nucleotides in the landing zone of each
	// potential R2L root.
	const int nPenalty = 150;
	const int endBonus = 150;
	const size_t qlen = q.length();
	// Calculate interval length
	int interval = rootIval_.f<int>((double)qlen);
	size_t sizeTarget = qlen - landing_ + 1;
	sizeTarget = (size_t)(ceilf((sizeTarget / (float)interval)));
	sizeTarget *= 4;
	// Set up initial score arrays
	for(int i = 0; i < 2; i++) {
		bool fw = (i == 0);
		scoresOrig_[i].resize(qlen);
		scores_[i].resize(qlen);
		for(size_t j = 0; j < qlen; j++) {
			size_t off5p = fw ? j : (qlen - j - 1);
			int c = q.getc(off5p, fw);
			int sc = q.getq(off5p) - ((c > 3) ? nPenalty : 0);
			scoresOrig_[i][j] = scores_[i][j] = sc;
		}
	}
	rootHeap_.clear();
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			continue;
		}
		int pri = 0;
		size_t revi = qlen;
		for(size_t i = 0; i < qlen; i++) {
			revi--;
			pri += scoresOrig_[fwi][i];
			if(i >= landing_) {
				pri -= scoresOrig_[fwi][i - landing_];
			}
			if(i >= landing_-1 && scoresOrig_[fwi][i] > 0) {
				rootHeap_.insert(DescentRoot(
					fw ? i : revi, // offset from 5' end
					false,         // left-to-right?
					fw,            // fw?
					landing_,      // landing length
					qlen,          // query length
					pri + ((revi == 0) ? endBonus : 0))); // root priority
				// Give priority boost for being flush with one end or the
				// other
			}
		}
		pri = 0;
		size_t i = qlen - revi;
		for(size_t revi = 0; revi < qlen; revi++) {
			i--;
			pri += scoresOrig_[fwi][i];
			if(revi >= landing_) {
				pri -= scoresOrig_[fwi][i + landing_];
			}
			if(revi >= landing_-1 && scoresOrig_[fwi][i] > 0) {
				rootHeap_.insert(DescentRoot(
					fw ? i : revi, // offset from 5' end
					true,          // left-to-right?
					fw,            // fw?
					landing_,      // landing length
					qlen,          // query length
					pri + ((i == 0) ? endBonus : 0))); // root priority
				// Give priority boost for being flush with one end or the
				// other
			}
		}
	}
	// Now that all the roots are in a heap, we select them one-by-one.
	// Each time we select a root beyond the first, we check to see if an
	// already-selected root's landing area overlaps.  If so, we take away
	// any benefit associated with the bases/qualities in the landing area
	// and then push it back onto the heap if that changes its priority.
	while(roots.size() < sizeTarget) {
		if(rootHeap_.empty()) {
			break;
		}
		DescentRoot r = rootHeap_.pop();
		const size_t off = r.fw ? r.off5p : (qlen - r.off5p - 1);
		int fwi = r.fw ? 0 : 1;
		// Re-calculate priority
		int pri = 0;
		if(r.l2r) {
			for(size_t i = 0; i < landing_; i++) {
				pri += scores_[fwi][off + i];
			}
		} else {
			for(size_t i = 0; i < landing_; i++) {
				pri += scores_[fwi][off - i];
			}
		}
		// Must take end bonus into account when re-calculating
		if((r.l2r && (off == 0)) || (!r.l2r && (off == qlen - 1))) {
			pri += endBonus;
		}
		if(pri == r.pri) {
			// Update the positions in this root's landing area
			if(r.l2r) {
				for(size_t i = 0; i < landing_; i++) {
					float frac = ((float)i / (float)landing_);
					scores_[fwi][off + i] = (int)(scores_[fwi][off + i] * frac);
				}
			} else {
				for(size_t i = 0; i < landing_; i++) {
					float frac = ((float)i / (float)landing_);
					scores_[fwi][off - i] = (int)(scores_[fwi][off - i] * frac);
				}
			}
			confs.expand();
			confs.back().cons.init(landing_, consExp_);
			roots.push_back(r);
		} else {
			// Re-insert the root, its priority now changed
			assert_gt(roots.size(), 0);
			r.pri = pri;
			rootHeap_.insert(r);
		}
	}
	assert(!roots.empty());
	//std::cerr << roots.size() << ", " << ncandidates << std::endl;
}

void IntervalRootSelector::select(
	const Read& q,
	const Read* qo,
	bool nofw,
	bool norc,
	EList<DescentConfig>& confs,
	EList<DescentRoot>& roots)
{
	// To specify a search root, we must specify an offset from the 5' end,
	// whether it is left-to-right, and whether it searchers over the read or
	// its reverse complement.
	
	// Note that it's not very sensible to pick a search root going
	// left-to-right but where its offset puts it very close to the
	// right-hand-side of the read.  I.e. that root will "bounce" almost
	// immediately and go in the other direction.
	
	// How to pick these roots?  One idea is to simply lay down roots every N
	// positions along the read and its reverse-complement.  That's what we do
	// here.
	
	// Calculate interval length for both mates
	int interval = rootIval_.f<int>((double)q.length());
	if(qo != NULL) {
		// Boost interval length by 20% for paired-end reads
		interval = (int)(interval * 1.2 + 0.5);
	}
	float pri = 0.0f;
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			continue;
		}
		// Put down left-to-right roots w/r/t forward and reverse-complement reads
		{
			bool first = true;
			size_t i = 0;
			while(first || (i + landing_ <= q.length())) {
				confs.expand();
				confs.back().cons.init(landing_, consExp_);
				roots.expand();
				roots.back().init(
					i,          // offset from 5' end
					true,       // left-to-right?
					fw,         // fw?
					1,          // landing
					q.length(), // query length
					pri);       // root priority
				i += interval;
				first = false;
			}
		}
		// Put down right-to-left roots w/r/t forward and reverse-complement reads
		{
			bool first = true;
			size_t i = 0;
			while(first || (i + landing_ <= q.length())) {
				confs.expand();
				confs.back().cons.init(landing_, consExp_);
				roots.expand();
				roots.back().init(
					q.length() - i - 1, // offset from 5' end
					false,              // left-to-right?
					fw,                 // fw?
					1,          // landing
					q.length(),         // query length
					pri);               // root priority
				i += interval;
				first = false;
			}
		}
	}
	//std::cerr << roots.size() << std::endl;
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
int AlignerDriver::go(
	const Scoring& sc,
	const Ebwt& ebwtFw,
	const Ebwt& ebwtBw,
	const BitPairReference& ref,
	DescentMetrics& met,
	WalkMetrics& wlm,
	PerReadMetrics& prm,
	RandomSource& rnd,
	AlnSinkWrap& sink)
{
	if(paired_) {
		// Paired-end - alternate between advancing dr1_ / dr2_ whenever a
		// new full alignment is discovered in the one currently being
		// advanced.  Whenever a new full alignment is found, check to see
		// if it pairs with a previously discovered alignment.
		bool first1 = rnd.nextBool();
		bool first = true;
		DescentStoppingConditions stopc1 = stop_;
		DescentStoppingConditions stopc2 = stop_;
		size_t totszIncr = (stop_.totsz + 7) / 8;
		stopc1.totsz = totszIncr;
		stopc2.totsz = totszIncr;
		while(stopc1.totsz <= stop_.totsz && stopc2.totsz <= stop_.totsz) {
			if(first && first1 && stopc1.totsz <= stop_.totsz) {
				dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
				stopc1.totsz += totszIncr;
			}
			if(stopc2.totsz <= stop_.totsz) {
				dr2_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
				stopc2.totsz += totszIncr;
			}
			first = false;
		}
	} else {
		// Unpaired
		size_t iter = 1;
		while(true) {
			int ret = dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
			if(ret == DESCENT_DRIVER_ALN) {
				cerr << iter << ". DESCENT_DRIVER_ALN" << endl;
			} else if(ret == DESCENT_DRIVER_MEM) {
				cerr << iter << ". DESCENT_DRIVER_MEM" << endl;
				break;
			} else if(ret == DESCENT_DRIVER_STRATA) {
				// DESCENT_DRIVER_STRATA is returned by DescentDriver.advance()
				// when it has finished with a "non-empty" stratum: a stratum
				// in which at least one alignment was found.  Here we report
				// the alignments in an arbitrary order.
				AlnRes res;
				// Initialize alignment selector with the DescentDriver's
				// alignment sink
				alsel_.init(
					dr1_.query(),
					dr1_.sink(),
					ebwtFw,
					ref,
					rnd,
					wlm);
				while(!alsel_.done() && !sink.state().doneWithMate(true)) {
					res.reset();
					bool ret2 = alsel_.next(
						dr1_,
						ebwtFw,
						ref,
						rnd,
						res,
						wlm,
						prm);
					if(ret2) {
						// Got an alignment
						assert(res.matchesRef(
							dr1_.query(),
							ref,
							tmp_rf_,
							tmp_rdseq_,
							tmp_qseq_,
							raw_refbuf_,
							raw_destU32_,
							raw_matches_));
						// Get reference interval involved in alignment
						Interval refival(res.refid(), 0, res.fw(), res.reflen());
						assert_gt(res.refExtent(), 0);
						// Does alignment falls off end of reference?
						if(gReportOverhangs &&
						   !refival.containsIgnoreOrient(res.refival()))
						{
							res.clipOutside(true, 0, res.reflen());
							if(res.refExtent() == 0) {
								continue;
							}
						}
						assert(gReportOverhangs ||
							   refival.containsIgnoreOrient(res.refival()));
						// Alignment fell entirely outside the reference?
						if(!refival.overlapsIgnoreOrient(res.refival())) {
							continue; // yes, fell outside
						}
						// Alignment redundant with one we've seen previously?
						if(red1_.overlap(res)) {
							continue; // yes, redundant
						}
						red1_.add(res); // so we find subsequent redundancies
						// Report an unpaired alignment
						assert(!sink.state().doneWithMate(true));
						assert(!sink.maxed());
						if(sink.report(0, &res, NULL)) {
							// Short-circuited because a limit, e.g. -k, -m or
							// -M, was exceeded
							return ALDRIVER_POLICY_FULFILLED;
						}
					}
				}
				dr1_.sink().advanceStratum();
			} else if(ret == DESCENT_DRIVER_BWOPS) {
				cerr << iter << ". DESCENT_DRIVER_BWOPS" << endl;
				break;
			} else if(ret == DESCENT_DRIVER_DONE) {
				cerr << iter << ". DESCENT_DRIVER_DONE" << endl;
				break;
			} else {
				assert(false);
			}
			iter++;
		}
	}
	return ALDRIVER_EXHAUSTED_CANDIDATES;
}
