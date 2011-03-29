/*
 * aligner_sw_driver.cpp
 *
 * Routines that drive the alignment process given a collection of seed hits.
 * This is generally done in a few stages: extendSeeds visits the set of
 * seed-hit BW elements in some order; for each element visited it resolves its
 * reference offset; once the reference offset is known, bounds for a dynamic
 * programming subproblem are established; if these bounds are distinct from
 * the bounds we've already tried, we solve the dynamic programming subproblem
 * and report the hit; if the AlnSinkWrap indicates that we can stop, we
 * return, otherwise we continue on to the next BW element.
 */

#include <iostream>
#include "aligner_cache.h"
#include "aligner_sw_driver.h"
#include "pe.h"
#include "dp_framer.h"

using namespace std;

/**
 * Try to convert a BW offset to the forward reference offset of the
 * leftmost position involved in the hit.  Return 0xffffffff if the
 * index doesn't contain the information needed to convert.
 */
uint32_t SwDriver::bwtOffToOff(
	const Ebwt& ebwt,
	bool fw,
	uint32_t bwoff,
	uint32_t hitlen,
	uint32_t& tidx,
	uint32_t& toff,
	uint32_t& tlen)
{
	assert_gt(hitlen, 0);
	uint32_t off = 0xffffffff;
	if((bwoff & ebwt.eh().offMask()) == bwoff) {
		// The index tells us the offset of this BW row directly
		uint32_t bwoffOff = bwoff >> ebwt.eh().offRate();
		assert_lt(bwoffOff, ebwt.eh().offsLen());
		off = ebwt.offs()[bwoffOff];
		assert_neq(0xffffffff, off);
		if(!fw) {
			assert_lt(off, ebwt.eh().len());
			off = ebwt.eh().len() - off - 1;
			assert_geq(off, hitlen-1);
			off -= (hitlen-1);
			assert_lt(off, ebwt.eh().len());
		}
	}
	if(off != 0xffffffff) {
		ebwt.joinedToTextOff(
			hitlen,
			off,
			tidx,
			toff,
			tlen);
	}
	return off;
}

/**
 * Given a collection of SeedHits for a single read, extend seed alignments
 * into full alignments.  Where possible, try to avoid redundant offset lookups
 * and dynamic programming wherever possible.  Optionally report alignments to
 * a AlnSinkWrap object as they are discovered.
 *
 * If 'reportImmediately' is true, returns true iff a call to msink->report()
 * returned true (indicating that the reporting policy is satisfied and we can
 * stop).  Otherwise, returns false.
 */
bool SwDriver::extendSeeds(
	const Read& rd,              // read to align
	bool color,                  // true -> read is colorspace
	SeedResults& sh,             // seed hits to extend into full alignments
	const Ebwt& ebwt,            // BWT
	const BitPairReference& ref, // Reference strings
	GroupWalk& gw,               // group walk left
	SwAligner& swa,              // dynamic programming aligner
	const SwParams& pa,          // parars1_meters for dyn prog aligner
	const Penalties& pen,        // penalties for edits
	int seedmms,                 // # mismatches allowed in seed
	int seedlen,                 // length of seed
	int seedival,                // interval between seeds
	int penceil,                 // maximum penalty allowed
	uint32_t maxposs,            // stop after examining this many positions (offset+orientation combos)
	uint32_t maxrows,            // stop examining a position after this many offsets are reported
	AlignmentCacheIface& sc,     // alignment cache for seed hits
	RandomSource& rnd,           // pseudo-random source
	WalkMetrics& wlm,            // group walk left metrics
	SwMetrics& swm,              // dynamic programming metrics
	ReportingMetrics& rpm,       // reporting metrics
	AlnSinkWrap* msink,        // AlnSink wrapper for multiseed-style aligner
	bool reportImmediately,      // whether to report hits immediately to msink
	EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
	EList<SwActionSink*>* swActionSinks)   // send action-list updates to these
{
	assert(!reportImmediately || msink != NULL);
	assert(!reportImmediately || msink->empty());
	assert(!reportImmediately || !msink->maxed());

	// Calculate the largest possible number of read and reference gaps
	// given 'penceil' and 'pen'
	int readGaps = pen.maxReadGaps(penceil);
	int refGaps  = pen.maxRefGaps(penceil);
	const size_t rdlen = rd.length();
	coords_.clear();   // ref coords tried so far
	coordsst_.clear(); // upstream coord for hits so far
	coordsen_.clear(); // downstream coord for hits so far

	DynProgFramer dpframe(!gReportOverhangs);

	// Iterate twice through levels seed hits from the lowest ranked
	// level to the highest ranked.  On the first iteration, look for
	// entries for which the offset is already known and try SWs.  On
	// the second iteration, resolve entries for which the offset is
	// unknown and try SWs.
	const size_t nonz = sh.nonzeroOffsets();
	const size_t poss = min<size_t>(nonz, maxposs);
	for(size_t i = 0; i < poss; i++) {
		bool fw = true;
		uint32_t offidx = 0, rdoff = 0, seedlen;
		QVal qv = sh.hitsByRank(i, offidx, rdoff, fw, seedlen);
		assert(qv.repOk(sc.current()));
		if(!fw) {
			// 'rdoff' and 'offidx' are with respect to the 5' end of
			// the read.  Here we convert rdoff to be with respect to
			// the upstream (3') end of ther read.
			rdoff = (uint32_t)(rdlen - rdoff - seedlen);
		}
		size_t rows = rdlen + (color ? 1 : 0);
		satups_.clear();
		sc.queryQval(qv, satups_);
		gw.initQval(ebwt, ref, qv, sc, rnd, maxrows, true, wlm);
		assert(gw.initialized());
		int advances = 0;
		ASSERT_ONLY(WalkResult lastwr);
		while(!gw.done()) {
			// Resolve next element offset
			WalkResult wr;
			gw.advancePos(false, 0, wr, wlm);
			assert(wr.elt != lastwr.elt);
			ASSERT_ONLY(lastwr = wr);
			advances++;
			ASSERT_ONLY(uint32_t off = wr.toff);
			assert_neq(0xffffffff, off);
			uint32_t tidx = 0, toff = 0, tlen = 0;
			ebwt.joinedToTextOff(
				wr.elt.len,
				wr.toff,
				tidx,
				toff,
				tlen);
			tlen += (color ? 1 : 0);
			if(tidx == 0xffffffff) {
				// The seed hit straddled a reference boundary so the seed hit
				// isn't valid
				continue;
			}
			// Now that we have a seed hit, there are many issues to solve
			// before we have a completely framed dynamic programming problem.
			// They include:
			//
			// 1. Setting reference offsets on either side of the seed hit,
			//    accounting for where the seed occurs in the read
			// 2. Adjusting the width of the banded dynamic programming problem
			//    and adjusting reference bounds to allow for gaps in the
			//    alignment
			// 3. Accounting for the edges of the reference, which can impact
			//    the width of the DP problem and reference bounds.
			// 4. Perhaps filtering the problem down to a smaller problem based
			//    on what DPs we've already solved for this read
			//
			// We do #1 here, since it is simple and we have all the seed-hit
			// information here.  #2 and #3 are handled in the DynProgFramer.
			
			// Find offset of alignment's upstream base assuming net gaps=0
			// between beginning of read and beginning of seed hit
			int64_t refoff = (int64_t)toff - rdoff;
			// TODO: need a more sophisticated filter here.  Really we care if
			// any of the start/end combos here have been covered by a previous
			// dynamic programming problem.
			Coord c(tidx, refoff, fw);
			if(!coords_.insert(c)) {
				// Already tried to find an alignment at these
				// coordinates
				swm.rshit++;
				continue;
			}
			size_t width = 0, trimup = 0, trimdn = 0;
			int64_t refl = 0, refr = 0;
			bool found = dpframe.frameSeedExtension(
				refoff,   // ref offset implied by seed hit assuming no gaps
				rows,     // length of read sequence used in DP table (so len
				          // of +1 nucleotide sequence for colorspace reads)
				tlen,     // length of reference
				readGaps, // max # of read gaps permitted in opp mate alignment
				refGaps,  // max # of ref gaps permitted in opp mate alignment
				width,    // out: calculated width stored here
				trimup,   // out: number of bases trimmed from upstream end
				trimdn,   // out: number of bases trimmed from downstream end
				refl,     // out: ref pos of upper LHS of parallelogram
				refr,     // out: ref pos of lower RHS of parallelogram
				st_,      // out: legal starting columns stored here
				en_);     // out: legal ending columns stored here
			if(!found) continue;
			assert_eq(width, st_.size());
			assert_eq(st_.size(), en_.size());
			res_.reset();
			assert(res_.empty());
			assert_neq(0xffffffff, tidx);
			// Given the boundaries defined by refl and refr, initilize the
			// SwAligner with the dynamic programming problem that aligns the
			// read to this reference stretch.
			swa.init(
				rd,        // read to align
				0,         // off of first char in 'rd' to consider
				rdlen,     // off of last char (excl) in 'rd' to consider
				fw,        // whether to align forward or revcomp read
				color,     // colorspace?
				tidx,      // reference aligned against
				refl,      // off of first character in ref to consider
				refr+1,    // off of last char (excl) in ref to consider
				ref,       // Reference strings
				tlen,      // length of reference sequence
				width,     // # bands to do (width of parallelogram)
				&st_,      // mask indicating which columns we can start in
				&en_,      // mask indicating which columns we can end in
				pa,        // dynamic programming parameters
				pen,       // penalty scheme
				penceil);  // penalty ceiling for valid alignments
			// Now fill the dynamic programming matrix and return true iff
			// there is at least one valid alignment
			found = swa.align(
				res_,
				rnd);
			assert( found ||  res_.empty());
			assert(!found || !res_.empty());
			swm.update(res_);
			if(found) {
				res_.swsucc++;
			} else {
				res_.swfail++;
				continue;
			}
			
			// User specified that alignments overhanging ends of reference
			// should be excluded...
			assert(gReportOverhangs || res_.alres.within(tidx, 0, fw, tlen));

			Coord st, en;
			res_.alres.getCoords(st, en);
			if(!coordsst_.insert(st) || !coordsen_.insert(en)) {
				// Redundant with an alignment we found already
				continue;
			}
			
			// Annotate the AlnRes object with some key parameters
			// that were used to obtain the alignment.
			res_.alres.setParams(
				seedmms,   // # mismatches allowed in seed
				seedlen,   // length of seed
				seedival,  // interval between seeds
				penceil);  // maximum penalty for valid alignment
			
			if(reportImmediately) {
				assert(msink != NULL);
				assert(res_.repOk());
				// Check that alignment accurately reflects the
				// reference characters aligned to
				assert(res_.alres.matchesRef(rd, ref));
				// Report an unpaired alignment
				assert(!msink->maxed());
				if(msink->report(0, &res_.alres, NULL)) {
					// Short-circuited because a limit, e.g. -k, -m or
					// -M, was exceeded
					return true;
				}
			}

			// At this point we know that we aren't bailing, and will continue to resolve seed hits.  

		} // while(!gw.done())
	}
	return false;
}

/**
 * Given a read, perform full dynamic programming against the entire
 * reference.  Optionally report alignments to a AlnSinkWrap object
 * as they are discovered.
 *
 * If 'reportImmediately' is true, returns true iff a call to
 * msink->report() returned true (indicating that the reporting
 * policy is satisfied and we can stop).  Otherwise, returns false.
 */
bool SwDriver::sw(
	const Read& rd,              // read to align
	bool color,                  // true -> read is colorspace
	const BitPairReference& ref, // Reference strings
	SwAligner& swa,              // dynamic programming aligner
	const SwParams& pa,          // parameters for dynamic programming aligner
	const Penalties& pen,        // penalties for edits
	int penceil,                 // maximum penalty allowed
	RandomSource& rnd,           // pseudo-random source
	SwMetrics& swm,              // dynamic programming metrics
	AlnSinkWrap* msink,          // HitSink for multiseed-style aligner
	bool reportImmediately,      // whether to report hits immediately to msink
	EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
	EList<SwActionSink*>* swActionSinks)   // send action-list updates to these
{
	assert(!reportImmediately || msink != NULL);
	return false;
}

/**
 * Given a collection of SeedHits for both mates in a read pair, extend seed
 * alignments into full alignments and then look for the opposite mate using
 * dynamic programming.  Where possible, try to avoid redundant offset lookups.
 * Optionally report alignments to a AlnSinkWrap object as they are discovered.
 *
 * If 'reportImmediately' is true, returns true iff a call to
 * msink->report() returned true (indicating that the reporting
 * policy is satisfied and we can stop).  Otherwise, returns false.
 *
 * REDUNDANT ALIGNMENTS
 *
 * In a paired-end context, it's not easy to pin down what a "redundant"
 * alignment is.  We break this down into cases:
 *
 * 1. A paired-end alignment is clearly redundant with another paired-end
 *    alignment if both the fragment extremes are identical.
 * 2. An anchor or unpaired alignment is clearly redundant with another anchor
 *    or unpaired alignment if either extreme is identical.
 *
 * On the other hand, it's also clear that:
 *
 * 1. A paired-end alignment is NOT redundant with another paired-end alignment
 *    if only *one* fragment extreme is identical.
 *
 * MIXING PAIRED AND UNPAIRED ALIGNMENTS
 *
 * There are distinct paired-end alignment modes for the cases where (a) the
 * user does or does not want to see unpaired alignments for individual mates
 * when there are no reportable paired-end alignments involving both mates, and
 * (b) the user does or does not want to see discordant paired-end alignments.
 * The modes have implications for this function and for the AlnSinkWrap, since
 * it affects when we're "done."  Also, whether the user has asked us to report
 * discordant alignments affects whether and how much searching for unpaired
 * alignments we must do (i.e. if there are no paired-end alignments, we must
 * at least do -m 1 for both mates).
 *
 * Mode 1: Just concordant paired-end.  Print only concordant paired-end
 * alignments.  As soon as any limits (-k/-m/-M) are reached, stop.
 *
 * Mode 2: Concordant and discordant paired-end.  If -k/-m/-M limits are
 * reached for paired-end alignments, stop.  Otherwise, if no paired-end
 * alignments are found, align both mates in an unpaired -m 1 fashion.  If
 * there is exactly one unpaired alignment for each mate, report the
 * combination as a discordant alignment.
 *
 * Mode 3: Concordant paired-end if possible, otherwise unpaired.  If -k/-M
 * limit is reached for paired-end alignmnts, stop.  If -m limit is reached for
 * paired-end alignments or no paired-end alignments are found, align both
 * mates in an unpaired fashion.  All the same settings governing validity and
 * reportability in paired-end mode apply here too (-k/-m/-M/etc).
 *
 * Mode 4: Concordant or discordant paired-end if possible, otherwise unpaired.
 * If -k/-M limit is reached for paired-end alignmnts, stop.  If -m limit is
 * reached for paired-end alignments or no paired-end alignments are found,
 * align both mates in an unpaired fashion.  If the -m limit was reached, there
 * is no need to search for a discordant alignment, and unapired alignment can
 * proceed as in Mode 3.  If no paired-end alignments were found, then unpaired
 * alignment proceeds as in Mode 3 but with this caveat: alignment must be at
 * least as thorough as dictated by -m 1 up until the point where
 *
 *Print paired-end alignments when there are reportable paired-end
 * alignments, otherwise report reportable unpaired alignments.  If -k limit is
 * reached for paired-end alignments, stop.  If -m/-M limit is reached for
 * paired-end alignments, stop searching for paired-end alignments and look
 * only for unpaired alignments.  If searching only for unpaired alignments,
 * respect -k/-m/-M limits separately for both mates.
 *
 * The return value from the AlnSinkWrap's report member function must be
 * specific enough to distinguish between:
 *
 * 1. Stop searching for paired-end alignments
 * 2. Stop searching for alignments for unpaired alignments for mate #1
 * 3. Stop searching for alignments for unpaired alignments for mate #2
 * 4. Stop searching for any alignments
 *
 * Note that in Mode 2, options affecting validity and reportability of
 * alignments apply .  E.g. if -m 1 is specified
 *
 * WORKFLOW
 *
 * Our general approach to finding paired and unpaired alignments here
 * is as follows:
 *
 * - For mate in mate1, mate2:
 *   - For each seed hit in mate:
 *     - Try to extend it into a full alignment; if we can't, continue
 *       to the next seed hit
 *     - Look for alignment for opposite mate; if we can't find one,
 *     - 
 *     - 
 *
 */
bool SwDriver::extendSeedsPaired(
	const Read& rd,              // mate to align as anchor
	const Read& ord,             // mate to align as opposite
	bool anchor1,                // true iff anchor mate is mate1
	bool color,                  // true -> reads are colorspace
	SeedResults& sh,             // seed hits for anchor
	const Ebwt& ebwt,            // BWT
	const BitPairReference& ref, // Reference strings
	GroupWalk& gw,               // group walk left
	SwAligner& swa,              // dynamic programming aligner
	const SwParams& pa,          // parameters for dynamic programming aligner
	const Penalties& pen,        // penalties for edits
	const PairedEndPolicy& pepol,// paired-end policy
	int seedmms,                 // # mismatches allowed in seed
	int seedlen,                 // length of seed
	int seedival,                // interval between seeds
	int penceil,                 // maximum penalty allowed for anchor
	int openceil,                // maximum penalty allowed for opposite
	uint32_t maxposs,            // stop after examining this many positions (offset+orientation combos)
	uint32_t maxrows,            // stop examining a position after this many offsets are reported
	AlignmentCacheIface& sc,     // alignment cache for seed hits
	RandomSource& rnd,           // pseudo-random source
	WalkMetrics& wlm,            // group walk left metrics
	SwMetrics& swm,              // dynamic programming metrics
	ReportingMetrics& rpm,       // reporting metrics
	AlnSinkWrap* msink,          // AlnSink wrapper for multiseed-style aligner
	bool swMateImmediately,      // whether to look for mate immediately
	bool reportImmediately,      // whether to report hits immediately to msink
	bool discord,                // look for discordant alignments?
	bool mixed,                  // look for unpaired as well as paired alns?
	EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
	EList<SwActionSink*>* swActionSinks)   // send action-list updates to these
{
	assert(!reportImmediately || msink != NULL);
	//assert(!reportImmediately || msink->empty());
	assert(!reportImmediately || !msink->maxed());
	assert(!msink->state().doneWithMate(anchor1));

	// Calculate the largest possible number of read and reference gaps
	// given 'penceil' and 'pen'
	int readGaps  = pen.maxReadGaps(penceil);
	int refGaps   = pen.maxRefGaps(penceil);
	int oreadGaps = pen.maxReadGaps(openceil);
	int orefGaps  = pen.maxRefGaps(openceil);

	const size_t rdlen  = rd.length();
	const size_t ordlen = ord.length();
	const size_t rows   = rdlen  + (color ? 1 : 0);
	const size_t orows  = ordlen + (color ? 1 : 0);
	coords_.clear();   // DP upstream char for seed hits
	coordsst_.clear(); // upstream coord for anchor hits so far
	coordsen_.clear(); // downstream coord for anchor hits so far
	
	DynProgFramer dpframe(!gReportOverhangs);

	// Iterate twice through levels seed hits from the lowest ranked
	// level to the highest ranked.  On the first iteration, look for
	// entries for which the offset is already known and try SWs.  On
	// the second iteration, resolve entries for which the offset is
	// unknown and try SWs.
	const size_t nonz = sh.nonzeroOffsets();
	const size_t poss = min<size_t>(nonz, maxposs);
	for(size_t i = 0; i < poss; i++) {
		bool fw = true;
		uint32_t offidx = 0, rdoff = 0, seedlen;
		QVal qv = sh.hitsByRank(i, offidx, rdoff, fw, seedlen);
		assert(qv.repOk(sc.current()));
		if(!fw) {
			// 'rdoff' and 'offidx' are with respect to the 5' end of
			// the read.  Here we convert rdoff to be with respect to
			// the upstream (3') end of ther read.
			rdoff = (uint32_t)(rdlen - rdoff - seedlen);
		}
		satups_.clear();
		sc.queryQval(qv, satups_);
		gw.initQval(ebwt, ref, qv, sc, rnd, maxrows, true, wlm);
		assert(gw.initialized());
		int advances = 0;
		ASSERT_ONLY(WalkResult lastwr);
		while(!gw.done()) {
			// Resolve next element offset
			WalkResult wr;
			gw.advancePos(false, 0, wr, wlm);
			assert(wr.elt != lastwr.elt);
			ASSERT_ONLY(lastwr = wr);
			advances++;
			ASSERT_ONLY(uint32_t joff = wr.toff);
			assert_neq(0xffffffff, joff);
			uint32_t tidx = 0, toff = 0, tlen = 0;
			ebwt.joinedToTextOff(
				wr.elt.len,
				wr.toff,
				tidx,
				toff,
				tlen);
			tlen += (color ? 1 : 0);
			if(tidx == 0xffffffff) {
				// The seed hit straddled a reference boundary so the seed hit
				// isn't valid
				continue;
			}
			// Now that we have a seed hit, there are many issues to solve
			// before we have a completely framed dynamic programming problem.
			// They include:
			//
			// 1. Setting reference offsets on either side of the seed hit,
			//    accounting for where the seed occurs in the read
			// 2. Adjusting the width of the banded dynamic programming problem
			//    and adjusting reference bounds to allow for gaps in the
			//    alignment
			// 3. Accounting for the edges of the reference, which can impact
			//    the width of the DP problem and reference bounds.
			// 4. Perhaps filtering the problem down to a smaller problem based
			//    on what DPs we've already solved for this read
			//
			// We do #1 here, since it is simple and we have all the seed-hit
			// information here.  #2 and #3 are handled in the DynProgFramer.
			
			// Find offset of alignment's upstream base assuming net gaps=0
			// between beginning of read and beginning of seed hit
			int64_t refoff = (int64_t)toff - rdoff;
			// TODO: need a more sophisticated filter here.  Really we care if
			// any of the start/end combos here have been covered by a previous
			// dynamic programming problem.
			Coord c(tidx, refoff, fw);
			if(!coords_.insert(c)) {
				// Already tried to find an alignment at these
				// coordinates
				swm.rshit++;
				continue;
			}
			size_t width = 0, trimup = 0, trimdn = 0;
			int64_t refl = 0, refr = 0;
			bool found = dpframe.frameSeedExtension(
				refoff,   // ref offset implied by seed hit assuming no gaps
				rows,     // length of read sequence used in DP table (so len
		                  // of +1 nucleotide sequence for colorspace reads)
				tlen,     // length of reference
				readGaps, // max # of read gaps permitted in opp mate alignment
				refGaps,  // max # of ref gaps permitted in opp mate alignment
				width,    // out: calculated width stored here
				trimup,   // out: number of bases trimmed from upstream end
				trimdn,   // out: number of bases trimmed from downstream end
				refl,     // out: ref pos of upper LHS of parallelogram
				refr,     // out: ref pos of lower RHS of parallelogram
				st_,      // out: legal starting columns stored here
				en_);     // out: legal ending columns stored here
			if(!found) {
				continue;
			}
			assert_eq(width, st_.size());
			assert_eq(st_.size(), en_.size());
			res_.reset();
			assert(res_.empty());
			assert_neq(0xffffffff, tidx);
			// Given the boundaries defined by refl and refr, initilize the
			// SwAligner with the dynamic programming problem that aligns the
			// read to this reference stretch.
			swa.init(
				rd,        // read to align
				0,         // off of first char in 'rd' to consider
				rdlen,     // off of last char (excl) in 'rd' to consider
				fw,        // whether to align forward or revcomp read
				color,     // colorspace?
				tidx,      // reference aligned against
				refl,      // off of first character in ref to consider
				refr+1,    // off of last char (excl) in ref to consider
				ref,       // Reference strings
				tlen,      // length of reference sequence
				width,     // # bands to do (width of parallelogram)
				&st_,      // mask indicating which columns we can start in
				&en_,      // mask indicating which columns we can end in
				pa,        // dynamic programming parameters
				pen,       // penalty scheme
				penceil);  // penalty ceiling for valid alignments
			// Now fill the dynamic programming matrix and return true iff
			// there is at least one valid alignment
			found = swa.align(
				res_,
				rnd);
			assert( found ||  res_.empty());
			assert(!found || !res_.empty());
			swm.update(res_);
			if(found) {
				res_.swsucc++;
			} else {
				res_.swfail++;
				continue;
			}

			// User specified that alignments overhanging ends of reference
			// should be excluded...
			assert(gReportOverhangs || res_.alres.within(tidx, 0, fw, tlen));

			Coord stAnchor, enAnchor;
			res_.alres.getCoords(stAnchor, enAnchor);
			if(!coordsst_.insert(stAnchor) || !coordsen_.insert(enAnchor)) {
				// Redundant with an alignment we found already
				continue;
			}
			
			// Annotate the AlnRes object with some key parameters
			// that were used to obtain the alignment.
			res_.alres.setParams(
				seedmms,   // # mismatches allowed in seed
				seedlen,   // length of seed
				seedival,  // interval between seeds
				penceil);  // maximum penalty for valid alignment
			
			bool foundMate = false;
			TRefOff off = res_.alres.refoff();
			if(found && swMateImmediately) {
				bool oleft = false, ofw = false;
				int64_t oll = 0, olr = 0, orl = 0, orr = 0;
				assert(!msink->state().done());
				if(!msink->state().doneConcordant()) {
					foundMate = pepol.otherMate(
						anchor1,             // anchor mate is mate #1?
						fw,                  // anchor aligned to Watson?
						off,                 // offset of anchor mate
						orows + oreadGaps,   // max # columns spanned by alignment
						tlen,                // reference length
						anchor1 ? rd.length()  : ord.length(), // mate #1 length
						anchor1 ? ord.length() : rd.length(),  // mate #2 length
						oleft,               // out: look left for opposite mate?
						oll,
						olr,
						orl,
						orr,
						ofw);
				} else {
					// We're no longer interested in finding additional
					// concordant paired-end alignments so we just report this
					// mate's alignment as an unpaired alignment (below)
				}
				size_t owidth = 0, otrimup = 0, otrimdn = 0;
				int64_t orefl = 0, orefr = 0;
				if(foundMate) {
					foundMate = dpframe.frameFindMate(
						!oleft,      // true iff anchor alignment is to the left
						oll,         // leftmost Watson off for LHS of opp aln
						olr,         // rightmost Watson off for LHS of opp aln
						orl,         // leftmost Watson off for RHS of opp aln
						orr,         // rightmost Watson off for RHS of opp aln
						orows,       // length of opposite mate
						tlen,        // length of reference sequence aligned to
						oreadGaps,   // max # of read gaps in opp mate aln
						orefGaps,    // max # of ref gaps in opp mate aln
						owidth,      // out: calculated width stored here
						otrimup,     // out: # bases trimmed from upstream end
						otrimdn,     // out: # bases trimmed from downstream end
						orefl,       // out: ref pos of upper LHS of parallelogram
						orefr,       // out: ref pos of lower RHS of parallelogram
						st_,         // out: legal starting columns stored here
						en_);        // out: legal ending columns stored here
				}
				if(foundMate) {
					ores_.reset();
					assert(ores_.empty());
					// Given the boundaries defined by refi and reff, initilize
					// the SwAligner with the dynamic programming problem that
					// aligns the read to this reference stretch.
					swa.init(
						ord,       // read to align
						0,         // off of first char in rd to consider
						ordlen,    // off of last char (excl) in rd to consider
						ofw,       // whether to align forward or revcomp read
						color,     // colorspace?
						tidx,      // reference aligned against
						oll,       // off of first character in rf to consider
						orr,       // off of last char (excl) in rf to consider
						ref,       // Reference strings
						tlen,      // length of reference sequence
						owidth,    // # bands to do (width of parallelogram)
						&st_,      // mask of which cols we can start in
						&en_,      // mask of which cols we can end in
						pa,        // dynamic programming parameters
						pen,       // penalty scheme
						openceil); // penalty ceiling for valid alignments
					// Now fill the dynamic programming matrix and return true
					// iff there is at least one valid alignment
					foundMate = swa.align(ores_, rnd);
				}
				if(foundMate) {
					assert_eq(ofw, ores_.alres.fw());
					// Annotate the AlnRes object with some key parameters
					// that were used to obtain the alignment.
					ores_.alres.setParams(
						seedmms,   // # mismatches allowed in seed
						seedlen,   // length of seed
						seedival,  // interval between seeds
						openceil); // maximum penalty for valid alignment
					if(!gReportOverhangs &&
					   !ores_.alres.within(tidx, 0, ofw, tlen))
					{
						foundMate = false;
					}
				}
				TRefId refid;
				TRefOff off1, off2;
				TRefOff fragoff;
				size_t len1, len2, fraglen;
				bool fw1, fw2;
				int pairCl;
				if(foundMate) {
					refid = res_.alres.refid();
					assert_eq(refid, ores_.alres.refid());
					off1 = anchor1 ? off : ores_.alres.refoff();
					off2 = anchor1 ? ores_.alres.refoff() : off;
					len1 = anchor1 ? res_.alres.extent() : ores_.alres.extent();
					len2 = anchor1 ? ores_.alres.extent() : res_.alres.extent();
					fw1  = anchor1 ? res_.alres.fw() : ores_.alres.fw();
					fw2  = anchor1 ? ores_.alres.fw() : res_.alres.fw();
					fragoff = min<TRefOff>(off1, off2);
					fraglen = max<TRefOff>(
						off1 - fragoff + len1,
						off2 - fragoff + len2);
					// Check that final mate alignments are consistent with
					// paired-end fragment constraints
					pairCl = pepol.peClassifyPair(
						off1,
						len1,
						fw1,
						off2,
						len2,
						fw2);
					foundMate = pairCl != PE_ALS_DISCORD;
				}
				if(foundMate) {
					// Check if this fragment is redundant with one found
					// previously, i.e. if the extents are the same and mate 1
					// has the same orientation.
					Interval ival(refid, fragoff, fw1, fraglen);
					foundMate = frags_.insert(ival);
				}
				// Check whether we've already seen these *mate* alignments
				bool seen1 = false, seen2 = false;
				if(foundMate) {
					// Remember that we saw each of these two mates
					Coord c1l(tidx, off1,            fw1);
					Coord c1r(tidx, off1 + len1 - 1, fw1);
					seen1 = !coords1seenup_.insert(c1l) ||
							!coords1seendn_.insert(c1r);
					Coord c2l(tidx, off2,            fw2);
					Coord c2r(tidx, off2 + len2 - 1, fw2);
					seen2 = !coords2seenup_.insert(c2l) ||
							!coords2seendn_.insert(c2r);
				}
				bool doneAnchor = anchor1 ?
					msink->state().doneUnpaired(true) :
					msink->state().doneUnpaired(false);
				if(reportImmediately) {
					if(foundMate) {
						// Report pair to the AlnSinkWrap
						assert(msink != NULL);
						assert(res_.repOk());
						assert(ores_.repOk());
						// Check that alignment accurately reflects the
						// reference characters aligned to
						assert(res_.alres.matchesRef(rd, ref));
						assert(ores_.alres.matchesRef(ord, ref));
						// Report an unpaired alignment
						assert(!msink->maxed());
						if(msink->report(
							0,
							anchor1 ? &res_.alres : &ores_.alres,
							anchor1 ? &ores_.alres : &res_.alres))
						{
							// Short-circuited because a limit, e.g. -k, -m or
							// -M, was exceeded
							return true;
						}
						if(mixed || discord) {
							// Report alignment for mate #1 as an unpaired
							// alignment
							bool seenAnchor   = anchor1 ? seen1 : seen2;
							bool seenOpposite = anchor1 ? seen2 : seen1;
							bool doneOpposite = anchor1 ?
								msink->state().doneUnpaired(false) :
								msink->state().doneUnpaired(true);
							if(!seenAnchor && !doneAnchor && msink->report(
								0,
								anchor1 ? &res_.alres : NULL,
								anchor1 ? NULL : &res_.alres))
							{
								return true; // Short-circuited
							}
							// Report alignment for mate #2 as an unpaired
							// alignment
							if(!seenOpposite && !doneOpposite && msink->report(
								0,
								anchor1 ? NULL : &ores_.alres,
								anchor1 ? &ores_.alres : NULL))
							{
								return true; // Short-circuited
							}
						}
						if(msink->state().doneWithMate(anchor1)) {
							// We're now done with the mate that we're
							// currently using as our anchor.  We're not with
							// the read overall.
							return false;
						}
					} else if(mixed || discord) {
						// Report unpaired hit for anchor
						assert(msink != NULL);
						assert(res_.repOk());
						// Check that alignment accurately reflects the
						// reference characters aligned to
						assert(res_.alres.matchesRef(rd, ref));
						// Report an unpaired alignment
						assert(!msink->maxed());
						bool seen = false;
						Coord cl(tidx, off,                           fw);
						Coord cr(tidx, off + res_.alres.extent() - 1, fw);
						if(anchor1) {
							seen = !coords1seenup_.insert(cl) ||
							       !coords1seendn_.insert(cr);
						} else {
							seen = !coords2seenup_.insert(cl) ||
							       !coords2seendn_.insert(cr);
						}
						if(!seen && !doneAnchor && msink->report(
							0,
							anchor1 ? &res_.alres : NULL,
							anchor1 ? NULL : &res_.alres))
						{
							return true; // Short-circuited
						}
						if(msink->state().doneWithMate(anchor1)) {
							// Done with this mate, but not with the read
							// overall
							return false;
						}
					}
				}
				
			} // if(found && swMateImmediately)
			
			// At this point we know that we aren't bailing, and will continue to resolve seed hits.  

		} // while(!gw.done())
	
	} // for(size_t i = 0; i < poss; i++)
	return false;
}
