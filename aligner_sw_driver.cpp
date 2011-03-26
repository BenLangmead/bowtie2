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
	int refGaps = pen.maxRefGaps(penceil);
	const size_t rdlen = rd.length();
	coords1_.clear();   // ref coords tried so far
	coords1st_.clear(); // upstream coord for hits so far
	coords1en_.clear(); // downstream coord for hits so far

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
			if(!coords1_.insert(c)) {
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
			if(!coords1st_.insert(st) || !coords1en_.insert(en)) {
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
	EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
	EList<SwActionSink*>* swActionSinks)   // send action-list updates to these
{
	assert(!reportImmediately || msink != NULL);
	//assert(!reportImmediately || msink->empty());
	assert(!reportImmediately || !msink->maxed());

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
	coords1_.clear();
	coords1st_.clear(); // upstream coord for mate 1 hits so far
	coords1en_.clear(); // downstream coord for mate 2 hits so far
	coords2_.clear();
	coords2st_.clear(); // upstream coord for mate 2 hits so far
	coords2en_.clear(); // downstream coord for mate 2 hits so far
	
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
		ESet<Coord> *coords, *coordsst, *coordsen;
		// Let coords, coordsst and coordsen point to coordinates for the
		// anchor mate
		if(anchor1) {
			coords   = &coords1_;
			coordsst = &coords1st_;
			coordsen = &coords1en_;
		} else {
			coords   = &coords2_;
			coordsst = &coords2st_;
			coordsen = &coords2en_;
		}
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
			if(!coords->insert(c)) {
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
			if(!coordsst->insert(stAnchor) || !coordsen->insert(enAnchor)) {
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
			if(found && swMateImmediately) {
				bool oleft = false, ofw = false;
				int64_t oll = 0, olr = 0, orl = 0, orr = 0;
				foundMate = pepol.otherMate(
					anchor1,             // anchor mate is mate #1?
					fw,                  // anchor aligned to Watson?
					res_.alres.refoff(), // offset of anchor mate
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
					off1 = anchor1 ? res_.alres.refoff() : ores_.alres.refoff();
					off2 = anchor1 ? ores_.alres.refoff() : res_.alres.refoff();
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
						off2,
						len1,
						len2,
						fw1,
						fw2);
					foundMate = pairCl != PE_ALS_DISCORD;
				}
				if(foundMate) {
					// Check if this fragment is redundant with one found
					// previously, i.e. if the extents are the same and mate 1
					// has the same orientation.
					Interval ival(refid, fragoff, fw1, fraglen);
					foundMate = coordspair_.insert(ival);
				}
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
					} else {
						// Report unpaired hit for anchor
						assert(msink != NULL);
						assert(res_.repOk());
						// Check that alignment accurately reflects the
						// reference characters aligned to
						assert(res_.alres.matchesRef(rd, ref));
						// Report an unpaired alignment
						assert(!msink->maxed());
						if(msink->report(
							0,
							anchor1 ? &res_.alres : NULL,
							anchor1 ? NULL : &res_.alres))
						{
							// Should not short-circuit on an unpaired
							// alignment when more paired alignments are
							// possible
							cerr << "Should not get true return value from report() for unpaired mate alignment" << endl;
							throw 1;
						}
					}
				}
			}
			
			// At this point we know that we aren't bailing, and will continue to resolve seed hits.  

		} // while(!gw.done())
	}
	return false;
}
