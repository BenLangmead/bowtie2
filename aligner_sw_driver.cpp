/*
 * aligner_sw_driver.cpp
 *
 * This file contains routines that drive the alignment process given
 * a collection of seed hits.  This is generally done in a few stages:
 * extendSeeds visits the set of seed-hit BW elements in some order;
 * for each element visited it resolves its reference offset; once the
 * reference offset is known, bounds for a Smith-Waterman subproblem
 * are established; if these bounds are distinct from the bounds we've
 * already tried, we solve the Smith-Waterman subproblem and report the
 * hit; if the AlnSinkWrap indicates that we can stop, we return,
 * otherwise we continue on to the next BW element.
 *
 * 
 */

#include <iostream>
#include "aligner_cache.h"
#include "aligner_sw_driver.h"
#include "pe.h"

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
 * Given a collection of SeedHits for a single read, extend seed
 * alignments into full alignments.  Where possible, try to avoid
 * redundant offset lookups and Smith-Watermans wherever possible.
 * Optionally report alignments to a AlnSinkWrap object as they
 * are discovered.
 *
 * If 'reportImmediately' is true, returns true iff a call to
 * msink->report() returned true (indicating that the reporting
 * policy is satisfied and we can stop).  Otherwise, returns false.
 */
bool SwDriver::extendSeeds(
	const Read& rd,              // read to align
	bool color,                  // true -> read is colorspace
	SeedResults& sh,             // seed hits to extend into full alignments
	const Ebwt& ebwt,            // BWT
	const BitPairReference& ref, // Reference strings
	GroupWalk& gw,               // group walk left
	SwAligner& swa,              // Smith-Waterman aligner
	const SwParams& pa,          // parars1_meters for Smith-Waterman aligner
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
	SwMetrics& swm,              // Smith-Waterman metrics
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
	int maxGaps = max(readGaps, refGaps);
	const size_t rdlen = rd.length();
	coords1fw_.clear();
	coords1rc_.clear();

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
			ASSERT_ONLY(uint32_t off = wr.toff);
			assert_neq(0xffffffff, off);
			uint32_t tidx = 0, toff = 0, tlen = 0;
			ebwt.joinedToTextOff(
				wr.elt.len,
				wr.toff,
				tidx,
				toff,
				tlen);
			if(tidx == 0xffffffff) {
				// The seed hit straddled a reference boundary
				continue;
			}
			// Hit's a hit at a known offset
			int64_t refi = (int64_t)toff;
			if(rdoff > refi) {
				// Traveling from the seed to the upstream edge of the
				// read takes us off the end of the reference
				// continue;
			}
			refi -= rdoff; // refi might be negative now
			if((int64_t)maxGaps > refi) {
				// The maximum number of gaps that we can add on the
				// left flank takes us off the end of the reference
				// continue;
			}
			refi -= maxGaps; // refi might be negative
			int64_t reff = toff + rdlen + maxGaps;
			if(reff > tlen) {
				// Traveling from the seed to the downstream edge of
				// the read *and* adding gaps to the right flank takes
				// us off the end of the reference
				// continue;
			}
			// Now we know what reference positions we're concerned
			// with.  If we have already investigated a window that
			// includes all of these reference characters, we can skip
			// this attempt.
			int64_t refiInner = max<int64_t>(refi, 0);
			Coord c(tidx, refiInner, fw);
			if((fw  && !coords1fw_.insert(c)) ||
			   (!fw && !coords1rc_.insert(c)))
			{
				// Already tried to find an alignment at these
				// coordinates
				swm.rshit++;
				continue;
			}
			res_.reset();
			assert(res_.empty());
			// Given the boundaries defined by refi and reff, align the
			// read to this reference stretch and, if a valid alignment
			// is found, store the result in res_.
			bool found = swa.alignToBitPairReference(
				rd,
				color,
				0,
				rdlen,
				fw,
				tidx,
				refi,
				reff,
				ref,
				tlen,
				pa,
				pen,
				penceil,
				res_,
				swCounterSinks,
				swActionSinks);
			assert_neq(0xffffffff, tidx);
			assert( found ||  res_.empty());
			assert(!found || !res_.empty());
			if(found) {
				res_.swsucc++;
			} else {
				res_.swfail++;
			}
			swm.update(res_);
			
			if(found) {
				// Annotate the AlnRes object with some key parameters
				// that were used to obtain the alignment.
				res_.alres.setParams(
					seedmms,   // # mismatches allowed in seed
					seedlen,   // length of seed
					seedival,  // interval between seeds
					penceil);  // maximum penalty for valid alignment
			}
			
			// If the user specified that alignments overhanging the
			// end of the reference be excluded, exclude them here.
			// TODO: exclude them in the traceback process so that we
			// don't prefer an invalid overhanging alignment over a
			// valid one with the same score during traceback.
			if(found && !gReportOverhangs) {
				if(!res_.alres.within(tidx, 0, fw, tlen)) {
					found = false;
				}
			}
			
			// Report this hit to an AlnSink
			if(found && reportImmediately) {
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
 * Given a read, perform full Smith-Waterman against the entire
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
	SwAligner& swa,              // Smith-Waterman aligner
	const SwParams& pa,          // parameters for Smith-Waterman aligner
	const Penalties& pen,        // penalties for edits
	int penceil,                 // maximum penalty allowed
	RandomSource& rnd,           // pseudo-random source
	SwMetrics& swm,              // Smith-Waterman metrics
	AlnSinkWrap* msink,        // HitSink for multiseed-style aligner
	bool reportImmediately,      // whether to report hits immediately to msink
	EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
	EList<SwActionSink*>* swActionSinks)   // send action-list updates to these
{
	assert(!reportImmediately || msink != NULL);
	return false;
}

/**
 * Given a collection of SeedHits for both mates in a read pair, extend
 * seed alignments into full alignments and then look for the opposite
 * mate using dynamic programming.  Where possible, try to avoid
 * redundant offset lookups.  Optionally report alignments to a
 * AlnSinkWrap object as they are discovered.
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
 * If 'reportImmediately' is true, returns true iff a call to
 * msink->report() returned true (indicating that the reporting
 * policy is satisfied and we can stop).  Otherwise, returns false.
 *
 */
bool SwDriver::extendSeedsPaired(
	const Read& rd1,             // mate1 to align
	const Read& rd2,             // mate2 to align
	bool color,                  // true -> reads are colorspace
	SeedResults& sh1,            // seed hits for mate1
	SeedResults& sh2,            // seed hits for mate2
	const Ebwt& ebwt,            // BWT
	const BitPairReference& ref, // Reference strings
	GroupWalk& gw,               // group walk left
	SwAligner& swa,              // Smith-Waterman aligner
	const SwParams& pa,          // parameters for Smith-Waterman aligner
	const Penalties& pen,        // penalties for edits
	const PairedEndPolicy& pepol,// paired-end policy
	int seedmms,                 // # mismatches allowed in seed
	int seedlen,                 // length of seed
	int seedival,                // interval between seeds
	int penceil1,                // maximum penalty allowed for rd1
	int penceil2,                // maximum penalty allowed for rd2
	uint32_t maxposs,            // stop after examining this many positions (offset+orientation combos)
	uint32_t maxrows,            // stop examining a position after this many offsets are reported
	AlignmentCacheIface& sc,     // alignment cache for seed hits
	RandomSource& rnd,           // pseudo-random source
	WalkMetrics& wlm,            // group walk left metrics
	SwMetrics& swm,              // Smith-Waterman metrics
	ReportingMetrics& rpm,       // reporting metrics
	AlnSinkWrap* msink,        // AlnSink wrapper for multiseed-style aligner
	bool swMateImmediately,      // whether to look for mate immediately
	bool reportImmediately,      // whether to report hits immediately to msink
	EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
	EList<SwActionSink*>* swActionSinks)   // send action-list updates to these
{
	assert(!reportImmediately || msink != NULL);
	assert(!reportImmediately || msink->empty());
	assert(!reportImmediately || !msink->maxed());

	// Calculate the largest possible number of read and reference gaps
	// given 'penceil' and 'pen'
	int readGaps1 = pen.maxReadGaps(penceil1);
	int refGaps1 = pen.maxRefGaps(penceil1);
	int maxGaps1 = max(readGaps1, refGaps1);
	int readGaps2 = pen.maxReadGaps(penceil2);
	int refGaps2 = pen.maxRefGaps(penceil2);
	int maxGaps2 = max(readGaps2, refGaps2);
	const size_t rd1len = rd1.length();
	const size_t rd2len = rd2.length();
	coords1fw_.clear();
	coords1rc_.clear();
	coords2fw_.clear();
	coords2rc_.clear();

	// Iterate twice through levels seed hits from the lowest ranked
	// level to the highest ranked.  On the first iteration, look for
	// entries for which the offset is already known and try SWs.  On
	// the second iteration, resolve entries for which the offset is
	// unknown and try SWs.
	const size_t nonz1 = sh1.nonzeroOffsets();
	const size_t poss1 = min<size_t>(nonz1, maxposs);
	const size_t nonz2 = sh2.nonzeroOffsets();
	const size_t poss2 = min<size_t>(nonz2, maxposs);
	for(size_t i = 0; i < poss1 + poss2; i++) {
		bool do2 = i >= poss1;
		bool fw = true;
		uint32_t offidx = 0, rdoff = 0, seedlen;
		QVal qv;
		size_t rdlen, ordlen;
		const Read *rd, *ord;
		ESet<Coord> *coordsfw, *coordsrc;
		int penceil, openceil;
		int maxGaps, omaxGaps;
		if(!do2) {
			qv = sh1.hitsByRank(i,       offidx, rdoff, fw, seedlen);
			rdlen = rd1len;
			ordlen = rd2len;
			rd = &rd1;
			ord = &rd2;
			coordsfw = &coords1fw_;
			coordsrc = &coords1rc_;
			penceil = penceil1;
			openceil = penceil2;
			maxGaps = maxGaps1;
			omaxGaps = maxGaps2;
		} else {
			qv = sh2.hitsByRank(i-poss1, offidx, rdoff, fw, seedlen);
			rdlen = rd2len;
			ordlen = rd1len;
			rd = &rd2;
			ord = &rd1;
			coordsfw = &coords2fw_;
			coordsrc = &coords2rc_;
			penceil = penceil2;
			openceil = penceil1;
			maxGaps = maxGaps2;
			omaxGaps = maxGaps1;
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
			if(tidx == 0xffffffff) {
				// The seed hit straddled a reference boundary
				continue;
			}
			// Hit's a hit at a known offset
			int64_t refi = (int64_t)toff;
			if(rdoff > refi) {
				// Traveling from the seed to the upstream edge of the
				// read takes us off the end of the reference
				// continue;
			}
			refi -= rdoff; // refi might be negative now
			if((int64_t)maxGaps > refi) {
				// The maximum number of gaps that we can add on the
				// left flank takes us off the end of the reference
				// continue;
			}
			refi -= maxGaps; // refi might be negative
			int64_t reff = toff + rdlen + maxGaps;
			if(reff > tlen) {
				// Traveling from the seed to the downstream edge of
				// the read *and* adding gaps to the right flank takes
				// us off the end of the reference
				// continue;
			}
			// Now we know what reference positions we're concerned
			// with.  If we have already investigated a window that
			// includes all of these reference characters, we can skip
			// this attempt.
			int64_t refiInner = max<int64_t>(refi, 0);
			Coord c(tidx, refiInner, fw);
			if((fw  && !coordsfw->insert(c)) ||
			   (!fw && !coordsrc->insert(c)))
			{
				// Already tried to find an alignment at these
				// coordinates
				swm.rshit++;
				continue;
			}
			res_.reset();
			assert(res_.empty());
			// Given the boundaries defined by refi and reff, align the
			// read to this reference stretch and, if a valid alignment
			// is found, store the result in res_.
			bool foundAnchor = swa.alignToBitPairReference(
				*rd,
				color,
				0,
				rdlen,
				fw,
				tidx,
				refi,
				reff,
				ref,
				tlen,
				pa,
				pen,
				penceil,
				res_,
				swCounterSinks,
				swActionSinks);
			assert_neq(0xffffffff, tidx);
			assert( foundAnchor ||  res_.empty());
			assert(!foundAnchor || !res_.empty());
			if(foundAnchor) {
				res_.swsucc++;
			} else {
				res_.swfail++;
			}
			swm.update(res_);
			
			if(foundAnchor) {
				// Annotate the AlnRes object with some key parameters
				// that were used to obtain the alignment.
				res_.alres.setParams(
					seedmms,   // # mismatches allowed in seed
					seedlen,   // length of seed
					seedival,  // interval between seeds
					penceil);  // maximum penalty for valid alignment
			}
			
			// If the user specified that alignments overhanging the
			// end of the reference be excluded, exclude them here.
			// TODO: exclude them in the traceback process so that we
			// don't prefer an invalid overhanging alignment over a
			// valid one with the same score during traceback.
			if(foundAnchor && !gReportOverhangs) {
				if(!res_.alres.within(tidx, 0, fw, tlen)) {
					foundAnchor = false;
				}
			}
			
			bool foundMate = false;
			if(foundAnchor && swMateImmediately) {
				bool oleft = false, ofw = false;
				int64_t oll = 0, olr = 0, orl = 0, orr = 0;
				foundMate = pepol.otherMate(
					!do2,
					fw,
					res_.alres.refoff(),
					tlen,
					rd1.length(),
					rd2.length(),
					oleft,
					oll,
					olr,
					orl,
					orr,
					ofw);
#if 0
				if(foundMate) {
					ores_.reset();
					assert(ores_.empty());
					foundMate = swa.alignToBitPairReference(
						*ord,
						color,
						0,
						ordlen,
						ofw,
						tidx,
						oll,
						olr,
						orl,
						orr,
						ref,
						tlen,
						pa,
						pen,
						openceil,
						ores_,
						swCounterSinks,
						swActionSinks);
					
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
				}
#endif
			}

			if(reportImmediately) {
				if(foundMate) {
					// Report pair to the AlnSinkWrap
					assert(msink != NULL);
					assert(res_.repOk());
					assert(ores_.repOk());
					// Check that alignment accurately reflects the
					// reference characters aligned to
					assert(res_.alres.matchesRef(*rd, ref));
					assert(ores_.alres.matchesRef(*ord, ref));
					// Report an unpaired alignment
					assert(!msink->maxed());
					if(msink->report(0, &res_.alres, &ores_.alres)) {
						// Short-circuited because a limit, e.g. -k, -m or
						// -M, was exceeded
						return true;
					}
				} else if(foundAnchor) {
					// Report unpaired hit for 
					assert(msink != NULL);
					assert(res_.repOk());
					// Check that alignment accurately reflects the
					// reference characters aligned to
					assert(res_.alres.matchesRef(*rd, ref));
					// Report an unpaired alignment
					assert(!msink->maxed());
					if(msink->report(0, &res_.alres, NULL)) {
						// Should not short-circuit on an unpaired
						// alignment when more paired alignments are
						// possible
						cerr << "Should not get true return value from report() for unpaired mate alignment" << endl;
						throw 1;
					}
				}
			}
			
			// At this point we know that we aren't bailing, and will continue to resolve seed hits.  

		} // while(!gw.done())
	}
	return false;
}
