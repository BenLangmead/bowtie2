/*
 *  aligner_sw_driver.h
 */

#ifndef ALIGNER_SW_DRIVER_H_
#define ALIGNER_SW_DRIVER_H_

#include <utility>
#include "ds.h"
#include "aligner_seed.h"
#include "aligner_sw.h"
#include "aligner_cache.h"
#include "reference.h"
#include "group_walk.h"
#include "ebwt.h"
#include "mem_ids.h"
#include "aln_sink.h"
#include "pe.h"

class SwDriver {

public:

	SwDriver() :
		satups_(SW_CAT),
		coords_(1024, SW_CAT),
		coordsst_(1024, SW_CAT),
		coordsen_(1024, SW_CAT),
		frags_(1024, SW_CAT),
		coords1seenup_(1024, SW_CAT),
		coords1seendn_(1024, SW_CAT),
		coords2seenup_(1024, SW_CAT),
		coords2seendn_(1024, SW_CAT),
		st_(1024, SW_CAT),
		en_(1024, SW_CAT),
		res_(), ores_()
	{ }

	/**
	 * Given a collection of SeedHits for a single read, extend seed alignments
	 * into full alignments.  Where possible, try to avoid redundant offset
	 * lookups and dynamic programming problems.  Optionally report alignments
	 * to a AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	bool extendSeeds(
		const Read& rd,              // read to align
		bool color,                  // true -> read is colorspace
		SeedResults& sh,             // seed hits to extend into full alignments
		const Ebwt& ebwt,            // BWT
		const BitPairReference& ref, // Reference strings
		GroupWalk& gw,               // group walk left
		SwAligner& swa,              // dynamic programming aligner
		const SwParams& pa,          // parameters for dynamic programming aligner
		const Penalties& pen,        // penalties for edits
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		int penceil,                 // maximum penalty allowed
		int nceil,                   // maximum # Ns permitted in ref portion
		uint32_t maxposs,            // stop after examining this many positions (offset+orientation combos)
		uint32_t maxrows,            // stop examining a position after this many offsets are reported
		AlignmentCacheIface& sc,     // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random source
		WalkMetrics& wlm,            // group walk left metrics
		SwMetrics& swm,              // dynamic programming metrics
		ReportingMetrics& rpm,       // reporting metrics
		AlnSinkWrap* mhs,            // HitSink for multiseed-style aligner
		bool reportImmediately,      // whether to report hits immediately to mhs
		EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
		EList<SwActionSink*>* swActionSinks);  // send action-list updates to these

	/**
	 * Given a read, perform full dynamic programming against the entire
	 * reference.  Optionally report alignments to a AlnSinkWrap object
	 * as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	bool sw(
		const Read& rd,              // read to align
		bool color,                  // true -> read is colorspace
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dynamic programming aligner
		const SwParams& pa,          // parameters for dynamic prog aligner
		const Penalties& pen,        // penalties for edits
		int penceil,                 // maximum penalty allowed
		RandomSource& rnd,           // pseudo-random source
		SwMetrics& swm,              // dynamic programming metrics
		AlnSinkWrap* mhs,            // HitSink for multiseed-style aligner
		bool reportImmediately,      // whether to report hits immediately to mhs
		EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
		EList<SwActionSink*>* swActionSinks);  // send action-list updates to these

	/**
	 * Given a collection of SeedHits for a read pair, extend seed
	 * alignments into full alignments and then look for the opposite
	 * mate using dynamic programming.  Where possible, try to avoid
	 * redundant offset lookups.  Optionally report alignments to a
	 * AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	bool extendSeedsPaired(
		const Read& rd,              // mate to align as anchor
		const Read& ord,             // mate to align as opposite
		bool anchor1,                // true iff anchor mate is mate1
		bool color,                  // true -> reads are colorspace
		SeedResults& sh,             // seed hits for anchor
		const Ebwt& ebwt,            // BWT
		const BitPairReference& ref, // Reference strings
		GroupWalk& gw,               // group walk left
		SwAligner& swa,              // dyn programming aligner for anchor
		SwAligner& swao,             // dyn programming aligner for opposite
		const SwParams& pa,          // parameters for dynamic programming aligner
		const Penalties& pen,        // penalties for edits
		const PairedEndPolicy& pepol,// paired-end policy
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		int penceil,                 // maximum penalty allowed for anchor
		int openceil,                // maximum penalty allowed for opposite
		int nceil,                   // max # Ns permitted in ref for anchor
		int onceil,                  // max # Ns permitted in ref for opposite
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
		EList<SwActionSink*>* swActionSinks);  // send action-list updates to these

	/**
	 * Prepare for a new read.  For paired-end reads, this means clearing state
	 * that would otherwise survive across calls to extendSeedsPaired.
	 */
	void nextRead(bool paired) {
		if(paired) {
			frags_.clear();
			coords1seenup_.clear();
			coords1seendn_.clear();
			coords2seenup_.clear();
			coords2seendn_.clear();
		}
	}

protected:

	/**
	 * Try to convert a BW offset to the forward reference offset of the
	 * leftmost position involved in the hit.  Return 0xffffffff if the
	 * index doesn't contain the information needed to convert.
	 */
	uint32_t bwtOffToOff(
		const Ebwt& ebwt,
		bool fw,
		uint32_t bwoff,
		uint32_t hitlen,
		uint32_t& tidx,
		uint32_t& toff,
		uint32_t& tlen);

	EList<SATuple> satups_;     // temporary holder for range lists
	ESet<Coord>    coords_;     // ref coords tried for mate 1 so far
	ESet<Coord>    coordsst_;   // upstream coord for mate 1 hits so far
	ESet<Coord>    coordsen_;   // downstream coord for mate 1 hits so far
	ESet<Interval> frags_;      // intervals for paired-end fragments
	ESet<Coord> coords1seenup_; // LHSs seen for mate1s - to avoid double-
	                            // reporting unpaired alignments
	ESet<Coord> coords1seendn_; // RHSs seen for mate1s - to avoid double-
	                            // reporting unpaired alignments
	ESet<Coord> coords2seenup_; // LHSs seen for mate2s - to avoid double-
	                            // reporting unpaired alignments
	ESet<Coord> coords2seendn_; // RHSs seen for mate2s - to avoid double-
	                            // reporting unpaired alignments
	EList<bool>    st_;         // temp holder for dyn prog starting mask
	EList<bool>    en_;         // temp holder for dyn prog ending mask
	SwResult       res_;        // temp holder for SW results
	SwResult       ores_;       // temp holder for SW results for opp mate
};

#endif /*ALIGNER_SW_DRIVER_H_*/
