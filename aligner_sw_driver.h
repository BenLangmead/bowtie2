/*
 *  aligner_sw_driver.h
 */

#ifndef ALIGNER_SW_DRIVER_H_
#define ALIGNER_SW_DRIVER_H_

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
		coords1fw_(1024, SW_CAT),
		coords1rc_(1024, SW_CAT),
		coords2fw_(1024, SW_CAT),
		coords2rc_(1024, SW_CAT),
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
		const Read& rd1,             // mate1 to align
		const Read& rd2,             // mate2 to align
		bool color,                  // true -> reads are colorspace
		SeedResults& sh1,            // seed hits for mate1
		SeedResults& sh2,            // seed hits for mate2
		const Ebwt& ebwt,            // BWT
		const BitPairReference& ref, // Reference strings
		GroupWalk& gw,               // group walk left
		SwAligner& swa,              // dynamic programming aligner
		const SwParams& pa,          // parameters for dynamic prog aligner
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
		SwMetrics& swm,              // dynamic programming metrics
		ReportingMetrics& rpm,       // reporting metrics
		AlnSinkWrap* msink,          // AlnSink wrapper for multiseed aligner
		bool swMateImmediately,      // whether to look for mate immediately
		bool reportImmediately,      // whether to report hits right to msink
		EList<SwCounterSink*>* swCounterSinks, // send counter updates to these
		EList<SwActionSink*>* swActionSinks);  // send action-list updates to these

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

	EList<SATuple> satups_;    // temporary holder for range lists
	ESet<Coord>    coords1fw_; // ref coords tried for fw mate 1 so far
	ESet<Coord>    coords1rc_; // ref coords tried for rc mate 1 so far
	ESet<Coord>    coords2fw_; // ref coords tried for fw mate 2 so far
	ESet<Coord>    coords2rc_; // ref coords tried for rc mate 2 so far
	EList<bool>    st_;        // temporary holder for dyn prog starting mask
	EList<bool>    en_;        // temporary holder for dyn prog ending mask
	SwResult       res_;       // temporary holder for SW results
	SwResult       ores_;      // temporary holder for SW results for opposite mate
};

#endif /*ALIGNER_SW_DRIVER_H_*/
