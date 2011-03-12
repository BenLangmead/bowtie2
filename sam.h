/*
 * sam.h
 *
 *  Created on: Sep 23, 2009
 *      Author: Ben Langmead
 */

#ifndef SAM_H_
#define SAM_H_

#include "refmap.h"
#include "annot.h"
#include "pat.h"
#include "random_source.h"
#include "ds.h"

class ReferenceMap;
class AnnotationMap;
class PatternSourcePerThread;

enum {
	SAM_FLAG_PAIRED = 1,
	SAM_FLAG_MAPPED_PAIRED = 2,
	SAM_FLAG_UNMAPPED = 4,
	SAM_FLAG_MATE_UNMAPPED = 8,
	SAM_FLAG_QUERY_STRAND = 16,
	SAM_FLAG_MATE_STRAND = 32,
	SAM_FLAG_FIRST_IN_PAIR = 64,
	SAM_FLAG_SECOND_IN_PAIR = 128,
	SAM_FLAG_NOT_PRIMARY = 256,
	SAM_FLAG_FAILS_CHECKS = 512,
	SAM_FLAG_DUPLICATE = 1024
};

/**
 * Sink that prints lines in SAM format:
 */
class SAMHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream VerboseHitSink (default)
	 */
	SAMHitSink(
		OutFileBuf *out,
		ReadSink *readSink,
		const EList<std::string>* refnames,
		int offBase,
		ReferenceMap *rmap,
		AnnotationMap *amap,
		bool fullRef,
		bool sampleMax,
		int defaultMapq) :
		HitSink(out, readSink, refnames),
		offBase_(offBase),
		defaultMapq_(defaultMapq),
		rmap_(rmap),
		amap_(amap),
		fullRef_(fullRef),
		sampleMax_(sampleMax)
	{ }

	/**
	 * Append a SAM alignment to the given output stream.
	 */
	static void append(ostream& ss,
	                   const Hit& h,
	                   int mapq,
	                   int xms,
	                   const EList<string>* refnames,
	                   ReferenceMap *rmap,
	                   AnnotationMap *amap,
	                   bool fullRef,
	                   int offBase);

	/**
	 * Append a SAM alignment for an aligned read to the given output
	 * stream.
	 */
	static void appendAligned(ostream& ss,
	                          const Hit& h,
	                          int mapq,
	                          int xms,
	                          const EList<string>* refnames,
	                          ReferenceMap *rmap,
	                          AnnotationMap *amap,
	                          bool fullRef,
	                          int offBase);

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(ostream& ss, const Hit& h) {
		SAMHitSink::append(ss, h, defaultMapq_, 0, refnames_, rmap_, amap_, fullRef_, offBase_);
	}

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(ostream& ss, const Hit& h, int mapq, int xms) {
		SAMHitSink::append(ss, h, mapq, xms, refnames_, rmap_, amap_, fullRef_, offBase_);
	}

	/**
	 * Write the SAM header lines.
	 */
	void appendHeaders(OutFileBuf& os,
	                   size_t numRefs,
	                   const EList<string>& refnames,
	                   bool nosq,
	                   ReferenceMap *rmap,
	                   const uint32_t* plen,
	                   bool fullRef,
	                   const char *cmdline,
	                   const char *rgline);

protected:

	/**
	 *
	 */
	void reportUnOrMax(PatternSourcePerThread& p,
	                   EList<Hit>* hs,
	                   bool un);

	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream.
	 */
	virtual void reportHit(const Hit& h) {
		reportHit(h, defaultMapq_, 0);
	}

	/**
	 * Report a SAM alignment with the given mapping quality and XM
	 * field.
	 */
	virtual void reportHit(const Hit& h, int mapq, int xms);

	/**
	 * Report a batch of SAM alignments (e.g. two mates that should be
	 * printed together) with the given mapping quality and XM field.
	 */
	virtual void reportHits(EList<Hit>& hs,
	                        size_t start,
	                        size_t end,
	                        int mapq,
	                        int xms);

	/**
	 * See sam.cpp
	 */
	virtual void reportMaxed(EList<Hit>& hs, PatternSourcePerThread& p);

	/**
	 * See sam.cpp
	 */
	virtual void reportUnaligned(PatternSourcePerThread& p) {
		reportUnOrMax(p, NULL, true);
	}

private:
	int  offBase_;        /// Add this to reference offsets before outputting.
	                      /// (An easy way to make things 1-based instead of
	                      /// 0-based)
	int  defaultMapq_;    /// Default mapping quality to report when one is
	                      /// not specified
	ReferenceMap *rmap_;  /// mapping to reference coordinate system.
	AnnotationMap *amap_; ///
	bool fullRef_;        /// print full reference name, not just up to whitespace
	bool sampleMax_;      /// user specified -M
};

#endif /* SAM_H_ */
