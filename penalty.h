/*
 *  penalty.h
 */

#ifndef PENALTY_H_
#define PENALTY_H_

#include "qual.h"

enum {
	COST_MODEL_ROUNDED_QUAL = 1,
	COST_MODEL_QUAL,
	COST_MODEL_CONSTANT
};

/**
 * How to penalize various types of sequence dissimilarity.
 */
class Penalties {

public:

	Penalties(
		int   mat,        // reward for a match
		int   mmcType,    // how to penalize mismatches
	    int   mmc,        // constant if mm pelanty is a constant
		int   sn,         // penalty for nuc mismatch in decoded colorspace als
		float ncConst,    // max # ref Ns allowed in alignment; const coeff
		float ncLinear,   // max # ref Ns allowed in alignment; linear coeff
	    int   nType,      // how to penalize Ns in the read
	    int   n,          // constant if N pelanty is a constant
		bool  ncat,       // whether to concatenate mates before N filtering
	    int   rdGpConst,  // constant coeff for cost of gap in the read
	    int   rfGpConst,  // constant coeff for cost of gap in the ref
	    int   rdGpLinear, // coeff of linear term for cost of gap in read
	    int   rfGpLinear) // coeff of linear term for cost of gap in ref
	{
		match        = mat;
		mmcostType   = mmcType;
		mmcost       = mmc;
		snp          = sn;
		nCeilConst   = ncConst;
		nCeilLinear  = ncLinear;
		npenType     = nType;
		npen         = n;
		ncatpair     = ncat;
		rdGapConst   = rdGpConst;
		rfGapConst   = rfGpConst;
		rdGapLinear  = rdGpLinear;
		rfGapLinear  = rfGpLinear;
		qualsMatter_ = mmcostType != COST_MODEL_CONSTANT;
		initPens(mmpens, mmcostType, mmcost);
		initPens(npens, npenType, npen);
		assert(repOk());
	}
	
	/**
	 * Set the N penalty.
	 */
	void setNPen(int nType, int n) {
		npenType     = nType;
		npen         = n;
		initPens(npens, npenType, npen);
	}
	
	/**
	 * Check that penalties are internally consistent.
	 */
	bool repOk() const {
		assert_geq(match, 0);
		assert_gt(snp, 0);
		assert_gt(rdGapConst, 0);
		assert_gt(rdGapLinear, 0);
		assert_gt(rfGapConst, 0);
		assert_gt(rfGapLinear, 0);
		return true;
	}
	
	/**
	 * Init an array that maps quality to penalty according to 'type'
	 * and 'cons'
	 */
	void initPens(
		int *pens,  // array to fill
		int type,   // penalty type; qual | rounded qual | constant
		int cons);  // constant for when penalty type is constant

	/**
	 * Return the penalty incurred by a mismatch at an alignment column
	 * with read character 'rdc' reference mask 'refm' and quality 'q'.
	 *
	 * qs should be clamped to 63 on the high end before this query.
	 */
	int mm(int rdc, int refm, int q) const {
		assert_range(0, 255, q);
		return (rdc > 3 || refm > 15) ? npens[q] : mmpens[q];
	}

	/**
	 * Return the penalty incurred by a mismatch at an alignment column
	 * with read character 'rdc' and quality 'q'.  We assume the
	 * reference character is non-N.
	 */
	int mm(int rdc, int q) const {
		assert_range(0, 255, q);
		return (rdc > 3) ? npens[q] : mmpens[q];
	}
	
	/**
	 * Return the marginal penalty incurred by a mismatch at a read
	 * position with quality 'q'.
	 */
	int mm(int q) const {
		assert_geq(q, 0);
		return q < 255 ? mmpens[q] : mmpens[255];
	}

	/**
	 * Return true iff the penalities are such that two reads with the
	 * same sequence but different qualities might yield different
	 * alignments.
	 */
	inline bool qualitiesMatter() const { return qualsMatter_; }
	
	/**
	 * Return the marginal penalty incurred by an N mismatch at a read
	 * position with quality 'q'.
	 */
	int n(int q) const {
		assert_geq(q, 0);
		return q < 255 ? npens[q] : npens[255];
	}

	
	/**
	 * Return the marginal penalty incurred by a gap in the read,
	 * given that this is the 'ext'th extension of the gap (0 = open,
	 * 1 = first, etc).
	 */
	int ins(int ext) const {
		assert_geq(ext, 0);
		if(ext == 0) return readGapOpen();
		return readGapExtend();
	}

	/**
	 * Return the marginal penalty incurred by a gap in the reference,
	 * given that this is the 'ext'th extension of the gap (0 = open,
	 * 1 = first, etc).
	 */
	int del(int ext) const {
		assert_geq(ext, 0);
		if(ext == 0) return refGapOpen();
		return refGapExtend();
	}

	/**
	 * Given the score floor for valid alignments and the length of the read,
	 * calculate the maximum possible number of read gaps that could occur in a
	 * valid alignment.
	 */
	int maxReadGaps(
		int64_t minsc,
		size_t rdlen) const;

	/**
	 * Given the score floor for valid alignments and the length of the read,
	 * calculate the maximum possible number of reference gaps that could occur
	 * in a valid alignment.
	 */
	int maxRefGaps(
		int64_t minsc,
		size_t rdlen) const;

	/**
	 * Given a read sequence, return true iff the read passes the N filter.
	 * The N filter rejects reads with more than the number of Ns calculated by
	 * taking nCeilConst + nCeilLinear * read length.
	 */
	bool nFilter(const BTDnaString& rd) const;

	/**
	 * Given a read sequence, return true iff the read passes the N filter.
	 * The N filter rejects reads with more than the number of Ns calculated by
	 * taking nCeilConst + nCeilLinear * read length.
	 *
	 * For paired-end reads, there is a	question of how to apply the filter.
	 * The filter could be applied to both mates separately, which might then
	 * prevent paired-end alignment.  Or the filter could be applied to the
	 * reads as though they're concatenated together.  The latter approach has
	 * pros and cons.  The pro is that we can use paired-end information to
	 * recover alignments for mates that would not have passed the N filter on
	 * their own.  The con is that we might not want to do that, since the
	 * non-N portion of the bad mate might contain particularly unreliable
	 * information.
	 */
	void nFilterPair(
		const BTDnaString* rd1, // mate 1
		const BTDnaString* rd2, // mate 2
		bool& filt1,            // true -> mate 1 rejected by filter
		bool& filt2)            // true -> mate 2 rejected by filter
		const;
	
	/**
	 * Given a sequence length, return the number of Ns that are
	 * allowed in the sequence.
	 */
	inline size_t nCeil(size_t rdlen) const {
		return (size_t)(nCeilConst + nCeilLinear * rdlen);
	}
	
	/**
	 * The penalty associated with opening a new read gap.
	 */
	inline int readGapOpen() const { 
		return rdGapConst + rdGapLinear;
	}

	/**
	 * The penalty associated with opening a new ref gap.
	 */
	inline int refGapOpen() const { 
		return rfGapConst + rfGapLinear;
	}

	/**
	 * The penalty associated with extending a read gap by one character.
	 */
	inline int readGapExtend() const { 
		return rdGapLinear;
	}

	/**
	 * The penalty associated with extending a ref gap by one character.
	 */
	inline int refGapExtend() const { 
		return rfGapLinear;
	}

	int   match;        // reward for a match
	int   mmcostType;   // based on qual? rounded? just a constant?
	int   mmcost;       // if mmcosttype=constant, this is the const
	int   snp;          // penalty for nuc mismatch in decoded colorspace als
	float nCeilConst;   // max # Ns involved in alignment, constant coeff
	float nCeilLinear;  // max # Ns involved in alignment, linear coeff
	int   npenType;     // N: based on qual? rounded? just a constant?
	int   npen;         // N: if mmcosttype=constant, this is the const
	bool  ncatpair;     // true -> do N filtering on concated pair
	int   rdGapConst;   // constant term coeffecient in extend cost
	int   rfGapConst;   // constant term coeffecient in extend cost
	int   rdGapLinear;  // linear term coeffecient in extend cost
	int   rfGapLinear;  // linear term coeffecient in extend cost
	int   mmpens[256];  // map from qualities to mm penalty
	int   npens[256];   // map from N qualities to penalty
	
	static Penalties bwaLike() {
		return Penalties(
			1,                       // reward for a match
			COST_MODEL_CONSTANT,     // how to penalize mismatches
			3,                       // constant if mm pelanty is a constant
			3,                       // penalty for nuc mm when decoding colors
			2.0f,                    // constant coeff for max Ns
			0.1f,                    // linear coeff for max Ns
			COST_MODEL_CONSTANT,     // how to penalize Ns in the read
			3,                       // constant if N pelanty is a constant
			false,                   // concatenate mates before N filtering?
			11,                      // constant coeff for gap in read
			11,                      // constant coeff for gap in ref
			4,                       // linear coeff for gap in read
			4);                      // linear coeff for gap in ref
	}

	static Penalties naLike() {
		return Penalties(
			0,                       // reward for a match
			COST_MODEL_ROUNDED_QUAL, // how to penalize mismatches
			0,                       // constant if mm pelanty is a constant
			30,                      // penalty for nuc mm when decoding colors
			2.0f,                    // constant coeff for max Ns
			0.1f,                    // linear coeff for max Ns
			COST_MODEL_ROUNDED_QUAL, // how to penalize Ns in the read
			0,                       // constant if N pelanty is a constant
			false,                   // concatenate mates before N filtering?
			40,                      // constant coeff for gap in read
			40,                      // constant coeff for gap in ref
			15,                      // linear coeff for gap in read
			15);                     // linear coeff for gap in ref
	}

protected:

	bool qualsMatter_;
};

#endif /*PENALTY_H_*/
