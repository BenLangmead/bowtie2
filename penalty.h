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
		int mmcType,    // how to penalize mismatches
	    int mmc,        // constant if mm pelanty is a constant
		int sn,         // penalty for nucleotide mismatch in decoded colorspace als
		float ncConst,  //
		float ncLinear, //
	    int nType,      // how to penalize Ns in the read
	    int n,          // constant if N pelanty is a constant
		bool ncat,      // whether to concatenate mates before N filtering
	    int rdOpen,     // cost of opening a gap in the read
	    int rfOpen,     // cost of opening a gap in the reference
	    int rdExConst,  // constant cost of extending a gap in the read
	    int rfExConst,  // constant cost of extending a gap in the reference
	    int rdExLinear, // coefficient of linear term for cost of gap extension in read
	    int rfExLinear) // coefficient of linear term for cost of gap extension in ref
	{
		mmcostType   = mmcType;
		mmcost       = mmc;
		snp          = sn;
		nCeilConst   = ncConst;
		nCeilLinear  = ncLinear;
		npenType     = nType;
		npen         = n;
		ncatpair     = ncat;
		readOpen     = rdOpen;
		refOpen      = rfOpen;
		readExConst  = rdExConst;
		refExConst   = rfExConst;
		readExLinear = rdExLinear;
		refExLinear  = rfExLinear;
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
		assert_gt(snp, 0);
		assert_gt(readOpen, 0);
		assert_gt(refOpen, 0);
		assert_gt(readExConst, 0);
		assert_gt(refExConst, 0);
		assert_geq(readExLinear, 0);
		assert_geq(refExLinear, 0);
		assert_geq(readOpen, readExConst + readExLinear);
		assert_geq(refOpen,  refExConst  + refExLinear);
		return true;
	}
	
	/**
	 * Init an array that maps quality to penalty according to 'type'
	 * and 'cons'
	 */
	void initPens(int *pens, int type, int cons) {
		qualsMatter_ = true;
		if(type == COST_MODEL_ROUNDED_QUAL) {
			for(int i = 0; i < 64; i++) {
				pens[i] = qualRounds[i];
			}
		} else if(type == COST_MODEL_QUAL) {
			for(int i = 0; i < 64; i++) {
				pens[i] = i;
			}
		} else if(type == COST_MODEL_CONSTANT) {
			qualsMatter_ = false;
			assert_gt(mmcost, 0);
			for(int i = 0; i < 64; i++) {
				pens[i] = cons;
			}
		} else {
			throw 1;
		}
		for(int i = 0; i < 64; i++) {
			assert_geq(pens[i], 0);
		}
	}

	/**
	 * Return the penalty incurred by a mismatch at an alignment column
	 * with read character 'rdc' reference mask 'refm' and quality 'q'.
	 *
	 * qs should be clamped to 63 on the high end before this query.
	 */
	int mm(int rdc, int refm, int q) const {
		assert_range(0, 63, q);
		return (rdc > 3 || refm > 15) ? npens[q] : mmpens[q];
	}

	/**
	 * Return the penalty incurred by a mismatch at an alignment column
	 * with read character 'rdc' and quality 'q'.  We assume the
	 * reference character is non-N.
	 */
	int mm(int rdc, int q) const {
		assert_range(0, 63, q);
		return (rdc > 3) ? npens[q] : mmpens[q];
	}
	
	/**
	 * Return the marginal penalty incurred by a mismatch at a read
	 * position with quality 'q'.
	 */
	int mm(int q) const {
		assert_geq(q, 0);
		return q < 64 ? mmpens[q] : mmpens[63];
	}

	/**
	 * Return true iff the penalities are such that two reads with the
	 * same sequence but different qualities might yield different
	 * alignments.
	 */
	bool qualitiesMatter() const {
		return qualsMatter_;
	}
	
	/**
	 * Return the marginal penalty incurred by an N mismatch at a read
	 * position with quality 'q'.
	 */
	int n(int q) const {
		assert_geq(q, 0);
		return q < 64 ? npens[q] : npens[63];
	}
	
	/**
	 * Return the marginal penalty incurred by a gap in the read,
	 * given that this is the 'ext'th extension of the gap (0 = open,
	 * 1 = first, etc).
	 */
	int ins(int ext) const {
		assert_geq(ext, 0);
		if(ext == 0) return readOpen;
		return std::max<int>(1, readExConst + (ext-1)*readExLinear);
	}

	/**
	 * Return the marginal penalty incurred by a gap in the reference,
	 * given that this is the 'ext'th extension of the gap (0 = open,
	 * 1 = first, etc).
	 */
	int del(int ext) const {
		assert_geq(ext, 0);
		if(ext == 0) return refOpen;
		return std::max<int>(1, refExConst + (ext-1)*refExLinear);
	}

	/**
	 * Calculate the maximum possible number of read gaps that could
	 * occur in a valid alignment given the penalties assoicated with
	 * read gaps.
	 */
	int maxReadGaps(int penceil) const {
		int max = 0;
		penceil -= readOpen;
		if(penceil >= 0) {
			max++;
		}
		int lin = 1;
		while(true) {
			int pen = (readExConst + readExLinear * lin++);
			penceil -= pen;
			if(penceil < 0) break;
			max++;
		}
		return max;
	}

	/**
	 * Calculate the maximum possible number of reference gaps that
	 * could occur in a valid alignment given the penalties assoicated
	 * with reference gaps.
	 */
	int maxRefGaps(int penceil) const {
		int max = 0;
		penceil -= refOpen;
		if(penceil >= 0) {
			max++;
		}
		int lin = 1;
		while(true) {
			int pen = (refExConst + refExLinear * lin++);
			penceil -= pen;
			if(penceil < 0) break;
			max++;
		}
		return max;
	}

	/**
	 * Given a read sequence, return true iff the read passes the N filter.
	 * The N filter rejects reads with more than the number of Ns calculated by
	 * taking nCeilConst + nCeilLinear * read length.
	 */
	bool nFilter(const BTDnaString& rd) const {
		float ns = 0.0f;
		size_t rdlen = rd.length();
		float maxns = nCeilConst + nCeilLinear * rdlen;
		if(maxns < 0.0f) maxns = 0.0f;
		assert_gt(rd.length(), 0);
		for(size_t i = 0; i < rdlen; i++) {
			if(rd[i] == 4) {
				ns += 1.0f;
				if(ns > maxns) return false; // doesn't pass
			}
		}
		return true; // passes
	}

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
		const
	{
		// Both fail to pass by default
		filt1 = filt2 = false;
		if(rd1 != NULL && rd2 != NULL && ncatpair) {
			size_t rdlen1 = rd1->length();
			size_t rdlen2 = rd2->length();
			float maxns = nCeilConst + nCeilLinear * (rdlen1 + rdlen2);
			size_t ns = 0;
			for(size_t i = 0; i < rdlen1; i++) {
				if((*rd1)[i] == 4) ns++;
				if((float)ns > maxns) {
					// doesn't pass
					return;
				}
			}
			for(size_t i = 0; i < rdlen2; i++) {
				if((*rd2)[i] == 4) ns++;
				if((float)ns > maxns) {
					// doesn't pass
					return;
				}
			}
			// Both pass
			filt1 = filt2 = true;
		} else {
			if(rd1 != NULL) filt1 = nFilter(*rd1);
			if(rd2 != NULL) filt2 = nFilter(*rd2);
		}
	}
	
	/**
	 * Given a sequence length, return the number of Ns that are
	 * allowed in the sequence.
	 */
	size_t nCeil(size_t rdlen) const {
		return (size_t)(nCeilConst + nCeilLinear * rdlen);
	}

	int   mmcostType;   // based on qual? rounded? just a constant?
	int   mmcost;       // if mmcosttype=constant, this is the const
	int   snp;          // penalty associated with a nucleotide mismatch in a decoded colorspace alignment
	float nCeilConst;   // max # Ns involved in alignment, constant coeff
	float nCeilLinear;  // max # Ns involved in alignment, linear coeff
	int   npenType;     // N: based on qual? rounded? just a constant?
	int   npen;         // N: if mmcosttype=constant, this is the const
	bool  ncatpair;     // true -> do N filtering on concated pair
	int   readOpen;     // read gap open penalty
	int   refOpen;      // reference gap open penalty
	int   readExConst;  // constant term coeffecient in extend cost
	int   refExConst;   // constant term coeffecient in extend cost
	int   readExLinear; // linear term coeffecient in extend cost
	int   refExLinear;  // linear term coeffecient in extend cost
	
	int   mmpens[64];   // map from qualities to mm penalty
	int   npens[64];    // map from N qualities to penalty
	
	static Penalties bwaLike() {
		return Penalties(
			COST_MODEL_CONSTANT, // how to penalize mismatches
			3,                   // constant if mm pelanty is a constant
			30,                  // penalty for nuc mm when decoding colors
			2.0f,                // constant coeff for max Ns
			0.1f,                // linear coeff for max Ns
			COST_MODEL_CONSTANT, // how to penalize Ns in the read
			3,                   // constant if N pelanty is a constant
			false,               // concatenate mates before N filtering?
			11,                  // cost of opening a gap in the read
			11,                  // cost of opening a gap in the reference
			4,                   // constant coeff for extending gap in read
			4,                   // constant coeff for extending gap in ref
			0,                   // linear coeff for extending gap in read
			0);                  // linear coeff for extending gap in ref
	}

	static Penalties naLike() {
		return Penalties(
			COST_MODEL_ROUNDED_QUAL, // how to penalize mismatches
			0,                       // constant if mm pelanty is a constant
			30,                      // penalty for nuc mm when decoding colors
			2.0f,                    // constant coeff for max Ns
			0.1f,                    // linear coeff for max Ns
			COST_MODEL_ROUNDED_QUAL, // how to penalize Ns in the read
			0,                       // constant if N pelanty is a constant
			false,                   // concatenate mates before N filtering?
			40,                      // cost of opening a gap in the read
			40,                      // cost of opening a gap in the reference
			15,                      // constant coeff for extending gap in read
			15,                      // constant coeff for extending gap in ref
			0,                       // linear coeff for extending gap in read
			0);                      // linear coeff for extending gap in ref
	}

protected:

	bool qualsMatter_;
};

#endif /*PENALTY_H_*/
