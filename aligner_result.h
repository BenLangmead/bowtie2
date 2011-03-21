/*
 * aligner_result.h
 */

#ifndef ALIGNER_RESULT_H_
#define ALIGNER_RESULT_H_

#include <utility>
#include "mem_ids.h"
#include "ref_coord.h"
#include "read.h"
#include "filebuf.h"
#include "ds.h"
#include "edit.h"

typedef int64_t TAlScore;

#define VALID_AL_SCORE(x) ((x).score_ > (x).scoreLim)

/**
 * A generic score object for an alignment.  Used for accounting during
 * SW and elsewhere.  Encapsulates the score, the number of N positions
 * and the number gaps in the alignment.
 *
 * The scale for 'score' is such that a perfect alignment score is 0
 * and a score with non-zero penalty is less than 0.  So differences
 * between scores work as expected, but interpreting an individual
 * score (larger is better) as a penalty (smaller is better) requires
 * taking the absolute value.
 */
class AlignmentScore {

public:

	/**
	 * Gapped scores are invalid until proven valid.
	 */
	AlignmentScore() : scoreLim(std::numeric_limits<TAlScore>::min()) {
		score_ = ns_ = gaps_ = 0;
		invalidate();
	}

	/**
	 * Return an invalid SwScore.
	 */
	static AlignmentScore INVALID() {
		AlignmentScore s;
		s.invalidate();
		return s;
	}

	/**
	 * Return true iff gapped score is valid (i.e., represents one or
	 * more paths leading to a valid partial alignment).
	 */
	//inline bool valid() const {
	//	return score_ > scoreLim;
	//}

	/**
	 * Make this score invalid (and therefore <= all other scores).
	 */
	inline void invalidate() {
		score_ = std::numeric_limits<TAlScore>::min();
	}
	
	/**
	 * Increment the number of gaps.  If currently invalid, this makes
	 * the score valid with gaps == 1.
	 */
	void incNs(int nceil) {
		if(++ns_ > nceil) {
			invalidate();
		}
	}

	/**
	 * Return true iff this score is > score o.
	 * Note: An "invalid" score is <= all other scores.
	 */
	bool operator>(const AlignmentScore& o) const {
		if(!VALID_AL_SCORE(o)) {
			if(!VALID_AL_SCORE(*this)) {
				// both invalid
				return false;
			} else {
				// I'm valid, other is invalid
				return true;
			}
		} else if(!VALID_AL_SCORE(*this)) {
			// I'm invalid, other is valid
			return false;
		}
		return score_ > o.score_;
	}

	/**
	 * Scores are equal iff they're bitwise equal.
	 */
	inline AlignmentScore& operator=(const AlignmentScore& o) {
		// Profiling shows many cache misses on following lines
		gaps_  = o.gaps_;
		ns_    = o.ns_;
		score_ = o.score_;
		return *this;
	}

	/**
	 * Scores are equal iff they're bitwise equal.
	 */
	inline bool operator==(const AlignmentScore& o) const {
		// Profiling shows cache misses on following line
		return VALID_AL_SCORE(*this) && VALID_AL_SCORE(o) && score_ == o.score_;
	}

	/**
	 * Return true iff this score is >= score o.
	 */
	bool operator>=(const AlignmentScore& o) const {
		return operator==(o) || operator>(o);
	}

	/**
	 * Return true iff this score is < score o.
	 */
	bool operator<(const AlignmentScore& o) const {
		return !operator>=(o);
	}

	/**
	 * Return true iff this score is <= score o.
	 */
	bool operator<=(const AlignmentScore& o) const {
		return !operator>(o);
	}

	/**
	 * Calculate difference between two SwScores.
	 */
	AlignmentScore operator-(const AlignmentScore& o) const {
		if(!VALID_AL_SCORE(*this)) return *this;
		AlignmentScore s; 
		s.gaps_ = gaps_ - o.gaps_;
		s.ns_ = ns_;
		s.score_ = score_ - o.score_;
		return s;
	}

	/**
	 * Calculate sum of two SwScores.
	 */
	AlignmentScore operator+(const AlignmentScore& o) const {
		if(!VALID_AL_SCORE(*this)) return *this;
		AlignmentScore s;
		s.gaps_ = gaps_ + o.gaps_;
		s.ns_ = ns_;
		s.score_ = score_ + o.score_;
		return s;
	}

	/**
	 * Add given SwScore into this one.
	 */
	AlignmentScore operator+=(const AlignmentScore& o) {
		if(VALID_AL_SCORE(*this)) {
			gaps_ += o.gaps_;
			score_ += o.score_;
		}
		return (*this);
	}

	/**
	 * Subtract given SwScore from this one.
	 */
	inline AlignmentScore operator-=(const AlignmentScore& o) {
		if(VALID_AL_SCORE(*this)) {
			gaps_ -= o.gaps_;
			score_ -= o.score_;
		}
		return (*this);
	}

	/**
	 * Calculate difference between two SwScores.
	 */
	inline AlignmentScore operator-(int o) const {
		if(!VALID_AL_SCORE(*this)) return *this;
		AlignmentScore s;
		s.gaps_ = gaps_;
		s.ns_ = ns_;
		s.score_ = score_ - o;
		return s;
	}
	
	TAlScore scoreLim;
	
	TAlScore score()   const { return  score_; }
	TAlScore penalty() const { return -score_; }
	TAlScore gaps()    const { return  gaps_;  }
	TAlScore ns()      const { return  ns_;    }

	// Score accumulated so far (penalties are subtracted starting at 0)
	TAlScore score_;
	
	// Ns accumulated so far.  An N opposite a non-gap counts as 1 N
	// (even if it's N-to-N)
	TAlScore ns_;
	
	// # gaps encountered so far, unless that number exceeds the
	// target, in which case the score becomes invalid and therefore <=
	// all other scores
	TAlScore gaps_;
};

// Forward declaration
class BitPairReference;

/**
 * Encapsulates an alignment result, including for colorspace
 * alignments.  The result comprises:
 *
 * 1. All the nucleotide edits for both mates ('ned').
 * 2. All the "edits" where an ambiguous reference char is resolved to
 *    an unambiguous char ('aed').
 * 3. All the color miscalls (if in colorspace) ('ced').
 * 4. The score for the alginment, including summary information about
 *    the number of gaps and Ns involved.
 * 5. The reference coordinates for the alignment, including strand.
 * 6. 
 *
 * Note that the AlnRes, together with the Read and an AlnSetSumm
 * (*and* the opposite mate's AlnRes and Read in the case of a paired-
 * end alignment), should contain enough information to print an entire
 * alignment record.
 */
class AlnRes {

public:

	AlnRes() : ned_(RES_CAT), aed_(RES_CAT), ced_(RES_CAT) { reset(); }

	/**
	 * Clear all contents.
	 */
	void reset() {
		ned_.clear();
		aed_.clear();
		ced_.clear();
		score_.invalidate();
		refcoord_.invalidate();
		extent_ = 0;
		seedmms_  = 0; // number of mismatches allowed in seed
		seedlen_  = 0; // length of seed
		seedival_ = 0; // interval between seeds
		penceil_  = 0; // penalty ceiling
		assert(!refcoord_.valid());
	}
	
	/**
	 * Reverse all edit lists.
	 */
	void reverseEdits() {
		ned_.reverse();
		aed_.reverse();
		ced_.reverse();
	}
	
	/**
	 * Invert positions of edits so that they're with respect to the
	 * other end of the read.  Caller must specify the read length.  In
	 * the case of a colorspace read, caller specifies the length in
	 * colors.
	 */
	void invertEdits(size_t rdlen) {
		if(color_) {
			Edit::invertPoss(ned_, rdlen+1);
			Edit::invertPoss(aed_, rdlen+1);
			Edit::invertPoss(ced_, rdlen);
		} else {
			Edit::invertPoss(ned_, rdlen);
			Edit::invertPoss(aed_, rdlen);
			assert(ced_.empty());
		}
	}
	
	/**
	 * Return true iff no result has been installed.
	 */
	bool empty() const {
		if(!VALID_AL_SCORE(score_)) {
			assert(ned_.empty());
			assert(aed_.empty());
			assert(ced_.empty());
			assert(!refcoord_.valid());
			return true;
		} else {
			//assert(refcoord_.valid());
			return false;
		}
	}

	/**
	 * Return the identifier for the reference that the alignment
	 * occurred in.
	 */
	TRefId refid() const { return refcoord_.ref(); }
	
	/**
	 * Return the 0-based offset of the alignment into the reference
	 * sequence it aligned to.
	 */
	TRefOff refoff() const { return refcoord_.off(); }
	
	/**
	 * Set the upstream-most reference offset involved in the
	 * alignment, and the extent of the alignment (w/r/t the
	 * reference)
	 */
	void setCoord(
		TRefId id,
		TRefOff off,
		bool fw,
		size_t extent)
	{
		refcoord_.init(id, off, fw);
		extent_ = extent;
	}
	
	/**
	 * Return true iff the reference chars involved in this alignment
	 * result are entirely within with given bounds.
	 */
	bool within(
		TRefId id,
		TRefOff off,
		bool fw,
		size_t extent) const
	{
		if(refcoord_.ref() == id &&
		   refcoord_.off() >= off &&
		   refcoord_.off() + extent_ <= off + extent &&
		   refcoord_.fw() == fw)
		{
			return true;
		}
		return false;
	}
	
	/**
	 * Set whether this is a colorspace alignment.
	 */
	void setColor(bool col) {
		color_ = col;
	}
	
	/**
	 * Set alignment score for this alignment.
	 */
	void setScore(AlignmentScore score) {
		score_ = score;
	}

	/**
	 * Set the upstream-most and downstream-most nucleotides.
	 */
	void setNucs(bool fw, int nup, int ndn) {
		assert(color_);
		nuc5p_ = fw ? nup : ndn;
		nuc3p_ = fw ? ndn : nup;
	}
	
	/**
	 * Return the reference coordinate where this alignment result
	 * lies.
	 */
	const Coord& refcoord() const {
		return refcoord_;
	}

	/**
	 * Return the reference coordinate where this alignment result
	 * lies.
	 */
	Coord& refcoord() {
		return refcoord_;
	}
	
	AlignmentScore     score()  const { return score_; }
	EList<Edit>&       ned()          { return ned_; }
	EList<Edit>&       aed()          { return aed_; }
	EList<Edit>&       ced()          { return ced_; }
	const EList<Edit>& ned()    const { return ned_; }
	const EList<Edit>& aed()    const { return aed_; }
	const EList<Edit>& ced()    const { return ced_; }
	size_t             extent() const { return extent_; }
	
	/**
	 * Print the sequence for the read that aligned using A, C, G and
	 * T.  This will simply print the read sequence (or its reverse
	 * complement) unless this is a colorspace read and printColors is
	 * false.  In that case, we print the decoded sequence rather than
	 * the original ones.
	 */
 	void printSeq(
		const Read& rd,
		const BTDnaString* dns,
		bool printColors,
		bool exEnds,
		OutFileBuf& o) const;

	/**
	 * Print the quality string for the read that aligned.  This will
	 * simply print the read qualities (or their reverse) unless this
	 * is a colorspace read and printColors is false.  In that case,
	 * we print the decoded qualities rather than the original ones.
	 */
 	void printQuals(
		const Read& rd,
		const BTString* dqs,
		bool printColors,
		bool exEnds,
		OutFileBuf& o) const;
	
	/**
	 * Check that alignment score is internally consistent.
	 */
	bool repOk() const {
		assert(refcoord_.repOk());
		assert(empty() || refcoord_.valid());
		assert(empty() || extent_ > 0);
		return true;
	}
	
	/**
	 * Assuming this AlnRes is an alignment for 'rd', check that the
	 * alignment and 'rd' are compatible with the corresponding
	 * reference sequence.
	 */
	bool matchesRef(
		const Read& rd,
		const BitPairReference& ref);
	
	/**
	 * Set information about the alignment parameters that led to this
	 * alignment.
	 */
	void setParams(
		int seedmms,
		int seedlen,
		int seedival,
		int penceil)
	{
		seedmms_ = seedmms;
		seedlen_ = seedlen;
		seedival_ = seedival;
		penceil_ = penceil;
	}
	
	// Accessors for alignment parameters
	int  seedmms()  const { return seedmms_;  }
	int  seedlen()  const { return seedlen_;  }
	int  seedival() const { return seedival_; }
	int  penceil()  const { return penceil_;  }
	bool color()    const { return color_;    }

	/**
	 * Get the decoded nucleotide sequence 
	 */
	void decodedNucsAndQuals(
		const Read& rd,       // read that led to alignment
		BTDnaString& ns,      // out: decoded nucleotides
		BTString& qs) const;  // out: decoded qualities

protected:

	AlignmentScore score_;    // best SW score found
	EList<Edit>    ned_;      // base edits
	EList<Edit>    aed_;      // ambiguous base resolutions
	EList<Edit>    ced_;      // color miscalls
	Coord          refcoord_; // reference coordinates (sequence index, offset, orientation)
	size_t         extent_;   // number of reference chars involved in alignment
	int            seedmms_;  // number of mismatches allowed in seed
	int            seedlen_;  // length of seed
	int            seedival_; // interval between seeds
	int            penceil_;  // penalty ceiling
	
	bool           color_;    // colorspace alignment?
	int            nuc5p_;    // 5'-most decoded nucleotide; this is chopped off if excluding ends
	int            nuc3p_;    // 3'-most decoded nucleotide; this is chopped off if excluding ends
};

typedef uint64_t TNumAlns;

/**
 * Encapsulates a concise summary of a set of alignment results for a
 * given read.  Referring to the fields of this object should provide
 * enough information to print output records for the read.
 */
class AlnSetSumm {

public:

	AlnSetSumm() {
		best_.invalidate();
		secbest_.invalidate();
		other_ = 0;
	}

	/**
	 * Given an unpaired read (in either rd1 or rd2) or a read pair
	 * (mate 1 in rd1, mate 2 in rd2).
	 */
	AlnSetSumm(
		const Read* rd1,
		const Read* rd2,
		const EList<AlnRes>* rs1,
		const EList<AlnRes>* rs2)
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);
		assert((rd1 == NULL) == (rs1 == NULL));
		assert((rd2 == NULL) == (rs2 == NULL));
		AlignmentScore best, secbest;
		best.invalidate();
		secbest.invalidate();
		bool paired = (rs1 != NULL && rs2 != NULL);
		if(paired) {
			// Paired alignments
			assert_eq(rs1->size(), rs2->size());
			for(size_t i = 0; i < rs1->size(); i++) {
				AlignmentScore sc = (*rs1)[i].score() + (*rs2)[i].score();
				if(sc > best) {
					secbest = best;
					best = sc;
					assert(VALID_AL_SCORE(best));
				} else if(sc > secbest) {
					secbest = sc;
					assert(VALID_AL_SCORE(best));
					assert(VALID_AL_SCORE(secbest));
				}
			}
		} else {
			// Unpaired alignments
			const EList<AlnRes>* rs = (rs1 != NULL ? rs1 : rs2);
			assert(rs != NULL);
			for(size_t i = 0; i < rs->size(); i++) {
				AlignmentScore sc = (*rs)[i].score();
				if(sc > best) {
					secbest = best;
					best = sc;
					assert(VALID_AL_SCORE(best));
				} else if(sc > secbest) {
					secbest = sc;
					assert(VALID_AL_SCORE(best));
					assert(VALID_AL_SCORE(secbest));
				}
			}
		}
		init(best, secbest, rs1->size()-1);
	}

	AlnSetSumm(
		AlignmentScore best,
		AlignmentScore secbest,
		TNumAlns other)
	{
		init(best, secbest, other);
	}
	
	void init(
		AlignmentScore best,
		AlignmentScore secbest,
		TNumAlns other)
	{
		best_    = best;
		secbest_ = secbest;
		other_   = other;
		assert(repOk());
	}
	
	/**
	 * Return true iff there is at least a best alignment
	 */
	bool empty() const {
		assert(repOk());
		return !VALID_AL_SCORE(best_);
	}
	
	/**
	 * Check that the summary is internally consistent.
	 */
	bool repOk() const {
		assert(other_ == 0 ||  VALID_AL_SCORE(secbest_));
		assert(other_ != 0 || !VALID_AL_SCORE(secbest_));
		return true;
	}
	
	AlignmentScore best()    const { return best_;    }
	AlignmentScore secbest() const { return secbest_; }
	TNumAlns       other()   const { return other_;   }
	
protected:
	
	AlignmentScore best_;    // best full-alignment score found for this read
	AlignmentScore secbest_; // second-best
	TNumAlns       other_;   // # more alignments within N points of second-best
};

#endif
