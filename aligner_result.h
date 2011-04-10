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

#define VALID_AL_SCORE(x)   ((x).score_ > std::numeric_limits<TAlScore>::min())
#define VALID_SCORE(x)      ((x) > std::numeric_limits<TAlScore>::min())
#define INVALIDATE_SCORE(x) ((x) = std::numeric_limits<TAlScore>::min())

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
class AlnScore {

public:

	/**
	 * Gapped scores are invalid until proven valid.
	 */
	AlnScore() {
		score_ = ns_ = gaps_ = 0;
		invalidate();
		assert(!valid());
	}

	/**
	 * Return an invalid SwScore.
	 */
	inline static AlnScore INVALID() {
		AlnScore s;
		s.invalidate();
		assert(!s.valid());
		return s;
	}
	
	/**
	 * Return true iff this score has a valid value.
	 */
	inline bool valid() const {
		return score_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Make this score invalid (and therefore <= all other scores).
	 */
	inline void invalidate() {
		score_ = std::numeric_limits<TAlScore>::min();
		assert(!valid());
	}
	
	/**
	 * Increment the number of gaps.  If currently invalid, this makes
	 * the score valid with gaps == 1.
	 */
	inline void incNs(int nceil) {
		if(++ns_ > nceil) {
			invalidate();
		}
		assert_lt(ns_, 0x7fffffff);
	}

	/**
	 * Return true iff this score is > score o.
	 * Note: An "invalid" score is <= all other scores.
	 */
	bool operator>(const AlnScore& o) const {
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
	inline AlnScore& operator=(const AlnScore& o) {
		// Profiling shows many cache misses on following lines
		gaps_  = o.gaps_;
		ns_    = o.ns_;
		score_ = o.score_;
		assert_lt(ns_, 0x7fffffff);
		return *this;
	}

	/**
	 * Scores are equal iff they're bitwise equal.
	 */
	inline bool operator==(const AlnScore& o) const {
		// Profiling shows cache misses on following line
		return VALID_AL_SCORE(*this) && VALID_AL_SCORE(o) && score_ == o.score_;
	}

	/**
	 * Return true iff the two scores are unequal.
	 */
	inline bool operator!=(const AlnScore& o) const {
		return !(*this == o);
	}

	/**
	 * Return true iff this score is >= score o.
	 */
	bool operator>=(const AlnScore& o) const {
		return operator==(o) || operator>(o);
	}

	/**
	 * Return true iff this score is < score o.
	 */
	bool operator<(const AlnScore& o) const {
		return !operator>=(o);
	}

	/**
	 * Return true iff this score is <= score o.
	 */
	bool operator<=(const AlnScore& o) const {
		return !operator>(o);
	}

	/**
	 * Calculate difference between two SwScores.
	 */
	AlnScore operator-(const AlnScore& o) const {
		if(!VALID_AL_SCORE(*this)) return *this;
		AlnScore s; 
		s.gaps_ = gaps_ - o.gaps_;
		s.ns_ = ns_;
		s.score_ = score_ - o.score_;
		assert_lt(s.ns_, 0x7fffffff);
		return s;
	}

	/**
	 * Calculate sum of two SwScores.
	 */
	AlnScore operator+(const AlnScore& o) const {
		if(!VALID_AL_SCORE(*this)) return *this;
		AlnScore s;
		s.gaps_ = gaps_ + o.gaps_;
		s.ns_ = ns_;
		s.score_ = score_ + o.score_;
		assert_lt(s.ns_, 0x7fffffff);
		return s;
	}

	/**
	 * Add given SwScore into this one.
	 */
	AlnScore operator+=(const AlnScore& o) {
		if(VALID_AL_SCORE(*this)) {
			gaps_ += o.gaps_;
			score_ += o.score_;
		}
		return (*this);
	}

	/**
	 * Subtract given SwScore from this one.
	 */
	inline AlnScore operator-=(const AlnScore& o) {
		if(VALID_AL_SCORE(*this)) {
			gaps_ -= o.gaps_;
			score_ -= o.score_;
		}
		return (*this);
	}

	/**
	 * Calculate difference between two SwScores.
	 */
	inline AlnScore operator-(int o) const {
		if(!VALID_AL_SCORE(*this)) return *this;
		AlnScore s;
		s.gaps_ = gaps_;
		s.ns_ = ns_;
		s.score_ = score_ - o;
		assert_lt(s.ns_, 0x7fffffff);
		return s;
	}
	
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

static inline ostream& operator<<(ostream& os, const AlnScore& o) {
	os << o.score();
	return os;
}

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

	AlnRes() :
		ned_(RES_CAT),
		aed_(RES_CAT),
		ced_(RES_CAT)
	{
		reset();
	}

	/**
	 * Clear all contents.
	 */
	void reset() {
		ned_.clear();
		aed_.clear();
		ced_.clear();
		score_.invalidate();
		refcoord_.invalidate();
		rdlen_ = 0;
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
	void invertEdits() {
		assert_gt(rdlen_, 0);
		if(color_) {
			Edit::invertPoss(ned_, rdlen_+1);
			Edit::invertPoss(aed_, rdlen_+1);
			Edit::invertPoss(ced_, rdlen_);
		} else {
			Edit::invertPoss(ned_, rdlen_);
			Edit::invertPoss(aed_, rdlen_);
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
	inline TRefId refid() const { return refcoord_.ref(); }
	
	/**
	 * Return the 0-based offset of the alignment into the reference
	 * sequence it aligned to.
	 */
	inline TRefOff refoff() const { return refcoord_.off(); }
	
	/**
	 * Set arguments to coordinates for the upstream-most and downstream-most
	 * reference positions involved in the alignment.
	 */
	inline void getCoords(
		Coord& st,  // out: install starting coordinate here
		Coord& en)  // out: install ending coordinate here
	{
		st.init(refcoord_);
		en.init(refcoord_);
		en.setOff(en.off() + extent() - 1);
	}
	
	/**
	 * Set the upstream-most reference offset involved in the
	 * alignment, and the extent of the alignment (w/r/t the
	 * reference)
	 */
	inline void setCoord(
		TRefId id,
		TRefOff off,
		bool fw,
		size_t extent)
	{
		refcoord_.init(id, off, fw);
		extent_ = extent;
	}
	
	/**
	 * Set the length of the original nucleotide read (not the decoded read).
	 */
	inline void setReadLength(size_t rdlen) {
		rdlen_ = rdlen;
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
	void setScore(AlnScore score) {
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
	
	/**
	 * Return true if this alignment is to the Watson strand.
	 */
	inline bool fw() const {
		return refcoord_.fw();
	}
	
	AlnScore     score()  const { return score_; }
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
	 * Print a stacked alignment with the reference on top, query on bottom,
	 * and lines connecting matched-up positions.
	 */
	void printStacked(
		const Read& rd,
		std::ostream& o) const
	{
		printStacked(refcoord_.fw() ? rd.patFw : rd.patRc, o);
	}

	/**
	 * Print a stacked alignment with the reference on bottom, query on top,
	 * and lines connecting matched-up positions.
	 */
	void printStacked(
		const BTDnaString& seq,
		std::ostream& o) const
	{
		Edit::printQAlign(o, seq, ned_);
		// Print reference offset below reference string
		o << "^" << std::endl;
		o << "(" << refcoord_.ref() << "," << refcoord_.off() << ")" << std::endl;
	}
	
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
	
	/**
	 * Return true iff this AlnRes and the given AlnRes overlap.  Two AlnRess
	 * overlap if they share a cell in the overall dynamic programming table:
	 * i.e. if there exists a read position s.t. that position in both reads
	 * matches up with the same reference character.  E.g., the following
	 * alignments (drawn schematically as paths through a dynamic programming
	 * table) are redundant:
	 *
	 *  a  b           a  b
	 *  \  \           \  \
	 *   \  \           \  \
	 *    \  \           \  \
	 *     ---\           \  \
	 *         \           ---\---
	 *       ---\              \  \
	 *        \  \              \  \
	 *         \  \              \  \
	 *          \  \              \  \
	 *          a  b              a  b
	 */
	bool overlap(AlnRes& res) {
		if(fw() != res.fw() || refid() != res.refid()) {
			// Must be same reference and same strand in order to overlap
			return false;
		}
		TRefOff my_left     = refoff();
		TRefOff other_left  = res.refoff();
		TRefOff my_right    = my_left    + extent();
		TRefOff other_right = other_left + res.extent();
		if(my_right < other_left || other_right < my_left) {
			// Their rectangular hulls don't overlap, so they can't overlap
			return false;
		}
		// Reference and strand are the same and hulls overlap.  Now go read
		// position by read position testing if any align identically with the
		// reference.
		if(!fw()) {
			invertEdits();
		}
		if(!res.fw()) {
			res.invertEdits();
		}
		size_t nedidx = 0, onedidx = 0;
		size_t nrow = rdlen_ + (color_ ? 1 : 0);
		bool olap = false;
		// For each row...
		for(size_t i = 0; i < nrow; i++) {
			size_t diff = 1;  // amount to shift to right for next round
			size_t odiff = 1; // amount to shift to right for next round
			// Unless there are insertions before the next position, we say
			// that there is one cell in this row involved in the alignment
			my_right = my_left + 1;
			other_right = other_left + 1;
			while(nedidx < ned_.size() && ned_[nedidx].pos == i) {
				if(ned_[nedidx].isDelete()) {
					// Next my_left will be in same column as this round
					diff = 0;
				}
				nedidx++;
			}
			while(onedidx < res.ned_.size() && res.ned_[onedidx].pos == i) {
				if(res.ned_[onedidx].isDelete()) {
					// Next my_left will be in same column as this round
					odiff = 0;
				}
				onedidx++;
			}
			if(i < nrow-1) {
				// See how many inserts there are before the next read
				// character
				size_t nedidx_next  = nedidx;
				size_t onedidx_next = onedidx;
				while(nedidx_next < ned_.size() &&
				      ned_[nedidx_next].pos == i+1)
				{
					if(ned_[nedidx_next].isInsert()) {
						my_right++;
					}
					nedidx_next++;
				}
				while(onedidx_next < res.ned_.size() &&
				      res.ned_[onedidx_next].pos == i+1)
				{
					if(res.ned_[onedidx_next].isInsert()) {
						other_right++;
					}
					onedidx_next++;
				}
			}
			// Contained?
			olap =
				(my_left >= other_left && my_right <= other_right) ||
				(other_left >= my_left && other_right <= my_right);
			// Overlapping but not contained?
			if(!olap) {
				olap =
					(my_left <= other_left && my_right > other_left) ||
					(other_left <= my_left && other_right > my_left);
			}
			if(olap) {
				break;
			}
			// How to do adjust my_left and my_right
			my_left = my_right + diff - 1;
			other_left = other_right + odiff - 1;
		}
		if(!fw()) {
			invertEdits();
		}
		if(!res.fw()) {
			res.invertEdits();
		}
		return olap;
	}
	
	/**
	 * Initialize new AlnRes.
	 */
	void init(
		size_t             rdlen,
		AlnScore           score,
		const EList<Edit>* ned,
		const EList<Edit>* aed,
		const EList<Edit>* ced,
		Coord              refcoord,
		bool               color,
		int                seedmms  = -1,
		int                seedlen  = -1,
		int                seedival = -1,
		int                penceil  = -1,
		int                nuc5p    = -1,
		int                nuc3p    = -1)
	{
		rdlen_ = rdlen;
		score_ = score;
		ned_.clear();
		aed_.clear();
		ced_.clear();
		if(ned != NULL) ned_ = *ned;
		if(aed != NULL) aed_ = *aed;
		if(ced != NULL) ced_ = *ced;
		refcoord_ = refcoord;
		color_ = color;
		seedmms_ = seedmms;
		seedlen_ = seedlen;
		seedival_ = seedival;
		penceil_ = penceil;
		nuc5p_ = nuc5p;
		nuc3p_ = nuc3p;
		extent_ = rdlen;
		for(size_t i = 0; i < ned_.size(); i++) {
			if(ned_[i].isDelete()) extent_--;
			if(ned_[i].isInsert()) extent_++;
		}
	}

protected:

	size_t      rdlen_;    // length of the original read
	AlnScore    score_;    // best SW score found
	EList<Edit> ned_;      // base edits
	EList<Edit> aed_;      // ambiguous base resolutions
	EList<Edit> ced_;      // color miscalls
	Coord       refcoord_; // ref coordinates (seq idx, offset, orient)
	size_t      extent_;   // number of ref chars involved in alignment
	int         seedmms_;  // number of mismatches allowed in seed
	int         seedlen_;  // length of seed
	int         seedival_; // interval between seeds
	int         penceil_;  // penalty ceiling
	bool        color_;    // colorspace alignment?
	int         nuc5p_;    // 5'-most decoded base; clipped if excluding end
	int         nuc3p_;    // 3'-most decoded base; clipped if excluding end
};

typedef uint64_t TNumAlns;

/**
 * Encapsulates a concise summary of a set of alignment results for a
 * given read.  Referring to the fields of this object should provide
 * enough information to print output records for the read.
 */

enum {
	// This alignment is one of a pair of alignments that form a concordant
	// alignment for a read
	ALN_FLAG_PAIR_CONCORD = 1,

	// This alignment is one of a pair of alignments that form a discordant
	// alignment for a read
	ALN_FLAG_PAIR_DISCORD,
	
	// This is an unpaired alignment but the read in question is a pair;
	// usually, this happens because the read had no reportable paired-end
	// alignments
	ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR,

	// This is an unpaired alignment of an unpaired read
	ALN_FLAG_PAIR_UNPAIRED
};

class AlnFlags {

public:

	AlnFlags(int pairing, bool maxed, bool maxedPair) {
		init(pairing, maxed, maxedPair);
	}

	/**
	 * Initialize given values for all settings.
	 */
	void init(int pairing, bool maxed, bool maxedPair) {
		assert_gt(pairing, 0);
		assert_leq(pairing, ALN_FLAG_PAIR_UNPAIRED);
		pairing_ = pairing;
		maxed_ = maxed;
		maxedPair_ = maxedPair;
	}

	/**
	 * Return true iff this alignment is from a paired-end read.
	 */
	bool partOfPair() const {
		assert_gt(pairing_, 0);
		return pairing_ < ALN_FLAG_PAIR_UNPAIRED;
	}
	
	/**
	 * Check that the flags are internally consistent.
	 */
	bool repOk() const {
		assert(partOfPair() || !maxedPair_);
		return true;
	}
	
	inline int  pairing()   const { return pairing_; }
	inline bool maxed()     const { return maxed_; }
	inline bool maxedPair() const { return maxedPair_; }

protected:

	// See ALN_FLAG_PAIR_* above
	int pairing_;

	// This alignment is sampled from among many alignments that, taken
	// together, cause this mate to align non-uniquely
	bool maxed_;
	
	// The paired-end read of which this mate is part has repetitive concordant
	// alignments
	bool maxedPair_;
};

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
		const EList<AlnRes>* rs2);

	AlnSetSumm(
		AlnScore best,
		AlnScore secbest,
		TNumAlns other)
	{
		init(best, secbest, other);
	}
	
	void init(
		AlnScore best,
		AlnScore secbest,
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
	
	AlnScore best()    const { return best_;    }
	AlnScore secbest() const { return secbest_; }
	TNumAlns       other()   const { return other_;   }
	
protected:
	
	AlnScore best_;    // best full-alignment score found for this read
	AlnScore secbest_; // second-best
	TNumAlns       other_;   // # more alignments within N points of second-best
};

#endif
