/*
 *  scoring.h
 */

#ifndef SCORING_H_
#define SCORING_H_

#include "qual.h"
#include "simple_func.h"

// Default type of bonus to added for matches
#define DEFAULT_MATCH_BONUS_TYPE COST_MODEL_CONSTANT
// When match bonus type is constant, use this constant
#define DEFAULT_MATCH_BONUS 0
// Same settings but different defaults for --local mode
#define DEFAULT_MATCH_BONUS_TYPE_LOCAL COST_MODEL_CONSTANT
#define DEFAULT_MATCH_BONUS_LOCAL 2

// Default type of penalty to assess against mismatches
#define DEFAULT_MM_PENALTY_TYPE COST_MODEL_CONSTANT
// When mismatch penalty type is constant, use this constant
#define DEFAULT_MM_PENALTY 6

// Default type of penalty to assess against mismatches
#define DEFAULT_N_PENALTY_TYPE COST_MODEL_CONSTANT
// When mismatch penalty type is constant, use this constant
#define DEFAULT_N_PENALTY 1

// Constant coefficient b in linear function f(x) = ax + b determining
// minimum valid score f when read length is x
#define DEFAULT_MIN_CONST (-0.6f)
// Linear coefficient a
#define DEFAULT_MIN_LINEAR (-0.6f)
// Different defaults for --local mode
#define DEFAULT_MIN_CONST_LOCAL (0.0f)
#define DEFAULT_MIN_LINEAR_LOCAL (0.66f)

// Constant coefficient b in linear function f(x) = ax + b determining
// what local-alignment score floor to impose when read length is x
#define DEFAULT_FLOOR_CONST (-std::numeric_limits<float>::max())
// Linear coefficient a
#define DEFAULT_FLOOR_LINEAR 0.0f
// Different defaults for --local mode
#define DEFAULT_FLOOR_CONST_LOCAL 0.0f
#define DEFAULT_FLOOR_LINEAR_LOCAL 0.0f

// Constant coefficient b in linear function f(x) = ax + b determining
// maximum permitted number of Ns f in a read before it is filtered &
// the maximum number of Ns in an alignment before it is considered
// invalid.
#define DEFAULT_N_CEIL_CONST 0.0f
// Linear coefficient a
#define DEFAULT_N_CEIL_LINEAR 0.15f

// Default for whether to concatenate mates before the N filter (as opposed to
// filting each mate separately)
#define DEFAULT_N_CAT_PAIR false

// Default penalty to assess against SNPs in colorspace alignments.
// Decoding must have occurred in order to distinguish SNPs from
// patterns of mismatches.
#define DEFAULT_SNP_PENALTY 6

// Default read gap penalties for when homopolymer calling is reliable	
#define DEFAULT_READ_GAP_CONST 5
#define DEFAULT_READ_GAP_LINEAR 3

// Default read gap penalties for when homopolymer calling is not reliable
#define DEFAULT_READ_GAP_CONST_BADHPOLY 3
#define DEFAULT_READ_GAP_LINEAR_BADHPOLY 1

// Default reference gap penalties for when homopolymer calling is reliable
#define DEFAULT_REF_GAP_CONST 5
#define DEFAULT_REF_GAP_LINEAR 3

// Default reference gap penalties for when homopolymer calling is not reliable
#define DEFAULT_REF_GAP_CONST_BADHPOLY 3
#define DEFAULT_REF_GAP_LINEAR_BADHPOLY 1

enum {
	COST_MODEL_ROUNDED_QUAL = 1,
	COST_MODEL_QUAL,
	COST_MODEL_CONSTANT
};

/**
 * How to penalize various types of sequence dissimilarity, and other settings
 * that govern how dynamic programming tables should be filled in and how to
 * backtrace to find solutions.
 */
class Scoring {

	/**
	 * Init an array that maps quality to penalty or bonus according to 'type'
	 * and 'cons'
	 */
	template<typename T>
	void initPens(
		T *pens,  // array to fill
		int type,   // penalty type; qual | rounded qual | constant
		int cons)   // constant for when penalty type is constant
	{
		if(type == COST_MODEL_ROUNDED_QUAL) {
			for(int i = 0; i < 256; i++) {
				pens[i] = (T)qualRounds[i];
			}
		} else if(type == COST_MODEL_QUAL) {
			for(int i = 0; i < 256; i++) {
				int ii = i;
				if(ii > 30) {
					ii = 30;
				}
				float frac = (float)ii / 30.0f;
				pens[i] = (T)(frac * cons);
				if(pens[i] == 0) {
					pens[i] = ((cons > 0) ? (T)1 : (T)-1);
				}
			}
		} else if(type == COST_MODEL_CONSTANT) {
			for(int i = 0; i < 256; i++) {
				pens[i] = (T)cons;
			}
		} else {
			throw 1;
		}
	}

public:

	Scoring(
		int   mat,          // reward for a match
		int   mmcType,      // how to penalize mismatches
	    int   mmc,          // constant if mm pelanty is a constant
		int   sn,           // penalty for nuc mismatch in decoded color alns
		const SimpleFunc& scoreMin_,   // minimum score for valid alignment; const coeff
		const SimpleFunc& scoreFloor_, // local-alignment score floor; const coeff
		const SimpleFunc& nCeil_,      // max # ref Ns allowed in alignment; const coeff
	    int   nType,        // how to penalize Ns in the read
	    int   n,            // constant if N pelanty is a constant
		bool  ncat,         // whether to concatenate mates before N filtering
	    int   rdGpConst,    // constant coeff for cost of gap in the read
	    int   rfGpConst,    // constant coeff for cost of gap in the ref
	    int   rdGpLinear,   // coeff of linear term for cost of gap in read
	    int   rfGpLinear,   // coeff of linear term for cost of gap in ref
		int     gapbar_,    // # rows at top/bot can only be entered diagonally
		int64_t rowlo_,     // min row idx to backtrace from; -1 = no limit
		bool    rowFirst_)  // sort results first by row then by score?
	{
		matchType    = COST_MODEL_CONSTANT;
		matchConst   = mat;
		mmcostType   = mmcType;
		mmcost       = mmc;
		snp          = sn;
		scoreMin     = scoreMin_;
		scoreFloor   = scoreFloor_;
		nCeil        = nCeil_;
		npenType     = nType;
		npen         = n;
		ncatpair     = ncat;
		rdGapConst   = rdGpConst;
		rfGapConst   = rfGpConst;
		rdGapLinear  = rdGpLinear;
		rfGapLinear  = rfGpLinear;
		qualsMatter_ = mmcostType != COST_MODEL_CONSTANT;
		gapbar       = gapbar_;
		rowlo        = rowlo_;
		rowFirst     = rowFirst_;
		monotone     = matchType == COST_MODEL_CONSTANT && matchConst == 0;
		initPens<int>(mmpens, mmcostType, mmcost);
		initPens<int>(npens, npenType, npen);
		initPens<float>(matchBonuses, matchType, matchConst);
		assert(repOk());
	}
	
	/**
	 * Set a constant match bonus.
	 */
	void setMatchBonus(int bonus) {
		matchType  = COST_MODEL_CONSTANT;
		matchConst = bonus;
		initPens<float>(matchBonuses, matchType, matchConst);
		assert(repOk());
	}
	
	/**
	 * Set the mismatch penalty.
	 */
	void setMmPen(int mmType, int c) {
		mmcostType = mmType;
		mmcost     = c;
		initPens<int>(mmpens, mmcostType, mmcost);
	}
	
	/**
	 * Set the N penalty.
	 */
	void setNPen(int nType, int n) {
		npenType     = nType;
		npen         = n;
		initPens<int>(npens, npenType, npen);
	}
	
	/**
	 * Check that scoring scheme is internally consistent.
	 */
	bool repOk() const {
		assert_geq(matchConst, 0);
		assert_gt(snp, 0);
		assert_gt(rdGapConst, 0);
		assert_gt(rdGapLinear, 0);
		assert_gt(rfGapConst, 0);
		assert_gt(rfGapLinear, 0);
		return true;
	}

	/**
	 * Return a linear function of x where 'cnst' is the constant coefficiant
	 * and 'lin' is the linear coefficient.
	 */
	static float linearFunc(int64_t x, float cnst, float lin) {
		return cnst + (lin * x);
	}

	/**
	 * Return the penalty incurred by a mismatch at an alignment column
	 * with read character 'rdc' reference mask 'refm' and quality 'q'.
	 *
	 * qs should be clamped to 63 on the high end before this query.
	 */
	inline int mm(int rdc, int refm, int q) const {
		assert_range(0, 255, q);
		return (rdc > 3 || refm > 15) ? npens[q] : mmpens[q];
	}
	
	/**
	 * Return the score of the given read character with the given quality
	 * aligning to the given reference mask.
	 */
	inline int score(int rdc, int refm, int q) const {
		assert_range(0, 255, q);
		if(rdc > 3 || refm > 15) {
			return -npens[q];
		}
		if((refm & (1 << rdc)) != 0) {
			return (int)matchBonuses[q];
		} else {
			return -mmpens[q];
		}
	}

	/**
	 * Return the penalty incurred by a mismatch at an alignment column
	 * with read character 'rdc' and quality 'q'.  We assume the
	 * reference character is non-N.
	 */
	inline int mm(int rdc, int q) const {
		assert_range(0, 255, q);
		return (rdc > 3) ? npens[q] : mmpens[q];
	}
	
	/**
	 * Return the marginal penalty incurred by a mismatch at a read
	 * position with quality 'q'.
	 */
	inline int mm(int q) const {
		assert_geq(q, 0);
		return q < 255 ? mmpens[q] : mmpens[255];
	}

	/**
	 * Return the marginal penalty incurred by a mismatch at a read
	 * position with quality 30.
	 */
	inline int64_t match() const {
		return match(30);
	}

	/**
	 * Return the marginal penalty incurred by a mismatch at a read
	 * position with quality 'q'.
	 */
	inline int64_t match(int q) const {
		assert_geq(q, 0);
		return (int64_t)((q < 255 ? matchBonuses[q] : matchBonuses[255]) + 0.5f);
	}
	
	/**
	 * Return the best score achievable by a read of length 'rdlen'.
	 */
	inline int64_t perfectScore(size_t rdlen) const {
		if(monotone) {
			return 0;
		} else {
			return rdlen * match(30);
		}
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
	inline int n(int q) const {
		assert_geq(q, 0);
		return q < 255 ? npens[q] : npens[255];
	}

	
	/**
	 * Return the marginal penalty incurred by a gap in the read,
	 * given that this is the 'ext'th extension of the gap (0 = open,
	 * 1 = first, etc).
	 */
	inline int ins(int ext) const {
		assert_geq(ext, 0);
		if(ext == 0) return readGapOpen();
		return readGapExtend();
	}

	/**
	 * Return the marginal penalty incurred by a gap in the reference,
	 * given that this is the 'ext'th extension of the gap (0 = open,
	 * 1 = first, etc).
	 */
	inline int del(int ext) const {
		assert_geq(ext, 0);
		if(ext == 0) return refGapOpen();
		return refGapExtend();
	}

	/**
	 * Return true iff a read of length 'rdlen' passes the score filter, i.e.,
	 * has enough characters to rise above the minimum score threshold.
	 */
	bool scoreFilter(
		int64_t minsc,
		size_t rdlen) const;

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

	int     matchType;    // how to reward matches
	int     matchConst;   // reward for a match
	int     mmcostType;   // based on qual? rounded? just a constant?
	int     mmcost;       // if mmcosttype=constant, this is the const
	int     snp;          // penalty for nuc mismatch in decoded colorspace als
	SimpleFunc scoreMin;  // minimum score for valid alignment, constant coeff
	SimpleFunc scoreFloor;// local-alignment score floor, constant coeff
	SimpleFunc nCeil;     // max # Ns involved in alignment, constant coeff
	int     npenType;     // N: based on qual? rounded? just a constant?
	int     npen;         // N: if mmcosttype=constant, this is the const
	bool    ncatpair;     // true -> do N filtering on concated pair
	int     rdGapConst;   // constant term coeffecient in extend cost
	int     rfGapConst;   // constant term coeffecient in extend cost
	int     rdGapLinear;  // linear term coeffecient in extend cost
	int     rfGapLinear;  // linear term coeffecient in extend cost
	int     gapbar;       // # rows at top/bot can only be entered diagonally
	int64_t rowlo;        // min row idx to backtrace from; -1 = no limit
	bool    rowFirst;     // sort results first by row then by score?
	bool    monotone;     // scores can only go down?
	float   matchBonuses[256]; // map from qualities to match bonus
	int     mmpens[256];  // map from qualities to mm penalty
	int     npens[256];   // map from N qualities to penalty
	
	static Scoring bwaSwLike() {
		const double DMAX = std::numeric_limits<double>::max();
		SimpleFunc scoreFloor(SIMPLE_FUNC_CONST, 0.0f, 0.0f, 0.0f, 0.0f);
		SimpleFunc scoreMin(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 37.0f, 0.3f);
		SimpleFunc nCeil(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
		return Scoring(
			1,                       // reward for a match
			COST_MODEL_CONSTANT,     // how to penalize mismatches
			3,                       // constant if mm pelanty is a constant
			3,                       // penalty for nuc mm when decoding colors
			scoreMin,                // score min: 37 + 0.3x
			scoreFloor,              // score floor: constant = 0
			nCeil,                   // n ceiling: 2 + 0.1x
			COST_MODEL_CONSTANT,     // how to penalize Ns in the read
			3,                       // constant if N pelanty is a constant
			false,                   // concatenate mates before N filtering?
			5,                       // constant coeff for gap in read
			5,                       // constant coeff for gap in ref
			2,                       // linear coeff for gap in read
			2,                       // linear coeff for gap in ref
			5,                       // 5 rows @ top/bot diagonal-entrance-only
			-1,                      // no restriction on row
			false);                  // score prioritized over row
	}

	static Scoring base1() {
		const double DMAX = std::numeric_limits<double>::max();
		SimpleFunc scoreFloor(SIMPLE_FUNC_CONST, 0.0f, 0.0f, 0.0f, 0.0f);
		SimpleFunc scoreMin(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 37.0f, 0.3f);
		SimpleFunc nCeil(SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
		return Scoring(
			1,                       // reward for a match
			COST_MODEL_CONSTANT,     // how to penalize mismatches
			3,                       // constant if mm pelanty is a constant
			3,                       // penalty for nuc mm when decoding colors
			scoreMin,                // score min: 37 + 0.3x
			scoreFloor,              // score floor: constant = 0
			nCeil,                   // n ceiling: 2 + 0.1x
			COST_MODEL_CONSTANT,     // how to penalize Ns in the read
			3,                       // constant if N pelanty is a constant
			false,                   // concatenate mates before N filtering?
			11,                      // constant coeff for gap in read
			11,                      // constant coeff for gap in ref
			4,                       // linear coeff for gap in read
			4,                       // linear coeff for gap in ref
			5,                       // 5 rows @ top/bot diagonal-entrance-only
			-1,                      // no restriction on row
			false);                  // score prioritized over row
	}

protected:

	bool qualsMatter_;
};

#endif /*SCORING_H_*/
