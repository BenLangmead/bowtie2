/*
 *  aligner_seed_policy.h
 */

#ifndef ALIGNER_SEED_POLICY_H_
#define ALIGNER_SEED_POLICY_H_

#include "penalty.h"

enum {
	SEED_IVAL_LINEAR = 1,
	SEED_IVAL_SQUARE_ROOT,
	SEED_IVAL_CUBE_ROOT
};

// Default type of penalty to assess against mismatches
#define DEFAULT_MM_PENALTY_TYPE COST_MODEL_CONSTANT
// When mismatch penalty type is constant, use this constant
#define DEFAULT_MM_PENALTY 30

// Default type of penalty to assess against mismatches
#define DEFAULT_N_PENALTY_TYPE COST_MODEL_CONSTANT
// When mismatch penalty type is constant, use this constant
#define DEFAULT_N_PENALTY 1

// Constant coefficient b in linear function f(x) = ax + b determining
// maximum allowed penalty f when read length is x
#define DEFAULT_CEIL_CONST 3.0f
// Linear coefficient a
#define DEFAULT_CEIL_LINEAR 3.0f

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

// Default penalty to asses against SNPs in colorspace alignments.
// Decoding must have occurred in order to distinguish SNPs from
// patterns of mismatches.
#define DEFAULT_SNP_PENALTY 30 

// Default read gap penalties for when homopolymer calling is reliable	
#define DEFAULT_READ_OPEN 40
#define DEFAULT_READ_EXTEND_CONST 15
#define DEFAULT_READ_EXTEND_LINEAR 0

// Default read gap penalties for when homopolymer calling is not reliable
#define DEFAULT_READ_OPEN_BADHPOLY 20
#define DEFAULT_READ_EXTEND_CONST_BADHPOLY 5
#define DEFAULT_READ_EXTEND_LINEAR_BADHPOLY 0

// Default reference gap penalties for when homopolymer calling is reliable
#define DEFAULT_REF_OPEN 40
#define DEFAULT_REF_EXTEND_CONST 15
#define DEFAULT_REF_EXTEND_LINEAR 0

// Default reference gap penalties for when homopolymer calling is not reliable
#define DEFAULT_REF_OPEN_BADHPOLY 20
#define DEFAULT_REF_EXTEND_CONST_BADHPOLY 5
#define DEFAULT_REF_EXTEND_LINEAR_BADHPOLY 0

#define DEFAULT_SEEDMMS 0
#define DEFAULT_SEEDLEN 22
#define DEFAULT_SEEDPERIOD -1

#define DEFAULT_IVAL SEED_IVAL_SQUARE_ROOT
#define DEFAULT_IVAL_A 1.0f
#define DEFAULT_IVAL_B 0.0f

/**
 * Encapsulates the set of all parameters that affect what the
 * SeedAligner does with reads.
 */
class SeedAlignmentPolicy {

public:

	/**
	 * Parse alignment policy when provided in this format:
	 * <lab>=<val>;<lab>=<val>;<lab>=<val>...
	 *
	 * And label=value possibilities are:
	 *
	 * Penalty for a mismatch
	 * ----------------------
	 *
	 * MMP={Cxx|Q|RQ} (default: MMP=C30)
	 *
	 *   Cxx = Each mismatch costs xx.  If MMP=Cxx is specified,
	 *         quality values are ignored when assessing penalities for
	 *         mismatches.
	 *   Q   = Each mismatch incurs a penalty equal to the mismatched
	 *         base's value.
	 *   R   = Each mismatch incurs a penalty equal to the mismatched
	 *         base's rounded quality value.  Qualities are rounded off
	 *         to the nearest 10, and qualities greater than 30 are
	 *         rounded to 30.
	 * 
	 * Penalty for position with N (in either read or reference)
	 * ---------------------------------------------------------
	 *
	 * NP={Cxx|Q|RQ} (default: NP=C1)
	 *
	 *   Cxx = Each alignment position with an N in either the read or
	 *         the reference costs xx.  If NP=Cxx is specified, quality
	 *         values are ignored when assessing penalities for Ns.
	 *   Q   = Each alignment position with an N in either the read or
	 *         the reference incurs a penalty equal to the read base's
	 *         quality value.
	 *   R   = Each alignment position with an N in either the read or
	 *         the reference incurs a penalty equal to the read base's
	 *         rounded quality value.  Qualities are rounded off to
	 *         the nearest 10, and qualities greater than 30 are
	 *         rounded to 30.
	 *
	 * Penalty for a read gap
	 * ----------------------
	 *
	 * RDG=xx,yy,zz (default: RDG=40,15,0)
	 *
	 *   xx    = Read gap open penalty.  Every contiguous read gap
	 *           incurs this penalty, which is added to any extension
	 *           penalties.  Must be > 0.
	 *   yy,zz = Read gap extension penalty.  Once the gap is opened,
	 *           each extension gets a
	 *           penalty = yy - zz * (gap length - 2).
	 *
	 * Penalty for a reference gap
	 * ---------------------------
	 *
	 * RFG=xx,yy,zz (default: RFG=40,15,0)
	 *
	 *   xx    = Reference gap open penalty.  Every contiguous
	 *           reference gap incurs this penalty, which is added to
	 *           any extension penalties.  Must be > 0.
	 *   yy,zz = Reference gap extension penalty.  Once the gap is
	 *           opened, each extension gets a
	 *           penalty = yy - zz * (gap length - 2).
	 *
	 * Read penalty ceiling
	 * --------------------
	 *
	 * CEIL=xx,yy (default: CEIL=3.0,3.0)
	 *
	 *   xx,yy = For a read of length N, the sum of all penalties
	 *           incurred by mismatches and gaps cannot exceed
	 *           ceiling = xx + (read length * yy).  If the ceiling is
	 *           exceeded, the alignment is considered invalid.
	 *
	 * N ceiling
	 * ---------
	 *
	 * NCEIL=xx,yy (default: NCEIL=0.0,0.15)
	 *
	 *   xx,yy = For a read of length N, the number of alignment
	 *           positions with an N in either the read or the
	 *           reference cannot exceed
	 *           ceiling = xx + (read length * yy).  If the ceiling is
	 *           exceeded, the alignment is considered invalid.
	 *
	 * Seeds
	 * -----
	 *
	 * SEED=mm,len,ival (default: SEED=0,22)
	 *
	 *   mm   = Maximum number of mismatches allowed within a seed.
	 *          Must be >= 0 and <= 2.  Note that 2-mismatch mode is
	 *          not fully sensitive; i.e. some 2-mismatch seed
	 *          alignments may be missed.
	 *   len  = Length of seed.
	 *   ival = Interval between seeds.  If not specified, seed
	 *          interval is determined by IVAL.
	 *
	 * Seed interval
	 * -------------
	 *
	 * IVAL={L|S|C},a,b (default: IVAL=S,1.0,0.0)
	 *
	 *   L  = let interval between seeds be a linear function of the
	 *        read length.  xx and yy are the constant and linear
	 *        coefficients respectively.  In other words, the interval
	 *        equals a * len + b, where len is the read length.
	 *        Intervals less than 1 are rounded up to 1.
	 *   S  = let interval between seeds be a function of the sqaure
	 *        root of the  read length.  xx and yy are the
	 *        coefficients.  In other words, the interval equals
	 *        a * sqrt(len) + b, where len is the read length.
	 *        Intervals less than 1 are rounded up to 1.
	 *   C  = Like S but uses cube root of length instead of square
	 *        root.
	 *
	 * Example 1:
	 *
	 *  SEED=1,10,5 and read sequence is TGCTATCGTACGATCGTAC:
	 *
	 *  The following seeds are extracted from the forward
	 *  representation of the read and aligned to the reference
	 *  allowing up to 1 mismatch:
	 *
	 *  Read:    TGCTATCGTACGATCGTACA
	 *
	 *  Seed 1+: TGCTATCGTA
	 *  Seed 2+:      TCGTACGATC
	 *  Seed 3+:           CGATCGTACA
	 *
	 *  ...and the following are extracted from the reverse-complement
	 *  representation of the read and align to the reference allowing
	 *  up to 1 mismatch:
	 *
	 *  Seed 1-: TACGATAGCA
	 *  Seed 2-:      GATCGTACGA
	 *  Seed 3-:           TGTACGATCG
	 *
	 * Example 2:
	 *
	 *  SEED=1,20,20 and read sequence is TGCTATCGTACGATC.  The seed
	 *  length is 20 but the read is only 15 characters long.  In this
	 *  case, Bowtie2 automatically shrinks the seed length to be equal
	 *  to the read length.
	 *
	 *  Read:    TGCTATCGTACGATC
	 *
	 *  Seed 1+: TGCTATCGTACGATC
	 *  Seed 1-: GATCGTACGATAGCA
	 *
	 * Example 3:
	 *
	 *  SEED=1,10,10 and read sequence is TGCTATCGTACGATC.  Only one seed
	 *  fits on the read; a second seed would overhang the end of the read
	 *  by 5 positions.  In this case, Bowtie2 extracts one seed.
	 *
	 *  Read:    TGCTATCGTACGATC
	 *
	 *  Seed 1+: TGCTATCGTA
	 *  Seed 1-: TACGATAGCA
	 */
	static void parseString(
		const std::string& s,
		bool noisyHpolymer,
		int& penMmcType,
		int& penMmc,
		int& penSnp,
		int& penNType,
		int& penN,
		int& penRdOpen,
		int& penRfOpen,
		int& penRdExConst,
		int& penRfExConst,
		int& penRdExLinear,
		int& penRfExLinear,
		float& costCeilConst,
		float& costCeilLinear,
		float& nCeilConst,
		float& nCeilLinear,
		int& multiseedMms,
		int& multiseedLen,
		int& multiseedPeriod,
		int& multiseedIvalType,
		float& multiseedIvalA,
		float& multiseedIvalB);
};

#endif /*ndef ALIGNER_SEED_POLICY_H_*/
