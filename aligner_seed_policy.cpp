/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include "ds.h"
#include "aligner_seed_policy.h"
#include "mem_ids.h"

using namespace std;

static int parseFuncType(const std::string& otype) {
	string type = otype;
	if(type == "C" || type == "Constant") {
		return SIMPLE_FUNC_CONST;
	} else if(type == "L" || type == "Linear") {
		return SIMPLE_FUNC_LINEAR;
	} else if(type == "S" || type == "Sqrt") {
		return SIMPLE_FUNC_SQRT;
	} else if(type == "G" || type == "Log") {
		return SIMPLE_FUNC_LOG;
	}
	std::cerr << "Error: Bad function type '" << otype
	          << "'.  Should be C (constant), L (linear), "
	          << "S (square root) or G (natural log)." << std::endl;
	throw 1;
}

#define PARSE_FUNC(fv) { \
	if(ctoks.size() >= 1) { \
		fv.setType(parseFuncType(ctoks[0])); \
	} \
	if(ctoks.size() >= 2) { \
		double co; \
		istringstream tmpss(ctoks[1]); \
		tmpss >> co; \
		fv.setConst(co); \
	} \
	if(ctoks.size() >= 3) { \
		double ce; \
		istringstream tmpss(ctoks[2]); \
		tmpss >> ce; \
		fv.setCoeff(ce); \
	} \
	if(ctoks.size() >= 4) { \
		double mn; \
		istringstream tmpss(ctoks[3]); \
		tmpss >> mn; \
		fv.setMin(mn); \
	} \
	if(ctoks.size() >= 5) { \
		double mx; \
		istringstream tmpss(ctoks[4]); \
		tmpss >> mx; \
		fv.setMin(mx); \
	} \
}

/**
 * Parse alignment policy when provided in this format:
 * <lab>=<val>;<lab>=<val>;<lab>=<val>...
 *
 * And label=value possibilities are:
 *
 * Bonus for a match
 * -----------------
 *
 * MA=xx (default: MA=0, or MA=2 if --local is set)
 *
 *    xx = Each position where equal read and reference characters match up
 *         in the alignment contriubtes this amount to the total score.
 *
 * Penalty for a mismatch
 * ----------------------
 *
 * MMP={Cxx|Q|RQ} (default: MMP=C6)
 *
 *   Cxx = Each mismatch costs xx.  If MMP=Cxx is specified, quality
 *         values are ignored when assessing penalities for mismatches.
 *   Q   = Each mismatch incurs a penalty equal to the mismatched base's
 *         value.
 *   R   = Each mismatch incurs a penalty equal to the mismatched base's
 *         rounded quality value.  Qualities are rounded off to the
 *         nearest 10, and qualities greater than 30 are rounded to 30.
 *
 * Penalty for position with N (in either read or reference)
 * ---------------------------------------------------------
 *
 * NP={Cxx|Q|RQ} (default: NP=C1)
 *
 *   Cxx = Each alignment position with an N in either the read or the
 *         reference costs xx.  If NP=Cxx is specified, quality values are
 *         ignored when assessing penalities for Ns.
 *   Q   = Each alignment position with an N in either the read or the
 *         reference incurs a penalty equal to the read base's quality
 *         value.
 *   R   = Each alignment position with an N in either the read or the
 *         reference incurs a penalty equal to the read base's rounded
 *         quality value.  Qualities are rounded off to the nearest 10,
 *         and qualities greater than 30 are rounded to 30.
 *
 * Penalty for a read gap
 * ----------------------
 *
 * RDG=xx,yy (default: RDG=5,3)
 *
 *   xx    = Read gap open penalty.
 *   yy    = Read gap extension penalty.
 *
 * Total cost incurred by a read gap = xx + (yy * gap length)
 *
 * Penalty for a reference gap
 * ---------------------------
 *
 * RFG=xx,yy (default: RFG=5,3)
 *
 *   xx    = Reference gap open penalty.
 *   yy    = Reference gap extension penalty.
 *
 * Total cost incurred by a reference gap = xx + (yy * gap length)
 *
 * Minimum score for valid alignment
 * ---------------------------------
 *
 * MIN=xx,yy (defaults: MIN=-0.6,-0.6, or MIN=0.0,0.66 if --local is set)
 *
 *   xx,yy = For a read of length N, the total score must be at least
 *           xx + (read length * yy) for the alignment to be valid.  The
 *           total score is the sum of all negative penalties (from
 *           mismatches and gaps) and all positive bonuses.  The minimum
 *           can be negative (and is by default in global alignment mode).
 *
 * Score floor for local alignment
 * -------------------------------
 *
 * FL=xx,yy (defaults: FL=-Infinity,0.0, or FL=0.0,0.0 if --local is set)
 *
 *   xx,yy = If a cell in the dynamic programming table has a score less
 *           than xx + (read length * yy), then no valid alignment can go
 *           through it.  Defaults are highly recommended.
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
 * IVAL={L|S|C},xx,yy (default: IVAL=S,1.0,0.0)
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
void SeedAlignmentPolicy::parseString(
	const       std::string& s,
	bool        local,
	bool        noisyHpolymer,
	bool        ignoreQuals,
	int&        bonusMatchType,
	int&        bonusMatch,
	int&        penMmcType,
	int&        penMmcMax,
	int&        penMmcMin,
	int&        penNType,
	int&        penN,
	int&        penRdExConst,
	int&        penRfExConst,
	int&        penRdExLinear,
	int&        penRfExLinear,
	SimpleFunc& costMin,
	SimpleFunc& nCeil,
	bool&       nCatPair,
	int&        multiseedMms,
	int&        multiseedLen,
	SimpleFunc& multiseedIval,
	size_t&     failStreak,
	size_t&     seedRounds)
{

	bonusMatchType    = local ? DEFAULT_MATCH_BONUS_TYPE_LOCAL : DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch        = local ? DEFAULT_MATCH_BONUS_LOCAL : DEFAULT_MATCH_BONUS;
	penMmcType        = ignoreQuals ? DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS :
	                                  DEFAULT_MM_PENALTY_TYPE;
	penMmcMax         = DEFAULT_MM_PENALTY_MAX;
	penMmcMin         = DEFAULT_MM_PENALTY_MIN;
	penNType          = DEFAULT_N_PENALTY_TYPE;
	penN              = DEFAULT_N_PENALTY;
	
	const double DMAX = std::numeric_limits<double>::max();
	costMin.init(
		local ? SIMPLE_FUNC_LOG : SIMPLE_FUNC_LINEAR,
		local ? DEFAULT_MIN_CONST_LOCAL  : DEFAULT_MIN_CONST,
		local ? DEFAULT_MIN_LINEAR_LOCAL : DEFAULT_MIN_LINEAR);
	nCeil.init(
		SIMPLE_FUNC_LINEAR, 0.0f, DMAX,
		DEFAULT_N_CEIL_CONST, DEFAULT_N_CEIL_LINEAR);
	multiseedIval.init(
		DEFAULT_IVAL, 1.0f, DMAX,
		DEFAULT_IVAL_B, DEFAULT_IVAL_A);
	nCatPair          = DEFAULT_N_CAT_PAIR;

	if(!noisyHpolymer) {
		penRdExConst  = DEFAULT_READ_GAP_CONST;
		penRdExLinear = DEFAULT_READ_GAP_LINEAR;
		penRfExConst  = DEFAULT_REF_GAP_CONST;
		penRfExLinear = DEFAULT_REF_GAP_LINEAR;
	} else {
		penRdExConst  = DEFAULT_READ_GAP_CONST_BADHPOLY;
		penRdExLinear = DEFAULT_READ_GAP_LINEAR_BADHPOLY;
		penRfExConst  = DEFAULT_REF_GAP_CONST_BADHPOLY;
		penRfExLinear = DEFAULT_REF_GAP_LINEAR_BADHPOLY;
	}
	
	multiseedMms      = DEFAULT_SEEDMMS;
	multiseedLen      = DEFAULT_SEEDLEN;
	
	EList<string> toks(MISC_CAT);
	string tok;
	istringstream ss(s);
	int setting = 0;
	// Get each ;-separated token
	while(getline(ss, tok, ';')) {
		setting++;
		EList<string> etoks(MISC_CAT);
		string etok;
		// Divide into tokens on either side of =
		istringstream ess(tok);
		while(getline(ess, etok, '=')) {
			etoks.push_back(etok);
		}
		// Must be exactly 1 =
		if(etoks.size() != 2) {
			cerr << "Error parsing alignment policy setting " << setting
			     << "; must be bisected by = sign" << endl
				 << "Policy: " << s << endl;
			assert(false); throw 1;
		}
		// LHS is tag, RHS value
		string tag = etoks[0], val = etoks[1];
		// Separate value into comma-separated tokens
		EList<string> ctoks(MISC_CAT);
		string ctok;
		istringstream css(val);
		while(getline(css, ctok, ',')) {
			ctoks.push_back(ctok);
		}
		if(ctoks.size() == 0) {
			cerr << "Error parsing alignment policy setting " << setting
			     << "; RHS must have at least 1 token" << endl
				 << "Policy: " << s << endl;
			assert(false); throw 1;
		}
		for(size_t i = 0; i < ctoks.size(); i++) {
			if(ctoks[i].length() == 0) {
				cerr << "Error parsing alignment policy setting " << setting
				     << "; token " << i+1 << " on RHS had length=0" << endl
					 << "Policy: " << s << endl;
				assert(false); throw 1;
			}
		}
		// Bonus for a match
		// MA=xx (default: MA=0, or MA=10 if --local is set)
		if(tag == "MA") {
			if(ctoks.size() != 1) {
				cerr << "Error parsing alignment policy setting " << setting
				     << "; RHS must have 1 token" << endl
					 << "Policy: " << s << endl;
				assert(false); throw 1;
			}
			string tmp = ctoks[0];
			istringstream tmpss(tmp);
			tmpss >> bonusMatch;
		}
		// Scoring for mismatches
		// MMP={Cxx|Q|RQ}
		//        Cxx = constant, where constant is integer xx
		//        Qxx = equal to quality, scaled
		//        R   = equal to maq-rounded quality value (rounded to nearest
		//              10, can't be greater than 30)
		else if(tag == "MMP") {
			if(ctoks.size() > 3) {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'"
				     << "; RHS must have at most 3 tokens" << endl
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
			if(ctoks[0][0] == 'C') {
				string tmp = ctoks[0].substr(1);
				// Parse constant penalty
				istringstream tmpss(tmp);
				tmpss >> penMmcMax;
				penMmcMin = penMmcMax;
				// Parse constant penalty
				penMmcType = COST_MODEL_CONSTANT;
			} else if(ctoks[0][0] == 'Q') {
				if(ctoks.size() >= 2) {
					string tmp = ctoks[1];
					istringstream tmpss(tmp);
					tmpss >> penMmcMax;
				} else {
					penMmcMax = DEFAULT_MM_PENALTY_MAX;
				}
				if(ctoks.size() >= 3) {
					string tmp = ctoks[2];
					istringstream tmpss(tmp);
					tmpss >> penMmcMin;
				} else {
					penMmcMin = DEFAULT_MM_PENALTY_MIN;
				}
				if(penMmcMin > penMmcMax) {
					cerr << "Error: Maximum mismatch penalty (" << penMmcMax
					     << ") is less than minimum penalty (" << penMmcMin
						 << endl;
					throw 1;
				}
				// Set type to =quality
				penMmcType = COST_MODEL_QUAL;
			} else if(ctoks[0][0] == 'R') {
				// Set type to=Maq-quality
				penMmcType = COST_MODEL_ROUNDED_QUAL;
			} else {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'"
				     << "; RHS must start with C, Q or R" << endl
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
		}
		// Scoring for mismatches where read char=N
		// NP={Cxx|Q|RQ}
		//        Cxx = constant, where constant is integer xx
		//        Q   = equal to quality
		//        R   = equal to maq-rounded quality value (rounded to nearest
		//              10, can't be greater than 30)
		else if(tag == "NP") {
			if(ctoks.size() != 1) {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'"
				     << "; RHS must have 1 token" << endl
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
			if(ctoks[0][0] == 'C') {
				string tmp = ctoks[0].substr(1);
				// Parse constant penalty
				istringstream tmpss(tmp);
				tmpss >> penN;
				// Parse constant penalty
				penNType = COST_MODEL_CONSTANT;
			} else if(ctoks[0][0] == 'Q') {
				// Set type to =quality
				penNType = COST_MODEL_QUAL;
			} else if(ctoks[0][0] == 'R') {
				// Set type to=Maq-quality
				penNType = COST_MODEL_ROUNDED_QUAL;
			} else {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'"
				     << "; RHS must start with C, Q or R" << endl
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
		}
		// Scoring for read gaps
		// RDG=xx,yy,zz
		//        xx = read gap open penalty
		//        yy = read gap extension penalty constant coefficient
		//             (defaults to open penalty)
		//        zz = read gap extension penalty linear coefficient
		//             (defaults to 0)
		else if(tag == "RDG") {
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> penRdExConst;
			} else {
				penRdExConst = noisyHpolymer ?
					DEFAULT_READ_GAP_CONST_BADHPOLY :
					DEFAULT_READ_GAP_CONST;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> penRdExLinear;
			} else {
				penRdExLinear = noisyHpolymer ?
					DEFAULT_READ_GAP_LINEAR_BADHPOLY :
					DEFAULT_READ_GAP_LINEAR;
			}
		}
		// Scoring for reference gaps
		// RFG=xx,yy,zz
		//        xx = ref gap open penalty
		//        yy = ref gap extension penalty constant coefficient
		//             (defaults to open penalty)
		//        zz = ref gap extension penalty linear coefficient
		//             (defaults to 0)
		else if(tag == "RFG") {
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> penRfExConst;
			} else {
				penRfExConst = noisyHpolymer ?
					DEFAULT_REF_GAP_CONST_BADHPOLY :
					DEFAULT_REF_GAP_CONST;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> penRfExLinear;
			} else {
				penRfExLinear = noisyHpolymer ?
					DEFAULT_REF_GAP_LINEAR_BADHPOLY :
					DEFAULT_REF_GAP_LINEAR;
			}
		}
		// Minimum score as a function of read length
		// MIN=xx,yy
		//        xx = constant coefficient
		//        yy = linear coefficient
		else if(tag == "MIN") {
			PARSE_FUNC(costMin);
		}
		// Per-read N ceiling as a function of read length
		// NCEIL=xx,yy
		//        xx = N ceiling constant coefficient
		//        yy = N ceiling linear coefficient (set to 0 if unspecified)
		else if(tag == "NCEIL") {
			PARSE_FUNC(nCeil);
		}
		/*
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
		 */
		else if(tag == "SEED") {
			if(ctoks.size() > 2) {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'; RHS must have 1 or 2 tokens, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> multiseedMms;
				if(multiseedMms > 1) {
					cerr << "Error: -N was set to " << multiseedMms << ", but cannot be set greater than 1" << endl;
					throw 1;
				}
				if(multiseedMms < 0) {
					cerr << "Error: -N was set to a number less than 0 (" << multiseedMms << ")" << endl;
					throw 1;
				}
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> multiseedLen;
			} else {
				multiseedLen = DEFAULT_SEEDLEN;
			}
		}
		else if(tag == "SEEDLEN") {
			if(ctoks.size() > 1) {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> multiseedLen;
			}
		}
		else if(tag == "DPS") {
			if(ctoks.size() > 1) {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> failStreak;
			}
		}
		else if(tag == "ROUNDS") {
			if(ctoks.size() > 1) {
				cerr << "Error parsing alignment policy setting "
				     << "'" << tag << "'; RHS must have 1 token, "
					 << "had " << ctoks.size() << ".  "
					 << "Policy: '" << s << "'" << endl;
				assert(false); throw 1;
			}
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> seedRounds;
			}
		}
		/*
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
		 */
		else if(tag == "IVAL") {
			PARSE_FUNC(multiseedIval);
		}
		else {
			// Unknown tag
			cerr << "Unexpected alignment policy setting "
				 << "'" << tag << "'" << endl
				 << "Policy: '" << s << "'" << endl;
			assert(false); throw 1;
		}
	}
}

#ifdef ALIGNER_SEED_POLICY_MAIN
int main() {

	int bonusMatchType;
	int bonusMatch;
	int penMmcType;
	int penMmc;
	int penNType;
	int penN;
	int penRdExConst;
	int penRfExConst;
	int penRdExLinear;
	int penRfExLinear;
	SimpleFunc costMin;
	SimpleFunc costFloor;
	SimpleFunc nCeil;
	bool nCatPair;
	int multiseedMms;
	int multiseedLen;
	SimpleFunc msIval;
	SimpleFunc posfrac;
	SimpleFunc rowmult;
	uint32_t mhits;

	{
		cout << "Case 1: Defaults 1 ... ";
		const char *pol = "";
		SeedAlignmentPolicy::parseString(
			string(pol),
			false,              // --local?
			false,              // noisy homopolymers a la 454?
			false,              // ignore qualities?
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		
		assert_eq(DEFAULT_MATCH_BONUS_TYPE,   bonusMatchType);
		assert_eq(DEFAULT_MATCH_BONUS,        bonusMatch);
		assert_eq(DEFAULT_MM_PENALTY_TYPE,    penMmcType);
		assert_eq(DEFAULT_MM_PENALTY_MAX,     penMmcMax);
		assert_eq(DEFAULT_MM_PENALTY_MIN,     penMmcMin);
		assert_eq(DEFAULT_N_PENALTY_TYPE,     penNType);
		assert_eq(DEFAULT_N_PENALTY,          penN);
		assert_eq(DEFAULT_MIN_CONST,          costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR,         costMin.getCoeff());
		assert_eq(DEFAULT_FLOOR_CONST,        costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR,       costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST,       nCeil.getConst());
		assert_eq(DEFAULT_N_CAT_PAIR,         nCatPair);

		assert_eq(DEFAULT_READ_GAP_CONST,     penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR,    penRdExLinear);
		assert_eq(DEFAULT_REF_GAP_CONST,      penRfExConst);
		assert_eq(DEFAULT_REF_GAP_LINEAR,     penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS,            multiseedMms);
		assert_eq(DEFAULT_SEEDLEN,            multiseedLen);
		assert_eq(DEFAULT_IVAL,               msIval.getType());
		assert_eq(DEFAULT_IVAL_A,             msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B,             msIval.getConst());
		
		cout << "PASSED" << endl;
	}

	{
		cout << "Case 2: Defaults 2 ... ";
		const char *pol = "";
		SeedAlignmentPolicy::parseString(
			string(pol),
			false,              // --local?
			true,               // noisy homopolymers a la 454?
			false,              // ignore qualities?
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		
		assert_eq(DEFAULT_MATCH_BONUS_TYPE,   bonusMatchType);
		assert_eq(DEFAULT_MATCH_BONUS,        bonusMatch);
		assert_eq(DEFAULT_MM_PENALTY_TYPE,    penMmcType);
		assert_eq(DEFAULT_MM_PENALTY_MAX,     penMmc);
		assert_eq(DEFAULT_MM_PENALTY_MIN,     penMmc);
		assert_eq(DEFAULT_N_PENALTY_TYPE,     penNType);
		assert_eq(DEFAULT_N_PENALTY,          penN);
		assert_eq(DEFAULT_MIN_CONST,          costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR,         costMin.getCoeff());
		assert_eq(DEFAULT_FLOOR_CONST,        costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR,       costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST,       nCeil.getConst());
		assert_eq(DEFAULT_N_CAT_PAIR,         nCatPair);

		assert_eq(DEFAULT_READ_GAP_CONST_BADHPOLY,  penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR_BADHPOLY, penRdExLinear);
		assert_eq(DEFAULT_REF_GAP_CONST_BADHPOLY,   penRfExConst);
		assert_eq(DEFAULT_REF_GAP_LINEAR_BADHPOLY,  penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS,            multiseedMms);
		assert_eq(DEFAULT_SEEDLEN,            multiseedLen);
		assert_eq(DEFAULT_IVAL,               msIval.getType());
		assert_eq(DEFAULT_IVAL_A,             msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B,             msIval.getConst());
		
		cout << "PASSED" << endl;
	}

	{
		cout << "Case 3: Defaults 3 ... ";
		const char *pol = "";
		SeedAlignmentPolicy::parseString(
			string(pol),
			true,               // --local?
			false,              // noisy homopolymers a la 454?
			false,              // ignore qualities?
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		
		assert_eq(DEFAULT_MATCH_BONUS_TYPE_LOCAL,   bonusMatchType);
		assert_eq(DEFAULT_MATCH_BONUS_LOCAL,        bonusMatch);
		assert_eq(DEFAULT_MM_PENALTY_TYPE,    penMmcType);
		assert_eq(DEFAULT_MM_PENALTY_MAX,     penMmcMax);
		assert_eq(DEFAULT_MM_PENALTY_MIN,     penMmcMin);
		assert_eq(DEFAULT_N_PENALTY_TYPE,     penNType);
		assert_eq(DEFAULT_N_PENALTY,          penN);
		assert_eq(DEFAULT_MIN_CONST_LOCAL,    costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR_LOCAL,   costMin.getCoeff());
		assert_eq(DEFAULT_FLOOR_CONST_LOCAL,  costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR_LOCAL, costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST,       nCeil.getConst());
		assert_eq(DEFAULT_N_CEIL_LINEAR,      nCeil.getCoeff());
		assert_eq(DEFAULT_N_CAT_PAIR,         nCatPair);

		assert_eq(DEFAULT_READ_GAP_CONST,     penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR,    penRdExLinear);
		assert_eq(DEFAULT_REF_GAP_CONST,      penRfExConst);
		assert_eq(DEFAULT_REF_GAP_LINEAR,     penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS,            multiseedMms);
		assert_eq(DEFAULT_SEEDLEN,            multiseedLen);
		assert_eq(DEFAULT_IVAL,               msIval.getType());
		assert_eq(DEFAULT_IVAL_A,             msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B,             msIval.getConst());

		cout << "PASSED" << endl;
	}

	{
		cout << "Case 4: Simple string 1 ... ";
		const char *pol = "MMP=C44;MA=4;RFG=24,12;FL=C,8;RDG=2;NP=C4;MIN=C,7";
		SeedAlignmentPolicy::parseString(
			string(pol),
			true,               // --local?
			false,              // noisy homopolymers a la 454?
			false,              // ignore qualities?
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmc,
			penNType,
			penN,
			penRdExConst,
			penRfExConst,
			penRdExLinear,
			penRfExLinear,
			costMin,
			costFloor,
			nCeil,
			nCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			mhits);
		
		assert_eq(COST_MODEL_CONSTANT,        bonusMatchType);
		assert_eq(4,                          bonusMatch);
		assert_eq(COST_MODEL_CONSTANT,        penMmcType);
		assert_eq(44,                         penMmc);
		assert_eq(COST_MODEL_CONSTANT,        penNType);
		assert_eq(4.0f,                       penN);
		assert_eq(7,                          costMin.getConst());
		assert_eq(DEFAULT_MIN_LINEAR_LOCAL,   costMin.getCoeff());
		assert_eq(8,                          costFloor.getConst());
		assert_eq(DEFAULT_FLOOR_LINEAR_LOCAL, costFloor.getCoeff());
		assert_eq(DEFAULT_N_CEIL_CONST,       nCeil.getConst());
		assert_eq(DEFAULT_N_CEIL_LINEAR,      nCeil.getCoeff());
		assert_eq(DEFAULT_N_CAT_PAIR,         nCatPair);

		assert_eq(2.0f,                       penRdExConst);
		assert_eq(DEFAULT_READ_GAP_LINEAR,    penRdExLinear);
		assert_eq(24.0f,                      penRfExConst);
		assert_eq(12.0f,                      penRfExLinear);
		assert_eq(DEFAULT_SEEDMMS,            multiseedMms);
		assert_eq(DEFAULT_SEEDLEN,            multiseedLen);
		assert_eq(DEFAULT_IVAL,               msIval.getType());
		assert_eq(DEFAULT_IVAL_A,             msIval.getCoeff());
		assert_eq(DEFAULT_IVAL_B,             msIval.getConst());

		cout << "PASSED" << endl;
	}
}
#endif /*def ALIGNER_SEED_POLICY_MAIN*/
