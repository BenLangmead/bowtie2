/*
 *  aligner_seed_policy
 *  bowtie-beta1
 *
 *  Created by Benjamin Langmead on 10/24/10.
 *  Copyright 2010 Johns Hopkins University. All rights reserved.
 *
 */

#include <string>
#include <iostream>
#include <sstream>
#include "ds.h"
#include "aligner_seed_policy.h"
#include "mem_ids.h"

using namespace std;

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
 * Penalty for a SNP in a colorspace alignment
 * -------------------------------------------
 *
 * SNP=xx
 *
 *    xx = Each nucleotide difference in the decoded alignment
 *         costs xx.
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
 * CEIL=xx,yy (default: CEIL=3.0,2.0)
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
void SeedAlignmentPolicy::parseString(
	const string& s,
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
	float& multiseedIvalB)
{
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
			cerr << "Error parsing alignment policy setting " << setting << "; must be bisected by = sign" << endl
				 << "Policy: " << s << endl;
			throw 1;
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
			cerr << "Error parsing alignment policy setting " << setting << "; RHS must have at least 1 token" << endl
				 << "Policy: " << s << endl;
			throw 1;
		}
		// In no case is >3 tokens OK
		if(ctoks.size() > 3) {
			cerr << "Error parsing alignment policy setting " << setting << "; RHS must have at most 3 tokens" << endl
				 << "Policy: " << s << endl;
			throw 1;
		}
		for(size_t i = 0; i < ctoks.size(); i++) {
			if(ctoks[i].length() == 0) {
				cerr << "Error parsing alignment policy setting " << setting << "; token " << i+1 << " on RHS had length=0" << endl
					 << "Policy: " << s << endl;
				throw 1;
			}
		}
		// Penalties for mismatches
		// MMP={Cxx|Q|RQ}
		//        Cxx = constant, where constant is integer xx
		//        Q   = equal to quality value
		//        R   = equal to maq-rounded quality value (rounded to
		//              nearest 10, can't be greater than 30).
		if(tag == "MMP") {
			if(ctoks.size() != 1) {
				cerr << "Error parsing alignment policy setting " << setting << "; RHS must have 1 token" << endl
					 << "Policy: " << s << endl;
				throw 1;
			}
			if(ctoks[0][0] == 'C') {
				string tmp = ctoks[0].substr(1);
				// Parse constant penalty
				istringstream tmpss(tmp);
				tmpss >> penMmc;
				// Set type to constant
				penMmcType = COST_MODEL_CONSTANT;
			}
			else if(ctoks[0][0] == 'Q') {
				// Set type to =quality
				penMmcType = COST_MODEL_QUAL;
			}
			else if(ctoks[0][0] == 'R') {
				// Set type to=Maq-quality
				penMmcType = COST_MODEL_ROUNDED_QUAL;
			}
		}
		// Penalties for SNPs in colorspace alignments
		// SNP=xx
		//        xx = penalty
		if(tag == "SNP") {
			if(ctoks.size() != 1) {
				cerr << "Error parsing alignment policy setting " << setting << "; RHS must have 1 token" << endl
					 << "Policy: " << s << endl;
				throw 1;
			}
			string tmp = ctoks[0];
			// Parse SNP penalty
			istringstream tmpss(tmp);
			tmpss >> penSnp;
		}
		// Penalties for mismatches where read char=N
		// NP={Cxx|Q|RQ}
		//        Cxx = constant, where constant is integer xx
		//        Q   = equal to quality
		//        R   = equal to maq-rounded quality value (rounded to nearest
		//              10, can't be greater than 30)
		else if(tag == "NP") {
			if(ctoks.size() != 1) {
				cerr << "Error parsing alignment policy setting " << setting << "; RHS must have 1 token" << endl
					 << "Policy: " << s << endl;
				throw 1;
			}
			if(ctoks[0][0] == 'C') {
				string tmp = ctoks[0].substr(1);
				// Parse constant penalty
				istringstream tmpss(tmp);
				tmpss >> penN;
				// Parse constant penalty
				penNType = COST_MODEL_CONSTANT;
			}
			else if(ctoks[0][0] == 'Q') {
				// Set type to =quality
				penNType = COST_MODEL_QUAL;
			}
			else if(ctoks[0][0] == 'R') {
				// Set type to=Maq-quality
				penNType = COST_MODEL_ROUNDED_QUAL;
			}
		}
		// Penalties for read gaps
		// RDG=xx,yy,zz
		//        xx = read gap open penalty
		//        yy = read gap extension penalty constant coefficient
		//             (defaults to open penalty)
		//        zz = read gap extension penalty linear coefficient
		//             (defaults to 0)
		else if(tag == "RDG") {
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> penRdOpen;
			} else {
				penRdOpen = noisyHpolymer ? DEFAULT_READ_OPEN_BADHPOLY : DEFAULT_READ_OPEN;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> penRdExConst;
			} else {
				penRdExConst = noisyHpolymer ? DEFAULT_READ_EXTEND_CONST_BADHPOLY : DEFAULT_READ_EXTEND_CONST;
			}
			if(ctoks.size() >= 3) {
				istringstream tmpss(ctoks[2]);
				tmpss >> penRdExLinear;
			} else {
				penRdExLinear = noisyHpolymer ? DEFAULT_READ_EXTEND_LINEAR_BADHPOLY : DEFAULT_READ_EXTEND_LINEAR;
			}
		}
		// Penalties for reference gaps
		// RFG=xx,yy,zz
		//        xx = ref gap open penalty
		//        yy = ref gap extension penalty constant coefficient
		//             (defaults to open penalty)
		//        zz = ref gap extension penalty linear coefficient
		//             (defaults to 0)
		else if(tag == "RFG") {
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> penRfOpen;
			} else {
				penRfOpen = noisyHpolymer ? DEFAULT_REF_OPEN_BADHPOLY : DEFAULT_REF_OPEN;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> penRfExConst;
			} else {
				penRfExConst = noisyHpolymer ? DEFAULT_REF_EXTEND_CONST_BADHPOLY : DEFAULT_REF_EXTEND_CONST;
			}
			if(ctoks.size() >= 3) {
				istringstream tmpss(ctoks[2]);
				tmpss >> penRfExLinear;
			} else {
				penRfExLinear = noisyHpolymer ? DEFAULT_REF_EXTEND_LINEAR_BADHPOLY : DEFAULT_REF_EXTEND_LINEAR;
			}
		}
		// Per-read penalty ceiling as a function of read length
		// CEIL=xx,yy
		//        xx = cost ceiling constant coefficient
		//        yy = cost ceiling linear coefficient (set to 0 if
		//             unspecified)
		else if(tag == "CEIL") {
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> costCeilConst;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> costCeilLinear;
			} else {
				costCeilLinear = DEFAULT_CEIL_LINEAR;
			}
		}
		// Per-read N ceiling as a function of read length
		// NCEIL=xx,yy
		//        xx = N ceiling constant coefficient
		//        yy = N ceiling linear coefficient (set to 0 if unspecified)
		else if(tag == "NCEIL") {
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> nCeilConst;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> nCeilLinear;
			} else {
				nCeilLinear = DEFAULT_N_CEIL_LINEAR;
			}
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
			if(ctoks.size() >= 1) {
				istringstream tmpss(ctoks[0]);
				tmpss >> multiseedMms;
			}
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> multiseedLen;
			} else {
				multiseedLen = DEFAULT_SEEDLEN;
			}
			if(ctoks.size() >= 3) {
				istringstream tmpss(ctoks[2]);
				tmpss >> multiseedPeriod;
			} else {
				multiseedPeriod = DEFAULT_SEEDPERIOD;
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
			if(ctoks.size() >= 1) {
				if(ctoks[0][0] == 'L') {
					multiseedIvalType = SEED_IVAL_LINEAR;
				} else if(ctoks[0][0] == 'S') {
					multiseedIvalType = SEED_IVAL_SQUARE_ROOT;
				} else if(ctoks[0][0] == 'C') {
					multiseedIvalType = SEED_IVAL_CUBE_ROOT;
				}
			}
			// A = Linear coefficient
			if(ctoks.size() >= 2) {
				istringstream tmpss(ctoks[1]);
				tmpss >> multiseedIvalA;
			} else {
				multiseedIvalA = 1.0f;
			}
			// B = Constant coefficient
			if(ctoks.size() >= 3) {
				istringstream tmpss(ctoks[2]);
				tmpss >> multiseedIvalB;
			} else {
				multiseedIvalB = 0.0f;
			}
		}
	}
}
