/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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

#include <limits>
// -- BTL remove --
//#include <stdlib.h>
//#include <sys/time.h>
// -- --
#include "aligner_sw.h"
#include "aligner_result.h"
#include "search_globals.h"
#include "scoring.h"
#include "mask.h"

/**
 * Initialize with a new read.
 */
void SwAligner::initRead(
	const BTDnaString& rdfw, // forward read sequence
	const BTDnaString& rdrc, // revcomp read sequence
	const BTString& qufw,    // forward read qualities
	const BTString& qurc,    // reverse read qualities
	size_t rdi,              // offset of first read char to align
	size_t rdf,              // offset of last read char to align
	const Scoring& sc)       // scoring scheme
{
	assert_gt(rdf, rdi);
	int nceil = sc.nCeil.f<int>((double)rdfw.length());
	rdfw_    = &rdfw;      // read sequence
	rdrc_    = &rdrc;      // read sequence
	qufw_    = &qufw;      // read qualities
	qurc_    = &qurc;      // read qualities
	rdi_     = rdi;        // offset of first read char to align
	rdf_     = rdf;        // offset of last read char to align
	sc_      = &sc;        // scoring scheme
	nceil_   = nceil;      // max # Ns allowed in ref portion of aln
	readSse16_ = false;    // true -> sse16 from now on for this read
	initedRead_ = true;
#ifndef NO_SSE
	sseU8fwBuilt_  = false;  // built fw query profile, 8-bit score
	sseU8rcBuilt_  = false;  // built rc query profile, 8-bit score
	sseI16fwBuilt_ = false;  // built fw query profile, 16-bit score
	sseI16rcBuilt_ = false;  // built rc query profile, 16-bit score
#endif
}

/**
 * Initialize with a new alignment problem.
 */
void SwAligner::initRef(
	bool fw,               // whether to forward or revcomp read is aligning
	TRefId refidx,         // id of reference aligned against
	const DPRect& rect,    // DP rectangle
	char *rf,              // reference sequence
	size_t rfi,            // offset of first reference char to align to
	size_t rff,            // offset of last reference char to align to
	TRefOff reflen,        // length of reference sequence
	const Scoring& sc,     // scoring scheme
	TAlScore minsc,        // minimum score
	bool enable8,          // use 8-bit SSE if possible?
	size_t cminlen,        // minimum length for using checkpointing scheme
	size_t cpow2,          // interval b/t checkpointed diags; 1 << this
	bool doTri,            // triangular mini-fills?
	bool extend)           // is this a seed extension?
{
	size_t readGaps = sc.maxReadGaps(minsc, rdfw_->length());
	size_t refGaps  = sc.maxRefGaps(minsc, rdfw_->length());
	assert_geq(readGaps, 0);
	assert_geq(refGaps, 0);
	assert_gt(rff, rfi);
	rdgap_       = readGaps;  // max # gaps in read
	rfgap_       = refGaps;   // max # gaps in reference
	state_       = STATE_INITED;
	fw_          = fw;       // orientation
	rd_          = fw ? rdfw_ : rdrc_; // read sequence
	qu_          = fw ? qufw_ : qurc_; // quality sequence
	refidx_      = refidx;   // id of reference aligned against
	rf_          = rf;       // reference sequence
	rfi_         = rfi;      // offset of first reference char to align to
	rff_         = rff;      // offset of last reference char to align to
	reflen_      = reflen;   // length of entire reference sequence
	rect_        = &rect;    // DP rectangle
	minsc_       = minsc;    // minimum score
	cural_       = 0;        // idx of next alignment to give out
	initedRef_   = true;     // indicate we've initialized the ref portion
	enable8_     = enable8;  // use 8-bit SSE if possible?
	extend_      = extend;   // true iff this is a seed extension
	cperMinlen_  = cminlen;  // reads shorter than this won't use checkpointer
	cperPerPow2_ = cpow2;    // interval b/t checkpointed diags; 1 << this
	cperEf_      = true;     // whether to checkpoint H, E, and F
	cperTri_     = doTri;    // triangular mini-fills?
	bter_.initRef(
		fw_ ? rdfw_->buf() : // in: read sequence
			  rdrc_->buf(), 
		fw_ ? qufw_->buf() : // in: quality sequence
			  qurc_->buf(),
		rd_->length(),       // in: read sequence length
		rf_ + rfi_,          // in: reference sequence
		rff_ - rfi_,         // in: in-rectangle reference sequence length
		reflen,              // in: total reference sequence length
		refidx_,             // in: reference id
		rfi_,                // in: reference offset
		fw_,                 // in: orientation
		rect_,               // in: DP rectangle
		&cper_,              // in: checkpointer
		*sc_,                // in: scoring scheme
		nceil_);             // in: N ceiling
}
	
/**
 * Given a read, an alignment orientation, a range of characters in a referece
 * sequence, and a bit-encoded version of the reference, set up and execute the
 * corresponding dynamic programming problem.
 *
 * The caller has already narrowed down the relevant portion of the reference
 * using, e.g., the location of a seed hit, or the range of possible fragment
 * lengths if we're searching for the opposite mate in a pair.
 */
void SwAligner::initRef(
	bool fw,               // whether to forward or revcomp read is aligning
	TRefId refidx,         // reference aligned against
	const DPRect& rect,    // DP rectangle
	const BitPairReference& refs, // Reference strings
	TRefOff reflen,        // length of reference sequence
	const Scoring& sc,     // scoring scheme
	TAlScore minsc,        // minimum score
	bool enable8,          // use 8-bit SSE if possible?
	size_t cminlen,        // minimum length for using checkpointing scheme
	size_t cpow2,          // interval b/t checkpointed diags; 1 << this
	bool doTri,            // triangular mini-fills?
	bool extend,           // true iff this is a seed extension
	size_t  upto,          // count the number of Ns up to this offset
	size_t& nsUpto)        // output: the number of Ns up to 'upto'
{
	TRefOff rfi = rect.refl;
	TRefOff rff = rect.refr + 1;
	assert_gt(rff, rfi);
	// Capture an extra reference character outside the rectangle so that we
	// can check matches in the next column over to the right
	rff++;
	// rflen = full length of the reference substring to consider, including
	// overhang off the boundaries of the reference sequence
	const size_t rflen = (size_t)(rff - rfi);
	// Figure the number of Ns we're going to add to either side
	size_t leftNs  =
		(rfi >= 0               ? 0 : (size_t)std::abs(static_cast<long>(rfi)));
	leftNs = min(leftNs, rflen);
	size_t rightNs =
		(rff <= (TRefOff)reflen ? 0 : (size_t)std::abs(static_cast<long>(rff - reflen)));
	rightNs = min(rightNs, rflen);
	// rflenInner = length of just the portion that doesn't overhang ref ends
	assert_geq(rflen, leftNs + rightNs);
	const size_t rflenInner = rflen - (leftNs + rightNs);
#ifndef NDEBUG
	bool haveRfbuf2 = false;
	EList<char> rfbuf2(rflen);
	// This is really slow, so only do it some of the time
	if((rand() % 10) == 0) {
		TRefOff rfii = rfi;
		for(size_t i = 0; i < rflen; i++) {
			if(rfii < 0 || (TRefOff)rfii >= reflen) {
				rfbuf2.push_back(4);
			} else {
				rfbuf2.push_back(refs.getBase(refidx, (uint32_t)rfii));
			}
			rfii++;
		}
		haveRfbuf2 = true;
	}
#endif
	// rfbuf_ = uint32_t list large enough to accommodate both the reference
	// sequence and any Ns we might add to either side.
	rfwbuf_.resize((rflen + 16) / 4);
	int offset = refs.getStretch(
		rfwbuf_.ptr(),               // buffer to store words in
		refidx,                      // which reference
		(rfi < 0) ? 0 : (size_t)rfi, // starting offset (can't be < 0)
		rflenInner                   // length to grab (exclude overhang)
		ASSERT_ONLY(, tmp_destU32_));// for BitPairReference::getStretch()
	assert_leq(offset, 16);
	rf_ = (char*)rfwbuf_.ptr() + offset;
	// Shift ref chars away from 0 so we can stick Ns at the beginning
	if(leftNs > 0) {
		// Slide everyone down
		for(size_t i = rflenInner; i > 0; i--) {
			rf_[i+leftNs-1] = rf_[i-1];
		}
		// Add Ns
		for(size_t i = 0; i < leftNs; i++) {
			rf_[i] = 4;
		}
	}
	if(rightNs > 0) {
		// Add Ns to the end
		for(size_t i = 0; i < rightNs; i++) {
			rf_[i + leftNs + rflenInner] = 4;
		}
	}
#ifndef NDEBUG
	// Sanity check reference characters
	for(size_t i = 0; i < rflen; i++) {
		assert(!haveRfbuf2 || rf_[i] == rfbuf2[i]);
		assert_range(0, 4, (int)rf_[i]);
	}
#endif
	// Count Ns and convert reference characters into A/C/G/T masks.  Ambiguous
	// nucleotides (IUPAC codes) have more than one mask bit set.  If a
	// reference scanner was provided, use it to opportunistically resolve seed
	// hits.
	nsUpto = 0;
	for(size_t i = 0; i < rflen; i++) {
		// rf_[i] gets mask version of refence char, with N=16
		if(i < upto && rf_[i] > 3) {
			nsUpto++;
		}
		rf_[i] = (1 << rf_[i]);
	}
	// Correct for having captured an extra reference character
	rff--;
	initRef(
		fw,          // whether to forward or revcomp read is aligning
		refidx,      // id of reference aligned against
		rect,        // DP rectangle
		rf_,         // reference sequence, wrapped up in BTString object
		0,           // use the whole thing
		(size_t)(rff - rfi), // ditto
		reflen,      // reference length
		sc,          // scoring scheme
		minsc,       // minimum score
		enable8,     // use 8-bit SSE if possible?
		cminlen,     // minimum length for using checkpointing scheme
		cpow2,       // interval b/t checkpointed diags; 1 << this
		doTri,       // triangular mini-fills?
		extend);     // true iff this is a seed extension
}

/**
 * Given a read, an alignment orientation, a range of characters in a referece
 * sequence, and a bit-encoded version of the reference, set up and execute the
 * corresponding ungapped alignment problem.  There can only be one solution.
 *
 * The caller has already narrowed down the relevant portion of the reference
 * using, e.g., the location of a seed hit, or the range of possible fragment
 * lengths if we're searching for the opposite mate in a pair.
 */
int SwAligner::ungappedAlign(
	const BTDnaString&      rd,     // read sequence (could be RC)
	const BTString&         qu,     // qual sequence (could be rev)
	const Coord&            coord,  // coordinate aligned to
	const BitPairReference& refs,   // Reference strings
	size_t                  reflen, // length of reference sequence
	const Scoring&          sc,     // scoring scheme
	bool                    ohang,  // allow overhang?
	TAlScore                minsc,  // minimum score
	SwResult&               res)    // put alignment result here
{
	const size_t len = rd.length();
	int nceil = sc.nCeil.f<int>((double)len);
	int ns = 0;
	TRefOff rfi = coord.off();
	TRefOff rff = rfi + (TRefOff)len;
	TRefId refidx = coord.ref();
	assert_gt(rff, rfi);
	// Figure the number of Ns we're going to add to either side
	size_t leftNs = 0;
	if(rfi < 0) {
		if(ohang) {
			leftNs = (size_t)(-rfi);
		} else {
			return 0;
		}
	}
	size_t rightNs = 0;
	if(rff > (TRefOff)reflen) {
		if(ohang) {
			rightNs = (size_t)(rff - (TRefOff)reflen);
		} else {
			return 0;
		}
	}
	if((leftNs + rightNs) > (size_t)nceil) {
		return 0;
	}
	// rflenInner = length of just the portion that doesn't overhang ref ends
	assert_geq(len, leftNs + rightNs);
	const size_t rflenInner = len - (leftNs + rightNs);
#ifndef NDEBUG
	bool haveRfbuf2 = false;
	EList<char> rfbuf2(len);
	// This is really slow, so only do it some of the time
	if((rand() % 10) == 0) {
		TRefOff rfii = rfi;
		for(size_t i = 0; i < len; i++) {
			if(rfii < 0 || (size_t)rfii >= reflen) {
				rfbuf2.push_back(4);
			} else {
				rfbuf2.push_back(refs.getBase(refidx, (uint32_t)rfii));
			}
			rfii++;
		}
		haveRfbuf2 = true;
	}
#endif
	// rfbuf_ = uint32_t list large enough to accommodate both the reference
	// sequence and any Ns we might add to either side.
	rfwbuf_.resize((len + 16) / 4);
	int offset = refs.getStretch(
		rfwbuf_.ptr(),               // buffer to store words in
		refidx,                      // which reference
		(rfi < 0) ? 0 : (size_t)rfi, // starting offset (can't be < 0)
		rflenInner                   // length to grab (exclude overhang)
		ASSERT_ONLY(, tmp_destU32_));// for BitPairReference::getStretch()
	assert_leq(offset, 16);
	rf_ = (char*)rfwbuf_.ptr() + offset;
	// Shift ref chars away from 0 so we can stick Ns at the beginning
	if(leftNs > 0) {
		// Slide everyone down
		for(size_t i = rflenInner; i > 0; i--) {
			rf_[i+leftNs-1] = rf_[i-1];
		}
		// Add Ns
		for(size_t i = 0; i < leftNs; i++) {
			rf_[i] = 4;
		}
	}
	if(rightNs > 0) {
		// Add Ns to the end
		for(size_t i = 0; i < rightNs; i++) {
			rf_[i + leftNs + rflenInner] = 4;
		}
	}
#ifndef NDEBUG
	// Sanity check reference characters
	for(size_t i = 0; i < len; i++) {
		assert(!haveRfbuf2 || rf_[i] == rfbuf2[i]);
		assert_range(0, 4, (int)rf_[i]);
	}
#endif
	// Count Ns and convert reference characters into A/C/G/T masks.  Ambiguous
	// nucleotides (IUPAC codes) have more than one mask bit set.  If a
	// reference scanner was provided, use it to opportunistically resolve seed
	// hits.
	TAlScore score = 0;
	res.alres.reset();
	size_t rowi = 0;
	size_t rowf = len-1;
	if(sc.monotone) {
		for(size_t i = 0; i < len; i++) {
			// rf_[i] gets mask version of refence char, with N=16
			assert_geq(qu[i], 33);
			score += sc.score(rd[i], (int)(1 << rf_[i]), qu[i] - 33, ns);
			assert_leq(score, 0);
			if(score < minsc || ns > nceil) {
				// Fell below threshold
				return 0;
			}
		}
		// Got a result!  Fill in the rest of the result object.
	} else {
		// Definitely ways to short-circuit this.  E.g. if diff between cur
		// score and minsc can't be met by matches.
		TAlScore floorsc = 0;
		TAlScore scoreMax = floorsc;
		size_t lastfloor = 0;
		rowi = MAX_SIZE_T;
		size_t sols = 0;
		for(size_t i = 0; i < len; i++) {
			score += sc.score(rd[i], (int)(1 << rf_[i]), qu[i] - 33, ns);
			if(score >= minsc && score >= scoreMax) {
				scoreMax = score;
				rowf = i;
				if(rowi != lastfloor) {
					rowi = lastfloor;
					sols++;
				}
			}
			if(score <= floorsc) {
				score = floorsc;
				lastfloor = i+1;
			}
		}
		if(ns > nceil || scoreMax < minsc) {
			// Too many Ns
			return 0;
		}
		if(sols > 1) {
			// >1 distinct solution in this diag; defer to DP aligner
			return -1;
		}
		score = scoreMax;
		// Got a result!  Fill in the rest of the result object.  
	}
	// Now fill in the edits
	res.alres.setScore(AlnScore(score, ns, 0));
	assert_geq(rowf, rowi);
	EList<Edit>& ned = res.alres.ned();
	size_t refns = 0;
	ASSERT_ONLY(BTDnaString refstr);
	for(size_t i = rowi; i <= rowf; i++) {
		ASSERT_ONLY(refstr.append((int)rf_[i]));
		if(rf_[i] > 3 || rd[i] != rf_[i]) {
			// Add edit
			Edit e((int)i,
			       mask2dna[1 << (int)rf_[i]],
			       "ACGTN"[(int)rd[i]],
			       EDIT_TYPE_MM);
			ned.push_back(e);
			if(rf_[i] > 3) {
				refns++;
			}
		}
	}
	assert(Edit::repOk(ned, rd));
	bool fw = coord.fw();
	assert_leq(rowf, len-1);
	size_t trimEnd = (len-1) - rowf;
	res.alres.setShape(
		coord.ref(),  // ref id
		coord.off()+rowi, // 0-based ref offset
		reflen,       // length of reference sequence aligned to
		fw,           // aligned to Watson?
		len,          // read length
		true,         // pretrim soft?
		0,            // pretrim 5' end
		0,            // pretrim 3' end
		true,         // alignment trim soft?
		fw ? rowi : trimEnd,  // alignment trim 5' end
		fw ? trimEnd : rowi); // alignment trim 3' end
	res.alres.setRefNs(refns);
	assert(res.repOk());
#ifndef NDEBUG
	BTDnaString editstr;
	Edit::toRef(rd, ned, editstr, true, rowi, trimEnd);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << res.alres.score().score() << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << rd << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	if(!fw) {
		// All edits are currently w/r/t upstream end; if read aligned to Crick
		// strand, invert them to be w/r/t 5' end instead.
		res.alres.invertEdits();
	}
	return 1;
}

/**
 * Align read 'rd' to reference using read & reference information given
 * last time init() was called.
 */
bool SwAligner::align(
	RandomSource& rnd, // source of pseudo-randoms
	TAlScore& best)    // best alignment score observed in DP matrix
{
	assert(initedRef() && initedRead());
	assert_eq(STATE_INITED, state_);
	state_ = STATE_ALIGNED;
	// Reset solutions lists
	btncand_.clear();
	btncanddone_.clear();
	btncanddoneSucc_ = btncanddoneFail_ = 0;
	best = std::numeric_limits<TAlScore>::min();
	sse8succ_ = sse16succ_ = false;
	int flag = 0;
	size_t rdlen = rdf_ - rdi_;
	bool checkpointed = rdlen >= cperMinlen_;
	bool gathered = false; // Did gathering happen along with alignment?
	if(sc_->monotone) {
		// End-to-end
		if(enable8_ && !readSse16_ && minsc_ >= -254) {
			// 8-bit end-to-end
			if(checkpointed) {
				best = alignGatherEE8(flag, false);
				if(flag == 0) {
					gathered = true;
				}
			} else {
				best = alignNucleotidesEnd2EndSseU8(flag, false);
#ifndef NDEBUG
				int flagtmp = 0;
				TAlScore besttmp = alignGatherEE8(flagtmp, true); // debug
				assert_eq(flagtmp, flag);
				assert_eq(besttmp, best);
#endif
			}
			sse8succ_ = (flag == 0);
#ifndef NDEBUG
			{
				int flag2 = 0;
				TAlScore best2 = alignNucleotidesEnd2EndSseI16(flag2, true);
				{
					int flagtmp = 0;
					TAlScore besttmp = alignGatherEE16(flagtmp, true);
					assert_eq(flagtmp, flag2);
					assert(flag2 != 0 || best2 == besttmp);
				}
				assert(flag < 0 || best == best2);
				sse16succ_ = (flag2 == 0);
			}
#endif /*ndef NDEBUG*/
		} else {
			// 16-bit end-to-end
			if(checkpointed) {
				best = alignGatherEE16(flag, false);
				if(flag == 0) {
					gathered = true;
				}
			} else {
				best = alignNucleotidesEnd2EndSseI16(flag, false);
#ifndef NDEBUG
				int flagtmp = 0;
				TAlScore besttmp = alignGatherEE16(flagtmp, true);
				assert_eq(flagtmp, flag);
				assert_eq(besttmp, best);
#endif
			}
			sse16succ_ = (flag == 0);
		}
	} else {
		// Local
		flag = -2;
		if(enable8_ && !readSse16_) {
			// 8-bit local
			if(checkpointed) {
				best = alignGatherLoc8(flag, false);
				if(flag == 0) {
					gathered = true;
				}
			} else {
				best = alignNucleotidesLocalSseU8(flag, false);
#ifndef NDEBUG
				int flagtmp = 0;
				TAlScore besttmp = alignGatherLoc8(flagtmp, true);
				assert_eq(flag, flagtmp);
				assert_eq(best, besttmp);
#endif
			}
		}
		if(flag == -2) {
			// 16-bit local
			flag = 0;
			if(checkpointed) {
				best = alignNucleotidesLocalSseI16(flag, false);
				best = alignGatherLoc16(flag, false);
				if(flag == 0) {
					gathered = true;
				}
			} else {
				best = alignNucleotidesLocalSseI16(flag, false);
#ifndef NDEBUG
				int flagtmp = 0;
				TAlScore besttmp = alignGatherLoc16(flagtmp, true);
				assert_eq(flag, flagtmp);
				assert_eq(best, besttmp);
#endif
			}
			sse16succ_ = (flag == 0);
		} else {
			sse8succ_ = (flag == 0);
#ifndef NDEBUG
			int flag2 = 0;
			TAlScore best2 = alignNucleotidesLocalSseI16(flag2, true);
			{
				int flagtmp = 0;
				TAlScore besttmp = alignGatherLoc16(flagtmp, true);
				assert_eq(flag2, flagtmp);
				assert(flag2 != 0 || best2 == besttmp);
			}
			assert(flag2 < 0 || best == best2);
			sse16succ_ = (flag2 == 0);
#endif /*ndef NDEBUG*/
		}
	}
#ifndef NDEBUG
	if(!checkpointed && (rand() & 15) == 0 && sse8succ_ && sse16succ_) {
		SSEData& d8  = fw_ ? sseU8fw_  : sseU8rc_;
		SSEData& d16 = fw_ ? sseI16fw_ : sseI16rc_;
		assert_eq(d8.mat_.nrow(), d16.mat_.nrow());
		assert_eq(d8.mat_.ncol(), d16.mat_.ncol());
		for(size_t i = 0; i < d8.mat_.nrow(); i++) {
			for(size_t j = 0; j < colstop_; j++) {
				int h8  = d8.mat_.helt(i, j);
				int h16 = d16.mat_.helt(i, j);
				int e8  = d8.mat_.eelt(i, j);
				int e16 = d16.mat_.eelt(i, j);
				int f8  = d8.mat_.felt(i, j);
				int f16 = d16.mat_.felt(i, j);
				TAlScore h8s  =
					(sc_->monotone ? (h8  - 0xff  ) : h8);
				TAlScore h16s =
					(sc_->monotone ? (h16 - 0x7fff) : (h16 + 0x8000));
				TAlScore e8s  =
					(sc_->monotone ? (e8  - 0xff  ) : e8);
				TAlScore e16s =
					(sc_->monotone ? (e16 - 0x7fff) : (e16 + 0x8000));
				TAlScore f8s  =
					(sc_->monotone ? (f8  - 0xff  ) : f8);
				TAlScore f16s =
					(sc_->monotone ? (f16 - 0x7fff) : (f16 + 0x8000));
				if(h8s < minsc_) {
					h8s = minsc_ - 1;
				}
				if(h16s < minsc_) {
					h16s = minsc_ - 1;
				}
				if(e8s < minsc_) {
					e8s = minsc_ - 1;
				}
				if(e16s < minsc_) {
					e16s = minsc_ - 1;
				}
				if(f8s < minsc_) {
					f8s = minsc_ - 1;
				}
				if(f16s < minsc_) {
					f16s = minsc_ - 1;
				}
				if((h8 != 0 || (int16_t)h16 != (int16_t)0x8000) && h8 > 0) {
					assert_eq(h8s, h16s);
				}
				if((e8 != 0 || (int16_t)e16 != (int16_t)0x8000) && e8 > 0) {
					assert_eq(e8s, e16s);
				}
				if((f8 != 0 || (int16_t)f16 != (int16_t)0x8000) && f8 > 0) {
					assert_eq(f8s, f16s);
				}
			}
		}
	}
#endif
	assert(repOk());
	cural_ = 0;
	if(best == MIN_I64 || best < minsc_) {
		return false;
	}
	if(!gathered) {
		// Look for solutions using SSE matrix
		assert(sse8succ_ || sse16succ_);
		if(sc_->monotone) {
			if(sse8succ_) {
				gatherCellsNucleotidesEnd2EndSseU8(best);
#ifndef NDEBUG
				if(sse16succ_) {
					cand_tmp_ = btncand_;
					gatherCellsNucleotidesEnd2EndSseI16(best);
					cand_tmp_.sort();
					btncand_.sort();
					assert(cand_tmp_ == btncand_);
				}
#endif /*ndef NDEBUG*/
			} else {
				gatherCellsNucleotidesEnd2EndSseI16(best);
			}
		} else {
			if(sse8succ_) {
				gatherCellsNucleotidesLocalSseU8(best);
#ifndef NDEBUG
				if(sse16succ_) {
					cand_tmp_ = btncand_;
					gatherCellsNucleotidesLocalSseI16(best);
					cand_tmp_.sort();
					btncand_.sort();
					assert(cand_tmp_ == btncand_);
				}
#endif /*ndef NDEBUG*/
			} else {
				gatherCellsNucleotidesLocalSseI16(best);
			}
		}
	}
	if(!btncand_.empty()) {
		btncand_.sort();
	}
	return !btncand_.empty();
}

/**
 * Populate the given SwResult with information about the "next best"
 * alignment if there is one.  If there isn't one, false is returned.  Note
 * that false might be returned even though a call to done() would have
 * returned false.
 */
bool SwAligner::nextAlignment(
	SwResult& res,
	TAlScore minsc,
	RandomSource& rnd)
{
	assert(initedRead() && initedRef());
	assert_eq(STATE_ALIGNED, state_);
	assert(repOk());
	if(done()) {
		res.reset();
		return false;
	}
	assert(!done());
	size_t off = 0, nbts = 0;
	assert_lt(cural_, btncand_.size());
	assert(res.repOk());
	// For each candidate cell that we should try to backtrack from...
	const size_t candsz = btncand_.size();
	size_t SQ = dpRows() >> 4;
	if(SQ == 0) SQ = 1;
	size_t rdlen = rdf_ - rdi_;
	bool checkpointed = rdlen >= cperMinlen_;
	while(cural_ < candsz) {
		// Doing 'continue' anywhere in here simply causes us to move on to the
		// next candidate
		if(btncand_[cural_].score < minsc) {
			btncand_[cural_].fate = BT_CAND_FATE_FILT_SCORE;
			nbtfiltsc_++; cural_++; continue;
		}
		nbts = 0;
		assert(sse8succ_ || sse16succ_);
		size_t row = btncand_[cural_].row;
		size_t col = btncand_[cural_].col;
		assert_lt(row, dpRows());
		assert_lt((TRefOff)col, rff_-rfi_);
		if(sse16succ_) {
			SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
			if(!checkpointed && d.mat_.reset_[row] && d.mat_.reportedThrough(row, col)) {
				// Skipping this candidate because a previous candidate already
				// moved through this cell
				btncand_[cural_].fate = BT_CAND_FATE_FILT_START;
				//cerr << "  skipped becuase starting cell was covered" << endl;
				nbtfiltst_++; cural_++; continue;
			}
		} else if(sse8succ_) {
			SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
			if(!checkpointed && d.mat_.reset_[row] && d.mat_.reportedThrough(row, col)) {
				// Skipping this candidate because a previous candidate already
				// moved through this cell
				btncand_[cural_].fate = BT_CAND_FATE_FILT_START;
				//cerr << "  skipped becuase starting cell was covered" << endl;
				nbtfiltst_++; cural_++; continue;
			}
		}
		if(sc_->monotone) {
			bool ret = false;
			if(sse8succ_) {
				uint32_t reseed = rnd.nextU32() + 1;
				rnd.init(reseed);
				res.reset();
				if(checkpointed) {
					size_t maxiter = MAX_SIZE_T;
					size_t niter = 0;
					ret = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res,      // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter,  // max # extensions to try
						niter,    // # extensions tried
						rnd);     // random gen, to choose among equal paths
				} else {
					ret = backtraceNucleotidesEnd2EndSseU8(
						btncand_[cural_].score, // in: expected score
						res,    // out: store results (edits and scores) here
						off,    // out: store diagonal projection of origin
						nbts,   // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
				}
#ifndef NDEBUG
				// if(...) statement here should check not whether the primary
				// alignment was checkpointed, but whether a checkpointed
				// alignment was done at all.
				if(!checkpointed) {
					SwResult res2;
					size_t maxiter2 = MAX_SIZE_T;
					size_t niter2 = 0;
					bool ret2 = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res2,     // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter2, // max # extensions to try
						niter2,   // # extensions tried
						rnd);     // random gen, to choose among equal paths
					// After the first alignment, there's no guarantee we'll
					// get the same answer from both backtrackers because of
					// differences in how they handle marking cells as
					// reported-through.
					assert(cural_ > 0 || !ret || ret == ret2);
					assert(cural_ > 0 || !ret || res.alres == res2.alres);
				}
				if(sse16succ_ && !checkpointed) {
					SwResult res2;
					size_t off2, nbts2 = 0;
					rnd.init(reseed);
					bool ret2 = backtraceNucleotidesEnd2EndSseI16(
						btncand_[cural_].score, // in: expected score
						res2,   // out: store results (edits and scores) here
						off2,   // out: store diagonal projection of origin
						nbts2,  // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
					assert_eq(ret, ret2);
					assert_eq(nbts, nbts2);
					assert(!ret || res2.alres.score() == res.alres.score());
#if 0
					if(!checkpointed && (rand() & 15) == 0) {
						// Check that same cells are reported through
						SSEData& d8  = fw_ ? sseU8fw_  : sseU8rc_;
						SSEData& d16 = fw_ ? sseI16fw_ : sseI16rc_;
						for(size_t i = d8.mat_.nrow(); i > 0; i--) {
							for(size_t j = 0; j < d8.mat_.ncol(); j++) {
								assert_eq(d8.mat_.reportedThrough(i-1, j),
										  d16.mat_.reportedThrough(i-1, j));
							}
						}
					}
#endif
				}
#endif
				rnd.init(reseed+1); // debug/release pseudo-randoms in lock step
			} else if(sse16succ_) {
				uint32_t reseed = rnd.nextU32() + 1;
				res.reset();
				if(checkpointed) {
					size_t maxiter = MAX_SIZE_T;
					size_t niter = 0;
					ret = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res,      // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter,  // max # extensions to try
						niter,    // # extensions tried
						rnd);     // random gen, to choose among equal paths
				} else {
					ret = backtraceNucleotidesEnd2EndSseI16(
						btncand_[cural_].score, // in: expected score
						res,    // out: store results (edits and scores) here
						off,    // out: store diagonal projection of origin
						nbts,   // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
				}
#ifndef NDEBUG
				// if(...) statement here should check not whether the primary
				// alignment was checkpointed, but whether a checkpointed
				// alignment was done at all.
				if(!checkpointed) {
					SwResult res2;
					size_t maxiter2 = MAX_SIZE_T;
					size_t niter2 = 0;
					bool ret2 = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res2,     // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter2, // max # extensions to try
						niter2,   // # extensions tried
						rnd);     // random gen, to choose among equal paths
					// After the first alignment, there's no guarantee we'll
					// get the same answer from both backtrackers because of
					// differences in how they handle marking cells as
					// reported-through.
					assert(cural_ > 0 || !ret || ret == ret2);
					assert(cural_ > 0 || !ret || res.alres == res2.alres);
				}
#endif
				rnd.init(reseed); // debug/release pseudo-randoms in lock step
			}
			if(ret) {
				btncand_[cural_].fate = BT_CAND_FATE_SUCCEEDED;
				break;
			} else {
				btncand_[cural_].fate = BT_CAND_FATE_FAILED;
			}
		} else {
			// Local alignment
			// Check if this solution is "dominated" by a prior one.
			// Domination is a heuristic designed to eliminate the vast
			// majority of valid-but-redundant candidates lying in the
			// "penumbra" of a high-scoring alignment.
			bool dom = false;
			{
				size_t donesz = btncanddone_.size();
				const size_t col = btncand_[cural_].col;
				const size_t row = btncand_[cural_].row;
				for(size_t i = 0; i < donesz; i++) {
					assert_gt(btncanddone_[i].fate, 0);
					size_t colhi = col, rowhi = row;
					size_t rowlo = btncanddone_[i].row;
					size_t collo = btncanddone_[i].col;
					if(colhi < collo) swap(colhi, collo);
					if(rowhi < rowlo) swap(rowhi, rowlo);
					if(colhi - collo <= SQ && rowhi - rowlo <= SQ) {
						// Skipping this candidate because it's "dominated" by
						// a previous candidate
						dom = true;
						break;
					}
				}
			}
			if(dom) {
				btncand_[cural_].fate = BT_CAND_FATE_FILT_DOMINATED;
				nbtfiltdo_++;
				cural_++;
				continue;
			}
			bool ret = false;
			if(sse8succ_) {
				uint32_t reseed = rnd.nextU32() + 1;
				res.reset();
				rnd.init(reseed);
				if(checkpointed) {
					size_t maxiter = MAX_SIZE_T;
					size_t niter = 0;
					ret = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res,      // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter,  // max # extensions to try
						niter,    // # extensions tried
						rnd);     // random gen, to choose among equal paths
				} else {
					ret = backtraceNucleotidesLocalSseU8(
						btncand_[cural_].score, // in: expected score
						res,    // out: store results (edits and scores) here
						off,    // out: store diagonal projection of origin
						nbts,   // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
				}
#ifndef NDEBUG
				// if(...) statement here should check not whether the primary
				// alignment was checkpointed, but whether a checkpointed
				// alignment was done at all.
				if(!checkpointed) {
					SwResult res2;
					size_t maxiter2 = MAX_SIZE_T;
					size_t niter2 = 0;
					bool ret2 = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res2,     // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter2, // max # extensions to try
						niter2,   // # extensions tried
						rnd);     // random gen, to choose among equal paths
					// After the first alignment, there's no guarantee we'll
					// get the same answer from both backtrackers because of
					// differences in how they handle marking cells as
					// reported-through.
					assert(cural_ > 0 || !ret || ret == ret2);
					assert(cural_ > 0 || !ret || res.alres == res2.alres);
				}
				if(!checkpointed && sse16succ_) {
					SwResult res2;
					size_t off2, nbts2 = 0;
					rnd.init(reseed); // same b/t backtrace calls
					bool ret2 = backtraceNucleotidesLocalSseI16(
						btncand_[cural_].score, // in: expected score
						res2,   // out: store results (edits and scores) here
						off2,   // out: store diagonal projection of origin
						nbts2,  // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
					assert_eq(ret, ret2);
					assert_eq(nbts, nbts2);
					assert(!ret || res2.alres.score() == res.alres.score());
#if 0
					if(!checkpointed && (rand() & 15) == 0) {
						// Check that same cells are reported through
						SSEData& d8  = fw_ ? sseU8fw_  : sseU8rc_;
						SSEData& d16 = fw_ ? sseI16fw_ : sseI16rc_;
						for(size_t i = d8.mat_.nrow(); i > 0; i--) {
							for(size_t j = 0; j < d8.mat_.ncol(); j++) {
								assert_eq(d8.mat_.reportedThrough(i-1, j),
										  d16.mat_.reportedThrough(i-1, j));
							}
						}
					}
#endif
				}
#endif
				rnd.init(reseed+1); // debug/release pseudo-randoms in lock step
			} else if(sse16succ_) {
				uint32_t reseed = rnd.nextU32() + 1;
				res.reset();
				if(checkpointed) {
					size_t maxiter = MAX_SIZE_T;
					size_t niter = 0;
					ret = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res,      // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter,  // max # extensions to try
						niter,    // # extensions tried
						rnd);     // random gen, to choose among equal paths
				} else {
					ret = backtraceNucleotidesLocalSseI16(
						btncand_[cural_].score, // in: expected score
						res,    // out: store results (edits and scores) here
						off,    // out: store diagonal projection of origin
						nbts,   // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
				}
#ifndef NDEBUG
				// if(...) statement here should check not whether the primary
				// alignment was checkpointed, but whether a checkpointed
				// alignment was done at all.
				if(!checkpointed) {
					SwResult res2;
					size_t maxiter2 = MAX_SIZE_T;
					size_t niter2 = 0;
					bool ret2 = backtrace(
						btncand_[cural_].score, // in: expected score
						true,     // in: use mini-fill?
						true,     // in: use checkpoints?
						res2,     // out: store results (edits and scores) here
						off,      // out: store diagonal projection of origin
						row,      // start in this rectangle row
						col,      // start in this rectangle column
						maxiter2, // max # extensions to try
						niter2,   // # extensions tried
						rnd);     // random gen, to choose among equal paths
					// After the first alignment, there's no guarantee we'll
					// get the same answer from both backtrackers because of
					// differences in how they handle marking cells as
					// reported-through.
					assert(cural_ > 0 || !ret || ret == ret2);
					assert(cural_ > 0 || !ret || res.alres == res2.alres);
				}
#endif
				rnd.init(reseed); // same b/t backtrace calls
			}
			if(ret) {
				btncand_[cural_].fate = BT_CAND_FATE_SUCCEEDED;
				btncanddone_.push_back(btncand_[cural_]);
				btncanddoneSucc_++;
				assert(res.repOk());
				break;
			} else {
				btncand_[cural_].fate = BT_CAND_FATE_FAILED;
				btncanddone_.push_back(btncand_[cural_]);
				btncanddoneFail_++;
			}
		}
		cural_++;
	} // while(cural_ < btncand_.size())
	if(cural_ == btncand_.size()) {
		assert(res.repOk());
		return false;
	}
	assert(!res.alres.empty());
	assert(res.repOk());
	if(!fw_) {
		// All edits are currently w/r/t upstream end; if read aligned
		// to Crick strand, we need to invert them so that they're
		// w/r/t the read's 5' end instead.
		res.alres.invertEdits();
	}
	cural_++;
	assert(res.repOk());
	return true;
}

#ifdef MAIN_ALIGNER_SW

#include <sstream>
#include <utility>
#include <getopt.h>
#include "scoring.h"
#include "aligner_seed_policy.h"

int  gGapBarrier;
int  gSnpPhred;
static int bonusMatchType;   // how to reward matches
static int bonusMatch;       // constant if match bonus is a constant
static int penMmcType;       // how to penalize mismatches
static int penMmc;           // constant if mm pelanty is a constant
static int penNType;         // how to penalize Ns in the read
static int penN;             // constant if N pelanty is a constant
static bool nPairCat;        // true -> concatenate mates before N filter
static int penRdExConst;     // constant coeff for cost of gap in read
static int penRfExConst;     // constant coeff for cost of gap in ref
static int penRdExLinear;    // linear coeff for cost of gap in read
static int penRfExLinear;    // linear coeff for cost of gap in ref
static float costMinConst;   // constant coeff for min score w/r/t read len
static float costMinLinear;  // linear coeff for min score w/r/t read len
static float costFloorConst; // constant coeff for score floor w/r/t read len
static float costFloorLinear;// linear coeff for score floor w/r/t read len
static float nCeilConst;     // constant coeff for N ceiling w/r/t read len
static float nCeilLinear;    // linear coeff for N ceiling w/r/t read len
static bool  nCatPair;       // concat mates before applying N filter?
static int multiseedMms;     // mismatches permitted in a multiseed seed
static int multiseedLen;     // length of multiseed seeds
static int multiseedIvalType;
static float multiseedIvalA;
static float multiseedIvalB;
static float posmin;
static float posfrac;
static float rowmult;

enum {
	ARG_TESTS = 256
};

static const char *short_opts = "s:m:r:d:i:";
static struct option long_opts[] = {
	{(char*)"snppen",       required_argument, 0, 's'},
	{(char*)"misspen",      required_argument, 0, 'm'},
	{(char*)"seed",         required_argument, 0, 'r'},
	{(char*)"align-policy", no_argument,       0, 'A'},
	{(char*)"test",         no_argument,       0, ARG_TESTS},
};

static void printUsage(ostream& os) {
	os << "Usage: aligner_sw <read-seq> <ref-nuc-seq> [options]*" << endl;
	os << "Options:" << endl;
	os << "  -s/--snppen <int>   penalty incurred by SNP; used for decoding"
	   << endl;
	os << "  -m/--misspen <int>  quality to use for read chars" << endl;
	os << "  -r/-seed <int>      seed for pseudo-random generator" << endl;
}

/**
 * Parse a T from a string 's'
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

static EList<bool> stbuf, enbuf;
static BTDnaString btread;
static BTString btqual;
static BTString btref;
static BTString btref2;

static BTDnaString readrc;
static BTString qualrc;

/**
 * Helper function for running a case consisting of a read (sequence
 * and quality), a reference string, and an offset that anchors the 0th
 * character of the read to a reference position.
 */
static void doTestCase(
	SwAligner&         al,
	const BTDnaString& read,
	const BTString&    qual,
	const BTString&    refin,
	TRefOff            off,
	EList<bool>       *en,
	const Scoring&     sc,
	TAlScore           minsc,
	SwResult&          res,
	bool               nsInclusive,
	bool               filterns,
	uint32_t           seed)
{
	RandomSource rnd(seed);
	btref2 = refin;
	assert_eq(read.length(), qual.length());
	size_t nrow = read.length();
	TRefOff rfi, rff;
	// Calculate the largest possible number of read and reference gaps given
	// 'minsc' and 'pens'
	size_t maxgaps;
	size_t padi, padf;
	{
		int readGaps = sc.maxReadGaps(minsc, read.length());
		int refGaps = sc.maxRefGaps(minsc, read.length());
		assert_geq(readGaps, 0);
		assert_geq(refGaps, 0);
		int maxGaps = max(readGaps, refGaps);
		padi = 2 * maxGaps;
		padf = maxGaps;
		maxgaps = (size_t)maxGaps;
	}
	size_t nceil = (size_t)sc.nCeil.f((double)read.length());
	size_t width = 1 + padi + padf;
	rfi = off;
	off = 0;
	// Pad the beginning of the reference with Ns if necessary
	if(rfi < padi) {
		size_t beginpad = (size_t)(padi - rfi);
		for(size_t i = 0; i < beginpad; i++) {
			btref2.insert('N', 0);
			off--;
		}
		rfi = 0;
	} else {
		rfi -= padi;
	}
	assert_geq(rfi, 0);
	// Pad the end of the reference with Ns if necessary
	while(rfi + nrow + padi + padf > btref2.length()) {
		btref2.append('N');
	}
	rff = rfi + nrow + padi + padf;
	// Convert reference string to masks
	for(size_t i = 0; i < btref2.length(); i++) {
		if(toupper(btref2[i]) == 'N' && !nsInclusive) {
			btref2.set(16, i);
		} else {
			int num = 0;
			int alts[] = {4, 4, 4, 4};
			decodeNuc(toupper(btref2[i]), num, alts);
			assert_leq(num, 4);
			assert_gt(num, 0);
			btref2.set(0, i);
			for(int j = 0; j < num; j++) {
				btref2.set(btref2[i] | (1 << alts[j]), i);
			}
		}
	}
	bool fw = true;
	uint32_t refidx = 0;
	size_t solwidth = width;
	if(maxgaps >= solwidth) {
		solwidth = 0;
	} else {
		solwidth -= maxgaps;
	}
	if(en == NULL) {
		enbuf.resize(solwidth);
		enbuf.fill(true);
		en = &enbuf;
	}
	assert_geq(rfi, 0);
	assert_gt(rff, rfi);
	readrc = read;
	qualrc = qual;
	al.initRead(
		read,          // read sequence
		readrc,
		qual,          // read qualities
		qualrc,
		0,             // offset of first character within 'read' to consider
		read.length(), // offset of last char (exclusive) in 'read' to consider
		floorsc);      // local-alignment score floor
	al.initRef(
		fw,            // 'read' is forward version of read?
		refidx,        // id of reference aligned to
		off,           // offset of upstream ref char aligned against
		btref2.wbuf(), // reference sequence (masks)
		rfi,           // offset of first char in 'ref' to consider
		rff,           // offset of last char (exclusive) in 'ref' to consider
		width,         // # bands to do (width of parallelogram)
		solwidth,      // # rightmost cols where solns can end
		sc,            // scoring scheme
		minsc,         // minimum score for valid alignment
		maxgaps,       // max of max # read gaps, ref gaps
		0,             // amount to truncate on left-hand side
		en);           // mask indicating which columns we can end in
	if(filterns) {
		al.filter((int)nceil);
	}
	al.align(rnd);
}

/**
 * Another interface for running a case.
 */
static void doTestCase2(
	SwAligner&         al,
	const char        *read,
	const char        *qual,
	const char        *refin,
	TRefOff            off,
	const Scoring&     sc,
	float              costMinConst,
	float              costMinLinear,
	SwResult&          res,
	bool               nsInclusive = false,
	bool               filterns = false,
	uint32_t           seed = 0)
{
	btread.install(read, true);
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		NULL,
		sc,  
		minsc,
		floorsc,
		res,
		nsInclusive,
		filterns,
		seed
	);
}

/**
 * Another interface for running a case.
 */
static void doTestCase3(
	SwAligner&         al,
	const char        *read,
	const char        *qual,
	const char        *refin,
	TRefOff            off,
	Scoring&           sc,
	float              costMinConst,
	float              costMinLinear,
	float              nCeilConst,
	float              nCeilLinear,
	SwResult&          res,
	bool               nsInclusive = false,
	bool               filterns = false,
	uint32_t           seed = 0)
{
	btread.install(read, true);
	// Calculate the penalty ceiling for the read
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	sc.nCeil.setType(SIMPLE_FUNC_LINEAR);
	sc.nCeil.setConst(costMinConst);
	sc.nCeil.setCoeff(costMinLinear);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		NULL,
		sc,  
		minsc,
		floorsc,
		res,
		nsInclusive,
		filterns,
		seed
	);
}

/**
 * Another interface for running a case.  Like doTestCase3 but caller specifies
 * st_ and en_ lists.
 */
static void doTestCase4(
	SwAligner&         al,
	const char        *read,
	const char        *qual,
	const char        *refin,
	TRefOff            off,
	EList<bool>&       en,
	Scoring&           sc,
	float              costMinConst,
	float              costMinLinear,
	float              nCeilConst,
	float              nCeilLinear,
	SwResult&          res,
	bool               nsInclusive = false,
	bool               filterns = false,
	uint32_t           seed = 0)
{
	btread.install(read, true);
	// Calculate the penalty ceiling for the read
	TAlScore minsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costMinConst,
		costMinLinear));
	TAlScore floorsc = (TAlScore)(Scoring::linearFunc(
		btread.length(),
		costFloorConst,
		costFloorLinear));
	btqual.install(qual);
	btref.install(refin);
	sc.nCeil.setType(SIMPLE_FUNC_LINEAR);
	sc.nCeil.setConst(costMinConst);
	sc.nCeil.setCoeff(costMinLinear);
	doTestCase(
		al,
		btread,
		btqual,
		btref,
		off,
		&en,
		sc,  
		minsc,
		floorsc,
		res,
		nsInclusive,
		filterns,
		seed
	);
}

/**
 * Do a set of unit tests.
 */
static void doTests() {
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	nPairCat        = DEFAULT_N_CAT_PAIR;
	penRdExConst    = DEFAULT_READ_GAP_CONST;
	penRfExConst    = DEFAULT_REF_GAP_CONST;
	penRdExLinear   = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear   = DEFAULT_REF_GAP_LINEAR;
	costMinConst    = DEFAULT_MIN_CONST;
	costMinLinear   = DEFAULT_MIN_LINEAR;
	costFloorConst  = DEFAULT_FLOOR_CONST;
	costFloorLinear = DEFAULT_FLOOR_LINEAR;
	nCeilConst      = 1.0f; // constant factor in N ceil w/r/t read len
	nCeilLinear     = 0.1f; // coeff of linear term in N ceil w/r/t read len
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	// Set up penalities
	Scoring sc(
		bonusMatch,
		penMmcType,    // how to penalize mismatches
		30,        // constant if mm pelanty is a constant
		30,        // penalty for decoded SNP
		costMinConst,  // constant factor in N ceiling w/r/t read length
		costMinLinear, // coeff of linear term in N ceiling w/r/t read length
		costFloorConst,  // constant factor in N ceiling w/r/t read length
		costFloorLinear, // coeff of linear term in N ceiling w/r/t read length
		nCeilConst,    // constant factor in N ceiling w/r/t read length
		nCeilLinear,   // coeff of linear term in N ceiling w/r/t read length
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nPairCat,      // true -> concatenate mates before N filtering
		25,  // constant coeff for cost of gap in read
		25,  // constant coeff for cost of gap in ref
		15, // linear coeff for cost of gap in read
		15, // linear coeff for cost of gap in ref
		1,             // # rows at top/bot can only be entered diagonally
		-1,            // min row idx to backtrace from; -1 = no limit
		false          // sort results first by row then by score?
	);
	// Set up alternative penalities
	Scoring sc2(
		bonusMatch,
		COST_MODEL_QUAL, // how to penalize mismatches
		30,          // constant if mm pelanty is a constant
		30,          // penalty for decoded SNP
		costMinConst,  // constant factor in N ceiling w/r/t read length
		costMinLinear, // coeff of linear term in N ceiling w/r/t read length
		costFloorConst,  // constant factor in N ceiling w/r/t read length
		costFloorLinear, // coeff of linear term in N ceiling w/r/t read length
		1.0f,            // constant factor in N ceiling w/r/t read length
		1.0f,            // coeff of linear term in N ceiling w/r/t read length
		penNType,        // how to penalize Ns in the read
		penN,            // constant if N pelanty is a constant
		nPairCat,        // true -> concatenate mates before N filtering
		25,    // constant coeff for cost of gap in read
		25,    // constant coeff for cost of gap in ref
		15,   // linear coeff for cost of gap in read
		15,   // linear coeff for cost of gap in ref
		1,               // # rows at top/bot can only be entered diagonally
		-1,              // min row idx to backtrace from; -1 = no limit
		false            // sort results first by row then by score?
	);
	SwResult res;
	
	//
	// Basic nucleotide-space tests
	//
	cerr << "Running tests..." << endl;
	int tests = 1;
	bool nIncl = false;
	bool nfilter = false;

	SwAligner al;
	RandomSource rnd(73);
	for(int i = 0; i < 3; i++) {
		cerr << "  Test " << tests++ << " (nuc space, offset "
		     << (i*4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		sc.rdGapLinear = 15;
		sc.rfGapLinear = 15;
	//        A           C           G           T           A           C           G           T
	//    H   E   F   H   E   F   H   E   F   H   E   F   H   E   F   H   E   F   H   E   F   H   E   F
	// A  0   lo  lo -30  lo  lo -30  lo  lo -30 lo lo 0 lo lo -30 lo lo-30 lo lo-30 lo lo
	// C -30  lo -55  0  -85 -85 -55 -55 -85
	// G -30  lo -70 -55 -85 -55  0 -100-100
	// T -30  lo -85 -60 -85 -70 -55-100 -55
	// A  0   lo -85 -55 -55 -85 -70 -70 -70
	// C -30  lo -55  0  -85-100 -55 -55 -85
	// G -30  lo -70 -55 -85 -55  0 -100-100
	// T -30  lo -85 -60 -85 -70 -55-100 -55
		doTestCase2(
			al,
			"ACGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc)...";
		sc.setMmPen(COST_MODEL_CONSTANT, 30);
		//sc.setMatchBonus(10);
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc, check qual 1)...";
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"ABCDEFGH",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc2,                // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		size_t lo, hi;
		if(i == 0) {
			lo = 0; hi = 1;
		} else if(i == 1) {
			lo = 1; hi = 2;
		} else {
			lo = 2; hi = 3;
		}
		for(size_t j = lo; j < hi; j++) {
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(j*4, res.alres.refoff());
			assert_eq(8, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), -36);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
		}
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc, check qual 2)...";
		doTestCase2(
			al,
			"ACGAACGT",         // read
			"ABCDEFGH",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc2,                // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc, check qual )...";
		assert(res.empty());
		doTestCase2(
			al,
			"TCGTACGT",         // read
			"ABCDEFGH",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc2,                // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -32);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm at the beginning, allowed by minsc)...";
		doTestCase2(
			al,
			"CCGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -30);
		assert_eq(res.alres.score().ns(), 0);
		assert_eq(1, res.alres.ned().size());
		assert_eq(0, res.alres.aed().size());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 n in read, allowed)...";
		doTestCase3(
			al,
			"ACGTNCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0f,               // allow 1 N
			0.0f,               // allow 1 N
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 n in read, allowed)...";
		doTestCase3(
			al,
			"ACGNNCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			2.0f,               // const coeff for N ceiling
			0.0f,               // linear coeff for N ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 n in read, 1 at beginning, allowed)...";
		doTestCase2(
			al,
			"NCGTNCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -2);
		assert_eq(res.alres.score().ns(), 2);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 n in ref, allowed)...";
		doTestCase2(
			al,
			"ACGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTNCGTACGTANGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), -1);
		assert_eq(res.alres.score().ns(), 1);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-10.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		// Read gap with equal read and ref gap penalties
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", read gap allowed by minsc)...";
		assert(res.empty());
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", read gap disallowed by minsc)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();

		cerr << "PASSED" << endl;
		// Ref gap with equal read and ref gap penalties
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", read gap disallowed by gap barrier)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 15;
		sc.gapbar = 4;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		sc.gapbar = 1;
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();

		cerr << "PASSED" << endl;
		// Ref gap with equal read and ref gap penalties
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap allowed by minsc, gapbar=3)...";
		sc.gapbar = 3;
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		sc.gapbar = 1;
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;

		// Ref gap with equal read and ref gap penalties
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap allowed by minsc, gapbar=4)...";
		sc.gapbar = 4;
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		sc.gapbar = 1;
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", ref gap disallowed by gap barrier)...";
		sc.gapbar = 5;
		doTestCase2(
			al,
			"ACGTAACGT",        // read
			"IIIIIIIII",        // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		sc.gapbar = 1;
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		
		// Read gap with one read gap and zero ref gaps allowed
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 read gap, ref gaps disallowed by minsc)...";
		sc.rfGapConst = 35;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 20;
		sc.rdGapLinear = 10;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -35);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", gaps disallowed by minsc)...";
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		doTestCase2(
			al,
			"ACGTCGT",          // read
			"IIIIIII",          // qual 
			"ACGTACGTACGTACGT", // ref 
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(res.empty());
		cerr << "PASSED" << endl;
		
		// Ref gap with one ref gap and zero read gaps allowed
		sc.rfGapConst = 25;
		sc.rdGapConst = 35;
		sc.rfGapLinear = 12;
		sc.rdGapLinear = 22;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, read gaps disallowed by minsc)...";
		assert(res.empty());
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -37);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
			<< ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;
		
		// Read gap with one read gap and two ref gaps allowed
		sc.rfGapConst = 20;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 15;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 read gap, 2 ref gaps allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		cerr << "PASSED" << endl;

		// Ref gap with one ref gap and two read gaps allowed
		sc.rfGapConst = 25;
		sc.rdGapConst = 11;  // if this were 10, we'd have ties
		sc.rfGapLinear = 15;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, 2 read gaps allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4) << ", gaps disallowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(res.empty());
		res.reset();
		assert(al.done());
		cerr << "PASSED" << endl;
		
		// Read gap with two read gaps and two ref gaps allowed
		sc.rfGapConst = 15;
		sc.rdGapConst = 15;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 ref gaps, 2 read gaps allowed by minsc)...";
		doTestCase3(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-40.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0,                // const coeff for N ceiling
			0.0,                // linear coeff for N ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			true);              // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		if(!res.empty()) {
			//al.printResultStacked(res, cerr); cerr << endl;
		}
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -25);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		// The following alignment is possible when i == 2:
		//   ACGTACGTACGTACGTN
		// A             x
		// C              x
		// G               x
		// T                x
		// C                x
		// G                x
		// T                 x
		assert(i == 2 || res.empty());
		res.reset();
		cerr << "PASSED" << endl;
		
		sc.rfGapConst = 10;
		sc.rdGapConst = 10;
		sc.rfGapLinear = 10;
		sc.rdGapLinear = 10;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, 1 read gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTCGT",
			"IIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -20);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		
		// Ref gap with two ref gaps and zero read gaps allowed
		sc.rfGapConst = 15;
		sc.rdGapConst = 15;
		sc.rfGapLinear = 5;
		sc.rdGapLinear = 5;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 2 ref gaps, 2 read gaps allowed by minsc)...";
		// Careful: it might be possible for the read to align with overhang
		// instead of with a gap
		doTestCase3(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-35.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			1.0f,               // needed to avoid overhang alignments
			0.0f,               // needed to avoid overhang alignments
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			true);              // filter Ns
		if(i == 0) {
			lo = 0; hi = 1;
		} else if(i == 1) {
			lo = 1; hi = 2;
		} else {
			lo = 2; hi = 3;
		}
		for(size_t j = lo; j < hi; j++) {
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			//al.printResultStacked(res, cerr); cerr << endl;
			assert(res.alres.refoff() == 0 ||
			       res.alres.refoff() == 4 ||
				   res.alres.refoff() == 8);
			assert_eq(8, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 1);
			assert_eq(res.alres.score().score(), -20);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
		}
		al.nextAlignment(res, rnd);
		//assert(res.empty());
		//res.reset();
		cerr << "PASSED" << endl;
		
		sc.rfGapConst = 25;
		sc.rdGapConst = 25;
		sc.rfGapLinear = 4;
		sc.rdGapLinear = 4;
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1 ref gap, 1 read gap allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTAACGT",
			"IIIIIIIII",
			"ACGTACGTACGTACGT",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 1);
		assert_eq(res.alres.score().score(), -29);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", short read)...";
		doTestCase2(
			al,
			"A",
			"I",
			"AAAAAAAAAAAA",
			i*4,                // off
			sc,                 // scoring scheme
			-30.0f,             // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++
			     << " (nuc space, offset 0, short read & ref)...";
			doTestCase2(
				al,
				"A",
				"I",
				"A",
				0,                  // off
				sc,                 // scoring scheme
				-30.0f,             // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", short read, many allowed gaps)...";
		doTestCase2(
			al,
			"A",
			"I",
			"AAAAAAAAAAAA",
			i*4,                // off
			sc,                 // scoring scheme
			-150.0f,            // const coeff for cost ceiling
			0.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 0);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++
			     << " (nuc space, offset 0, short read & ref, "
				 << "many allowed gaps)...";
			doTestCase2(
				al,
				"A",
				"I",
				"A",
				0,                  // off
				sc,                 // scoring scheme
				-150.0f,            // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 0);
			assert_eq(res.alres.score().ns(), 0);
			res.reset();
			cerr << "PASSED" << endl;
		}
	}

	// A test case where a valid alignment with a worse score should be
	// accepted over a valid alignment with a better score but too many
	// Ns
	cerr << "  Test " << tests++ << " (N ceiling 1)...";
	sc.mmcostType = COST_MODEL_CONSTANT;
	sc.mmcost = 10;
	sc.snp = 30;
	sc.nCeilConst  = 0.0f;
	sc.nCeilLinear = 0.0f;
	sc.rfGapConst  = 10;
	sc.rdGapLinear = 10;
	sc.rfGapConst  = 10;
	sc.rfGapLinear = 10;
	sc.setNPen(COST_MODEL_CONSTANT, 2);
	sc.gapbar = 1;
	// No Ns allowed, so this hit should be filtered
	doTestCase2(
		al,
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		0.0f,       // linear coeff for cost ceiling
		res,        // result
		false,      // ns are in inclusive
		true,       // nfilter
		0);
	al.nextAlignment(res, rnd);
	assert(res.empty());
	cerr << "PASSED" << endl;
	res.reset();

	// 1 N allowed, so this hit should stand
	cerr << "  Test " << tests++ << " (N ceiling 2)...";
	doTestCase3(
		al,
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		0.0f,       // linear coeff for cost ceiling
		1.0f,       // constant coefficient for # Ns allowed
		0.0f,       // linear coefficient for # Ns allowed
		res,        // result
		false,      // ns are in inclusive
		false,      // nfilter - NOTE: FILTER OFF
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(0,  res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1,  res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();

	// 1 N allowed, but we set st_ such that this hit should not stand
	for(size_t i = 0; i < 2; i++) {
		cerr << "  Test " << tests++ << " (N ceiling 2 with st_ override)...";
		EList<bool> en;
		en.resize(3); en.fill(true);
		if(i == 1) {
			en[1] = false;
		}
		sc.rfGapConst  = 10;
		sc.rdGapLinear = 10;
		sc.rfGapConst  = 10;
		sc.rfGapLinear = 10;
		doTestCase4(
			al,
			"ACGTACGT", // read seq
			"IIIIIIII", // read quals
			"NCGTACGT", // ref seq
			0,          // offset
			en,         // rectangle columns where solution can end
			sc,         // scoring scheme
			-25.0f,     // const coeff for cost ceiling
			0.0f,       // linear coeff for cost ceiling
			1.0f,       // constant coefficient for # Ns allowed
			0.0f,       // linear coefficient for # Ns allowed
			res,        // result
			false,      // ns are in inclusive
			false,      // nfilter - NOTE: FILTER OFF
			0);
		al.nextAlignment(res, rnd);
		if(i > 0) {
			assert(res.empty());
		} else {
			assert(!res.empty());
		}
		cerr << "PASSED" << endl;
		res.reset();
	}

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (N ceiling 3)...";
	sc.nCeilConst = 1.0f;
	sc.nCeilLinear = 0.0f;
	doTestCase2(
		al,
		"ACGTACGT", // read seq
		"IIIIIIII", // read quals
		"NCGTACGT", // ref seq
		0,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		0.0f,       // linear coeff for cost ceiling
		res,        // result
		false,      // ns are in inclusive
		true,       // nfilter - NOTE: FILTER ON
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(0,  res.alres.score().gaps());
	assert_eq(-2, res.alres.score().score());
	assert_eq(1,  res.alres.score().ns());
	cerr << "PASSED" << endl;
	res.reset();

	// No Ns allowed, so this hit should be filtered
	cerr << "  Test " << tests++ << " (redundant alignment elimination 1)...";
	sc.nCeilConst = 1.0f;
	sc.nCeilLinear = 0.0f;
	sc.rfGapConst  = 25;
	sc.rdGapLinear = 15;
	sc.rfGapConst  = 25;
	sc.rfGapLinear = 15;
	doTestCase2(
		al,
		//                   1         2         3         4
		//         01234567890123456789012345678901234567890123456
		          "AGGCTATGCCTCTGACGCGATATCGGCGCCCACTTCAGAGCTAACCG",
		          "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
		  "TTTTTTTTAGGCTATGCCTCTGACGCGATATCGGCGCCCACTTCAGAGCTAACCGTTTTTTT",
		// 01234567890123456789012345678901234567890123456789012345678901
		//           1         2         3         4         5         6
		8,          // offset
		sc,         // scoring scheme
		-25.0f,     // const coeff for cost ceiling
		-5.0f,      // linear coeff for cost ceiling
		res,        // result
		false,      // ns are in inclusive
		true,       // nfilter - NOTE: FILTER ON
		0);
	al.nextAlignment(res, rnd);
	assert(!res.empty());
	assert_eq(8, res.alres.refoff());
	assert_eq(47, res.alres.refExtent());
	assert_eq(0, res.alres.score().gaps());
	assert_eq(0, res.alres.score().score());
	assert_eq(0, res.alres.score().ns());
	res.reset();
	al.nextAlignment(res, rnd);
	assert(res.empty());
	assert(al.done());
	cerr << "PASSED" << endl;
	res.reset();
	
}

/**
 * Do a set of unit tests for local alignment.
 */
static void doLocalTests() {
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS_LOCAL;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	nPairCat        = DEFAULT_N_CAT_PAIR;
	penRdExConst    = DEFAULT_READ_GAP_CONST;
	penRfExConst    = DEFAULT_REF_GAP_CONST;
	penRdExLinear   = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear   = DEFAULT_REF_GAP_LINEAR;
	costMinConst    = DEFAULT_MIN_CONST_LOCAL;
	costMinLinear   = DEFAULT_MIN_LINEAR_LOCAL;
	costFloorConst  = DEFAULT_FLOOR_CONST_LOCAL;
	costFloorLinear = DEFAULT_FLOOR_LINEAR_LOCAL;
	nCeilConst      = 1.0f; // constant factor in N ceil w/r/t read len
	nCeilLinear     = 0.1f; // coeff of linear term in N ceil w/r/t read len
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	// Set up penalities
	Scoring sc(
		10,
		penMmcType,    // how to penalize mismatches
		30,            // constant if mm pelanty is a constant
		penSnp,        // penalty for decoded SNP
		costMinConst,  // constant factor in N ceiling w/r/t read length
		costMinLinear, // coeff of linear term in N ceiling w/r/t read length
		costFloorConst,  // constant factor in N ceiling w/r/t read length
		costFloorLinear, // coeff of linear term in N ceiling w/r/t read length
		nCeilConst,    // constant factor in N ceiling w/r/t read length
		nCeilLinear,   // coeff of linear term in N ceiling w/r/t read length
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nPairCat,      // true -> concatenate mates before N filtering
		25,            // constant coeff for cost of gap in read
		25,            // constant coeff for cost of gap in ref
		15,            // linear coeff for cost of gap in read
		15,            // linear coeff for cost of gap in ref
		1,             // # rows at top/bot can only be entered diagonally
		-1,            // min row idx to backtrace from; -1 = no limit
		false          // sort results first by row then by score?
	);
	SwResult res;
	
	//
	// Basic nucleotide-space tests
	//
	cerr << "Running local tests..." << endl;
	int tests = 1;
	bool nIncl = false;
	bool nfilter = false;

	SwAligner al;
	RandomSource rnd(73);
	for(int i = 0; i < 3; i++) {
		cerr << "  Test " << tests++ << " (short nuc space, offset "
		     << (i*4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGT",             // read
			"IIII",             // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			8.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(4, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 40);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;

		//     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
		//     A C G T A C G T A C G T A C G T
		// 0 C
		// 1 C   x
		// 2 G     x
		// 3 T       x

		cerr << "  Test " << tests++ << " (short nuc space, offset "
		     << (i*4) << ", 1mm)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"CCGT",             // read
			"IIII",             // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			7.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4+1, res.alres.refoff());
		assert_eq(3, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 30);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (short nuc space, offset "
		     << (i*4) << ", 1mm)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGA",             // read
			"IIII",             // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			7.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(3, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 30);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;

		if(i == 0) {
			cerr << "  Test " << tests++ << " (short nuc space, offset "
				 << (i*4) << ", 1mm, match bonus=20)...";
			sc.rdGapConst = 40;
			sc.rfGapConst = 40;
			sc.setMatchBonus(20);
			doTestCase2(
				al,
				"TTGT",             // read
				"IIII",             // qual
				"TTGA",             // ref in
				i*4,                // off
				sc,                 // scoring scheme
				25.0f,               // const coeff for cost ceiling
				0.0f,               // linear coeff for cost ceiling
				res,                // result
				nIncl,              // Ns inclusive (not mismatches)
				nfilter);           // filter Ns
			assert(!al.done());
			al.nextAlignment(res, rnd);
			assert(!res.empty());
			assert_eq(i*4, res.alres.refoff());
			assert_eq(3, res.alres.refExtent());
			assert_eq(res.alres.score().gaps(), 0);
			assert_eq(res.alres.score().score(), 60);
			assert_eq(res.alres.score().ns(), 0);
			assert(res.alres.ned().empty());
			assert(res.alres.aed().empty());
			res.reset();
			al.nextAlignment(res, rnd);
			assert(res.empty());
			assert(al.done());
			res.reset();
			sc.setMatchBonus(10);
			cerr << "PASSED" << endl;
		}

		cerr << "  Test " << tests++ << " (nuc space, offset "
		     << (i*4) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGTACGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			8.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 80);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		assert(res.empty());
		assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;
		
		cerr << "  Test " << tests++ << " (long nuc space, offset "
		     << (i*8) << ", exact)...";
		sc.rdGapConst = 40;
		sc.rfGapConst = 40;
		doTestCase2(
			al,
			"ACGTACGTACGTACGTACGTA", // read
			"IIIIIIIIIIIIIIIIIIIII",  // qual
			"ACGTACGTACGTACGTACGTACGTACGTACGTACGTA", // ref in
		//   ACGTACGTACGTACGTACGT
		//           ACGTACGTACGTACGTACGT
		//                   ACGTACGTACGTACGTACGT
			i*8,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			8.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*8, res.alres.refoff());
		assert_eq(21, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 210);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		//assert(res.empty());
		//assert(al.done());
		res.reset();
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (nuc space, offset " << (i*4)
		     << ", 1mm allowed by minsc)...";
		doTestCase2(
			al,
			"ACGTTCGT",         // read
			"IIIIIIII",         // qual
			"ACGTACGTACGTACGT", // ref in
			i*4,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			5.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*4, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 40);
		assert_eq(res.alres.score().ns(), 0);
		res.reset();
		al.nextAlignment(res, rnd);
		//assert(res.empty());
		//assert(al.done());
		cerr << "PASSED" << endl;

		cerr << "  Test " << tests++ << " (long nuc space, offset "
		     << (i*8) << ", 6mm allowed by minsc)...";
		sc.rdGapConst = 50;
		sc.rfGapConst = 50;
		sc.rdGapLinear = 45;
		sc.rfGapLinear = 45;
		doTestCase2(
			al,
			"ACGTACGATGCATCGTACGTA", // read
			"IIIIIIIIIIIIIIIIIIIII",  // qual
			"ACGTACGTACGTACGTACGTACGTACGTACGTACGTA", // ref in
		//   ACGTACGTACGTACGTACGT
		//           ACGTACGTACGTACGTACGT
		//                   ACGTACGTACGTACGTACGT
			i*8,                // off
			sc,                 // scoring scheme
			0.0f,               // const coeff for cost ceiling
			1.0f,               // linear coeff for cost ceiling
			res,                // result
			nIncl,              // Ns inclusive (not mismatches)
			nfilter);           // filter Ns
		assert(!al.done());
		al.nextAlignment(res, rnd);
		assert(!res.empty());
		assert_eq(i*8 + 13, res.alres.refoff());
		assert_eq(8, res.alres.refExtent());
		assert_eq(res.alres.score().gaps(), 0);
		assert_eq(res.alres.score().score(), 80);
		assert_eq(res.alres.score().ns(), 0);
		assert(res.alres.ned().empty());
		assert(res.alres.aed().empty());
		res.reset();
		al.nextAlignment(res, rnd);
		res.reset();
		cerr << "PASSED" << endl;
	}
}

int main(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	unsigned seed = 0;
	gGapBarrier = 1;
	gSnpPhred = 30;
	bool nsInclusive = false;
	bool nfilter = false;
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	penRdExConst    = DEFAULT_READ_GAP_CONST;
	penRfExConst    = DEFAULT_REF_GAP_CONST;
	penRdExLinear   = DEFAULT_READ_GAP_LINEAR;
	penRfExLinear   = DEFAULT_REF_GAP_LINEAR;
	costMinConst    = DEFAULT_MIN_CONST;
	costMinLinear   = DEFAULT_MIN_LINEAR;
	costFloorConst  = DEFAULT_FLOOR_CONST;
	costFloorLinear = DEFAULT_FLOOR_LINEAR;
	nCeilConst      = 1.0f; // constant factor in N ceiling w/r/t read length
	nCeilLinear     = 1.0f; // coeff of linear term in N ceiling w/r/t read length
	nCatPair        = false;
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	multiseedIvalType = DEFAULT_IVAL;
	multiseedIvalA    = DEFAULT_IVAL_A;
	multiseedIvalB    = DEFAULT_IVAL_B;
	mhits           = 1;
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 's': gSnpPhred  = parse<int>(optarg); break;
			case 'r': seed       = parse<unsigned>(optarg); break;
			case ARG_TESTS: {
				doTests();
				cout << "PASSED end-to-ends" << endl;
				doLocalTests();
				cout << "PASSED locals" << endl;
				return 0;
			}
			case 'A': {
				bool localAlign = false;
				bool noisyHpolymer = false;
				bool ignoreQuals = false;
				SeedAlignmentPolicy::parseString(
					optarg,
					localAlign,
					noisyHpolymer,
					ignoreQuals,
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
					costMinConst,
					costMinLinear,
					costFloorConst,
					costFloorLinear,
					nCeilConst,
					nCeilLinear,
					nCatPair,
					multiseedMms,
					multiseedLen,
					multiseedIvalType,
					multiseedIvalA,
					multiseedIvalB,
					posmin);
				break;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	srand(seed);
	if(argc - optind < 4) {
		cerr << "Need at least 4 arguments" << endl;
		printUsage(cerr);
		exit(1);
	}
	BTDnaString read;
	BTString ref, qual;
	// Get read
	read.installChars(argv[optind]);
	// Get qualities
	qual.install(argv[optind+1]);
	assert_eq(read.length(), qual.length());
	// Get reference
	ref.install(argv[optind+2]);
	// Get reference offset
	size_t off = parse<size_t>(argv[optind+3]);
	// Set up penalities
	Scoring sc(
		false,         // local alignment?
		false,         // bad homopolymer?
		bonusMatchType,
		bonusMatch,
		penMmcType,    // how to penalize mismatches
		penMmc,        // constant if mm pelanty is a constant
		costMinConst,
		costMinLinear,
		costFloorConst,
		costFloorLinear,
		nCeilConst,    // N ceiling constant coefficient
		nCeilLinear,   // N ceiling linear coefficient
		penNType,      // how to penalize Ns in the read
		penN,          // constant if N pelanty is a constant
		nCatPair,      // true -> concatenate mates before N filtering
		penRdExConst,  // constant cost of extending a gap in the read
		penRfExConst,  // constant cost of extending a gap in the reference
		penRdExLinear, // coeff of linear term for cost of gap extension in read
		penRfExLinear  // coeff of linear term for cost of gap extension in ref
	);
	// Calculate the penalty ceiling for the read
	TAlScore minsc = Scoring::linearFunc(
		read.length(),
		costMinConst,
		costMinLinear);
	TAlScore floorsc = Scoring::linearFunc(
		read.length(),
		costFloorConst,
		costFloorLinear);
	SwResult res;
	SwAligner al;
	doTestCase(
		al,
		read,
		qual,
		ref,
		off,
		NULL,
		sc,  
		minsc,
		res,
		nsInclusive,
		nfilter,
		seed);
}
#endif /*MAIN_ALIGNER_SW*/
