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
	if(dpLog_ != NULL) {
		if(!firstRead_) {
			(*dpLog_) << '\n';
		}
		(*dpLog_) << rdfw.toZBuf() << '\t' << qufw.toZBuf();
	}
	firstRead_ = false;
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
	// Record the reference sequence in the log
	if(dpLog_ != NULL) {
		(*dpLog_) << '\t';
		(*dpLog_) << refidx_ << ',';
		(*dpLog_) << reflen_ << ',';
		(*dpLog_) << minsc_ << ',';
		(*dpLog_) << (fw ? '+' : '-') << ',';
		rect_->write(*dpLog_);
		(*dpLog_) << ',';
		for(TRefOff i = rfi_; i < rff_; i++) {
			(*dpLog_) << mask2dna[(int)rf[i]];
		}
	}
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
				rfbuf2.push_back(refs.getBase(refidx, (size_t)rfii));
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
 * Given a read, an alignment orientation, a range of characters in a
 * referece sequence, and a bit-encoded version of the reference, set up
 * and execute the corresponding ungapped alignment problem.  There can
 * only be one solution.
 *
 * The caller has already narrowed down the relevant portion of the
 * reference.
 *
 * Does not handle the case where we'd like to scan a large section of the
 * reference for an ungapped alignment, e.g., if we're searching for the
 * opposite mate after finding an alignment for the anchor mate.
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
				rfbuf2.push_back(refs.getBase(refidx, (size_t)rfii));
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
	res.alres.setScore(AlnScore(score,
								(int)(rd.length() - ned.size()),
								(int)ned.size(), ns, 0));
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
		if(dpLog_ != NULL) {
			(*dpLog_) << ",0,0";
		}
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
	if(dpLog_ != NULL) {
		(*dpLog_) << ",1," << best;
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
					// TODO: I find that sometimes there is disagreement here
					// where the alignments are in the same place with
					// identical scores, but one is more soft-trimmed than the other
					//assert(cural_ > 0 || !ret || res.alres == res2.alres);
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
					// Both backtraces start from the same cell, so they agree
					// on success, score and where the alignment ends.  They
					// can break a tie differently where it starts, e.g. one
					// soft-clips a prefix that adds nothing to the score
					// (matches whose bonus a gap cancels) while the other
					// keeps it, so the alignments themselves need not be
					// equal.
					assert(cural_ > 0 || !ret || ret == ret2);
					assert(cural_ > 0 || !ret || res.alres.score() == res2.alres.score());
					assert(cural_ > 0 || !ret ||
					       res.alres.refoff() + res.alres.refExtent() ==
					       res2.alres.refoff() + res2.alres.refExtent());
					assert(cural_ > 0 || !ret ||
					       res.alres.trimmedRight(true) == res2.alres.trimmedRight(true));
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
