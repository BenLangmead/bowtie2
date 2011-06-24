/*
 * aligner_result.cpp
 */

#include <iostream>
#include "reference.h"
#include "aligner_result.h"
#include "read.h"
#include "edit.h"
#include "sstring.h"
#include "ds.h"
#include "util.h"

using namespace std;

/**
 * Clear all contents.
 */
void AlnRes::reset() {
	ned_.clear();
	aed_.clear();
	ced_.clear();
	score_.invalidate();
	refcoord_.invalidate();
	refival_.invalidate();
	shapeSet_     = false;
	rdlen_        = 0;
	rdrows_       = 0;
	rdextent_     = 0;
	rdexrows_     = 0;
	rfextent_     = 0;
	refns_        = 0;
	type_         = ALN_RES_TYPE_UNPAIRED;
	fraglen_      = -1;
	// Trimming of the nucleotide read or the decoded nucleotides of a
	// colorspace read
	trimSoft_     = false;
	trim5p_       = 0;
	trim3p_       = 0;
	pretrimSoft_  = true;
	pretrim5p_    = 0;
	pretrim3p_    = 0;
	// Trimming of the colorspace read
	cTrimSoft_    = false;
	cTrim5p_      = 0;
	cTrim3p_      = 0;
	cPretrimSoft_ = true;
	cPretrim5p_   = 0;
	cPretrim3p_   = 0;
	seedmms_      = 0; // number of mismatches allowed in seed
	seedlen_      = 0; // length of seed
	seedival_     = 0; // interval between seeds
	minsc_        = 0; // minimum score
	floorsc_      = 0; // score floor
	assert(!refcoord_.valid());
	assert(!refival_.valid());
}

/**
 * Set the upstream-most reference offset involved in the alignment, and
 * the extent of the alignment (w/r/t the reference)
 */
void AlnRes::setShape(
	TRefId  id,          // id of reference aligned to
	TRefOff off,         // offset of first aligned char into ref seq
	bool    fw,          // aligned to Watson strand?
	size_t  rdlen,       // length of read after hard trimming, before soft
	bool    color,       // colorspace alignment?
	bool    pretrimSoft, // whether trimming prior to alignment was soft
	size_t  pretrim5p,   // # poss trimmed form 5p end before alignment
	size_t  pretrim3p,   // # poss trimmed form 3p end before alignment
	bool    trimSoft,    // whether local-alignment trimming was soft
	size_t  trim5p,      // # poss trimmed form 5p end during alignment
	size_t  trim3p)      // # poss trimmed form 3p end during alignment
{
	rdlen_       = rdlen;
	rdrows_      = rdlen + (color ? 1 : 0);
	color_       = color;
	refcoord_.init(id, off, fw);
	if(color) {
		// In effect, all that trimming becomes hard trimming for the
		// decoded nucleotide-space sequence.
		pretrimSoft_  = false;
		pretrim5p_    = pretrim5p;
		pretrim3p_    = pretrim3p;
		trimSoft_     = false;
		trim5p_       = trim5p;
		trim3p_       = trim3p;
		// All trimming so far has been to colorspace sequence
		cPretrimSoft_ = pretrimSoft;
		cPretrim5p_   = pretrim5p;
		cPretrim3p_   = pretrim3p;
		cTrimSoft_    = trimSoft;
		cTrim5p_      = trim5p;
		cTrim3p_      = trim3p;
	} else {
		pretrimSoft_  = pretrimSoft;
		pretrim5p_    = pretrim5p;
		pretrim3p_    = pretrim3p;
		trimSoft_     = trimSoft;
		trim5p_       = trim5p;
		trim3p_       = trim3p;
		cTrimSoft_    = false;
		cTrim5p_      = 0;
		cTrim3p_      = 0;
		cPretrimSoft_ = true;
		cPretrim5p_   = 0;
		cPretrim3p_   = 0;
	}
	// Propagate trimming to the edits.  We assume that the pos fields of the
	// edits are set w/r/t to the rows of the dynamic programming table, and
	// haven't taken trimming into account yet.
	//
	// TODO: The division of labor between the aligner and the AlnRes is not
	// clean.  Perhaps the trimming and *all* of its side-effects should be
	// handled by the aligner.
	size_t trimBeg = fw ? trim5p : trim3p;
	if(trimBeg > 0) {
		for(size_t i = 0; i < ned_.size(); i++) {
			// Shift by trim5p, since edits are w/r/t 5p end
			assert_geq(ned_[i].pos, trimBeg);
			ned_[i].pos -= trimBeg;
		}
	}
	// Length after all soft trimming and any hard trimming that occurred
	// during alignment
	rdextent_ = rdlen;
	if(pretrimSoft_) {
		rdextent_ -= (pretrim5p + pretrim3p); // soft trim
	}
	rdextent_ -= (trim5p + trim3p); // soft or hard trim from alignment
	assert_gt(rdextent_, 0);
	rdexrows_ = rdextent_ + (color ? 1 : 0);
	calcRefExtent();
	refival_.init(id, off, fw, rfextent_);
	shapeSet_ = true;
}
	
/**
 * Get the decoded nucleotide sequence.  The alignment must have been
 * colorspace.
 */
void AlnRes::decodedNucsAndQuals(
	const Read& rd,        // read that led to alignment
	BTDnaString& ns,       // out: decoded nucleotides
	BTString& qs) const    // out: decoded qualities
{
	assert(shapeSet_);
	assert(color_);
	// Walk along the colors
	bool fw = refcoord_.fw();
	// How long is the color-to-color alignment after all trimming?
	const size_t rdext = readExtent();
	assert_gt(rdext, 0);
	// How long is the decoded string?  (Note that if user has requested that
	// ends be excluded, that happens later on.)  
	size_t len = rdext+1;
	ns.resize(rdext);
	qs.resize(rdext);
	for(size_t i = 0; i < rdext; i++) {
		ns.set(rd.patFw[trim5p_ + i], i); // set to original colors at first
		qs.set(rd.qual[trim5p_ + i],  i); // set to original qualities at first
	}
	if(!fw) {
		// Reverse the read to make it upstream-to-downstream.  It's in
		// colorspace, so no need to complement.
		ns.reverse();
		qs.reverse();
		// Flip edit,making them upstream-to-downstream.
		Edit::invertPoss(const_cast<EList<Edit>& >(ced_), rdextent_);
	}
	ns.resize(len); ns.set(4, len-1);
	qs.resize(len);
	
	// Convert 3' and 5' nucleotides to upstream and downstream
	// nucleotides
	int nup = fw ? nuc5p_ : nuc3p_;
	ASSERT_ONLY(int ndn = fw ? nuc3p_ : nuc5p_);
	
	// Note: the nucleotides in the ned and aed lists are already
	// w/r/t to the Watson strand so there's no need to complement
	// them
	int lastn = nup;
	size_t cedidx = 0;
	int c = ns[0], q = qs[0]-33;
	int lastq = 0;
	for(size_t i = 0; i < rdext; i++) {
		// If it was a miscall, get the true subject color
		if(cedidx < ced_.size() && ced_[cedidx].pos == i) {
			assert_neq("ACGTN"[c], ced_[cedidx].chr);
			assert_eq ("ACGTN"[c], ced_[cedidx].qchr);
			c = ced_[cedidx].chr;
			c = asc2dnaOrCol[c];
			q = -q;
			cedidx++;
		}
		// Determine next nucleotide by combining previous nucleotide and
		// current color
		int n = nuccol2nuc[lastn][c];
		c = ns[i+1];
		ns.set(n, i+1);
		
		int dq = max(q + lastq, 0);
		dq = min(dq, 127);
		lastq = q;
		q = qs[i+1]-33;
		qs.set(dq+33, i);
		lastn = n;
	}
	ns.set(nup, 0);
	int dq = max(lastq, 0)+33;
	qs.set(min(dq, 127) , len-1);
	assert_eq(ndn, ns[rdext]);
	assert_eq(cedidx, ced_.size());
	if(!fw) {
		// Need to re-flip edits to make them 5'-to-3' again.
		Edit::invertPoss(const_cast<EList<Edit>& >(ced_), rdextent_);
	}
}

/**
 * Initialize new AlnRes.
 */
void AlnRes::init(
	size_t             rdlen,           // # chars after hard trimming
	AlnScore           score,           // alignment score
	const EList<Edit>* ned,             // nucleotide edits
	const EList<Edit>* aed,             // ambiguous base resolutions
	const EList<Edit>* ced,             // color edits
	Coord              refcoord,        // leftmost ref pos of 1st al char
	bool               color,           // colorspace?
	int                seedmms,         // # seed mms allowed
	int                seedlen,         // seed length
	int                seedival,        // space between seeds
	int64_t            minsc,           // minimum score for valid aln
	int64_t            floorsc,         // local-alignment floor
	int                nuc5p,
	int                nuc3p,
	bool               pretrimSoft,
	size_t             pretrim5p,       // trimming prior to alignment
	size_t             pretrim3p,       // trimming prior to alignment
	bool               trimSoft,
	size_t             trim5p,          // trimming from alignment
	size_t             trim3p,          // trimming from alignment
	bool               cPretrimSoft,
	size_t             cPretrim5p,      // trimming prior to alignment
	size_t             cPretrim3p,      // trimming prior to alignment
	bool               cTrimSoft,
	size_t             cTrim5p,         // trimming from alignment
	size_t             cTrim3p)         // trimming from alignment
{
	rdlen_  = rdlen;
	rdrows_ = rdlen + (color ? 1 : 0);
	score_  = score;
	ned_.clear();
	aed_.clear();
	ced_.clear();
	if(ned != NULL) ned_ = *ned;
	if(aed != NULL) aed_ = *aed;
	if(ced != NULL) ced_ = *ced;
	refcoord_     = refcoord;
	color_        = color;
	seedmms_      = seedmms;
	seedlen_      = seedlen;
	seedival_     = seedival;
	minsc_        = minsc;
	floorsc_      = floorsc;
	nuc5p_        = nuc5p;
	nuc3p_        = nuc3p;
	pretrimSoft_  = pretrimSoft;
	pretrim5p_    = pretrim5p;
	pretrim3p_    = pretrim3p;
	trimSoft_     = trimSoft;
	trim5p_       = trim5p;
	trim3p_       = trim3p;
	cPretrimSoft_ = cPretrimSoft;
	cPretrim5p_   = cPretrim5p;
	cPretrim3p_   = cPretrim3p;
	cTrimSoft_    = cTrimSoft;
	cTrim5p_      = cTrim5p;
	cTrim3p_      = cTrim3p;
	rdextent_     = rdlen;      // # read characters after any hard trimming
	if(pretrimSoft) {
		rdextent_ -= (pretrim5p + pretrim3p);
	}
	if(trimSoft) {
		rdextent_ -= (trim5p + trim3p);
	}
	rdexrows_ = rdextent_ + (color ? 1 : 0);
	calcRefExtent();
	shapeSet_ = true;
}

/**
 * Shift all edits left-to-right along the Watson strand by the given amount.
 */
void AlnRes::clipLeft(TRefOff amt) {
#if 0
	ned_;          // base edits
	aed_;          // ambiguous base resolutions
	ced_;          // color miscalls
	refcoord_.adjustOff(amt); // shift
	// We have to determine the shape of the portion of the alignment that
	// we're clipping
	size_t      rdextent_;     // number of read chars involved in alignment
	size_t      rdexrows_;     // number of read rows involved in alignment
	size_t      rfextent_;     // number of ref chars involved in alignment

	int         nuc5p_;        // 5'-most decoded base; clipped if excluding end
	int         nuc3p_;        // 3'-most decoded base; clipped if excluding end
	size_t      refns_;        // # of reference Ns overlapped
	
	// Nucleotide-sequence trimming
	bool        trimSoft_;     // trimming by local alignment is soft?
	size_t      trim5p_;       // # bases trimmed from 5p end by local alignment
	size_t      trim3p_;       // # bases trimmed from 3p end by local alignment

	// Colorspace-sequence trimming; only relevant in colorspace
	bool        cPretrimSoft_; // trimming prior to alignment is soft?
	size_t      cPretrim5p_;   // # bases trimmed from 5p end prior to alignment
	size_t      cPretrim3p_;   // # bases trimmed from 3p end prior to alignment
	bool        cTrimSoft_;    // trimming by local alignment is soft?
	size_t      cTrim5p_;      // # bases trimmed from 5p end by local alignment
	size_t      cTrim3p_;      // # bases trimmed from 3p end by local alignment
#endif
}

/**
 * Shift all edits left-to-right along the Watson strand by the given amount.
 */
void AlnRes::clipRight(TRefOff amt) {
}

/**
 * Clip away portions of the alignment that are outside the given bounds.
 * Clipping is soft if soft == true, hard otherwise.  Assuming for now that
 * there isn't any other clipping.
 */
void AlnRes::clipOutside(bool soft, TRefOff refi, TRefOff reff) {
	// Overhang on LHS
	TRefOff left = refcoord_.off();
	if(left < refi) {
		clipLeft(refi - refcoord_.off());
	}
	// Overhang on RHS
	TRefOff right = refcoord_.off() + refNucExtent();
	if(right > reff) {
		clipRight(right - reff);
	}
}

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
 *
 * We iterate over each read position that hasn't been hard-trimmed, but
 * only overlaps at positions that have also not been soft-trimmed are
 * considered.
 */
bool AlnRes::overlap(AlnRes& res) {
	if(fw() != res.fw() || refid() != res.refid()) {
		// Must be same reference and same strand in order to overlap
		return false;
	}
	TRefOff my_left     = refoff();     // my leftmost aligned char
	TRefOff other_left  = res.refoff(); // other leftmost aligned char
	TRefOff my_right    = my_left    + refExtent();
	TRefOff other_right = other_left + res.refExtent();
	if(my_right < other_left || other_right < my_left) {
		// The rectangular hulls of the two alignments don't overlap, so
		// they can't overlap at any cell
		return false;
	}
	// Reference and strand are the same and hulls overlap.  Now go read
	// position by read position testing if any align identically with the
	// reference.
	
	// Edits are ordered and indexed from 5' to 3' to start with.  We
	// reorder them to go from left to right along the Watson strand.
	if(!fw()) {
		invertEdits();
	}
	if(!res.fw()) {
		res.invertEdits();
	}
	size_t nedidx = 0, onedidx = 0;
	bool olap = false;
	// For each row, going left to right along Watson reference strand...
	for(size_t i = 0; i < rdexrows_; i++) {
		size_t fivep = i;
		if(!fw()) fivep = rdexrows_ - i - 1;
		size_t diff = 1;  // amount to shift to right for next round
		size_t odiff = 1; // amount to shift to right for next round
		// Unless there are insertions before the next position, we say
		// that there is one cell in this row involved in the alignment
		my_right = my_left + 1;
		other_right = other_left + 1;
		while(nedidx < ned_.size() && ned_[nedidx].pos == i) {
			if(ned_[nedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				diff = 0;
			}
			nedidx++;
		}
		while(onedidx < res.ned_.size() && res.ned_[onedidx].pos == i) {
			if(res.ned_[onedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				odiff = 0;
			}
			onedidx++;
		}
		if(i < rdexrows_ - 1) {
			// See how many inserts there are before the next read
			// character
			size_t nedidx_next  = nedidx;
			size_t onedidx_next = onedidx;
			while(nedidx_next < ned_.size() &&
				  ned_[nedidx_next].pos == i+1)
			{
				if(ned_[nedidx_next].isReadGap()) {
					my_right++;
				}
				nedidx_next++;
			}
			while(onedidx_next < res.ned_.size() &&
				  res.ned_[onedidx_next].pos == i+1)
			{
				if(res.ned_[onedidx_next].isReadGap()) {
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

#ifndef NDEBUG

/**
 * Assuming this AlnRes is an alignment for 'rd', check that the alignment and
 * 'rd' are compatible with the corresponding reference sequence.
 */
bool AlnRes::matchesRef(
	const Read& rd,
	const BitPairReference& ref,
	BTDnaString& rf,
	BTDnaString& rdseq,
	BTString& qseq,
	SStringExpandable<char>& raw_refbuf,
	SStringExpandable<uint32_t>& destU32,
	EList<bool>& matches)
{
	assert(!empty());
	assert(repOk());
	assert(refcoord_.valid());
	bool fw = refcoord_.fw();
	// Adjust reference string length according to edits
	size_t refallen = refNucExtent();
	raw_refbuf.resize(refallen + 16);
	raw_refbuf.clear();
	int nsOnLeft = 0;
	if(refcoord_.off() < 0) {
		nsOnLeft = -((int)refcoord_.off());
	}
	int off = ref.getStretch(
		reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
		refcoord_.ref(),
		max<TRefOff>(refcoord_.off(), 0),
		refallen,
		destU32);
	assert_leq(off, 16);
	char *refbuf = raw_refbuf.wbuf() + off;
	size_t trimBeg = 0, trimEnd = 0;
	if(trimSoft_) {
		trimBeg += (fw ? trim5p_ : trim3p_);
		trimEnd += (fw ? trim3p_ : trim5p_);
	}
	if(pretrimSoft_) {
		trimBeg += (fw ? pretrim5p_ : pretrim3p_);
		trimEnd += (fw ? pretrim3p_ : pretrim5p_);
	}
	// All decoded-nucleotide trimming for colorspace alignments is hard
	// trimming
	assert(!color_ || trimBeg == 0);
	assert(!color_ || trimEnd == 0);
	rf.clear();
	rdseq.clear();
	if(rd.color) {
		// Decode the nucleotide sequence from the alignment
		qseq.clear();
		decodedNucsAndQuals(rd, rdseq, qseq);
		assert_eq(rdexrows_, rdseq.length());
		assert_eq(rdseq.length(), qseq.length());
	} else {
		rdseq = rd.patFw;
		if(!fw) {
			rdseq.reverseComp(false);
		}
		assert_eq(rdrows_, rdseq.length());
	}
	if(!fw) {
		// Invert the nucleotide edits so that they go from upstream to
		// downstream on the Watson strand
		Edit::invertPoss(ned_, rdexrows_);
	}
	// rdseq is the nucleotide sequence (decoded in the case of a
	// colorspace read) from upstream to downstream on the Watson
	// strand.  ned_ are the nucleotide edits from upstream to
	// downstream.  rf contains the reference characters.
	Edit::toRef(rdseq, ned_, rf, trimBeg, trimEnd);
	if(!fw) {
		// Re-invert the nucleotide edits so that they go from 5' to 3'
		Edit::invertPoss(ned_, rdexrows_);
	}
	assert_eq(refallen, rf.length());
	matches.clear();
	bool matchesOverall = true;
	matches.resize(refallen);
	matches.fill(true);
	for(size_t i = 0; i < refallen; i++) {
		if((int)i < nsOnLeft) {
			if((int)rf[i] != 4) {
				matches[i] = false;
				matchesOverall = false;
			}
		} else {
			if((int)rf[i] != (int)refbuf[i-nsOnLeft]) {
				matches[i] = false;
				matchesOverall = false;
			}
		}
	}
	if(!matchesOverall) {
		// Print a friendly message showing the difference between the
		// reference sequence obtained with Edit::toRef and the actual
		// reference sequence
		cerr << endl;
		Edit::printQAlignNoCheck(
			cerr,
			"    ",
			rdseq,
			ned_);
		cerr << "    ";
		for(size_t i = 0; i < refallen; i++) {
			cerr << (matches[i] ? " " : "*");
		}
		cerr << endl;
		Edit::printQAlign(
			cerr,
			"    ",
			rdseq,
			ned_);
		cerr << endl;
	}
	return matchesOverall;
}

#endif /*ndef NDEBUG*/

#define COPY_BUF() { \
	char *bufc = buf; \
	while(*bufc != '\0') { \
		*occ = *bufc; \
		occ++; \
		bufc++; \
	} \
}

/**
 * Print a CIGAR-string representation of the alignment.  In the
 * CIGAR-string representation, edit operations are printed in an order
 * that corresponds to moving left-to-right along the Watson strand.  The
 * operators indicate one of: match, mismatch, read gap, and reference gap.
 * With each operator is an associated run length (printed prior to the
 * operator) indicating how many times in a row that feature occurs.
 */
void AlnRes::printCigar(
	bool printColors,     // print CIGAR for colorspace alignment?
	bool exEnds,          // exclude ends in CIGAR?
	bool distinguishMm,   // use =/X instead of just M
	EList<char>& op,      // stick CIGAR operations here
	EList<size_t>& run,   // stick CIGAR run lengths here
	OutFileBuf *o,        // write to this buf if o != NULL
	char *oc) const       // write to this buf if oc != NULL
{
	char *occ = oc;
	op.clear();
	run.clear();
	// Any hard or soft clipping on the Beginning?
	size_t trimHardBeg = 0, trimSoftBeg = 0;
	size_t trimHardEnd = 0, trimSoftEnd = 0;
	if(pretrimSoft_) {
		trimSoftBeg = fw() ? pretrim5p_ : pretrim3p_;
		trimSoftEnd = fw() ? pretrim3p_ : pretrim5p_;
	} else {
		trimHardBeg = fw() ? pretrim5p_ : pretrim3p_;
		trimHardEnd = fw() ? pretrim3p_ : pretrim5p_;
	}
	if(trimSoft_) {
		trimSoftBeg += fw() ? trim5p_ : trim3p_;
		trimSoftEnd += fw() ? trim3p_ : trim5p_;
	} else {
		trimHardBeg += fw() ? trim5p_ : trim3p_;
		trimHardEnd += fw() ? trim3p_ : trim5p_;
	}
	// Print hard clipping
	if(trimHardBeg > 0) {
		op.push_back('H');
		run.push_back(trimHardBeg);
	}
	// Print soft clipping
	if(trimSoftBeg > 0) {
		op.push_back('S');
		run.push_back(trimSoftBeg);
	}
	// Go through edits from front to back
	if(!fw()) {
		const_cast<AlnRes*>(this)->invertEdits();
	}
	const EList<Edit>& ed = (printColors ? ced_ : ned_);
	size_t last = 0;
	for(size_t i = 0; i < ed.size(); i++) {
		if(ed[i].isMismatch() && !distinguishMm) {
			// If we're not distinguishing matches from mismatches, ignore
			// mismatches here
			continue;
		}
		// Print previous run of matches
		if(ed[i].pos > last) {
			// There's a run of matches prior to this edit.  
			op.push_back(distinguishMm ? '=' : 'M');
			run.push_back(ed[i].pos - last);
		}
		last = ed[i].pos;
		// Print edit
		if(ed[i].isMismatch()) {
			size_t len = 1;
			last++;
			// Mismatches at successive positions?
			while(i+1 < ed.size()) {
				if(ed[i+1].isMismatch() && ed[i+1].pos == ed[i].pos+1) {
					len++;  // increment length of mismatch run
					i++;    // move to next edit
					last++; // adjust beginning of next run
				} else {
					break;
				}
			}
			op.push_back('X');
			run.push_back(len);
		} else if(ed[i].isRefGap()) {
			size_t len = 1;
			last++;
			// Deletes at successive positions?
			while(i+1 < ed.size()) {
				if(ed[i+1].isRefGap() && ed[i+1].pos == ed[i].pos+1) {
					len++;  // increment length of deletion run
					i++;    // move to next edit
					last++; // adjust beginning of next run
				} else {
					break;
				}
			}
			op.push_back('I');
			run.push_back(len);
		} else if(ed[i].isReadGap()) {
			size_t len = 1;
			// Deletes at successive positions?
			while(i+1 < ed.size()) {
				if(ed[i+1].isReadGap() && ed[i+1].pos == ed[i].pos) {
					len++;  // increment length of deletion run
					i++;    // move to next edit
				} else {
					break;
				}
			}
			op.push_back('D');
			run.push_back(len);
		}
	}
	size_t end = printColors ? rdrows_ : rdexrows_;
	if(last < end) {
		// There's a run of matches prior to the end.
		op.push_back(distinguishMm ? '=' : 'M');
		run.push_back(end - last);
	}
	if(!fw()) {
		const_cast<AlnRes*>(this)->invertEdits();
	}
	// Print soft clipping
	if(trimSoftEnd) {
		op.push_back('S');
		run.push_back(trimSoftEnd);
	}
	// Print hard clipping
	if(trimHardEnd) {
		op.push_back('H');
		run.push_back(trimHardEnd);
	}
	// Write to the output file buffer and/or string buffer.
	assert_eq(op.size(), run.size());
	if(o != NULL || oc != NULL) {
		char buf[128];
		bool printed = false;
		for(size_t i = 0; i < op.size(); i++) {
			bool first = (i == 0);
			bool last  = (i == op.size()-1);
			size_t r = run[i];
			if(first && exEnds && r > 0) {
				r--;
			}
			if(last && exEnds && r > 0) {
				r--;
			}
			if(r > 0) {
				itoa10<size_t>(r, buf);
				printed = true;
				if(o != NULL) {
					o->writeChars(buf);
					o->write(op[i]);
				}
				if(oc != NULL) {
					COPY_BUF();
					*occ = op[i];
					occ++;
				}
			}
		}
		assert(printed);
		if(oc != NULL) {
			*occ = '\0';
		}
	}
}

/**
 * Print a MD:Z:-string representation of the alignment, a la BWA.  In this
 * representation runs of either matches or reference gaps are represented
 * by a single number indicating the length of the run.  Mismatches are
 * indicated by the DNA character that occurs in the reference part of the
 * mismatch.  Read gaps are indicated by a carat (^) followed by the string
 * of reference characters that occur in the gap.  If a mismatch follows a read
 * gap, the read gap string (e.g. "^AAG") and the mismatch string (e.g. "T")
 * are separated by a "0" (e.g. "^AAAG0T").  Also, if a mismatch occurs at
 * either end, the end is capped with a "0".
 */
void AlnRes::printMD(
	bool printColors,     // print CIGAR for colorspace alignment?
	bool exEnds,          // exclude ends nucleotides for decoded nucs?
	EList<char>& op,      // stick operations here
	EList<char>& ch,      // stick reference characters here
	EList<size_t>& run,   // stick run lengths here
	OutFileBuf* o,        // write to this buf if o != NULL
	char* oc) const       // write to this buf if oc != NULL
{
	char *occ = oc;
	op.clear();
	ch.clear();
	run.clear();
	// Go through edits from front to back
	if(!fw()) {
		const_cast<AlnRes*>(this)->invertEdits();
	}
	const EList<Edit>& ed = (printColors ? ced_ : ned_);
	size_t last = 0;
	for(size_t i = 0; i < ed.size(); i++) {
		// Ignore ref gaps
		if(ed[i].isRefGap()) {
			// Ref gaps take up rows in the DP table, but don't count toward
			// run length
			last++;
			continue;
		}
		// Print previous run of matches
		if(ed[i].pos > last) {
			// There's a run of matches prior to this edit.  
			op.push_back('=');
			ch.push_back('-');
			run.push_back(ed[i].pos - last);
		}
		last = ed[i].pos;		// Print edit
		if(ed[i].isMismatch()) {
			last++;
			op.push_back('X');
			ch.push_back(ed[i].chr);
			assert_neq('-', ed[i].chr);
			run.push_back(1);
		} else if(ed[i].isReadGap()) {
			op.push_back('G');
			ch.push_back(ed[i].chr);
			assert_neq('-', ed[i].chr);
			run.push_back(1);
		}
	}
	size_t end = printColors ? rdrows_ : rdexrows_;
	if(last < end) {
		// There's a run of matches prior to the end.
		op.push_back('=');
		ch.push_back('-');
		run.push_back(end - last);
	}
	if(!fw()) {
		const_cast<AlnRes*>(this)->invertEdits();
	}
	// Write to the output file buffer and/or string buffer.
	assert_eq(op.size(), run.size());
	assert_eq(op.size(), ch.size());
	if(o != NULL || oc != NULL) {
		char buf[128];
		bool mm_last = false;
		bool rdgap_last = false;
		bool first_print = true;
		for(size_t i = 0; i < op.size(); i++) {
			bool first = (i == 0);
			bool last  = (i == op.size()-1);
			size_t r = run[i];
			if(first && exEnds && r > 0) {
				r--;
			}
			if(last && exEnds && r > 0) {
				r--;
			}
			if(r > 0) {
				if(op[i] == '=') {
					itoa10<size_t>(r, buf);
					if(o != NULL) {
						o->writeChars(buf);
					}
					if(oc != NULL) {
						COPY_BUF();
					}
					first_print = false;
					mm_last = false;
					rdgap_last = false;
				} else if(op[i] == 'X') {
					if(o != NULL) {
						if(rdgap_last || first_print) {
							o->write('0');
						}
						o->write(ch[i]);
					}
					if(oc != NULL) {
						if(rdgap_last || first_print) {
							*occ = '0';
							occ++;
						}
						*occ = ch[i];
						occ++;
					}
					first_print = false;
					mm_last = true;
					rdgap_last = false;
				} else if(op[i] == 'G') {
					if(o != NULL) {
						if(first_print) {
							o->write('0');
						}
						if(!rdgap_last) {
							o->write('^');
						}
						o->write(ch[i]);
					}
					if(oc != NULL) {
						if(first_print) {
							*occ = '0';
							occ++;
						}
						if(!rdgap_last) {
							*occ = '^';
							occ++;
						}
						*occ = ch[i];
						occ++;
					}
					first_print = false;
					mm_last = false;
					rdgap_last = true;
				}
			} // if r > 0
		} // for loop over ops
		if(mm_last || rdgap_last) {
			if(o != NULL) {
				o->write('0');
			}
			if(oc != NULL) {
				*occ = '0';
				occ++;
			}
		}
		if(oc != NULL) {
			*occ = '\0';
		}
	}
}

/**
 * Print the sequence for the read that aligned using A, C, G and
 * T.  This will simply print the read sequence (or its reverse
 * complement) unless this is a colorspace read and printColors is
 * false.  In that case, we print the decoded sequence rather than
 * the original ones.
 */
void AlnRes::printSeq(
	const Read& rd,         // read
	const BTDnaString* dns, // already-decoded nucleotides
	bool printColors,       // print colors instead of decoded nucleotides?
	bool exEnds,            // exclude ends when printing decoded nucleotides?
	OutFileBuf& o) const    // output stream to write to
{
	assert(!rd.patFw.empty());
	bool fw = refcoord_.fw();
	assert(!printColors || rd.color);
	ASSERT_ONLY(size_t written = 0);
	if(!rd.color || printColors) {
		// Should not have had hard clipping during alignment
		assert(trimSoft_ || (trim3p_ + trim5p_ == 0));
		// Print nucleotides or colors
		size_t len = rd.patFw.length();
		for(size_t i = 0; i < len; i++) {
			int c;
			if(fw) {
				c = rd.patFw[i];
			} else {
				// Reverse-complement
				c = rd.patFw[len-i-1];
				if(c < 4) c = c ^ 3;
			}
			assert_range(0, 4, c);
			o.write("ACGTN"[c]);
			ASSERT_ONLY(written++);
		}
#ifndef NDEBUG
		for(size_t i = 0; i < ced_.size(); i++) {
			if(ced_[i].isReadGap()) {
				assert_leq(ced_[i].pos, written);
			} else {
				assert_lt(ced_[i].pos, written);
			}
		}
#endif
	} else {
		// Print decoded nucleotides
		assert(dns != NULL);
		size_t len = dns->length();
		size_t st = 0;
		size_t en = len;
		if(exEnds) {
			st++; en--;
		}
		for(size_t i = st; i < en; i++) {
			int c = dns->get(i);
			assert_range(0, 3, c);
			o.write("ACGT"[c]);
			ASSERT_ONLY(written++);
		}
#ifndef NDEBUG
		for(size_t i = 0; i < ned_.size(); i++) {
			if(ned_[i].isReadGap()) {
				assert_leq(ned_[i].pos, dns->length());
			} else {
				assert_lt(ned_[i].pos, dns->length());
			}
		}
#endif
	}
}

/**
 * Print the quality string for the read that aligned.  This will simply print
 * the read qualities (or their reverse) unless this is a colorspace read and
 * printColors is false.  In that case, we print the decoded qualities rather
 * than the original ones.
 */
void AlnRes::printQuals(
	const Read& rd,         // read
	const BTString* dqs,    // already-decoded qualities
	bool printColors,       // true -> print colors instead of decoded nucleotides for colorspace alignment
	bool exEnds,            // true -> exclude ends when printing decoded nucleotides
	OutFileBuf& o) const    // output stream to write to
{
	bool fw = refcoord_.fw();
	assert(!printColors || rd.color);
	if(!rd.color || printColors) {
		size_t len = rd.qual.length();
		// Print original qualities from upstream to downstream Watson
		for(size_t i = 0; i < len; i++) {
			int c = (fw ? rd.qual[i] : rd.qual[len-i-1]);
			o.write(c);
		}
	} else {
		assert(dqs != NULL);
		size_t len = dqs->length();
		// Print decoded qualities from upstream to downstream Watson
		if(!exEnds) {
			// Print upstream-most quality
			o.write(dqs->get(0));
		}
		for(size_t i = 1; i < len-1; i++) {
			o.write(dqs->get(i));
		}
		if(!exEnds) {
			// Print downstream-most quality
			o.write(dqs->get(len-1));
		}
	}
}

/**
 * Add all of the cells involved in the given alignment to the database.
 */
void RedundantAlns::add(const AlnRes& res) {
	assert(!cells_.empty());
	TRefOff left = res.refoff(), right;
	const size_t len = res.readExtentRows();
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
	const EList<Edit>& ned = res.ned();
	size_t nedidx = 0;
	assert_leq(len, cells_.size());
	// For each row...
	for(size_t i = 0; i < len; i++) {
		size_t diff = 1;  // amount to shift to right for next round
		right = left + 1;
		while(nedidx < ned.size() && ned[nedidx].pos == i) {
			if(ned[nedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				diff = 0;
			}
			nedidx++;
		}
		if(i < len - 1) {
			// See how many inserts there are before the next read
			// character
			size_t nedidx_next = nedidx;
			while(nedidx_next < ned.size() && ned[nedidx_next].pos == i+1)
			{
				if(ned[nedidx_next].isReadGap()) {
					right++;
				}
				nedidx_next++;
			}
		}
		for(TRefOff j = left; j < right; j++) {
			// Add to db
			RedundantCell c(res.refid(), res.fw(), j, i);
			ASSERT_ONLY(bool ret =) cells_[i].insert(c);
			assert(ret);
		}
		left = right + diff - 1;
	}
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
}

/**
 * Return true iff the given alignment has at least one cell that overlaps
 * one of the cells in the database.
 */
bool RedundantAlns::overlap(const AlnRes& res) {
	assert(!cells_.empty());
	TRefOff left = res.refoff(), right;
	const size_t len = res.readExtentRows();
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
	const EList<Edit>& ned = res.ned();
	size_t nedidx = 0;
	// For each row...
	bool olap = false;
	assert_leq(len, cells_.size());
	for(size_t i = 0; i < len; i++) {
		size_t diff = 1;  // amount to shift to right for next round
		right = left + 1;
		while(nedidx < ned.size() && ned[nedidx].pos == i) {
			if(ned[nedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				diff = 0;
			}
			nedidx++;
		}
		if(i < len - 1) {
			// See how many inserts there are before the next read
			// character
			size_t nedidx_next = nedidx;
			while(nedidx_next < ned.size() && ned[nedidx_next].pos == i+1)
			{
				if(ned[nedidx_next].isReadGap()) {
					right++;
				}
				nedidx_next++;
			}
		}
		for(TRefOff j = left; j < right; j++) {
			// Add to db
			RedundantCell c(res.refid(), res.fw(), j, i);
			if(cells_[i].contains(c)) {
				olap = true;
				break;
			}
		}
		if(olap) {
			break;
		}
		left = right + diff - 1;
	}
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
	return olap;
}

/**
 * Given an unpaired read (in either rd1 or rd2) or a read pair
 * (mate 1 in rd1, mate 2 in rd2).
 */
AlnSetSumm::AlnSetSumm(
	const Read* rd1,
	const Read* rd2,
	const EList<AlnRes>* rs1,
	const EList<AlnRes>* rs2)
{
	assert(rd1 != NULL || rd2 != NULL);
	AlnScore best, secbest;
	best.invalidate();
	secbest.invalidate();
	bool paired = (rs1 != NULL && rs2 != NULL);
	bool noResult = (rs1 == NULL && rs2 == NULL);
	size_t sz;
	if(paired) {
		// Paired alignments
		assert_eq(rs1->size(), rs2->size());
		sz = rs1->size();
		assert_gt(sz, 0);
		for(size_t i = 0; i < rs1->size(); i++) {
			AlnScore sc = (*rs1)[i].score() + (*rs2)[i].score();
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
		init(best, secbest, sz-1);
		assert(!empty());
	} else if(!noResult) {
		// Unpaired alignments
		const EList<AlnRes>* rs = (rs1 != NULL ? rs1 : rs2);
		assert(rs != NULL);
		sz = rs->size();
		assert_gt(sz, 0);
		for(size_t i = 0; i < rs->size(); i++) {
			AlnScore sc = (*rs)[i].score();
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
		init(best, secbest, sz-1);
		assert(!empty());
	} else {
		// No result - leave best and secbest as invalid
		init(best, secbest, 0);
		assert(empty());
	}
}

/**
 * Print out string representation of these flags.
 */
void AlnFlags::printYM(OutFileBuf& o) const {
	o.writeChars("YM:i:");
	o.write(maxed() ? '1' : '0');
}

/**
 * Print out string representation of these flags.
 */
void AlnFlags::printYP(OutFileBuf& o) const {
	o.writeChars("YP:i:");
	o.write(maxedPair() ? '1' : '0');
}

/**
 * Print out string representation of these flags.
 */
void AlnFlags::printYT(OutFileBuf& o) const {
	o.writeChars("YT:Z:");
	if(pairing() == ALN_FLAG_PAIR_CONCORD) {
		o.writeChars("CP");
	} else if(pairing() == ALN_FLAG_PAIR_DISCORD) {
		o.writeChars("DP");
	} else if(pairing() == ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR) {
		o.writeChars("UP");
	} else if(pairing() == ALN_FLAG_PAIR_UNPAIRED) {
		o.writeChars("UU");
	} else { throw 1; }
}

#ifdef ALIGNER_RESULT_MAIN

#include "mem_ids.h"

int main() {
	EList<char> op;
	EList<char> ch;
	EList<size_t> run;
	{
		// On top of each other, same length
		cerr << "Test case 1, simple overlap 1 ... ";
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		assert(res1.overlap(res2));
		
		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "10M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "10="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "10M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "10="));

		char buf3[1024];
		res1.printMD(false, false, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "10"));
		res1.printMD(false, true, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "8"));

		char buf4[1024];
		res2.printMD(false, false, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "10"));
		res2.printMD(false, true, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "8"));

		cerr << "PASSED" << endl;
	}

	{
		// On top of each other, different lengths
		cerr << "Test case 2, simple overlap 2 ... ";
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		AlnRes res2;
		res2.init(
			11,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		assert(res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(11);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "10M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "10="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "11M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "11="));

		char buf3[1024];
		res1.printMD(false, false, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "10"));
		res1.printMD(false, true, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "8"));

		char buf4[1024];
		res2.printMD(false, false, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "11"));
		res2.printMD(false, true, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "9"));

		cerr << "PASSED" << endl;
	}

	{
		// Different references
		cerr << "Test case 3, simple overlap 3 ... ";
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 1, true),
			false);
		AlnRes res2;
		res2.init(
			11,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(11);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		cerr << "PASSED" << endl;
	}

	{
		// Different references
		cerr << "Test case 4, simple overlap 4 ... ";
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(1, 0, true),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		cerr << "PASSED" << endl;
	}

	{
		// Different strands
		cerr << "Test case 5, simple overlap 5 ... ";
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, true),
			false);
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 0, false),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		cerr << "PASSED" << endl;
	}

	{
		// Different strands
		cerr << "Test case 6, simple overlap 6 ... ";
		EList<Edit> ned1(RES_CAT);
		ned1.expand();
		// 1 step to the right in the middle of the alignment
		ned1.back().init(5, 'A' /*chr*/, '-' /*qchr*/, EDIT_TYPE_READ_GAP);
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, false),
			false);
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 6, false),
			false);
		assert(res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M1D5M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=1D5="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "10M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "10="));

		char buf3[1024];
		res1.printMD(false, false, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "5^A5"));
		res1.printMD(false, true, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "4^A4"));

		char buf4[1024];
		res2.printMD(false, false, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "10"));
		res2.printMD(false, true, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "8"));

		cerr << "PASSED" << endl;
	}

	{
		// Different strands
		cerr << "Test case 7, simple overlap 7 ... ";
		EList<Edit> ned1(RES_CAT);
		// 3 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, 'A', '-', EDIT_TYPE_READ_GAP));
		ned1.push_back(Edit(5, 'C', '-', EDIT_TYPE_READ_GAP));
		ned1.push_back(Edit(5, 'G', '-', EDIT_TYPE_READ_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, false),
			false);
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			NULL,
			NULL,
			NULL,
			Coord(0, 6, false),
			false);
		assert(res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M3D5M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=3D5="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "10M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "10="));

		char buf3[1024];
		res1.printMD(false, false, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "5^GCA5"));
		res1.printMD(false, true, op, ch, run, NULL, buf3);
		assert_eq(0, strcmp(buf3, "4^GCA4"));

		char buf4[1024];
		res2.printMD(false, false, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "10"));
		res2.printMD(false, true, op, ch, run, NULL, buf4);
		assert_eq(0, strcmp(buf4, "8"));

		cerr << "PASSED" << endl;
	}

	{
		// Both with horizontal movements; overlap
		cerr << "Test case 8, simple overlap 8 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, 'A', '-', EDIT_TYPE_READ_GAP));
		ned1.push_back(Edit(5, 'C', '-', EDIT_TYPE_READ_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, false),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(5, 'A', '-', EDIT_TYPE_READ_GAP));
		ned2.push_back(Edit(5, 'C', '-', EDIT_TYPE_READ_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 6, false),
			false);
		assert(res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M2D5M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=2D5="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "5M2D5M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "5=2D5="));

		cerr << "PASSED" << endl;
	}

	{
		// Both with horizontal movements; no overlap
		cerr << "Test case 9, simple overlap 9 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(6, 'A', '-', EDIT_TYPE_READ_GAP));
		ned1.push_back(Edit(6, 'C', '-', EDIT_TYPE_READ_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, true),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(5, 'A', '-', EDIT_TYPE_READ_GAP));
		ned2.push_back(Edit(5, 'C', '-', EDIT_TYPE_READ_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 6, true),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "6M2D4M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "6=2D4="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "5M2D5M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "5=2D5="));

		cerr << "PASSED" << endl;
	}

	{
		// Both with horizontal movements; no overlap.  Reverse strand.
		cerr << "Test case 10, simple overlap 10 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, 'A', '-', EDIT_TYPE_READ_GAP));
		ned1.push_back(Edit(5, 'C', '-', EDIT_TYPE_READ_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, false),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(6, 'A', '-', EDIT_TYPE_READ_GAP));
		ned2.push_back(Edit(6, 'C', '-', EDIT_TYPE_READ_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 6, false),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M2D5M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=2D5="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "4M2D6M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "4=2D6="));

		cerr << "PASSED" << endl;
	}

	{
		// Both with vertical movements; no overlap
		cerr << "Test case 11, simple overlap 11 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, '-', 'A', EDIT_TYPE_REF_GAP));
		ned1.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, true),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(6, '-', 'A', EDIT_TYPE_REF_GAP));
		ned2.push_back(Edit(7, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 6, true),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M2I3M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=2I3="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "6M2I2M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "6=2I2="));

		cerr << "PASSED" << endl;
	}

	{
		// Both with vertical movements; no overlap
		cerr << "Test case 12, simple overlap 12 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, '-', 'A', EDIT_TYPE_REF_GAP));
		ned1.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, true),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(5, '-', 'A', EDIT_TYPE_REF_GAP));
		ned2.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 6, true),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M2I3M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=2I3="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "5M2I3M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "5=2I3="));

		cerr << "PASSED" << endl;
	}

	{
		// Both with vertical movements; overlap
		cerr << "Test case 13, simple overlap 13 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, '-', 'A', EDIT_TYPE_REF_GAP));
		ned1.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, true),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(4, '-', 'A', EDIT_TYPE_REF_GAP));
		ned2.push_back(Edit(5, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 6, true),
			false);
		assert(res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(ra.overlap(res2));

		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M2I3M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=2I3="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "4M2I4M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "4=2I4="));

		cerr << "PASSED" << endl;
	}

	{
		// Not even close
		cerr << "Test case 14, simple overlap 14 ... ";
		EList<Edit> ned1(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned1.push_back(Edit(5, '-', 'A', EDIT_TYPE_REF_GAP));
		ned1.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res1;
		res1.init(
			10,
			AlnScore(),
			&ned1,
			NULL,
			NULL,
			Coord(0, 5, true),
			false);
		EList<Edit> ned2(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned2.push_back(Edit(4, '-', 'A', EDIT_TYPE_REF_GAP));
		ned2.push_back(Edit(5, '-', 'C', EDIT_TYPE_REF_GAP));
		AlnRes res2;
		res2.init(
			10,
			AlnScore(),
			&ned2,
			NULL,
			NULL,
			Coord(0, 400, true),
			false);
		assert(!res1.overlap(res2));

		// Try again, but using the redundant-alignment database
		RedundantAlns ra;
		ra.reset();
		ra.init(10);
		ra.add(res1);
		assert(ra.overlap(res1));
		assert(!ra.overlap(res2));
		
		char buf1[1024];
		res1.printCigar(false, false, false, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5M2I3M"));
		res1.printCigar(false, false, true, op, run, NULL, buf1);
		assert_eq(0, strcmp(buf1, "5=2I3="));

		char buf2[1024];
		res2.printCigar(false, false, false, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "4M2I4M"));
		res2.printCigar(false, false, true, op, run, NULL, buf2);
		assert_eq(0, strcmp(buf2, "4=2I4="));

		cerr << "PASSED" << endl;
	}

	{
		cerr << "Test case 15, CIGAR string with mismatches ... ";
		EList<Edit> ned(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned.push_back(Edit(0, 'C', 'A', EDIT_TYPE_MM));
		ned.push_back(Edit(4, '-', 'C', EDIT_TYPE_REF_GAP));
		ned.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		ned.push_back(Edit(7, '-', 'C', EDIT_TYPE_REF_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(10, '-', 'A', EDIT_TYPE_MM));
		AlnRes res; res.init(
			11,
			AlnScore(),
			&ned,
			NULL,
			NULL,
			Coord(0, 44, true),
			false);
		char buf[1024];
		res.printCigar(false, false, false, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "4M1I1M2I1M4D2M"));
		res.printCigar(false, false, true, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "1X3=1I1=2I1=4D1=1X"));
		cerr << "PASSED" << endl;
	}

	{
		cerr << "Test case 16, Single colorspace 1 ... ";
		EList<Edit> ned(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned.push_back(Edit(0, 'C', 'A', EDIT_TYPE_MM));
		ned.push_back(Edit(4, '-', 'C', EDIT_TYPE_REF_GAP));
		ned.push_back(Edit(6, '-', 'C', EDIT_TYPE_REF_GAP));
		ned.push_back(Edit(7, '-', 'C', EDIT_TYPE_REF_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(9, '-', 'A', EDIT_TYPE_READ_GAP));
		ned.push_back(Edit(10, '-', 'A', EDIT_TYPE_MM));
		AlnRes res; res.init(
			11,
			AlnScore(),
			&ned,
			NULL,
			NULL,
			Coord(0, 44, true),
			false);
		char buf[1024];
		// Don't exclude ends
		res.printCigar(false, false, false, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "4M1I1M2I1M4D2M"));
		res.printCigar(false, false, true, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "1X3=1I1=2I1=4D1=1X"));
		// Exclude ends
		res.printCigar(false, true, false, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "3M1I1M2I1M4D1M"));
		res.printCigar(false, true, true, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "3=1I1=2I1=4D1="));
		cerr << "PASSED" << endl;
	}

	{
		cerr << "Test case 17, Overhang ... ";
		EList<Edit> ned(RES_CAT);
		// 2 steps to the right in the middle of the alignment
		ned.push_back(Edit(0, 'N', 'A', EDIT_TYPE_MM));
		ned.push_back(Edit(5, 'C', 'A', EDIT_TYPE_MM));
		AlnRes res; res.init(
			10,
			AlnScore(),
			&ned,
			NULL,
			NULL,
			Coord(0, -1, true),
			false);
		
		char buf[1024];
		res.printCigar(false, false, false, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "10M"));
		res.printCigar(false, false, true, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "1X4=1X4="));
		res.printMD(false, false, op, ch, run, NULL, buf);
		assert_eq(0, strcmp(buf, "0N4C4"));
		
		#if 0
		AlnRes res2(res);
		// Now soft-clip away the overhang
		res2.clipOutside(
			true,  // soft clip
			0,     // ref begins
			40);   // ref ends (excl)
		res2.printCigar(false, false, false, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "1S9M"));
		res2.printCigar(false, false, true, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "4=1X4="));
		res2.printMD(false, false, op, ch, run, NULL, buf);
		assert_eq(0, strcmp(buf, "4C4"));

		AlnRes res3 = res;
		// Now hard-clip away the overhang
		res3.clipOutside(
			false, // hard clip
			0,     // ref begins
			40);   // ref ends (excl)
		res3.printCigar(false, false, false, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "9M"));
		res3.printCigar(false, false, true, op, run, NULL, buf);
		assert_eq(0, strcmp(buf, "4=1X4="));
		res3.printMD(false, false, op, ch, run, NULL, buf);
		assert_eq(0, strcmp(buf, "4C4"));
		#endif

		cerr << "PASSED" << endl;
	}
}

#endif /*def ALIGNER_RESULT_MAIN*/
