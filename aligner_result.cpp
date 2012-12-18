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

#include <iostream>
#include "reference.h"
#include "aligner_result.h"
#include "read.h"
#include "edit.h"
#include "sstring.h"
#include "ds.h"
#include "util.h"
#include "alphabet.h"

using namespace std;

/**
 * Clear all contents.
 */
void AlnRes::reset() {
	ned_.clear();
	aed_.clear();
	score_.invalidate();
	refcoord_.invalidate();
	refival_.invalidate();
	shapeSet_     = false;
	rdlen_        = 0;
	reflen_       = 0;
	rdrows_       = 0;
	rdextent_     = 0;
	rdexrows_     = 0;
	rfextent_     = 0;
	refns_        = 0;
	type_         = ALN_RES_TYPE_UNPAIRED;
	fraglen_      = -1;
	trimSoft_     = false;
	trim5p_       = 0;
	trim3p_       = 0;
	pretrimSoft_  = true;
	pretrim5p_    = 0;
	pretrim3p_    = 0;
	seedmms_      = 0; // number of mismatches allowed in seed
	seedlen_      = 0; // length of seed
	seedival_     = 0; // interval between seeds
	minsc_        = 0; // minimum score
	nuc5p_        = 0;
	nuc3p_        = 0;
	fraglenSet_   = false;
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
	TRefOff reflen,      // length of reference sequence aligned to
	bool    fw,          // aligned to Watson strand?
	size_t  rdlen,       // length of read after hard trimming, before soft
	bool    pretrimSoft, // whether trimming prior to alignment was soft
	size_t  pretrim5p,   // # poss trimmed form 5p end before alignment
	size_t  pretrim3p,   // # poss trimmed form 3p end before alignment
	bool    trimSoft,    // whether local-alignment trimming was soft
	size_t  trim5p,      // # poss trimmed form 5p end during alignment
	size_t  trim3p)      // # poss trimmed form 3p end during alignment
{
	rdlen_       = rdlen;
	rdrows_      = rdlen;
	refcoord_.init(id, off, fw);
	pretrimSoft_  = pretrimSoft;
	pretrim5p_    = pretrim5p;
	pretrim3p_    = pretrim3p;
	trimSoft_     = trimSoft;
	trim5p_       = trim5p;
	trim3p_       = trim3p;
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
			ned_[i].pos -= (uint32_t)trimBeg;
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
	rdexrows_ = rdextent_;
	calcRefExtent();
	refival_.init(id, off, fw, rfextent_);
	reflen_ = reflen;
	shapeSet_ = true;
}

/**
 * Initialize new AlnRes.
 */
void AlnRes::init(
	size_t             rdlen,           // # chars after hard trimming
	AlnScore           score,           // alignment score
	const EList<Edit>* ned,             // nucleotide edits
	size_t             ned_i,           // first position to copy
	size_t             ned_n,           // # positions to copy
	const EList<Edit>* aed,             // ambiguous base resolutions
	size_t             aed_i,           // first position to copy
	size_t             aed_n,           // # positions to copy
	Coord              refcoord,        // leftmost ref pos of 1st al char
	TRefOff            reflen,          // length of ref aligned to
	int                seedmms,         // # seed mms allowed
	int                seedlen,         // seed length
	int                seedival,        // space between seeds
	int64_t            minsc,           // minimum score for valid aln
	int                nuc5p,
	int                nuc3p,
	bool               pretrimSoft,
	size_t             pretrim5p,       // trimming prior to alignment
	size_t             pretrim3p,       // trimming prior to alignment
	bool               trimSoft,
	size_t             trim5p,          // trimming from alignment
	size_t             trim3p)          // trimming from alignment
{
	rdlen_  = rdlen;
	rdrows_ = rdlen;
	score_  = score;
	ned_.clear();
	aed_.clear();
	if(ned != NULL) {
		for(size_t i = ned_i; i < ned_i + ned_n; i++) {
			ned_.push_back((*ned)[i]);
		}
	}
	if(aed != NULL) {
		for(size_t i = aed_i; i < aed_i + aed_n; i++) {
			aed_.push_back((*aed)[i]);
		}
	}
	refcoord_     = refcoord;
	reflen_       = reflen;
	seedmms_      = seedmms;
	seedlen_      = seedlen;
	seedival_     = seedival;
	minsc_        = minsc;
	nuc5p_        = nuc5p;
	nuc3p_        = nuc3p;
	pretrimSoft_  = pretrimSoft;
	pretrim5p_    = pretrim5p;
	pretrim3p_    = pretrim3p;
	trimSoft_     = trimSoft;
	trim5p_       = trim5p;
	trim3p_       = trim3p;
	rdextent_     = rdlen;      // # read characters after any hard trimming
	if(pretrimSoft) {
		rdextent_ -= (pretrim5p + pretrim3p);
	}
	if(trimSoft) {
		rdextent_ -= (trim5p + trim3p);
	}
	rdexrows_ = rdextent_;
	calcRefExtent();
	setShape(
		refcoord.ref(), // id of reference aligned to
		refcoord.off(), // offset of first aligned char into ref seq
		reflen,         // length of reference sequence aligned to
		refcoord.fw(),  // aligned to Watson strand?
		rdlen,          // length of read after hard trimming, before soft
		pretrimSoft,    // whether trimming prior to alignment was soft
		pretrim5p,      // # poss trimmed form 5p end before alignment
		pretrim3p,      // # poss trimmed form 3p end before alignment
		trimSoft,       // whether local-alignment trimming was soft
		trim5p,         // # poss trimmed form 5p end during alignment
		trim3p);        // # poss trimmed form 3p end during alignment
	shapeSet_ = true;
}

/**
 * Clip given number of characters from the Watson-upstream end of the
 * alignment.
 */
void AlnRes::clipLeft(size_t rd_amt, size_t rf_amt) {
	assert_geq(rd_amt, 0);
	assert_geq(rf_amt, 0);
	assert_leq(rd_amt, rdexrows_);
	assert_leq(rf_amt, rfextent_);
	assert(trimSoft_);
	if(fw()) {
		trim5p_ += rd_amt;
		Edit::clipLo(ned_, rdexrows_, rd_amt);
		Edit::clipLo(aed_, rdexrows_, rd_amt);
	} else {
		trim3p_ += rd_amt;
		Edit::clipHi(ned_, rdexrows_, rd_amt);
		Edit::clipHi(aed_, rdexrows_, rd_amt);
	}
	rdexrows_ -= rd_amt;
	rdextent_ -= rd_amt;
	rfextent_ -= rf_amt;
	refcoord_.adjustOff(rf_amt);
	refival_.adjustOff(rf_amt);
	// Adjust refns_?
}

/**
 * Clip given number of characters from the Watson-downstream end of the
 * alignment.
 */
void AlnRes::clipRight(size_t rd_amt, size_t rf_amt) {
	assert_geq(rd_amt, 0);
	assert_geq(rf_amt, 0);
	assert_leq(rd_amt, rdexrows_);
	assert_leq(rf_amt, rfextent_);
	assert(trimSoft_);
	if(fw()) {
		trim3p_ += rd_amt;
		Edit::clipHi(ned_, rdexrows_, rd_amt);
		Edit::clipHi(aed_, rdexrows_, rd_amt);
	} else {
		trim5p_ += rd_amt;
		Edit::clipLo(ned_, rdexrows_, rd_amt);
		Edit::clipLo(aed_, rdexrows_, rd_amt);
	}
	rdexrows_ -= rd_amt;
	rdextent_ -= rd_amt;
	rfextent_ -= rf_amt;
	// Adjust refns_?
}

/**
 * Clip away portions of the alignment that are outside the given bounds.
 * Clipping is soft if soft == true, hard otherwise.  Assuming for now that
 * there isn't any other clipping.
 *
 * Note that all clipping is expressed in terms of read positions.  So if there
 * are reference gaps in the overhanging portion, we must 
 */
void AlnRes::clipOutside(bool soft, TRefOff refi, TRefOff reff) {
	// Overhang on LHS
	TRefOff left = refcoord_.off();
	if(left < refi) {
		size_t rf_amt = (size_t)(refi - left);
		size_t rf_i = rf_amt;
		size_t nedsz = ned_.size();
		if(!fw()) {
			Edit::invertPoss(ned_, rdexrows_, false);
		}
		for(size_t i = 0; i < nedsz; i++) {
			assert_lt(ned_[i].pos, rdexrows_);
			if(ned_[i].pos > rf_i) break;
			if(ned_[i].isRefGap()) rf_i++;
		}
		if(!fw()) {
			Edit::invertPoss(ned_, rdexrows_, false);
		}
		clipLeft(rf_i, rf_amt);
	}
	// Overhang on RHS
	TRefOff right = refcoord_.off() + refNucExtent();
	if(right > reff) {
		size_t rf_amt = (size_t)(right - reff);
		size_t rf_i = rf_amt;
		size_t nedsz = ned_.size();
		if(fw()) {
			Edit::invertPoss(ned_, rdexrows_, false);
		}
		for(size_t i = 0; i < nedsz; i++) {
			assert_lt(ned_[i].pos, rdexrows_);
			if(ned_[i].pos > rf_i) break;
			if(ned_[i].isRefGap()) rf_i++;
		}
		if(fw()) {
			Edit::invertPoss(ned_, rdexrows_, false);
		}
		clipRight(rf_i, rf_amt);
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
	assert_geq(refcoord_.ref(), 0);
	int off = ref.getStretch(
		reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
		(size_t)refcoord_.ref(),
		(size_t)max<TRefOff>(refcoord_.off(), 0),
		refallen,
		destU32);
	assert_leq(off, 16);
	char *refbuf = raw_refbuf.wbuf() + off;
	size_t trim5 = 0, trim3 = 0;
	if(trimSoft_) {
		trim5 += trim5p_;
		trim3 += trim3p_;
	}
	if(pretrimSoft_) {
		trim5 += pretrim5p_;
		trim3 += pretrim3p_;
	}
	rf.clear();
	rdseq.clear();
	rdseq = rd.patFw;
	if(!fw) {
		rdseq.reverseComp(false);
	}
	assert_eq(rdrows_, rdseq.length());
	// rdseq is the nucleotide sequence from upstream to downstream on the
	// Watson strand.  ned_ are the nucleotide edits from upstream to
	// downstream.  rf contains the reference characters.
	assert(Edit::repOk(ned_, rdseq, fw, trim5, trim3));
	Edit::toRef(rdseq, ned_, rf, fw, trim5, trim3);
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
		cerr << "    ";
		for(size_t i = 0; i < refallen-nsOnLeft; i++) {
			cerr << "ACGTN"[(int)refbuf[i]];
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
 * Initialized the stacked alignment with respect to a read string, a list of
 * edits (expressed left-to-right), and integers indicating how much hard and
 * soft trimming has occurred on either end of the read.
 *
 * s: read sequence
 * ed: all relevant edits, including ambiguous nucleotides
 * trimLS: # bases soft-trimmed from LHS
 * trimLH: # bases hard-trimmed from LHS
 * trimRS: # bases soft-trimmed from RHS
 * trimRH: # bases hard-trimmed from RHS
 */
void StackedAln::init(
	const BTDnaString& s,
	const EList<Edit>& ed,
	size_t trimLS,
	size_t trimLH,
	size_t trimRS,
	size_t trimRH)
{
	trimLS_ = trimLS;
	trimLH_ = trimLH;
	trimRS_ = trimRS;
	trimRH_ = trimRH;
	ASSERT_ONLY(size_t ln_postsoft = s.length() - trimLS - trimRS);
	stackRef_.clear();
	stackRel_.clear();
	stackRead_.clear();
	size_t rdoff = trimLS;
	for(size_t i = 0; i < ed.size(); i++) {
		assert_lt(ed[i].pos, ln_postsoft);
		size_t pos = ed[i].pos + trimLS;
		while(rdoff < pos) {
			int c = s[rdoff++];
			assert_range(0, 4, c);
			stackRef_.push_back("ACGTN"[c]);
			stackRel_.push_back('=');
			stackRead_.push_back("ACGTN"[c]);
		}
		if(ed[i].isMismatch()) {
			int c = s[rdoff++];
			assert_range(0, 4, c);
			assert_eq(c, asc2dna[(int)ed[i].qchr]);
			assert_neq(c, asc2dna[(int)ed[i].chr]);
			stackRef_.push_back(ed[i].chr);
			stackRel_.push_back('X');
			stackRead_.push_back("ACGTN"[c]);
		} else if(ed[i].isRefGap()) {
			int c = s[rdoff++];
			assert_range(0, 4, c);
			assert_eq(c, asc2dna[(int)ed[i].qchr]);
			stackRef_.push_back('-');
			stackRel_.push_back('I');
			stackRead_.push_back("ACGTN"[c]);
		} else if(ed[i].isReadGap()) {
			stackRef_.push_back(ed[i].chr);
			stackRel_.push_back('D');
			stackRead_.push_back('-');
		}
	}
	while(rdoff < s.length() - trimRS) {
		int c = s[rdoff++];
		assert_range(0, 4, c);
		stackRef_.push_back("ACGTN"[c]);
		stackRel_.push_back('=');
		stackRead_.push_back("ACGTN"[c]);
	}
	inited_ = true;
}

/**
 * Left-align all the gaps.  If this changes the alignment and the CIGAR or
 * MD:Z strings have already been calculated, this renders them invalid.
 *
 * We left-align gaps with in the following way: for each gap, we check
 * whether the character opposite the rightmost gap character is the same
 * as the character opposite the character just to the left of the gap.  If
 * this is the case, we can slide the gap to the left and make the
 * rightmost position previously covered by the gap into a non-gap.
 *
 * This scheme allows us to push the gap past a mismatch.  BWA does seem to
 * allow this.  It's not clear that Bowtie 2 should, since moving the
 * mismatch could cause a mismatch with one base quality to be replaced
 * with a mismatch with a different base quality.
 */
void StackedAln::leftAlign(bool pastMms) {
	assert(inited_);
	bool changed = false;
	size_t ln = stackRef_.size();
	// Scan left-to-right
	for(size_t i = 0; i < ln; i++) {
		int rel = stackRel_[i];
		if(rel != '=' && rel != 'X') {
			// Neither a match nor a mismatch - must be a gap
			assert(rel == 'I' || rel == 'D');
			size_t glen = 1;
			// Scan further right to measure length of gap
			for(size_t j = i+1; j < ln; j++) {
				if(rel != (int)stackRel_[j]) break;
				glen++;
			}
			// We've identified a gap of type 'rel' (D = deletion or read
			// gap, I = insertion or ref gap) with length 'glen'.  Now we
			// can try to slide it to the left repeatedly.
			size_t l = i - 1;
			size_t r = l + glen;
			EList<char>& gp  = ((rel == 'I') ? stackRef_ : stackRead_);
			const EList<char>& ngp = ((rel == 'I') ? stackRead_ : stackRef_);
			while(l > 0 && ngp[l] == ngp[r]) {
				assert(stackRel_[l] == '=' || stackRel_[l] == 'X');
				assert(stackRel_[r] == 'D' || stackRel_[r] == 'I');
				if(!pastMms && stackRel_[l] == 'X') {
					break;
				}
				swap(gp[l], gp[r]);
				swap(stackRel_[l], stackRel_[r]);
				assert_neq('-', gp[r]);
				assert_eq('-', gp[l]);
				l--; r--;
				changed = true;
			}
			i += (glen-1);
		}
	}
	if(changed) {
		cigCalc_ = mdzCalc_ = false;
	}
}

/**
 * Build the CIGAR list, if it hasn't already built.  Returns true iff it
 * was built for the first time.
 */
bool StackedAln::buildCigar(bool xeq) {
	assert(inited_);
	if(cigCalc_) {
		return false; // already done
	}
	cigOp_.clear();
	cigRun_.clear();
	if(trimLS_ > 0) {
		cigOp_.push_back('S');
		cigRun_.push_back(trimLS_);
	}
	size_t ln = stackRef_.size();
	for(size_t i = 0; i < ln; i++) {
		char op = stackRel_[i];
		if(!xeq && (op == 'X' || op == '=')) {
			op = 'M';
		}
		size_t run = 1;
		for(; i + run < ln; run++) {
			char op2 = stackRel_[i + run];
			if(!xeq && (op2 == 'X' || op2 == '=')) {
				op2 = 'M';
			}
			if(op2 != op) {
				break;
			}
		}
		i += (run-1);
		cigOp_.push_back(op);
		cigRun_.push_back(run);
	}
	if(trimRS_ > 0) {
		cigOp_.push_back('S');
		cigRun_.push_back(trimRS_);
	}
	cigCalc_ = true;
	return true;
}

/**
 * Build the CIGAR list, if it hasn't already built.  Returns true iff it
 * was built for the first time.
 */
bool StackedAln::buildMdz() {
	assert(inited_);
	if(mdzCalc_) {
		return false; // already done
	}
	mdzOp_.clear();
	mdzChr_.clear();
	mdzRun_.clear();
	size_t ln = stackRef_.size();
	for(size_t i = 0; i < ln; i++) {
		char op = stackRel_[i];
		if(op == '=') {
			size_t run = 1;
			size_t ninserts = 0;
			// Skip over matches and insertions (ref gaps)
			for(; i+run < ln; run++) {
				if(stackRel_[i + run] == '=') {
					// do nothing
				} else if(stackRel_[i + run] == 'I') {
					ninserts++;
				} else {
					break;
				}
			}
			i += (run - 1);
			mdzOp_.push_back('='); // = X or G
			mdzChr_.push_back('-');
			mdzRun_.push_back(run - ninserts);
		} else if(op == 'X') {
			assert_neq(stackRef_[i], stackRead_[i]);
			mdzOp_.push_back('X'); // = X or G
			mdzChr_.push_back(stackRef_[i]);
			mdzRun_.push_back(1);
		} else if(op == 'D') {
			assert_neq('-', stackRef_[i]);
			mdzOp_.push_back('G'); // = X or G
			mdzChr_.push_back(stackRef_[i]);
			mdzRun_.push_back(1);
		}
	}
	mdzCalc_ = true;
	return true;
}

/**
 * Write a CIGAR representation of the alignment to the given string and/or
 * char buffer.
 */
void StackedAln::writeCigar(
	BTString* o,      // if non-NULL, string to append to
	char* occ) const  // if non-NULL, character string to append to
{
	const EList<char>& op = cigOp_;
	const EList<size_t>& run = cigRun_;
	assert_eq(op.size(), run.size());
	if(o != NULL || occ != NULL) {
		char buf[128];
		ASSERT_ONLY(bool printed = false);
		for(size_t i = 0; i < op.size(); i++) {
			size_t r = run[i];
			if(r > 0) {
				itoa10<size_t>(r, buf);
				ASSERT_ONLY(printed = true);
				if(o != NULL) {
					o->append(buf);
					o->append(op[i]);
				}
				if(occ != NULL) {
					COPY_BUF();
					*occ = op[i];
					occ++;
				}
			}
		}
		assert(printed);
		if(occ != NULL) {
			*occ = '\0';
		}
	}
}

/**
 * Write an MD:Z representation of the alignment to the given string and/or
 * char buffer.
 */
void StackedAln::writeMdz(BTString* o, char* occ) const {
	char buf[128];
	bool mm_last = false;
	bool rdgap_last = false;
	bool first_print = true;
	const EList<char>& op = mdzOp_;
	const EList<char>& ch = mdzChr_;
	const EList<size_t>& run = mdzRun_;
	for(size_t i = 0; i < op.size(); i++) {
		size_t r = run[i];
		if(r > 0) {
			if(op[i] == '=') {
				// Write run length
				itoa10<size_t>(r, buf);
				if(o != NULL)  { o->append(buf); }
				if(occ != NULL) { COPY_BUF(); }
				first_print = false;
				mm_last = false;
				rdgap_last = false;
			} else if(op[i] == 'X') {
				if(o != NULL) {
					if(rdgap_last || mm_last || first_print) {
						o->append('0');
					}
					o->append(ch[i]);
				}
				if(occ != NULL) {
					if(rdgap_last || mm_last || first_print) {
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
					if(mm_last || first_print) {
						o->append('0');
					}
					if(!rdgap_last) {
						o->append('^');
					}
					o->append(ch[i]);
				}
				if(occ != NULL) {
					if(mm_last || first_print) {
						*occ = '0'; occ++;
					}
					if(!rdgap_last) {
						*occ = '^'; occ++;
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
		if(o   != NULL) { o->append('0'); }
		if(occ != NULL) { *occ = '0'; occ++; }
	}
	if(occ != NULL) { *occ = '\0'; }
}

/**
 * Print the sequence for the read that aligned using A, C, G and
 * T.  This will simply print the read sequence (or its reverse
 * complement).
 */
void AlnRes::printSeq(
	const Read& rd,         // read
	const BTDnaString* dns, // already-decoded nucleotides
	BTString& o) const      // buffer to write to
{
	assert(!rd.patFw.empty());
	ASSERT_ONLY(size_t written = 0);
	// Print decoded nucleotides
	assert(dns != NULL);
	size_t len = dns->length();
	size_t st = 0;
	size_t en = len;
	for(size_t i = st; i < en; i++) {
		int c = dns->get(i);
		assert_range(0, 3, c);
		o.append("ACGT"[c]);
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

/**
 * Print the quality string for the read that aligned.  This will simply print
 * the read qualities (or their reverse).
 */
void AlnRes::printQuals(
	const Read& rd,         // read
	const BTString* dqs,    // already-decoded qualities
	BTString& o) const      // output stream to write to
{
	assert(dqs != NULL);
	size_t len = dqs->length();
	// Print decoded qualities from upstream to downstream Watson
	for(size_t i = 1; i < len-1; i++) {
		o.append(dqs->get(i));
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
 * Given all the paired and unpaired results involving mates #1 and #2,
 * calculate best and second-best scores for both mates.  These are
 * used for future MAPQ calculations.
 */
void AlnSetSumm::init(
	const Read* rd1,
	const Read* rd2,
	const EList<AlnRes>* rs1,
	const EList<AlnRes>* rs2,
	const EList<AlnRes>* rs1u,
	const EList<AlnRes>* rs2u,
	bool exhausted1,
	bool exhausted2,
	TRefId orefid,
	TRefOff orefoff)
{
	assert(rd1 != NULL || rd2 != NULL);
	assert((rs1 == NULL) == (rs2 == NULL));
	AlnScore best[2], secbest[2], bestPaired, secbestPaired;
	size_t szs[2];
	best[0].invalidate();    secbest[0].invalidate();
	best[1].invalidate();    secbest[1].invalidate();
	bestPaired.invalidate(); secbestPaired.invalidate();
	bool paired = (rs1 != NULL && rs2 != NULL);
	szs[0] = szs[1] = 0;
	if(paired) {
		// Paired alignments
		assert_eq(rs1->size(), rs2->size());
		szs[0] = szs[1] = rs1->size();
		assert_gt(szs[0], 0);
		for(size_t i = 0; i < rs1->size(); i++) {
			AlnScore sc = (*rs1)[i].score() + (*rs2)[i].score();
			if(sc > bestPaired) {
				secbestPaired = bestPaired;
				bestPaired = sc;
				assert(VALID_AL_SCORE(bestPaired));
			} else if(sc > secbestPaired) {
				secbestPaired = sc;
				assert(VALID_AL_SCORE(bestPaired));
				assert(VALID_AL_SCORE(secbestPaired));
			}
		}
	}
	for(int j = 0; j < 2; j++) {
		const EList<AlnRes>* rs = (j == 0 ? rs1u : rs2u);
		if(rs == NULL) {
			continue;
		}
		assert(rs != NULL);
		szs[j] = rs->size();
		//assert_gt(szs[j], 0);
		for(size_t i = 0; i < rs->size(); i++) {
			AlnScore sc = (*rs)[i].score();
			if(sc > best[j]) {
				secbest[j] = best[j];
				best[j] = sc;
				assert(VALID_AL_SCORE(best[j]));
			} else if(sc > secbest[j]) {
				secbest[j] = sc;
				assert(VALID_AL_SCORE(best[j]));
				assert(VALID_AL_SCORE(secbest[j]));
			}
		}
	}
	if(szs[0] > 0 || szs[1] > 0) {
		init(
			best[0],
			secbest[0],
			best[1],
			secbest[1],
			bestPaired,
			secbestPaired,
			(szs[0] == 0) ? 0 : (szs[0] - 1),
			(szs[1] == 0) ? 0 : (szs[1] - 1),
			paired,
			exhausted1,
			exhausted2,
			orefid,
			orefoff);
	} else {
		reset();
		orefid_ = orefid;
		orefoff_ = orefoff;
	}
}

/**
 * Print out string representation of YF:i flag for indicating whether and
 * why the mate was filtered.
 */
bool AlnFlags::printYF(BTString& o, bool first) const {
	const char *flag = "";
	if     (!lenfilt_) flag = "LN";
	else if(!nfilt_  ) flag = "NS";
	else if(!scfilt_ ) flag = "SC";
	else if(!qcfilt_ ) flag = "QC";
	if(flag > 0) {
		if(!first) o.append('\t');
		o.append("YF:Z:");
		o.append(flag);
		return false;
	}
	return true;
}


/**
 * Print out string representation of YM:i flag for indicating with the
 * mate per se aligned repetitively.
 */
void AlnFlags::printYM(BTString& o) const {
	o.append("YM:i:");
	o.append(maxed() ? '1' : '0');
}

/**
 * Print out string representation of YM:i flag for indicating with the
 * pair containing the mate aligned repetitively.
 */
void AlnFlags::printYP(BTString& o) const {
	o.append("YP:i:");
	o.append(maxedPair() ? '1' : '0');
}

/**
 * Print out string representation of these flags.
 */
void AlnFlags::printYT(BTString& o) const {
	o.append("YT:Z:");
	if(alignedConcordant()) {
		o.append("CP");
	} else if(alignedDiscordant()) {
		o.append("DP");
	} else if(alignedUnpairedMate()) {
		o.append("UP");
	} else if(alignedUnpaired()) {
		o.append("UU");
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
