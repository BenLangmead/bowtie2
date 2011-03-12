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

using namespace std;

/**
 * Get the decoded nucleotide sequence 
 */
void AlnRes::decodedNucsAndQuals(
	const Read& rd,        // read that led to alignment
	BTDnaString& ns,       // out: decoded nucleotides
	BTString& qs) const    // out: decoded qualities
{
	// Walk along the colors
	bool fw = refcoord_.fw();
	assert(color_);
	const size_t rdlen = rd.length();
	const size_t len = rdlen+1;
	ns.resize(len);
	qs.resize(len);
	ns = rd.patFw;
	qs = rd.qual;
	if(!fw) {
		// Reverse the read to make it upstream-to-downstream.  Recall
		// that it's in colorspace, so no need to complement.
		ns.reverse();
		qs.reverse();
		// Need to flip edits around to make them
		// upstream-to-downstream.
		Edit::invertPoss(const_cast<EList<Edit>& >(ced_), rdlen);
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
	for(size_t i = 0; i < rdlen; i++) {
		// If it was determined to have been miscalled, get the
		// decoded subject color
		if(cedidx < ced_.size() && ced_[cedidx].pos == i) {
			assert_neq("ACGTN"[c], ced_[cedidx].chr);
			assert_eq ("ACGTN"[c], ced_[cedidx].qchr);
			c = ced_[cedidx].chr;
			c = asc2dnaOrCol[c];
			q = -q;
			cedidx++;
		}
		// Determine next nucleotide by combining the previous
		// nucleotide and the current color
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
	qs.set(min(dq, 127) , rdlen);
	assert_eq(ndn, ns[rdlen]);
	assert_eq(cedidx, ced_.size());
	
	if(!fw) {
		// Need to re-flip edits to make them 5'-to-3' again.
		Edit::invertPoss(const_cast<EList<Edit>& >(ced_), rdlen);
	}
}

#ifndef NDEBUG

/**
 * Assuming this AlnRes is an alignment for 'rd', check that the
 * alignment and 'rd' are compatible with the corresponding
 * reference sequence.
 */
bool AlnRes::matchesRef(
	const Read& rd,
	const BitPairReference& ref)
{
	assert(!empty());
	assert(repOk());
	assert(refcoord_.valid());
	bool fw = refcoord_.fw();
	const size_t rdlen = rd.length();
	size_t rflen = rdlen + (rd.color ? 1 : 0);
	// Adjust reference string length according to edits
	const size_t nedlen = ned_.size();
	for(size_t i = 0; i < nedlen; i++) {
		if(ned_[i].isInsert()) {
			rflen++;
		} else if(ned_[i].isDelete()) {
			rflen--;
		}
	}
	char *raw_refbuf = new char[rflen + 16];
	int nsOnLeft = 0;
	if(refcoord_.off() < 0) {
		nsOnLeft = -refcoord_.off();
	}
	int off = ref.getStretch(
		reinterpret_cast<uint32_t*>(raw_refbuf),
		refcoord_.ref(),
		max<TRefOff>(refcoord_.off(), 0),
		rflen);
	char *refbuf = raw_refbuf + off;
	
	BTDnaString rf;
	BTDnaString rdseq;
	BTString qseq;
	if(rd.color) {
		// Decode the nucleotide sequence from the alignment
		decodedNucsAndQuals(rd, rdseq, qseq);
	} else {
		// Get the nucleotide sequence from the read
		rdseq = rd.patFw;
		if(!fw) rdseq.reverseComp(false);
	}
	if(!fw) {
		// Invert the nucleotide edits so that they go from upstream to
		// downstream on the Watson strand
		Edit::invertPoss(ned_, rdlen + (color() ? 1 : 0));
	}
	// rdseq is the nucleotide sequence (decoded in the case of a
	// colorspace read) from upstream to downstream on the Watson
	// strand.  ned_ are the nucleotide edits from upstream to
	// downstream.  rf contains the reference characters.
	Edit::toRef(rdseq, ned_, rf);
	if(!fw) {
		// Re-invert the nucleotide edits so that they go from 5' to 3'
		Edit::invertPoss(ned_, rdlen + (color() ? 1 : 0));
	}
	assert_eq(rflen, rf.length());
	EList<bool> matches;
	bool matchesOverall = true;
	matches.resize(rflen);
	matches.fill(true);
	for(size_t i = 0; i < rflen; i++) {
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
		for(size_t i = 0; i < rflen; i++) {
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
	delete[] raw_refbuf;
	return matchesOverall;
}

#endif /*ndef NDEBUG*/

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
	bool printColors,       // true -> print colors instead of decoded nucleotides for colorspace alignment
	bool exEnds,            // true -> exclude ends when printing decoded nucleotides
	OutFileBuf& o) const    // output stream to write to
{
	assert(!rd.patFw.empty());
	bool fw = refcoord_.fw();
	assert(!printColors || rd.color);
	if(!rd.color || printColors) {
		// Print nucleotides or colors
		size_t len = rd.patFw.length();
		for(size_t i = 0; i < len; i++) {
			int c;
			if(fw) {
				c = rd.patFw[i];
			} else {
				c = rd.patFw[len-i-1];
				if(c < 4) c = c ^ 3;
			}
			assert_range(0, 4, c);
			o.write("ACGTN"[c]);
		}
	} else {
		// Print decoded nucleotides
		assert(dns != NULL);
		size_t len = dns->length();
		assert_eq(rd.patFw.length(), len - 1);
		size_t st = 0;
		size_t en = len;
		if(exEnds) {
			st++; en--;
		}
		for(size_t i = st; i < en; i++) {
			int c = dns->get(i);
			assert_range(0, 3, c);
			o.write("ACGT"[c]);
		}
	}
}

/**
 * Print the quality string for the read that aligned.  This will
 * simply print the read qualities (or their reverse) unless this
 * is a colorspace read and printColors is false.  In that case,
 * we print the decoded qualities rather than the original ones.
 */
void AlnRes::printQuals(
	const Read& rd,         // read
	const BTString* dqs,    // already-decoded qualities
	bool printColors,       // true -> print colors instead of decoded nucleotides for colorspace alignment
	bool exEnds,            // true -> exclude ends when printing decoded nucleotides
	OutFileBuf& o) const    // output stream to write to
{
	bool fw = refcoord_.fw();
	size_t len = rd.qual.length();
	assert(!printColors || rd.color);
	if(!rd.color || printColors) {
		// Print original qualities from upstream to downstream Watson
		for(size_t i = 0; i < len; i++) {
			int c = (fw ? rd.qual[i] : rd.qual[len-i-1]);
			o.write(c);
		}
	} else {
		assert(dqs != NULL);
		assert_eq(len+1, dqs->length());
		// Print decoded qualities from upstream to downstream Watson
		if(!exEnds) {
			// Print upstream-most quality
			o.write(dqs->get(0));
		}
		for(size_t i = 0; i < len-1; i++) {
			o.write(dqs->get(i+1));
		}
		if(!exEnds) {
			// Print downstream-most quality
			o.write(dqs->get(len));
		}
	}
}
