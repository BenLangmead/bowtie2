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
#include "edit.h"

using namespace std;

/**
 * Print a single edit to a std::ostream.  Format is
 * (pos):(ref chr)>(read chr).  Where 'pos' is an offset from the 5'
 * end of the read, and the ref and read chrs are expressed w/r/t the
 * Watson strand.
 */
ostream& operator<< (ostream& os, const Edit& e) {
	os << e.pos << ":" << (char)e.chr << ">" << (char)e.qchr;
	return os;
}

/**
 * Print a list of edits to a std::ostream, separated by commas.
 */
void Edit::print(ostream& os, const EList<Edit>& edits, char delim) {
	for(size_t i = 0; i < edits.size(); i++) {
		os << edits[i];
		if(i < edits.size()-1) os << delim;
	}
}

/**
 * Flip all the edits.pos fields so that they're with respect to
 * the other end of the read (of length 'sz').
 */
void Edit::invertPoss(
	EList<Edit>& edits,
	size_t sz,
	size_t ei,
	size_t en,
	bool sort)
{
	// Invert elements
	size_t ii = 0;
	for(size_t i = ei; i < ei + en/2; i++) {
		Edit tmp = edits[i];
		edits[i] = edits[ei + en - ii - 1];
		edits[ei + en - ii - 1] = tmp;
		ii++;
	}
	for(size_t i = ei; i < ei + en; i++) {
		assert(edits[i].pos < sz ||
			   (edits[i].isReadGap() && edits[i].pos == sz));
		// Adjust pos
		edits[i].pos =
			(uint32_t)(sz - edits[i].pos - (edits[i].isReadGap() ? 0 : 1));
		// Adjust pos2
		if(edits[i].isReadGap()) {
			int64_t pos2diff = (int64_t)(uint64_t)edits[i].pos2 - (int64_t)((uint64_t)std::numeric_limits<uint32_t>::max() >> 1);
			int64_t pos2new = (int64_t)(uint64_t)edits[i].pos2 - 2*pos2diff;
			assert(pos2diff == 0 || (uint32_t)pos2new != (std::numeric_limits<uint32_t>::max() >> 1));
			edits[i].pos2 = (uint32_t)pos2new;
		}
	}
	if(sort) {
		// Edits might not necessarily be in same order after inversion
		edits.sortPortion(ei, en);
#ifndef NDEBUG
		for(size_t i = ei + 1; i < ei + en; i++) {
			assert_geq(edits[i].pos, edits[i-1].pos);
		}
#endif
	}
}

/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlign(
	std::ostream& os,
	const BTDnaString& read,
	const EList<Edit>& edits)
{
	printQAlign(os, "", read, edits);
}

/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlignNoCheck(
	std::ostream& os,
	const BTDnaString& read,
	const EList<Edit>& edits)
{
	printQAlignNoCheck(os, "", read, edits);
}

/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlign(
	std::ostream& os,
	const char *prefix,
	const BTDnaString& read,
	const EList<Edit>& edits)
{
	size_t eidx = 0;
	os << prefix;
	// Print read
	for(size_t i = 0; i < read.length(); i++) {
		bool del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				os << '-';
			} else if(edits[eidx].isRefGap()) {
				del = true;
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				os << read.toChar(i);
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				os << (char)edits[eidx].qchr;
			}
			eidx++;
		}
		if(!del && !mm) os << read.toChar(i);
	}
	os << endl;
	os << prefix;
	eidx = 0;
	// Print match bars
	for(size_t i = 0; i < read.length(); i++) {
		bool del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				os << ' ';
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << ' ';
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				os << ' ';
			}
			eidx++;
		}
		if(!del && !mm) os << '|';
	}
	os << endl;
	os << prefix;
	eidx = 0;
	// Print reference
	for(size_t i = 0; i < read.length(); i++) {
		bool del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				os << (char)edits[eidx].chr;
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << '-';
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				os << (char)edits[eidx].chr;
			}
			eidx++;
		}
		if(!del && !mm) os << read.toChar(i);
	}
	os << endl;
}

/**
 * For now, we pretend that the alignment is in the forward orientation
 * and that the Edits are listed from left- to right-hand side.
 */
void Edit::printQAlignNoCheck(
	std::ostream& os,
	const char *prefix,
	const BTDnaString& read,
	const EList<Edit>& edits)
{
	size_t eidx = 0;
	os << prefix;
	// Print read
	for(size_t i = 0; i < read.length(); i++) {
		bool del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				os << '-';
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << read.toChar(i);
			} else {
				mm = true;
				os << (char)edits[eidx].qchr;
			}
			eidx++;
		}
		if(!del && !mm) os << read.toChar(i);
	}
	os << endl;
	os << prefix;
	eidx = 0;
	// Print match bars
	for(size_t i = 0; i < read.length(); i++) {
		bool del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				os << ' ';
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << ' ';
			} else {
				mm = true;
				os << ' ';
			}
			eidx++;
		}
		if(!del && !mm) os << '|';
	}
	os << endl;
	os << prefix;
	eidx = 0;
	// Print reference
	for(size_t i = 0; i < read.length(); i++) {
		bool del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				os << (char)edits[eidx].chr;
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << '-';
			} else {
				mm = true;
				os << (char)edits[eidx].chr;
			}
			eidx++;
		}
		if(!del && !mm) os << read.toChar(i);
	}
	os << endl;
}

/**
 * Sort the edits in the provided list.
 */
void Edit::sort(EList<Edit>& edits) {
	edits.sort(); // simple!
}

/**
 * Given a read string and some edits, generate and append the corresponding
 * reference string to 'ref'.  If read aligned to the Watson strand, the caller
 * should pass the original read sequence and original edits.  If a read
 * aligned to the Crick strand, the caller should pass the reverse complement
 * of the read and a version of the edits list that has had Edit:invertPoss
 * called on it to cause edits to be listed in 3'-to-5' order.
 */
void Edit::toRef(
	const BTDnaString& read,
	const EList<Edit>& edits,
	BTDnaString& ref,
	bool fw,
	size_t trim5,
	size_t trim3)
{
	// edits should be sorted
	size_t eidx = 0;
	// Print reference
	const size_t rdlen = read.length();
	size_t trimBeg = fw ? trim5 : trim3;
	size_t trimEnd = fw ? trim3 : trim5;
	assert(Edit::repOk(edits, read, fw, trim5, trim3));
	if(!fw) {
		invertPoss(const_cast<EList<Edit>&>(edits), read.length()-trimBeg-trimEnd, false);
	}
	for(size_t i = 0; i < rdlen; i++) {
		ASSERT_ONLY(int c = read[i]);
		assert_range(0, 4, c);
		bool del = false, mm = false;
		bool append = i >= trimBeg && rdlen - i -1 >= trimEnd;
		bool appendIns = i >= trimBeg && rdlen - i >= trimEnd;
		while(eidx < edits.size() && edits[eidx].pos+trimBeg == i) {
			if(edits[eidx].isReadGap()) {
				// Inserted characters come before the position's
				// character
				if(appendIns) {
					ref.appendChar((char)edits[eidx].chr);
				}
			} else if(edits[eidx].isRefGap()) {
				assert_eq("ACGTN"[c], edits[eidx].qchr);
				del = true;
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				assert(edits[eidx].qchr != edits[eidx].chr || edits[eidx].qchr == 'N');
				assert_eq("ACGTN"[c], edits[eidx].qchr);
				if(append) {
					ref.appendChar((char)edits[eidx].chr);
				}
			}
			eidx++;
		}
		if(!del && !mm) {
			if(append) {
				ref.append(read[i]);
			}
		}
	}
	if(trimEnd == 0) {
		while(eidx < edits.size()) {
			assert_gt(rdlen, edits[eidx].pos);
			if(edits[eidx].isReadGap()) {
				ref.appendChar((char)edits[eidx].chr);
			}
			eidx++;
		}
	}
	if(!fw) {
		invertPoss(const_cast<EList<Edit>&>(edits), read.length()-trimBeg-trimEnd, false);
	}
}

#ifndef NDEBUG
/**
 * Check that the edit is internally consistent.
 */
bool Edit::repOk() const {
		assert(inited());
	// Ref and read characters cannot be the same unless they're both Ns
	assert(qchr != chr || qchr == 'N');
	// Type must match characters
	assert(isRefGap() ||  chr != '-');
	assert(isReadGap() || qchr != '-');
	assert(!isMismatch() || (qchr != '-' && chr != '-'));
	return true;
}

/**
 * Given a list of edits and a DNA string representing the query
 * sequence, check that the edits are consistent with respect to the
 * query.
 */
bool Edit::repOk(
	const EList<Edit>& edits,
	const BTDnaString& s,
	bool fw,
	size_t trimBeg,
	size_t trimEnd)
{
	if(!fw) {
		invertPoss(const_cast<EList<Edit>&>(edits), s.length()-trimBeg-trimEnd, false);
		swap(trimBeg, trimEnd);
	}
	for(size_t i = 0; i < edits.size(); i++) {
		const Edit& e = edits[i];
		size_t pos = e.pos;
		if(i > 0) {
			assert_geq(pos, edits[i-1].pos);
		}
		bool del = false, mm = false;
		while(i < edits.size() && edits[i].pos == pos) {
			const Edit& ee = edits[i];
			assert_lt(ee.pos, s.length());
			if(ee.qchr != '-') {
				assert(ee.isRefGap() || ee.isMismatch());
				assert_eq((int)ee.qchr, s.toChar(ee.pos+trimBeg));
			}
			if(ee.isMismatch()) {
				assert(!mm);
				mm = true;
				assert(!del);
			} else if(ee.isReadGap()) {
				assert(!mm);
			} else if(ee.isRefGap()) {
				assert(!mm);
				assert(!del);
				del = true;
			}
			i++;
		}
	}
	if(!fw) {
		invertPoss(const_cast<EList<Edit>&>(edits), s.length()-trimBeg-trimEnd, false);
	}
	return true;
}
#endif

/**
 * Merge second argument into the first.  Assume both are sorted to
 * begin with.
 */
void Edit::merge(EList<Edit>& dst, const EList<Edit>& src) {
	size_t di = 0, si = 0;
	while(di < dst.size()) {
		if(src[si].pos < dst[di].pos) {
			dst.insert(src[si], di);
			si++; di++;
		} else if(src[si].pos == dst[di].pos) {
			// There can be two inserts at a given position, but we
			// can't merge them because there's no way to know their
			// order
			assert(src[si].isReadGap() != dst[di].isReadGap());
			if(src[si].isReadGap()) {
				dst.insert(src[si], di);
				si++; di++;
			} else if(dst[di].isReadGap()) {
				di++;
			}
		}
	}
	while(si < src.size()) dst.push_back(src[si++]);
}

/**
 * Clip off some of the low-numbered positions.
 */
void Edit::clipLo(EList<Edit>& ed, size_t len, size_t amt) {
	size_t nrm = 0;
	for(size_t i = 0; i < ed.size(); i++) {
		assert_lt(ed[i].pos, len);
		if(ed[i].pos < amt) {
			nrm++;
		} else {
			// Shift everyone else up
			ed[i].pos -= (uint32_t)amt;
		}
	}
	ed.erase(0, nrm);
}

/**
 * Clip off some of the high-numbered positions.
 */
void Edit::clipHi(EList<Edit>& ed, size_t len, size_t amt) {
	assert_leq(amt, len);
	size_t max = len - amt;
	size_t nrm = 0;
	for(size_t i = 0; i < ed.size(); i++) {
		size_t ii = ed.size() - i - 1;
		assert_lt(ed[ii].pos, len);
		if(ed[ii].pos > max) {
			nrm++;
		} else if(ed[ii].pos == max && !ed[ii].isReadGap()) {
			nrm++;
		} else {
			break;
		}
	}
	ed.resize(ed.size() - nrm);
}
