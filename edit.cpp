/*
 * edit.cpp
 *
 *  Created on: Jul 14, 2009
 *      Author: Ben Langmead
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
void Edit::invertPoss(EList<Edit>& edits, size_t sz) {
	// Invert elements
	edits.reverse();
	// Invert all the .pos's
	for(size_t i = 0; i < edits.size(); i++) {
		assert_lt(edits[i].pos, sz);
		edits[i].pos = sz - edits[i].pos - (edits[i].type == EDIT_TYPE_INS ? 0 : 1);
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
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << '-';
			} else if(edits[eidx].isDelete()) {
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
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << ' ';
			} else if(edits[eidx].isDelete()) {
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
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << (char)edits[eidx].chr;
			} else if(edits[eidx].isDelete()) {
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
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << '-';
			} else if(edits[eidx].isDelete()) {
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
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << ' ';
			} else if(edits[eidx].isDelete()) {
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
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				os << (char)edits[eidx].chr;
			} else if(edits[eidx].isDelete()) {
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
 * Given a read string and some edits, generate and append the
 * corresponding reference string to 'ref'.  If read aligned to
 * the Watson strand, the caller should pass the original read
 * sequence and original edits.  If a read aligned to the Crick strand,
 * the caller should pass the reverse complement of the read and a
 * version of the edits list that has had Edit:invertPoss called on it
 * to cause edits to be listed in 3'-to-5' order.
 */
void Edit::toRef(
	const BTDnaString& read,
	const EList<Edit>& edits,
	BTDnaString& ref)
{
	// edits should be sorted
	size_t eidx = 0;
	// Print reference
	const size_t rdlen = read.length();
	for(size_t i = 0; i < rdlen; i++) {
		ASSERT_ONLY(int c = read[i]);
		assert_range(0, 4, c);
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isInsert()) {
				ins = true;
				// Inserted characters come before the position's
				// character
				ref.appendChar((char)edits[eidx].chr);
			} else if(edits[eidx].isDelete()) {
				assert_eq("ACGTN"[c], edits[eidx].qchr);
				del = true;
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				assert(edits[eidx].qchr != edits[eidx].chr || edits[eidx].qchr == 'N');
				assert_eq("ACGTN"[c], edits[eidx].qchr);
				ref.appendChar((char)edits[eidx].chr);
			}
			eidx++;
		}
		if(!del && !mm) {
			ref.append(read[i]);
		}
	}
}

/**
 * Check that the edit is internally consistent.
 */
bool Edit::repOk() const {
	// Ref and read characters cannot be the same unless they're both Ns
	assert(qchr != chr || qchr == 'N');
	// Type must match characters
	assert(isDelete() ||  chr != '-');
	assert(isInsert() || qchr != '-');
	assert(!isMismatch() || (qchr != '-' && chr != '-'));
	return true;
}

/**
 * Given a list of edits and a DNA string representing the query
 * sequence, check that the edits are consistent with respect to the
 * query.
 */
bool Edit::repOk(const EList<Edit>& edits,
                 const BTDnaString& s)
{
	for(size_t i = 0; i < edits.size(); i++) {
		const Edit& e = edits[i];
		size_t pos = e.pos;
		if(i > 0) {
			assert_geq(pos, edits[i-1].pos);
		}
		bool del = false, ins = false, mm = false;
		while(i < edits.size() && edits[i].pos == pos) {
			const Edit& ee = edits[i];
			assert_lt(ee.pos, s.length());
			if(ee.qchr != '-') {
				assert(ee.isDelete() || ee.isMismatch());
				assert_eq((int)ee.qchr, s.toChar(ee.pos));
			}
			if(ee.isMismatch()) {
				assert(!mm);
				mm = true;
				assert(!del);
			} else if(ee.isInsert()) {
				ins = true;
				assert(!mm);
			} else if(ee.isDelete()) {
				assert(!mm);
				assert(!del);
				del = true;
			}
			i++;
		}
	}
	return true;
}

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
			assert(src[si].isInsert() != dst[di].isInsert());
			if(src[si].isInsert()) {
				dst.insert(src[si], di);
				si++; di++;
			} else if(dst[di].isInsert()) {
				di++;
			}
		}
	}
	while(si < src.size()) dst.push_back(src[si++]);
}
