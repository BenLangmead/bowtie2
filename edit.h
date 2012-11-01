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

#ifndef EDIT_H_
#define EDIT_H_

#include <iostream>
#include <stdint.h>
#include <limits>
#include "assert_helpers.h"
#include "filebuf.h"
#include "sstring.h"
#include "ds.h"

/**
 * 3 types of edits; mismatch (substitution), insertion in the
 * reference, deletion in the reference.
 */
enum {
	EDIT_TYPE_READ_GAP = 1,
	EDIT_TYPE_REF_GAP,
	EDIT_TYPE_MM,
	EDIT_TYPE_SNP
};

/**
 * Encapsulates an edit between the read sequence and the reference sequence.
 * We obey a few conventions when populating its fields.  The fields are:
 *
 * 	uint8_t  chr;  // reference character involved (for subst and ins)
 *  uint8_t  qchr; // read character involved (for subst and del)
 *  uint8_t  type; // 1 -> mm, 2 -> SNP, 3 -> ins, 4 -> del
 *  uint32_t pos;  // position w/r/t search root
 *
 * One convention is that pos is always an offset w/r/t the 5' end of the read.
 *
 * Another is that chr and qchr are expressed in terms of the nucleotides on
 * the forward version of the read.  So if we're aligning the reverse
 * complement of the read, and an A in the reverse complement mismatches a C in
 * the reference, chr should be G and qchr should be T.
 */
struct Edit {

	Edit() { reset(); }

	Edit(
		uint32_t po,
		int ch,
		int qc,
		int ty,
		bool chrs = true)
	{
		init(po, ch, qc, ty, chrs);
	}
	
    /**
     * Reset Edit to uninitialized state.
     */
	void reset() {
		pos = pos2 = std::numeric_limits<uint32_t>::max();
		chr = qchr = type = 0;
	}
	
    /**
     * Return true iff the Edit is initialized.
     */
	bool inited() const {
		return pos != std::numeric_limits<uint32_t>::max();
	}
	
    /**
     * Initialize a new Edit.
     */
	void init(
		uint32_t po,
		int ch,
		int qc,
		int ty,
		bool chrs = true)
	{
		chr = ch;
		qchr = qc;
		type = ty;
		pos = po;
		pos2 = std::numeric_limits<uint32_t>::max() >> 1;
		if(!chrs) {
			assert_range(0, 4, (int)chr);
			assert_range(0, 4, (int)qchr);
			chr = "ACGTN"[chr];
			qchr = "ACGTN"[qchr];
		}
		assert_in(chr, "ACMGRSVTWYHKDBN-");
		assert_in(qchr, "ACGTN-");
		assert(chr != qchr || chr == 'N');
		assert(inited());
	}
	
	/**
	 * Return true iff one part of the edit or the other has an 'N'.
	 */
	bool hasN() const {
		assert(inited());
		return chr == 'N' || qchr == 'N';
	}

	/**
	 * Edit less-than overload.
	 */
	int operator< (const Edit &rhs) const {
		assert(inited());
		if(pos  < rhs.pos) return 1;
		if(pos  > rhs.pos) return 0;
		if(pos2 < rhs.pos2) return 1;
		if(pos2 > rhs.pos2) return 0;
		if(type < rhs.type) return 1;
		if(type > rhs.type) return 0;
		if(chr  < rhs.chr) return 1;
		if(chr  > rhs.chr) return 0;
		return (qchr < rhs.qchr)? 1 : 0;
	}

	/**
	 * Edit equals overload.
	 */
	int operator== (const Edit &rhs) const {
		assert(inited());
		return(pos  == rhs.pos &&
			   pos2 == rhs.pos2 &&
			   chr  == rhs.chr &&
			   qchr == rhs.qchr &&
			   type == rhs.type);
	}

	/**
	 * Return true iff this Edit is an initialized insertion.
	 */
	bool isReadGap() const {
		assert(inited());
		return type == EDIT_TYPE_READ_GAP;
	}

	/**
	 * Return true iff this Edit is an initialized deletion.
	 */
	bool isRefGap() const {
		assert(inited());
		return type == EDIT_TYPE_REF_GAP;
	}

	/**
	 * Return true if this Edit is either an initialized deletion or an
	 * initialized insertion.
	 */
	bool isGap() const {
		assert(inited());
		return (type == EDIT_TYPE_REF_GAP || type == EDIT_TYPE_READ_GAP);
	}
	
	/**
	 * Return the number of gaps in the given edit list.
	 */
	static size_t numGaps(const EList<Edit>& es) {
		size_t gaps = 0;
		for(size_t i = 0; i < es.size(); i++) {
			if(es[i].isGap()) gaps++;
		}
		return gaps;
	}

	/**
	 * Return true iff this Edit is an initialized mismatch.
	 */
	bool isMismatch() const {
		assert(inited());
		return type == EDIT_TYPE_MM;
	}

	/**
	 * Sort the edits in the provided list.
	 */
	static void sort(EList<Edit>& edits);

	/**
	 * Flip all the edits.pos fields so that they're with respect to
	 * the other end of the read (of length 'sz').
	 */
	static void invertPoss(EList<Edit>& edits, size_t sz, size_t ei, size_t en);

	/**
	 * Flip all the edits.pos fields so that they're with respect to
	 * the other end of the read (of length 'sz').
	 */
	static void invertPoss(EList<Edit>& edits, size_t sz) {
		invertPoss(edits, sz, 0, edits.size());
	}
	
	/**
	 * Clip off some of the low-numbered positions.
	 */
	static void clipLo(EList<Edit>& edits, size_t len, size_t amt);

	/**
	 * Clip off some of the high-numbered positions.
	 */
	static void clipHi(EList<Edit>& edits, size_t len, size_t amt);

	/**
	 * Given a read string and some edits, generate and append the
	 * corresponding reference string to 'ref'.
	 */
	static void toRef(
		const BTDnaString& read,
		const EList<Edit>& edits,
		BTDnaString& ref,
		bool fw = true,
		size_t trim5 = 0,
		size_t trim3 = 0);

	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.
	 */
	static void printQAlign(
		std::ostream& os,
		const BTDnaString& read,
		const EList<Edit>& edits);

	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.  Add 'prefix'
	 * before each line of output.
	 */
	static void printQAlign(
		std::ostream& os,
		const char *prefix,
		const BTDnaString& read,
		const EList<Edit>& edits);

	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.
	 */
	static void printQAlignNoCheck(
		std::ostream& os,
		const BTDnaString& read,
		const EList<Edit>& edits);

	/**
	 * Given a string and its edits with respect to some other string,
	 * print the alignment between the strings with the strings stacked
	 * vertically, with vertical bars denoting matches.  Add 'prefix'
	 * before each line of output.
	 */
	static void printQAlignNoCheck(
		std::ostream& os,
		const char *prefix,
		const BTDnaString& read,
		const EList<Edit>& edits);

	bool repOk() const;

	/**
	 * Given a list of edits and a DNA string representing the query
	 * sequence, check that the edits are consistent with respect to the
	 * query.
	 */
	static bool repOk(
		const EList<Edit>& edits,
		const BTDnaString& s,
		bool fw = true,
		size_t trim5 = 0,
		size_t trim3 = 0);

	uint8_t  chr;  // reference character involved (for subst and ins)
	uint8_t  qchr; // read character involved (for subst and del)
	uint8_t  type; // 1 -> mm, 2 -> SNP, 3 -> ins, 4 -> del
	uint32_t pos;  // position w/r/t search root
	uint32_t pos2; // Second int to take into account when sorting.  Useful for
	               // sorting read gap edits that are all part of the same long
				   // gap.

	friend std::ostream& operator<< (std::ostream& os, const Edit& e);

	/**
	 * Print a comma-separated list of Edits to given output stream.
	 */
	static void print(
		std::ostream& os,
		const EList<Edit>& edits,
		char delim = '\t');

	/**
	 * Merge second argument into the first.  Assume both are sorted to
	 * begin with.
	 */
	static void merge(EList<Edit>& dst, const EList<Edit>& src);
};

#endif /* EDIT_H_ */
