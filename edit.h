/*
 * edit.h
 *
 *  Created on: Jul 31, 2009
 *      Author: Ben Langmead
 */

#ifndef EDIT_H_
#define EDIT_H_

#include <iostream>
#include <stdint.h>
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
 * Encapsulates an edit between the read sequence and the reference
 * sequence.
 */
struct Edit {

	Edit() : pos(0xffffffff) { }

	Edit(
		uint32_t po,
		int ch,
		int qc,
		int ty,
		bool chrs = true)
	{
		init(po, ch, qc, ty, chrs);
	}
	
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
		if(!chrs) {
			assert_range(0, 4, (int)chr);
			assert_range(0, 4, (int)qchr);
			chr = "ACGTN"[chr];
			qchr = "ACGTN"[qchr];
		}
		assert_in(chr, "ACMGRSVTWYHKDBN-");
		assert_in(qchr, "ACGTN-");
		assert(chr != qchr || chr == 'N');
	}

	/**
	 * Edit less-than overload.
	 */
	int operator< (const Edit &rhs) const {
		if(pos  < rhs.pos) return 1;
		if(pos  > rhs.pos) return 0;
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
		return(pos  == rhs.pos &&
			   chr  == rhs.chr &&
			   qchr == rhs.qchr &&
			   type == rhs.type);
	}

	/**
	 * Return true iff this Edit is initialized.
	 */
	bool initialized() const {
		return pos != 0xffffffff;
	}

	/**
	 * Return true iff this Edit is an initialized insertion.
	 */
	bool isReadGap() const {
		return initialized() && type == EDIT_TYPE_READ_GAP;
	}

	/**
	 * Return true iff this Edit is an initialized deletion.
	 */
	bool isRefGap() const {
		return initialized() && type == EDIT_TYPE_REF_GAP;
	}

	/**
	 * Return true if this Edit is either an initialized deletion or an
	 * initialized insertion.
	 */
	bool isGap() const {
		return initialized() && (type == EDIT_TYPE_REF_GAP || type == EDIT_TYPE_READ_GAP);
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
		return initialized() && type == EDIT_TYPE_MM;
	}

	/**
	 * Sort the edits in the provided list.
	 */
	static void sort(EList<Edit>& edits);

	/**
	 * Flip all the edits.pos fields so that they're with respect to
	 * the other end of the read (of length 'sz').
	 */
	static void invertPoss(EList<Edit>& edits, size_t sz);
	
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
		size_t trimBeg = 0,
		size_t trimEnd = 0);

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

	static bool repOk(
		const EList<Edit>& edits,
		const BTDnaString& s);

	uint8_t  chr;  // reference character involved (for subst and ins)
	uint8_t  qchr; // read character involved (for subst and del)
	uint8_t  type; // 1 -> mm, 2 -> SNP, 3 -> ins, 4 -> del
	uint32_t pos;  // position w/r/t search root

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
