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

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "bt2_idx.h"

using namespace std;

#ifdef BOWTIE_64BIT_INDEX

const std::string gEbwt_ext("bt2l");

#else

const std::string gEbwt_ext("bt2");

#endif  // BOWTIE_64BIT_INDEX

string gLastIOErrMsg;

///////////////////////////////////////////////////////////////////////
//
// Functions for searching Ebwts
// (But most of them are defined in the header)
//
///////////////////////////////////////////////////////////////////////

/**
 * Take an offset into the joined text and translate it into the
 * reference of the index it falls on, the offset into the reference,
 * and the length of the reference.  Use a binary search through the
 * sorted list of reference fragment ranges t
 */
void Ebwt::joinedToTextOff(
	TIndexOffU qlen, 
	TIndexOffU off,
	TIndexOffU& tidx,
	TIndexOffU& textoff,
	TIndexOffU& tlen,
	bool rejectStraddle,
	bool& straddled) const
{
	assert(rstarts() != NULL); // must have loaded rstarts
	TIndexOffU top = 0;
	TIndexOffU bot = _nFrag; // 1 greater than largest addressable element
	TIndexOffU elt = OFF_MASK;
	// Begin binary search
	while(true) {
		ASSERT_ONLY(TIndexOffU oldelt = elt);
		elt = top + ((bot - top) >> 1);
		assert_neq(oldelt, elt); // must have made progress
		TIndexOffU lower = rstarts()[elt*3];
		TIndexOffU upper;
		if(elt == _nFrag-1) {
			upper = _eh._len;
		} else {
			upper = rstarts()[((elt+1)*3)];
		}
		assert_gt(upper, lower);
		TIndexOffU fraglen = upper - lower;
		if(lower <= off) {
			if(upper > off) { // not last element, but it's within
				// off is in this range; check if it falls off
				if(off + qlen > upper) {
					straddled = true;
					if(rejectStraddle) {
						// it falls off; signal no-go and return
						tidx = OFF_MASK;
						assert_lt(elt, _nFrag-1);
						return;
					}
				}
				// This is the correct text idx whether the index is
				// forward or reverse
				tidx = rstarts()[(elt*3)+1];
				assert_lt(tidx, this->_nPat);
				assert_leq(fraglen, this->plen()[tidx]);
				// it doesn't fall off; now calculate textoff.
				// Initially it's the number of characters that precede
				// the alignment in the fragment
				TIndexOffU fragoff = off - rstarts()[(elt*3)];
				if(!this->fw_) {
					fragoff = fraglen - fragoff - 1;
					fragoff -= (qlen-1);
				}
				// Add the alignment's offset into the fragment
				// ('fragoff') to the fragment's offset within the text
				textoff = fragoff + rstarts()[(elt*3)+2];
				assert_lt(textoff, this->plen()[tidx]);
				break; // done with binary search
			} else {
				// 'off' belongs somewhere in the region between elt
				// and bot
				top = elt;
			}
		} else {
			// 'off' belongs somewhere in the region between top and
			// elt
			bot = elt;
		}
		// continue with binary search
	}
	tlen = this->plen()[tidx];
}

/**
 * Walk 'steps' steps to the left and return the row arrived at.  If we
 * walk through the dollar sign, return max value.
 */
TIndexOffU Ebwt::walkLeft(TIndexOffU row, TIndexOffU steps) const {
	assert(offs() != NULL);
	assert_neq(OFF_MASK, row);
	SideLocus l;
	if(steps > 0) l.initFromRow(row, _eh, ebwt());
	while(steps > 0) {
		if(row == _zOff) return OFF_MASK;
		TIndexOffU newrow = this->mapLF(l ASSERT_ONLY(, false));
		assert_neq(OFF_MASK, newrow);
		assert_neq(newrow, row);
		row = newrow;
		steps--;
		if(steps > 0) l.initFromRow(row, _eh, ebwt());
	}
	return row;
}

/**
 * Resolve the reference offset of the BW element 'elt'.
 */
TIndexOffU Ebwt::getOffset(TIndexOffU row) const {
	assert(offs() != NULL);
	assert_neq(OFF_MASK, row);
	if(row == _zOff) return 0;
	if((row & _eh._offMask) == row) return this->offs()[row >> _eh._offRate];
	TIndexOffU jumps = 0;
	SideLocus l;
	l.initFromRow(row, _eh, ebwt());
	while(true) {
		TIndexOffU newrow = this->mapLF(l ASSERT_ONLY(, false));
		jumps++;
		assert_neq(OFF_MASK, newrow);
		assert_neq(newrow, row);
		row = newrow;
		if(row == _zOff) {
			return jumps;
		} else if((row & _eh._offMask) == row) {
			return jumps + this->offs()[row >> _eh._offRate];
		}
		l.initFromRow(row, _eh, ebwt());
	}
}

/**
 * Resolve the reference offset of the BW element 'elt' such that
 * the offset returned is at the right-hand side of the forward
 * reference substring involved in the hit.
 */
TIndexOffU Ebwt::getOffset(
	TIndexOffU elt,
	bool fw,
	TIndexOffU hitlen) const
{
	TIndexOffU off = getOffset(elt);
	assert_neq(OFF_MASK, off);
	if(!fw) {
		assert_lt(off, _eh._len);
		off = _eh._len - off - 1;
		assert_geq(off, hitlen-1);
		off -= (hitlen-1);
		assert_lt(off, _eh._len);
	}
	return off;
}

/**
 * Returns true iff the index contains the given string (exactly).  The given
 * string must contain only unambiguous characters.  TODO: support ambiguous
 * characters in 'str'.
 */
bool Ebwt::contains(
	const BTDnaString& str,
	TIndexOffU *otop,
	TIndexOffU *obot) const
{
	assert(isInMemory());
	SideLocus tloc, bloc;
	if(str.empty()) {
		if(otop != NULL && obot != NULL) *otop = *obot = 0;
		return true;
	}
	int c = str[str.length()-1];
	assert_range(0, 4, c);
	TIndexOffU top = 0, bot = 0;
	if(c < 4) {
		top = fchr()[c];
		bot = fchr()[c+1];
	} else {
		bool set = false;
		for(int i = 0; i < 4; i++) {
			if(fchr()[c] < fchr()[c+1]) {
				if(set) {
					return false;
				} else {
					set = true;
					top = fchr()[c];
					bot = fchr()[c+1];
				}
			}
		}
	}
	assert_geq(bot, top);
	tloc.initFromRow(top, eh(), ebwt());
	bloc.initFromRow(bot, eh(), ebwt());
	ASSERT_ONLY(TIndexOffU lastDiff = bot - top);
	for(TIndexOff i = (TIndexOff)str.length()-2; i >= 0; i--) {
		c = str[i];
		assert_range(0, 4, c);
		if(c <= 3) {
			top = mapLF(tloc, c);
			bot = mapLF(bloc, c);
		} else {
			TIndexOffU sz = bot - top;
			int c1 = mapLF1(top, tloc ASSERT_ONLY(, false));
			bot = mapLF(bloc, c1);
			assert_leq(bot - top, sz);
			if(bot - top < sz) {
				// Encountered an N and could not proceed through it because
				// there was more than one possible nucleotide we could replace
				// it with
				return false;
			}
		}
		assert_geq(bot, top);
		assert_leq(bot-top, lastDiff);
		ASSERT_ONLY(lastDiff = bot-top);
		if(i > 0) {
			tloc.initFromRow(top, eh(), ebwt());
			bloc.initFromRow(bot, eh(), ebwt());
		}
	}
	if(otop != NULL && obot != NULL) {
		*otop = top; *obot = bot;
	}
	return bot > top;
}

/**
 * Try to find the Bowtie index specified by the user.  First try the
 * exact path given by the user.  Then try the user-provided string
 * appended onto the path of the "indexes" subdirectory below this
 * executable, then try the provided string appended onto
 * "$BOWTIE2_INDEXES/".
 */
string adjustEbwtBase(const string& cmdline,
					  const string& ebwtFileBase,
					  bool verbose = false)
{
	string str = ebwtFileBase;
	ifstream in;
	if(verbose) cout << "Trying " << str.c_str() << endl;
	in.open((str + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		if(verbose) cout << "  didn't work" << endl;
		in.close();
		if(getenv("BOWTIE2_INDEXES") != NULL) {
			str = string(getenv("BOWTIE2_INDEXES")) + "/" + ebwtFileBase;
			if(verbose) cout << "Trying " << str.c_str() << endl;
			in.open((str + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
			if(!in.is_open()) {
				if(verbose) cout << "  didn't work" << endl;
				in.close();
			} else {
				if(verbose) cout << "  worked" << endl;
			}
		}
	}
	if(!in.is_open()) {
		cerr << "Could not locate a Bowtie index corresponding to basename \"" << ebwtFileBase.c_str() << "\"" << endl;
		throw 1;
	}
	return str;
}
