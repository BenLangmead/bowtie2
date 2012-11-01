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

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "bt2_idx.h"

///////////////////////////////////////////////////////////////////////
//
// Functions for printing and sanity-checking Ebwts
//
///////////////////////////////////////////////////////////////////////

/**
 * Check that the ebwt array is internally consistent up to (and not
 * including) the given side index by re-counting the chars and
 * comparing against the embedded occ[] arrays.
 */
void Ebwt::sanityCheckUpToSide(int upToSide) const {
	assert(isInMemory());
	uint32_t occ[] = {0, 0, 0, 0};
	ASSERT_ONLY(uint32_t occ_save[] = {0, 0, 0, 0});
	uint32_t cur = 0; // byte pointer
	const EbwtParams& eh = this->_eh;
	bool fw = false;
	while(cur < (upToSide * eh._sideSz)) {
		assert_leq(cur + eh._sideSz, eh._ebwtTotLen);
		for(uint32_t i = 0; i < eh._sideBwtSz; i++) {
			uint8_t by = this->ebwt()[cur + (fw ? i : eh._sideBwtSz-i-1)];
			for(int j = 0; j < 4; j++) {
				// Unpack from lowest to highest bit pair
				int twoBit = unpack_2b_from_8b(by, fw ? j : 3-j);
				occ[twoBit]++;
			}
			assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) % 4);
		}
		assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) % eh._sideBwtLen);
		// Finished forward bucket; check saved [A], [C], [G] and [T]
		// against the uint32_ts encoded here
		ASSERT_ONLY(const uint32_t *u32ebwt = reinterpret_cast<const uint32_t*>(&ebwt()[cur + eh._sideBwtSz]));
		ASSERT_ONLY(uint32_t as = u32ebwt[0]);
		ASSERT_ONLY(uint32_t cs = u32ebwt[1]);
		ASSERT_ONLY(uint32_t gs = u32ebwt[2]);
		ASSERT_ONLY(uint32_t ts = u32ebwt[3]);
		assert(as == occ_save[0] || as == occ_save[0]-1);
		assert_eq(cs, occ_save[1]);
		assert_eq(gs, occ_save[2]);
		assert_eq(ts, occ_save[3]);
#ifndef NDEBUG
		occ_save[0] = occ[0];
		occ_save[1] = occ[1];
		occ_save[2] = occ[2];
		occ_save[3] = occ[3];
#endif
		cur += eh._sideSz;
	}
}

/**
 * Sanity-check various pieces of the Ebwt
 */
void Ebwt::sanityCheckAll(int reverse) const {
	const EbwtParams& eh = this->_eh;
	assert(isInMemory());
	// Check ftab
	for(uint32_t i = 1; i < eh._ftabLen; i++) {
		assert_geq(this->ftabHi(i), this->ftabLo(i-1));
		assert_geq(this->ftabLo(i), this->ftabHi(i-1));
		assert_leq(this->ftabHi(i), eh._bwtLen+1);
	}
	assert_eq(this->ftabHi(eh._ftabLen-1), eh._bwtLen);
	
	// Check offs
	int seenLen = (eh._bwtLen + 31) >> 5;
	uint32_t *seen;
	try {
		seen = new uint32_t[seenLen]; // bitvector marking seen offsets
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating seen[] at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	memset(seen, 0, 4 * seenLen);
	uint32_t offsLen = eh._offsLen;
	for(uint32_t i = 0; i < offsLen; i++) {
		assert_lt(this->offs()[i], eh._bwtLen);
		int w = this->offs()[i] >> 5;
		int r = this->offs()[i] & 31;
		assert_eq(0, (seen[w] >> r) & 1); // shouldn't have been seen before
		seen[w] |= (1 << r);
	}
	delete[] seen;
	
	// Check nPat
	assert_gt(this->_nPat, 0);
	
	// Check plen, flen
	for(uint32_t i = 0; i < this->_nPat; i++) {
		assert_geq(this->plen()[i], 0);
	}
	
	// Check rstarts
	if(this->rstarts() != NULL) {
		for(uint32_t i = 0; i < this->_nFrag-1; i++) {
			assert_gt(this->rstarts()[(i+1)*3], this->rstarts()[i*3]);
			if(reverse == REF_READ_REVERSE) {
				assert(this->rstarts()[(i*3)+1] >= this->rstarts()[((i+1)*3)+1]);
			} else {
				assert(this->rstarts()[(i*3)+1] <= this->rstarts()[((i+1)*3)+1]);
			}
		}
	}
	
	// Check ebwt
	sanityCheckUpToSide(eh._numSides);
	VMSG_NL("Ebwt::sanityCheck passed");
}

/**
 * Transform this Ebwt into the original string in linear time by using
 * the LF mapping to walk backwards starting at the row correpsonding
 * to the end of the string.  The result is written to s.  The Ebwt
 * must be in memory.
 */
void Ebwt::restore(SString<char>& s) const {
	assert(isInMemory());
	s.resize(this->_eh._len);
	uint32_t jumps = 0;
	uint32_t i = this->_eh._len; // should point to final SA elt (starting with '$')
	SideLocus l(i, this->_eh, this->ebwt());
	while(i != _zOff) {
		assert_lt(jumps, this->_eh._len);
		//if(_verbose) cout << "restore: i: " << i << endl;
		// Not a marked row; go back a char in the original string
		uint32_t newi = mapLF(l);
		assert_neq(newi, i);
		s[this->_eh._len - jumps - 1] = rowL(l);
		i = newi;
		l.initFromRow(i, this->_eh, this->ebwt());
		jumps++;
	}
	assert_eq(jumps, this->_eh._len);
}

/**
 * Check that this Ebwt, when restored via restore(), matches up with
 * the given array of reference sequences.  For sanity checking.
 */
void Ebwt::checkOrigs(
	const EList<SString<char> >& os,
	bool color,
	bool mirror) const
{
	SString<char> rest;
	restore(rest);
	uint32_t restOff = 0;
	size_t i = 0, j = 0;
	if(mirror) {
		// TODO: FIXME
		return;
	}
	while(i < os.size()) {
		size_t olen = os[i].length();
		int lastorig = -1;
		for(; j < olen; j++) {
			size_t joff = j;
			if(mirror) joff = olen - j - 1;
			if((int)os[i][joff] == 4) {
				// Skip over Ns
				lastorig = -1;
				if(!mirror) {
					while(j < olen && (int)os[i][j] == 4) j++;
				} else {
					while(j < olen && (int)os[i][olen-j-1] == 4) j++;
				}
				j--;
				continue;
			}
			if(lastorig == -1 && color) {
				lastorig = os[i][joff];
				continue;
			}
			if(color) {
				assert_neq(-1, lastorig);
				assert_eq(dinuc2color[(int)os[i][joff]][lastorig], rest[restOff]);
			} else {
				assert_eq(os[i][joff], rest[restOff]);
			}
			lastorig = (int)os[i][joff];
			restOff++;
		}
		if(j == os[i].length()) {
			// Moved to next sequence
			i++;
			j = 0;
		} else {
			// Just jumped over a gap
		}
	}
}
