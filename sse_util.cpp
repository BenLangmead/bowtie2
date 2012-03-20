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

#include "sse_util.h"
#include "aligner_swsse.h"

/**
 * Given a column of filled-in cells, save the checkpointed cells in cs_.
 */
void Checkpointer::commitCol(
	__m128i *pvH,
	__m128i *pvE,
	__m128i *pvF,
	size_t coli)
{
	assert(doCheckpoints());
	if(coli < locol_) {
		locol_ = coli;
	}
	if(coli > hicol_) {
		hicol_ = coli;
	}
	size_t rc = coli;
	size_t rc_mod = rc & lomask_; // where does 0th row lie in period?
	assert_lt(rc_mod, per_);
	size_t diff1 = 0, diff2 = 0; // offset to 1-diag, 2-diag
	if(rc_mod == per_ - 1) {
		diff1 = per_ - 1;
		// diff2 stays at 0
	} else {
		assert_leq(rc_mod, per_ - 2);
		diff1 = per_ - 2 - rc_mod;
		diff2 = diff1 + 1;
	}
	assert_lt(diff1, per_);
	assert_lt(diff2, per_);
	// Commit all the diag1 scores
	for(size_t i = 0; i < nrow_; i += per_) {
		for(size_t j = 0; j < 2; j++) {
			size_t row = (j == 0 ? (i + diff1) : (i + diff2));
			EList<int16_t>& diags  = (j == 0 ? diag1s_  : diag2s_);
			EList<_CpQuad>& qdiags = (j == 0 ? qdiag1s_ : qdiag2s_);
			if(row < nrow_) {
				size_t rc_cur = rc + row;
				// Calc diag index
				size_t idx = rc_cur >> perpow2_;
				size_t elt = row / iter_;
				size_t word = row % iter_;
				if(is8_) {
					if(ef_) {
						assert(pvE != NULL);
						assert(pvF != NULL);
						assert_lt(idx * nrow_, qdiags.size());
						assert_lt(elt, 16);
						int16_t h_sc = ((uint8_t*)(pvH + (word * ROWSTRIDE_2COL)))[elt];
						int16_t e_sc = ((uint8_t*)(pvE + (word * ROWSTRIDE_2COL)))[elt];
						int16_t f_sc = ((uint8_t*)(pvF + (word * ROWSTRIDE_2COL)))[elt];
						if(!local_) {
							h_sc -= 0xff;
							e_sc -= 0xff;
							f_sc -= 0xff;
						}
						assert_leq(h_sc, perf_);
						assert_leq(e_sc, perf_);
						assert_leq(f_sc, perf_);
						// Calc diag offset
						const size_t delt = idx * nrow_ + row;
						qdiags[delt].sc[0] = h_sc;
						qdiags[delt].sc[1] = e_sc;
						qdiags[delt].sc[2] = f_sc;
					} else {
						assert_lt(idx * nrow_, diags.size());
						assert_lt(elt, 16);
						int sc = ((uint8_t*)(pvH + (word * ROWSTRIDE_2COL)))[elt];
						if(!local_) {
							sc -= 0xff;
						}
						assert_leq(sc, perf_);
						// Calc diag offset
						const size_t delt = idx * nrow_ + row;
						diags[delt] = (int16_t)sc;
					}
				} else {
					if(ef_) {
						assert(pvE != NULL);
						assert(pvF != NULL);
						assert_lt(idx * nrow_, qdiags.size());
						assert_lt(elt, 8);
						int16_t h_sc = ((int16_t*)(pvH + (word * ROWSTRIDE_2COL)))[elt];
						int16_t e_sc = ((int16_t*)(pvE + (word * ROWSTRIDE_2COL)))[elt];
						int16_t f_sc = ((int16_t*)(pvF + (word * ROWSTRIDE_2COL)))[elt];
						if(!local_) {
							if(h_sc != std::numeric_limits<int16_t>::min()) {
								h_sc -= 0x7fff;
							}
							if(e_sc != std::numeric_limits<int16_t>::min()) {
								e_sc -= 0x7fff;
							}
							if(f_sc != std::numeric_limits<int16_t>::min()) {
								f_sc -= 0x7fff;
							}
						} else {
							h_sc += 0x8000;
							e_sc += 0x8000;
							f_sc += 0x8000;
							assert_geq(h_sc, 0);
							assert_geq(e_sc, 0);
							assert_geq(f_sc, 0);
						}
						assert_leq(h_sc, perf_);
						assert_leq(e_sc, perf_);
						assert_leq(f_sc, perf_);
						// Calc diag offset
						const size_t delt = idx * nrow_ + row;
						qdiags[delt].sc[0] = h_sc;
						qdiags[delt].sc[1] = e_sc;
						qdiags[delt].sc[2] = f_sc;
					} else {
						assert_lt(idx * nrow_, diags.size());
						assert_lt(elt, 8);
						int sc = ((int16_t*)(pvH + (word * ROWSTRIDE_2COL)))[elt];
						if(!local_) {
							if(sc == std::numeric_limits<int16_t>::min()) {
								sc = std::numeric_limits<int>::min();
							} else {
								sc -= 0x7fff;
							}
						} else {
							sc += 0x8000;
							assert_geq(sc, 0);
						}
						assert_leq(sc, perf_);
						// Calc diag offset
						const size_t delt = idx * nrow_ + row;
						diags[delt] = (int16_t)sc;
					}
				}
			}
		}
	}
}
