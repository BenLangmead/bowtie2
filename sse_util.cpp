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
#include "limit.h"

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
	int64_t row = -rc_mod-1;
	int64_t row_mod = row;
	int64_t row_div = 0;

	// Initialize idx, elt and word
	size_t idx = rc >> perpow2_;
	size_t idxrow = idx * nrow_;
	assert_eq(4, ROWSTRIDE_2COL);

	// Commit all the diag1 scores
	bool done = false;
	while(true) {
		row += (per_ - 2);
		row_mod += (per_ - 2);
		for(size_t j = 0; j < 2; j++) {
			row++;
			row_mod++;
			if(row >= 0 && (size_t)row < nrow_) {
				// Update row divided by iter_ and mod iter_
				while(row_mod >= (int64_t)iter_) {
					row_mod -= (int64_t)iter_;
					row_div++;
				}
				size_t delt = idxrow + row;
				if(is8_) {
					size_t vecoff = (row_mod << 6) + row_div;
					if(ef_) {
						assert_lt(row_div, 16);
						int16_t h_sc = ((uint8_t*)pvH)[vecoff];
						int16_t e_sc = ((uint8_t*)pvE)[vecoff];
						int16_t f_sc = ((uint8_t*)pvF)[vecoff];
						if(!local_) {
							if(h_sc == 0) h_sc = MIN_I16;
							else h_sc -= 0xff;
							if(e_sc == 0) e_sc = MIN_I16;
							else e_sc -= 0xff;
							if(f_sc == 0) f_sc = MIN_I16;
							else f_sc -= 0xff;
						}
						assert_leq(h_sc, perf_);
						assert_leq(e_sc, perf_);
						assert_leq(f_sc, perf_);
						_CpQuad *qdiags = ((j == 0) ? qdiag1s_.ptr() : qdiag2s_.ptr());
						qdiags[delt].sc[0] = h_sc;
						qdiags[delt].sc[1] = e_sc;
						qdiags[delt].sc[2] = f_sc;
					} else {
						assert_lt(row_div, 16);
						int sc = ((uint16_t*)pvH)[vecoff];
						if(!local_) sc -= 0xff;
						assert_leq(sc, perf_);
						EList<int16_t>& diags = (j == 0 ? diag1s_  : diag2s_);
						diags[delt] = (int16_t)sc;
					}
				} else {
					size_t vecoff = (row_mod << 5) + row_div;
					if(ef_) {
						assert_lt(row_div, 8);
						int16_t h_sc = ((int16_t*)pvH)[vecoff];
						int16_t e_sc = ((int16_t*)pvE)[vecoff];
						int16_t f_sc = ((int16_t*)pvF)[vecoff];
						if(!local_) {
							if(h_sc != MIN_I16) h_sc -= 0x7fff;
							if(e_sc != MIN_I16) e_sc -= 0x7fff;
							if(f_sc != MIN_I16) f_sc -= 0x7fff;
						} else {
							h_sc += 0x8000; assert_geq(h_sc, 0);
							e_sc += 0x8000; assert_geq(e_sc, 0);
							f_sc += 0x8000; assert_geq(f_sc, 0);
						}
						assert_leq(h_sc, perf_);
						assert_leq(e_sc, perf_);
						assert_leq(f_sc, perf_);
						_CpQuad *qdiags = ((j == 0) ? qdiag1s_.ptr() : qdiag2s_.ptr());
						qdiags[delt].sc[0] = h_sc;
						qdiags[delt].sc[1] = e_sc;
						qdiags[delt].sc[2] = f_sc;
					} else {
						assert_lt(row_div, 16);
						int sc = ((uint16_t*)pvH)[vecoff];
						if(!local_) {
							if(sc == MIN_I16) {
								sc = MIN_I;
							} else sc -= 0x7fff;
						} else {
							sc += 0x8000; assert_geq(sc, 0);
						}
						assert_leq(sc, perf_);
						EList<int16_t>& diags = (j == 0 ? diag1s_  : diag2s_);
						diags[delt] = (int16_t)sc;
					}
				}
			} // if(row >= 0 && row < nrow_)
			else if(row >= 0 && (size_t)row >= nrow_) {
				done = true;
				break;
			}
		} // end of loop over anti-diags
		if(done) {
			break;
		}
		idx++;
		idxrow += nrow_;
	} // end of loop over strips
}
