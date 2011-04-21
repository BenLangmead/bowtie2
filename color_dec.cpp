/*
 * color_dec.cpp
 *
 *  Created on: October 15, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <string>
#include <limits>
#include <stdlib.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "color_dec.h"
#include "color.h"
#include "qual.h"
#include "mask.h"

using namespace std;

const TScore maxPenalty = 256;

static inline TScore penaltyToScore(TScore pen) {
	return maxPenalty - pen;
}

#define QUAL(i) mmPenalty(false, phredcToPhredq(qual[i]))

/**
 * Given the dynamic programming table, trace backwards from the last
 * column and populate the 's' and 'cmm' strings accordingly.  Whenever
 * there are multiple equally good ways of backtracking, choose one at
 * random.
 */
void ColorspaceDecoderUngapped::backtrack(
	const BTDnaString& read, // read
	size_t readi, // offset of first read char
	size_t /*readf*/, // offset of last read char
	BTString& ref, // reference
	size_t refi, // offset of first reference char
	size_t reff, // offset of last reference char
	BTDnaString& decoded, // dest for decoded bases
	EList<Edit>& nedits, // base edits
	EList<Edit>& aedits, // base edits
	EList<Edit>& cedits, // color miscalls
	EList<Edit>& ccedits, // all color edits (due both to base edits and miscalls)
	RandomSource& rand) // color edits
{
	assert(nedits.empty());
	assert(aedits.empty());
	assert(cedits.empty());
	assert(ccedits.empty());
	const size_t len = reff-refi;
	TScore min = MAX_SCORE;
	int bests = 0;
	// Determine best base in final column of table_
	for(int i = 0; i < 4; i++) {
		// Install minimum and backtrack info
		int m = table_[len-1].best[i];
		if(m < min) {
			min = m;
			bests = (1 << i);
		} else if(m == min) {
			bests |= (1 << i);
		}
	}
	decoded.resize(len);
	// i <- position of rightmost nucleotide
	int lenm = (int)len-1;
	// to <- rightmost nucleotide
	int to = randFromMask(rand, bests);
	assert_geq(to, 0);
	assert_lt(to, 4);
	while(true) {
		bests = table_[lenm].mask[to]; // get next best mask
		decoded.set(to, lenm--); // install best nucleotide
		if(lenm < 0) break; // done
		assert_gt(bests, 0);
		assert_lt(bests, 16);
		to = randFromMask(rand, bests); // select
	}
	// Determine what reference nucleotides were matched against
	for(size_t i = 0; i < len; i++) {
		assert_lt(i, reff);
		assert_geq(decoded[i], 0);
		assert_lt(decoded[i], 4);
		assert_gt(mask2popcnt[(int)ref[refi+i]], 0);
		assert_lt((int)ref[refi+i], 16);
		char refc = mask2dna[(int)ref[refi+i]];
		char qc = "ACGTN"[(int)decoded[i]];
		if(!matches(decoded[i], ref[refi+i])) {
			nedits.push_back(Edit((int)i, refc, qc, EDIT_TYPE_MM));
			assert(nedits.back().repOk());
		} else if(mask2popcnt[(int)ref[refi+i]] > 1) {
			aedits.push_back(Edit((int)i, refc, qc, EDIT_TYPE_MM));
			assert(aedits.back().repOk());
		}
	}
	for(size_t i = 0; i < len-1; i++) {
		int c1 = (int)read[readi+i]; // actual
		int c2ref = dnamasks2colormask[(int)ref[refi+i]][(int)ref[refi+i+1]];
		int c2sub = dinuc2color[(int)decoded[i]][(int)decoded[i+1]]; // decoded
		assert_leq(c1, 4); assert_geq(c1, 0);
		if(c1 != c2sub || c1 == 4) {
			// Read != decoding-induced
			char subc = "ACGT"[c2sub];
			char qc = "ACGTN"[c1];
			cedits.push_back(Edit((int)i, subc, qc, EDIT_TYPE_MM));
		}
		// The last condition in this if statement is important for the
		// case where a color matches w/r/t one of the ambiguous colors
		// that arises from the reference base being ambiguous, but
		// isn't compatible with the one we chose in our decoding.
		if(c1 == 4 || (c2ref & (1 << c1)) == 0) {
			// Read != reference-induced
			char refc = mask2dna[c2ref];
			char qc = "ACGTN"[c1];
			ccedits.push_back(Edit((int)i, refc, qc, EDIT_TYPE_MM));
			assert(ccedits.back().repOk());
		}
	}
	// done
}

/**
 * Decode the colorspace read 'read' as aligned against the reference
 * string 'ref', assuming that it's a hit without gaps (only matches
 * and mismatches).
 */
TScore ColorspaceDecoderUngapped::decode(
	const BTDnaString& read,
	const BTString& qual,
	size_t readi,
	size_t readf,
	BTString& ref,
	size_t refi,
	size_t reff,
	int snpPhred,
	BTDnaString& decoded,
	EList<Edit>& nedits,
	EList<Edit>& aedits,
	EList<Edit>& cedits,
	EList<Edit>& ccedits,
	RandomSource& rand)
{
	assert_lt(refi, reff);
	assert_lt(readi, readf);
	assert_eq(reff-refi-1, readf-readi);

	// The first column of the table just considers the first
	// nucleotide and whether it matches the ref nucleotide.
	table_.resize(1);
	for(int to = 0; to < 4; to++) {
		if(matches(to, ref[refi])) {
			// The assigned subject nucleotide matches the reference;
			// no penalty
			table_[0].best[to] = 0;  // best
			table_[0].mask[to] = 15; // mask
		} else {
			// The assigned subject nucleotide does not match the
			// reference nucleotide, so we add a SNP penalty
			table_[0].best[to] = snpPhred; // best
			table_[0].mask[to] = 15;     // mask
		}
		assert_gt(table_[0].mask[to], 0);
	}

	// Successive columns examine successive alignment positions
	TScore omin = MAX_SCORE;
	TScore lastOmin = MAX_SCORE;
	int t = 0;
	for(size_t c = readi; c < readf; c++) {
		const int readc = (int)read[c];
		assert_leq(readc, 4);
		assert_geq(readc, 0);
		lastOmin = omin;
		omin = MAX_SCORE;
		// t <- index of column in dynamic programming table
		t = (int)(c - readi + 1);
		int q = QUAL(c);
		assert_gt(t, 0);
		table_.resize(t+1);
		const int refc = ref[refi + t];
		int from[] = { table_[t-1].best[0], table_[t-1].best[1],
		               table_[t-1].best[2], table_[t-1].best[3] };
		// For each downstream nucleotide
		for(int to = 0; to < 4; to++) {
			// For each upstream nucleotide
			TScore min = MAX_SCORE;
			const int goodfrom = nuccol2nuc[to][readc];
			// Reward the preceding position
			if(goodfrom < 4) from[goodfrom] -= q;
			min = from[0];
			table_[t].mask[to] = 1;
			if(from[1] < min) {
				min = from[1];
				table_[t].mask[to] = 2;
			} else if(from[1] == min) {
				table_[t].mask[to] |= 2;
			}
			if(from[2] < min) {
				min = from[2];
				table_[t].mask[to] = 4;
			} else if(from[2] == min) {
				table_[t].mask[to] |= 4;
			}
			if(from[3] < min) {
				min = from[3];
				table_[t].mask[to] = 8;
			} else if(from[3] == min) {
				table_[t].mask[to] |= 8;
			}
			assert_gt(table_[t].mask[to], 0);
			min += q;
			if(!matches(to, refc)) {
				min += snpPhred;
			}
			table_[t].best[to] = min;
			if(min < omin) omin = min;
			if(goodfrom < 4) from[goodfrom] += q;
		}
	}

	t++;
	assert_eq(t, (int)(reff - refi));
	// Install the best backward path into ns, cmm, nmm
	backtrack(read, readi, readi + t - 1, ref, refi, refi + t,
	          decoded, nedits, aedits, cedits, ccedits, rand);
	return omin;
}

ostream& operator<<(ostream& os, const GappedScore& o) {
	os << o.gaps << ':' << o.score;
	return os;
}

/**
 * Given the dynamic programming table, trace backwards from the lower
 * right-hand corner and populate 'decoded' 'nedits', 'cedits' and
 * 'ccedits' accordingly.
 *
 * Initially, this routine was designed to avoid having to "backtrack"
 * or otherwise maintain per-path state besides what's already
 * calculated in the decode function.  This proved inadequate for some
 * gapped situtations because the decode function doesn't guarantee
 * that the steps marked as "best" adhere globally to the insAllow and
 * delAllow constraints.  The solution was to introduce a backtracking
 * facility to enforce the additional constraint that paths cannot
 * violoate insAllow and delAllow.
 */
GappedScore ColorspaceDecoderGapped::backtrack(
	const BTDnaString& read, // read
	const BTString& qual,
	size_t readi, // offset of first read char
	size_t readf, // offset of last read char
	BTString& ref, // reference
	size_t refi, // offset of first reference char
	size_t reff, // offset of last reference char
	int snpPhred,
	GappedScore decodeScore,
	int readGaps,
	int refGaps,      // # of reference gaps in
	int readOpenPen,  // penalty for opening a new gap in the read
	int readExtendPen,// penalty for extending a gap in the read
	int refOpenPen,   // penalty for opening a new gap in the reference
	int refExtendPen, // penalty for extending a gap in the reference
	int gapBarrier,   // # bases on either side of alignment that must be gap-free
	BTDnaString& decoded, // dest for decoded bases
	EList<Edit>& nedits, // base edits
	EList<Edit>& aedits, // ambiguous nucleotides resolved
	EList<Edit>& cedits, // color miscalls
	EList<Edit>& ccedits, // all color edits (due both to base edits and miscalls)
	bool constrainNumGaps,
	int lastC,
	RandomSource& rand)
{
	assert(cedits.empty());
	assert(nedits.empty());
	assert(aedits.empty());
	int row = (int)(readf - readi);
	assert_eq(row, (int)table_.size()-1);
	int col = (int)(reff - refi - 1);
	decoded.resize(readf-readi+1);
	int curC = lastC;
	bool refExtend = false;
	bool readExtend = false;
	GappedScore score;
	score.gaps = score.score = 0;
	return backtrackRecursive(
			read,
			qual,
			readi,
			readf,
			ref,
			refi,
			reff,
			snpPhred,
			decodeScore,
			readGaps,
			refGaps,
			readOpenPen,
			readExtendPen,
			refOpenPen,
			refExtendPen,
			gapBarrier,
			decoded,
			nedits,
			aedits,
			cedits,
			ccedits,
			constrainNumGaps,
			curC,
			rand,
			row,
			col,
			refExtend,
			readExtend,
			score,
			0);
}

/**
 * Given the dynamic programming table, trace backwards from the lower
 * right-hand corner and populate 'decoded' 'nedits', 'cedits' and
 * 'ccedits' accordingly.
 *
 * Tricky aspects to the recursion:
 *
 *  
 */
GappedScore ColorspaceDecoderGapped::backtrackRecursive(
	const BTDnaString& read, // read
	const BTString& qual,
	size_t readi, // offset of first read char
	size_t readf, // offset of last read char
	BTString& ref, // reference
	size_t refi, // offset of first reference char
	size_t reff, // offset of last reference char
	int snpPhred,
	GappedScore decodeScore,
	int readGaps,
	int refGaps,      // # of reference gaps in
	int readOpenPen,  // penalty for opening a new gap in the read
	int readExtendPen,// penalty for extending a gap in the read
	int refOpenPen,   // penalty for opening a new gap in the reference
	int refExtendPen, // penalty for extending a gap in the reference
	int gapBarrier,   // # bases on either side of alignment that must be gap-free
	BTDnaString& decoded, // dest for decoded bases
	EList<Edit>& nedits, // base edits
	EList<Edit>& aedits, // ambiguous nucleotides resolved
	EList<Edit>& cedits, // color miscalls
	EList<Edit>& ccedits, // all color edits (due both to base edits and miscalls)
	bool constrainNumGaps,
	int curC,
	RandomSource& rand,
	int row,
	int col,
	bool refExtend,
	bool readExtend,
	GappedScore score, // the score so far
	int recurDepth)
{
	// col and row are initialized to readf-readi and reff-refi-1 and
	// passed in by the non-recursive caller.
	assert_range(0, (int)table_.size()-1, row);
	assert_eq(readf-readi+1, decoded.length());
	int lastC = curC;
	while(col > 0 || row > 0) {
		assert_range(0, (int)table_.size()-1, row);
		int tabcol = col - row + refGaps;
		if(row < refGaps) tabcol -= (refGaps - row);
		assert_range(0, (int)table_[row].size()-1, tabcol);
		decoded.set(curC, row);
		assert_range(0, 3, curC);
		ASSERT_ONLY(GappedScore scoreThisRound = score);
		ASSERT_ONLY(GappedScore bestThisRound = table_[row][tabcol].best[curC]);
		assert_gt(table_[row][tabcol].mask[curC].numPossible(), 0);
		pair<int, int> cur =
			table_[row][tabcol].mask[curC].randBacktrack(rand);
		int altsLeft = table_[row][tabcol].mask[curC].numPossible();
		bool succeeded = true;
		bool oldRefExtend = refExtend;
		bool oldReadExtend = readExtend;
		size_t neditsSz  = nedits.size();
		size_t aeditsSz  = aedits.size();
		size_t ceditsSz  = cedits.size();
		size_t cceditsSz = ccedits.size();
		GappedScore lastScore = score;
		switch(cur.first) {
			case GAPPED_BT_DIAG: {
				refExtend = readExtend = false;
				assert_gt(row, 0); assert_gt(col, 0);
				assert_neq(-1, cur.second);
				// Check for nucleotide mismatch at source ("source" =
				// lower-right cell)
				if(!matches(curC, ref[refi + col])) {
					score.score -= snpPhred;
					int refM = (int)ref[refi + col];
					assert_range(1, 16, refM);
					assert_eq(1, alts5[refM]);
					int refC = firsts5[refM];
					assert_range(0, 3, refC);
					nedits.push_back(Edit(row, "ACGT"[refC], "ACGTN"[curC], EDIT_TYPE_MM));
					assert(nedits.back().repOk());
					assert_eq(decoded.toChar(row), "ACGTN"[curC]);
				}
				row--; col--;
				// Check for color mismatch
				int readC = read[row];
				int decC = dinuc2color[cur.second][(int)decoded[row+1]];
				int refCmask = dnamasks2colormask[(int)ref[refi+col]][(int)ref[refi+col+1]];
				assert_lt(decC, 4);
				if(decC != readC) {
					score.score -= QUAL(row);
					cedits.push_back(Edit(row, "ACGT"[decC], "ACGTN"[readC], EDIT_TYPE_MM));
					assert(cedits.back().repOk());
				}
				if(readC == 4 || (refCmask & (1 << readC)) == 0) {
					ccedits.push_back(Edit(row, mask2dna[refCmask], "ACGTN"[readC], EDIT_TYPE_MM));
					assert(ccedits.back().repOk());
				}
				assert_range(0, (int)table_.size()-1, row);
				ASSERT_ONLY(tabcol = col - row + refGaps);
				ASSERT_ONLY(if(row < refGaps) tabcol -= (refGaps - row));
				assert_range(0, (int)table_[row].size()-1, tabcol);
				assert(score.valid());
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= table_[row][tabcol].best[cur.second]);
				// Make sure that both changed in the same way
				assert_eq(scoreThisRound, bestThisRound);
				lastC = curC;
				curC = cur.second;
				break;
			}
			case GAPPED_BT_REF_OPEN: {
				if(delAllow_[row] == 0) {
					// This delete isn't allowed
					succeeded = false;
					break;
				}
				refExtend = true; readExtend = false;
				assert_gt(row, 0);
				assert_neq(-1, cur.second);
				nedits.push_back (Edit(row, '-', "ACGTN"[curC], EDIT_TYPE_REF_GAP));
				assert(nedits.back().repOk());
				assert_eq(decoded.toChar(row), "ACGTN"[curC]);
				assert_geq(row, gapBarrier);
				assert_geq((int)(readf-readi-row-1), gapBarrier-1);
				if(constrainNumGaps) {
					delAllow_.set(delAllow_[row]-1, row);
				}
				row--;
				// Check for color miscall
				int decC = dinuc2color[cur.second][(int)decoded[row+1]];
				assert_lt(decC, 4);
				if(decC != read[row]) {
					score.score -= QUAL(row);
					cedits.push_back(Edit(row, "ACGT"[decC], "ACGTN"[(int)read[row]], EDIT_TYPE_MM));
					assert(cedits.back().repOk());
				}
				ccedits.push_back(Edit(row, '-', "ACGTN"[(int)read[row]], EDIT_TYPE_REF_GAP));
				assert(ccedits.back().repOk());
				score.score -= refOpenPen;
				score.incGaps();
				assert_leq(score.gaps, readGaps + refGaps);
				assert_range(0, (int)table_.size()-1, row);
				ASSERT_ONLY(tabcol = col - row + refGaps);
				ASSERT_ONLY(if(row < refGaps) tabcol -= (refGaps - row));
				assert_range(0, (int)table_[row].size()-1, tabcol);
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= table_[row][tabcol].best[cur.second]);
				assert_eq(scoreThisRound, bestThisRound);
				lastC = curC;
				curC = cur.second;
				break;
			}
			case GAPPED_BT_REF_EXTEND: {
				if(delAllow_[row] == 0) {
					// This delete isn't allowed
					succeeded = false;
					break;
				}
				refExtend = true; readExtend = false;
				assert_gt(row, 1);
				assert_neq(-1, cur.second);
				nedits.push_back(Edit(row, '-', "ACGTN"[curC], EDIT_TYPE_REF_GAP));
				assert(nedits.back().repOk());
				assert_eq(decoded.toChar(row), "ACGTN"[curC]);
				assert_geq(row, gapBarrier);
				assert_geq((int)(readf-readi-row-1), gapBarrier-1);
				if(constrainNumGaps) {
					delAllow_.set(delAllow_[row]-1, row);
				}
				row--;
				// Check for color mismatch
				int decC = dinuc2color[cur.second][(int)decoded[row+1]];
				if(decC != read[row]) {
					score.score -= QUAL(row);
					cedits.push_back(Edit(row, "ACGT"[decC], "ACGTN"[(int)read[row]], EDIT_TYPE_MM));
					assert(cedits.back().repOk());
				}
				ccedits.push_back(Edit(row, '-', "ACGTN"[(int)read[row]], EDIT_TYPE_REF_GAP));
				assert(ccedits.back().repOk());
				score.score -= refExtendPen;
				score.incGaps();
				assert_leq(score.gaps, readGaps + refGaps);
				assert_range(0, (int)table_.size()-1, row);
				ASSERT_ONLY(tabcol = col - row + refGaps);
				ASSERT_ONLY(if(row < refGaps) tabcol -= (refGaps - row));
				assert_range(0, (int)table_[row].size()-1, tabcol);
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= table_[row][tabcol].best[cur.second]);
				assert_eq(scoreThisRound, bestThisRound);
				lastC = curC;
				curC = cur.second;
				break;
			}
			case GAPPED_BT_READ_OPEN: {
				if(insAllow_[row] == 0) {
					// This insert isn't allowed
					succeeded = false;
					break;
				}
				refExtend = false; readExtend = true;
				assert_gt(col, 0);
				assert_eq(-1, cur.second);
				nedits.push_back(Edit(row+1, "ACGT"[firsts5[(int)ref[refi+col]]], '-', EDIT_TYPE_READ_GAP));
				assert(nedits.back().repOk());
				assert_geq(row, gapBarrier);
				assert_geq((int)(readf-readi-row-1), gapBarrier-1);
				col--;
				assert_gt(insAllow_[row], 0);
				if(constrainNumGaps) {
					insAllow_.set(insAllow_[row]-1, row);
				}
				score.score -= readOpenPen;
				score.incGaps();
				// Get the reference color opposite the gap so that we
				// can add it to ccedits.
				int refCmask = dnamasks2colormask[(int)ref[refi+col]][(int)ref[refi+col+1]];
				char refc = mask2dna[refCmask];
				ccedits.push_back(Edit(row, refc, '-', EDIT_TYPE_READ_GAP));
				assert(ccedits.back().repOk());
				assert_leq(score.gaps, readGaps + refGaps);
				assert_range(0, (int)table_.size()-1, row);
				ASSERT_ONLY(tabcol = col - row + refGaps);
				ASSERT_ONLY(if(row < refGaps) tabcol -= (refGaps - row));
				assert_range(0, (int)table_[row].size()-1, tabcol);
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= table_[row][tabcol].best[curC]);
				assert_eq(scoreThisRound, bestThisRound);
				lastC = curC;
				break;
			}
			case GAPPED_BT_READ_EXTEND: {
				if(insAllow_[row] == 0) {
					// This insert isn't allowed
					succeeded = false;
					break;
				}
				refExtend = false; readExtend = true;
				assert_gt(col, 1);
				assert_eq(-1, cur.second);
				nedits.push_back(Edit(row+1, "ACGT"[firsts5[(int)ref[refi+col]]], '-', EDIT_TYPE_READ_GAP));
				assert(nedits.back().repOk());
				assert_geq(row, gapBarrier);
				assert_geq((int)(readf-readi-row-1), gapBarrier-1);
				col--;
				assert_gt(insAllow_[row], 0);
				if(constrainNumGaps) {
					insAllow_.set(insAllow_[row]-1, row);
				}
				score.score -= readExtendPen;
				score.incGaps();
				// Get the reference color opposite the gap so that we
				// can add it to ccedits.
				int refCmask = dnamasks2colormask[(int)ref[refi+col]][(int)ref[refi+col+1]];
				char refc = mask2dna[refCmask];
				ccedits.push_back(Edit(row, refc, '-', EDIT_TYPE_READ_GAP));
				assert(ccedits.back().repOk());
				assert_leq(score.gaps, readGaps + refGaps);
				assert_range(0, (int)table_.size()-1, row);
				ASSERT_ONLY(tabcol = col - row + refGaps);
				ASSERT_ONLY(if(row < refGaps) tabcol -= (refGaps - row));
				assert_range(0, (int)table_[row].size()-1, tabcol);
				ASSERT_ONLY(scoreThisRound = score - scoreThisRound);
				ASSERT_ONLY(bestThisRound -= table_[row][tabcol].best[curC]);
				assert_eq(scoreThisRound, bestThisRound);
				lastC = curC;
				break;
			}
			default: throw 1;
		}
		if(succeeded && altsLeft > 0) {
			score = backtrackRecursive(
				read,
				qual,
				readi,
				readf,
				ref,
				refi,
				reff,
				snpPhred,
				decodeScore,
				readGaps,
				refGaps,
				readOpenPen,
				readExtendPen,
				refOpenPen,
				refExtendPen,
				gapBarrier,
				decoded,
				nedits,
				aedits,
				cedits,
				ccedits,
				constrainNumGaps,
				curC,
				rand,
				row,
				col,
				refExtend,
				readExtend,
				score,
				recurDepth+1);
			if(score.gaps != 0) {
				// Success
				return score;
			}
			succeeded = false;
		}
		if(!succeeded) {
			if(altsLeft == 0) {
				score.gaps = score.score = 0;
				return score;
			}
			// Restore curC, score, row, col, nedits, aedits, cedits,
			// ccedits, refExtend, readExtend, insAllow_ and delAllow_
			// to how they were before the recursive call and continue
			nedits.resize(neditsSz);
			aedits.resize(aeditsSz);
			cedits.resize(ceditsSz);
			ccedits.resize(cceditsSz);
			switch(cur.first) {
				case GAPPED_BT_DIAG:
					row++; col++;
					break;
				case GAPPED_BT_REF_OPEN:
					if(constrainNumGaps) delAllow_.set(delAllow_[row]+1, row);
					row++;
					break;
				case GAPPED_BT_REF_EXTEND:
					if(constrainNumGaps) delAllow_.set(delAllow_[row]+1, row);
					row++;
					break;
				case GAPPED_BT_READ_OPEN:
					if(constrainNumGaps) insAllow_.set(insAllow_[row]+1, row);
					col++;
					break;
				case GAPPED_BT_READ_EXTEND:
					if(constrainNumGaps) insAllow_.set(insAllow_[row]+1, row);
					col++;
					break;
				default:
					break;
			}
			// Restore score
			score = lastScore;
			refExtend = oldRefExtend;
			readExtend = oldReadExtend;
			curC = lastC;
		}
	}
	assert_eq(0, row);
	assert_eq(0, col);
	// Last thing to do: factor in SNP penalty from upper-left-hand
	// cell, if applicable
	if(!matches(curC, ref[refi])) {
		score.score -= snpPhred;
		int refM = (int)ref[refi];
		assert_range(1, 16, refM);
		assert_eq(1, alts5[refM]);
		int refC = firsts5[refM];
		assert_range(0, 3, refC);
		nedits.push_back(Edit(0, "ACGT"[refC], "ACGTN"[curC], EDIT_TYPE_MM));
		assert(nedits.back().repOk());
	}
	decoded.set(curC, 0);
	nedits.reverse();
	aedits.reverse();
	cedits.reverse();
	ccedits.reverse();
	assert(Edit::repOk(nedits, decoded));
#ifndef NDEBUG
	BTDnaString refstr;
	for(size_t i = refi; i < reff; i++) {
		refstr.append(firsts5[(int)ref[i]]);
	}
	BTDnaString editstr;
	Edit::toRef(decoded, nedits, editstr);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           edits: ";
		Edit::print(cerr, nedits);
		cerr << endl;
		cerr << "    decoded nucs: " << decoded << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
	
	// done
	assert_eq(score, decodeScore);
	assert_eq(score.gaps, readGaps + refGaps);
	assert_gt(score.gaps, 0);
	return score;
}

/**
 * Update a GappedCell's best[] and mask[] arrays with respect to its
 * neighbor on the left.
 */
void GappedCell::updateHoriz(
	const GappedCell& lc,
	int totGaps,
	int readOpenPen,
	int readExtendPen)
{
	for(int to = 0; to < 4; to++) {
		if(lc.mask[to].empty()) {
			// There's no way to get to the cell to our left within
			// budget if the read character assigned to this row is
			// 'to', so try next character assignment.
			continue;
		}
		GappedScore leftBest = lc.best[to];
		if(!leftBest.valid()) continue;
		// *Don't* penalize for a nucleotide mismatch because we must
		// have already done that in a previous vertical or diagonal
		// step.
		if(lc.mask[to].readExtendPossible()) {
			// Read gap extension possible?
			GappedScore ex = leftBest;
			ex.incGaps();
			assert(ex.valid());
			if(ex.gaps <= totGaps) {
				ex.score -= readExtendPen;
				if(ex >= best[to]) {
					if(ex > best[to]) {
						mask[to].clear();
						best[to] = ex;
					}
					mask[to].rdex = 1;
				}
			} else {
				// # gaps exceeded; don't consider
			}
		}
		if(lc.mask[to].readOpenPossible()) {
			// Read gap open possible?
			GappedScore ex = leftBest;
			ex.incGaps();
			assert(ex.valid());
			if(ex.gaps <= totGaps) {
				ex.score -= readOpenPen;
				if(ex >= best[to]) {
					if(ex > best[to]) {
						mask[to].clear();
						best[to] = ex;
					}
					mask[to].rdop = 1;
				}
			} else {
				// # gaps exceeded; don't consider
			}
		}
	}
}

/**
 * Update a GappedCell's best[] and mask[] arrays with respect to its
 * neighbor on the left.
 *
 * A subtlety here is that, when traveling vertically through the
 * dynamic programming matrix, we are saying that each row from which
 * there is a vertical downward path is aligned against a gap in the
 * reference.  Thus, those rows should not be penalized as SNPs if the
 * chosen character
 */
void GappedCell::updateVert(const GappedCell& uc,
                            int totGaps,
                            int prevColor, // color b/t this row, one above
                            int prevQual,  // quality of color
                            int refOpenPen,
                            int refExtendPen)
{
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		// No SNP penalty because destination read char aligns to a
		// gap in the reference.
		GappedScore from[] = { uc.best[0], uc.best[1],
		                       uc.best[2], uc.best[3] };
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][prevColor];
		// Reward the 'from' that corroborates the color
		if(goodfrom < 4) from[goodfrom].score += prevQual;
		for(int fr = 0; fr < 4; fr++) {
			if(!from[fr].valid()) continue;
			from[fr].incGaps();
			if(from[fr].gaps > totGaps) {
				continue;
			}
			from[fr].score -= prevQual;
			if(uc.mask[fr].refExtendPossible()) {
				// Extend is possible
				from[fr].score -= refExtendPen;
				if(from[fr] >= best[to]) {
					if(from[fr] > best[to]) {
						best[to] = from[fr];
						mask[to].clear();
					}
					mask[to].rfex |= (1 << fr);
				}
				// put it back
				from[fr].score += refExtendPen;
			}
			if(uc.mask[fr].refOpenPossible()){
				// Open is possible
				from[fr].score -= refOpenPen;
				if(from[fr] >= best[to]) {
					if(from[fr] > best[to]) {
						best[to] = from[fr];
						mask[to].clear();
					}
					mask[to].rfop |= (1 << fr);
				}
			}
		}
	}
}

/**
 * Update a GappedCell's best[] and mask[] arrays with respect to its
 * neighbor up and to the left.  SNPs are charged
 */
void GappedCell::updateDiag(const GappedCell& uc,
                            int snpPhred,
                            int refMask,
                            int prevColor, // color b/t this row, one above
                            int prevQual)  // quality of color
{
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...
		int add = (matches(to, refMask) ? 0 : snpPhred);
		GappedScore from[] = { uc.best[0] - add, uc.best[1] - add,
		                       uc.best[2] - add, uc.best[3] - add };
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][prevColor];
		// Reward the preceding position
		if(goodfrom < 4) {
			from[goodfrom].score += prevQual;
		}
		for(int fr = 0; fr < 4; fr++) {
			if(!from[fr].valid()) continue;
			from[fr].score -= prevQual;
			if(from[fr] >= best[to]) {
				if(from[fr] > best[to]) {
					best[to] = from[fr];
					mask[to].clear();
					assert_eq(0, mask[to].diag);
				}
				mask[to].diag |= (1 << fr);
			}
		}
	}
}

/**
 * Update a GappedCell's best[] and mask[] arrays with respect to its
 * neighbors both above and above-and-to-the-left.  Assumes best[] and
 * mask[] have been initialized properly.
 */
void GappedCell::updateDiagVert(
	const GappedCell& dc,
	const GappedCell& uc,
	int totGaps,
	int snpPhred,
	int refMask,
	int prevColor, // color b/t this row, one above
	int prevQual,  // quality of color
	int refOpenPen,
	int refExtendPen)
{
	for(int to = 0; to < 4; to++) {
		// Assuming that the read character in this row is 'to'...

		// Penalize a value for 'to' that doesn't match the reference
		// character for this column
		int add = (matches(to, refMask) ? 0 : snpPhred);
		// Calculate which 'from' character isn't going to
		// cause a color mismatch
		const int goodfrom = nuccol2nuc[to][prevColor];

		//
		// Diagonal case
		//
		{
			GappedScore from[] = {
				dc.best[0] - add, dc.best[1] - add,
				dc.best[2] - add, dc.best[3] - add };
			// Reward the preceding position
			if(goodfrom < 4) {
				from[goodfrom].score += prevQual;
			}
			for(int fr = 0; fr < 4; fr++) {
				if(!from[fr].valid()) continue;
				from[fr].score -= prevQual;
				if(from[fr] >= best[to]) {
					if(from[fr] > best[to]) {
						best[to] = from[fr];
						mask[to].clear();
						assert_eq(0, mask[to].diag);
					}
					mask[to].diag |= (1 << fr);
				}
			}
		}

		//
		// Vertical case
		//
		{
			// No SNP penalty because destination read char aligns to a
			// gap in the reference
			GappedScore from[] = {
				uc.best[0], uc.best[1], uc.best[2], uc.best[3] };
			// Reward the source that doesn't violate the read color
			if(goodfrom < 4) {
				from[goodfrom].score += prevQual;
			}
			for(int fr = 0; fr < 4; fr++) {
				if(!from[fr].valid()) continue;
				from[fr].incGaps();
				if(from[fr].gaps > totGaps) {
					//
					continue;
				}
				from[fr].score -= prevQual;
				if(uc.mask[fr].refExtendPossible()) {
					// Extend is possible
					from[fr].score -= refExtendPen;
					if(from[fr] >= best[to]) {
						if(from[fr] > best[to]) {
							best[to] = from[fr];
							mask[to].clear();
						}
						mask[to].rfex |= (1 << fr);
					}
					// put it back
					from[fr].score += refExtendPen;
				}
				if(uc.mask[fr].refOpenPossible()){
					// Open is possible
					from[fr].score -= refOpenPen;
					if(from[fr] >= best[to]) {
						if(from[fr] > best[to]) {
							best[to] = from[fr];
							mask[to].clear();
						}
						mask[to].rfop |= (1 << fr);
					}
				}
			}
		}
	}
}

/**
 * Decode the colorspace read 'read' as aligned against the reference
 * string 'ref', assuming that it's a hit without gaps (only matches
 * and mismatches).
 */
TScore ColorspaceDecoderGapped::decode(
	const BTDnaString& read,
	const BTString& qual,
	size_t readi,
	size_t readf,
	BTString& ref, // this will be overwritten
	size_t refi,
	size_t reff,
	int snpPhred,
	int maxCost,
	int readGaps,
	int refGaps,      // # of reference gaps in
	const EList<Edit>& edits, // colorspace-to-colorspace edits
	int readOpenPen,  // penalty for opening a new gap in the read
	int readExtendPen,// penalty for extending a gap in the read
	int refOpenPen,   // penalty for opening a new gap in the reference
	int refExtendPen, // penalty for extending a gap in the reference
	int gapBarrier,   // # bases on either side of alignment that must be gap-free
	BTDnaString& decoded,
	EList<Edit>& nedits, // base edits
	EList<Edit>& aedits, // ambiguous base resolutions
	EList<Edit>& cedits, // color miscalls
	EList<Edit>& ccedits, // all color edits (due both to base edits and miscalls)
	RandomSource& rand)
{
	assert_lt(refi, reff);
	assert_lt(readi, readf);
	assert(readGaps > 0 || refGaps > 0);
	assert_eq((reff-refi-readGaps)-1, (readf-refGaps)-readi);
	assert_gt(edits.size(), 0);

	int totGaps = readGaps + refGaps;

	insAllow_.resize(read.length()+1);
	delAllow_.resize(read.length()+1);
	insAllow_.fill(0);
	delAllow_.fill(0);
	ASSERT_ONLY(int ins = 0, dels = 0);
	for(size_t i = 0; i < edits.size(); i++) {
		size_t pos = edits[i].pos;
		int dleft = (int)pos;
		int dright = (int)(read.length() - pos - 1);
		assert_leq(pos, read.length());
		if(edits[i].isReadGap()) {
			if(dleft >= gapBarrier && dright >= gapBarrier-1) {
				insAllow_.set(insAllow_[pos]+1, pos);
			}
			if(dleft+1 >= gapBarrier && (dright-1) >= gapBarrier-1) {
				insAllow_.set(insAllow_[pos+1]+1, pos+1);
			}
			ASSERT_ONLY(ins++);
		} else if(edits[i].isRefGap()) {
			if(dleft >= gapBarrier && dright >= gapBarrier) {
				delAllow_.set(delAllow_[pos]+1, pos);
			}
			if(dleft+1 >= gapBarrier && (dright-1) >= gapBarrier) {
				delAllow_.set(delAllow_[pos+1]+1, pos+1);
			}
			ASSERT_ONLY(dels++);
		}
	}
	assert_eq(ins, readGaps);
	assert_eq(dels, refGaps);

	//
	// Initialize the first row
	//
	table_.resize(1); // add first row to row list
	int wlo = 0;
	int whi = min(readGaps, (int)(reff-refi)-1);
	table_[0].resize(whi-wlo+1); // add columns to first row
	// Go from leftmost to rightmost column for the first row
	for(int to = 0; to < 4; to++) {
		table_[0][0].mask[to].clear();
		table_[0][0].best[to].gaps = 0;
		if(matches(to, ref[refi])) {
			// The assigned subject nucleotide matches the reference;
			// no penalty
			table_[0][0].best[to].score = 0;
			table_[0][0].mask[to].diag = 0xf;
		} else if(snpPhred <= maxCost) {
			// The assigned subject nucleotide does not match the
			// reference nucleotide, so we add a SNP penalty
			table_[0][0].best[to].score = -snpPhred;
			table_[0][0].mask[to].diag = 0xf;
		}
		// We're in the corner, so penalty for getting here via any
		// kind of gap is infinity.  Also because we're in a corner,
		// there's no opportunity yet to violate a color.
	}
	// Calculate starting values for the rest of the columns in the
	// first row.
	for(int col = 1; col <= whi; col++) {
		table_[0][col].clear();
		if(col <= insAllow_[0]) {
			table_[0][col].updateHoriz(
				table_[0][col-1], totGaps, readOpenPen, readExtendPen);
		}
	}

	//
	// Calculate all subsequent rows
	//

	// Do rest of table
	int lastwlo = wlo;
	int lastwhi = whi;
	for(int row = 1; row <= (int)(readf-readi); row++) {
		table_.expand(); // add another row
		// Index of leftmost column to visit, absolute (not relative to table)
		wlo = max(row - refGaps, 0);
		// Index of rightmost column to visit, absolute (not relative to table)
		whi = row + readGaps;
		bool rightCurtailed = false;
		if(whi > (int)(reff-refi-1)) {
			whi = (int)(reff-refi-1);
			rightCurtailed = true;
		}
		assert_leq(wlo, whi);
		assert_lt(whi, (int)(reff-refi));
		table_.back().resize(whi-wlo+1); // add enough space for columns

		assert_range(1, (int)qual.length(), row);
		assert_range(1, (int)read.length(), row);

		int q = QUAL(row-1); // quality for color b/t this row, previous
		int c = read[row-1]; // color b/t this row, previous

		bool vertOk  = delAllow_[row] > 0;
		bool horizOk = insAllow_[row] > 0;

		// Go from leftmost to rightmost column for this row
		ASSERT_ONLY(bool validInRow = false);
		int cols = whi - wlo + 1;
		for(int col = wlo; col <= whi; col++) {
			int r = ref[refi + col];
			GappedCell& cur = table_[row][col-wlo];
			cur.clear();
			// Can we do a diagonal update?
			if(col > lastwlo) {
				// Can do diagonal, what about vertical?
				if(vertOk && (col < whi || rightCurtailed)) {
					// Can do vertical
					assert_lt(col-wlo, (int)table_[row-1].size());
					GappedCell& dg = table_[row-1][col-wlo];
					assert_lt(col-wlo+1, (int)table_[row-1].size());
					GappedCell& up = table_[row-1][col-wlo+1];
					cur.updateDiagVert(dg, up, totGaps, snpPhred, r, c,
					                   q, refOpenPen, refExtendPen);
				} else {
					// Can't do vertical, just do diagonal
					int prevcols = lastwhi - lastwlo + 1;
					int prevcol = col-wlo;
					if(prevcols < cols) {
						prevcol -= (cols - prevcols);
					}
					assert_geq(prevcol, 0);
					assert_lt(prevcol, (int)table_[row-1].size());
					GappedCell& dg = table_[row-1][prevcol];
					cur.updateDiag(dg, snpPhred, r, c, q);
				}
			} else if(vertOk) {
				// No, can't do diagonal
				// Must be able to do vertical
				GappedCell& up = table_[row-1][0];
				cur.updateVert(up, totGaps, c, q, refOpenPen, refExtendPen);
			}
			if(horizOk && col > wlo) {
				// Can do horizontal
				GappedCell& lf = table_[row][col-wlo-1];
				cur.updateHoriz(lf, totGaps, readOpenPen, readExtendPen);
			}
			// 'cur' is now initialized
			ASSERT_ONLY(if(cur.valid()) validInRow = true);
		}
		assert(validInRow);
		lastwlo = wlo;
		lastwhi = whi;
	}

	// All table cells are now initialized
	const GappedCell& lr = table_.back().back();
	int lastReadC = 0;
	GappedScore bestScore = GappedScore::INVALID();
	for(int i = 0; i < 4; i++) {
		if(lr.best[i] >= bestScore) {
			if(lr.best[i] > bestScore) {
				lastReadC = 0;
				bestScore = lr.best[i];
				assert(bestScore.valid());
				assert_eq(bestScore.gaps, refGaps + readGaps);
			}
			lastReadC |= (1 << i);
		}
	}
	assert_neq(0, lastReadC);
	assert(bestScore.valid());
	assert_eq(bestScore.gaps, refGaps + readGaps);
#ifndef NDEBUG
	int upperBoundScore = 0;
	if(refGaps > 0) {
		upperBoundScore += refOpenPen;
		if(refGaps > 1) {
			upperBoundScore += ((refGaps-1)*refExtendPen);
		}
	}
	if(readGaps > 0) {
		upperBoundScore += readOpenPen;
		if(readGaps > 1) {
			upperBoundScore += ((readGaps-1)*readExtendPen);
		}
	}
	assert_leq(bestScore.score, -upperBoundScore);
#endif
	assert_eq(bestScore.gaps, refGaps + readGaps);

	GappedScore actualScore =
	backtrack(read, qual, readi, readf, ref, refi, reff, snpPhred,
	          bestScore, readGaps, refGaps, readOpenPen, readExtendPen,
	          refOpenPen, refExtendPen, gapBarrier, decoded, nedits,
	          aedits, cedits, ccedits, false, randFromMask(rand, lastReadC),
	          rand);
	assert_eq(actualScore, bestScore);
	return -actualScore.score;
}

TScore ColorspaceDecoder::decode(
	const BTDnaString& read, // read colors
	const BTString& qual, // read qualities
	size_t readi, // offset of first character within 'read' to consider
	size_t readf, // offset of last char (exclusive) in 'read' to consider
	BTString& ref,    // reference sequence, as masks
	size_t refi,  // offset of first character within 'ref' to consider
	size_t reff,  // offset of last char (exclusive) in 'ref' to consider
	int snpPen,   // cost of a SNP when decoding
	int maxCost,  // maximum permitted cost of decoded sequence
	int readGaps, // number of read gaps in the color-to-color alignment
	int refGaps,  // number of ref gaps in the color-to-color alignment
	const EList<Edit>& edits, // colorspace-to-colorspace edits
	int readOpenPen,  // penalty for opening a new gap in the read
	int readExtendPen,// penalty for extending a gap in the read
	int refOpenPen,   // penalty for opening a new gap in the reference
	int refExtendPen, // penalty for extending a gap in the reference
	int gapBarrier,   // no gaps permitted within this many chars of ends
	bool exEnds,   // true -> nucleotide alignment is 1 char shorter than color read
	bool maqRound, // true -> use Maq-like rounding
	BTDnaString& dseq, // decoded sequence installed here
	BTString& dqual,   // decoded qualities installed here
	EList<Edit>& nedits, // decoded nucleotide edits installed here
	EList<Edit>& aedits, // resolved ambiguous nucleotides installed here
	EList<Edit>& cedits, // decoded color miscalls installed here
	EList<Edit>& ccedits, // decoded color edits installed here
	RandomSource& rand) // pseudo-random generator
{
	dseq.clear();
	dqual.clear();
	assert(nedits.empty());
	assert(aedits.empty());
	assert(cedits.empty());
	assert(ccedits.empty());
	TScore score;
#ifndef NDEBUG
	for(size_t i = refi; i < reff; i++) {
		assert_range(0, 15, (int)ref[i]);
	}
#endif
	if(readGaps == 0 && refGaps == 0) {
		score = ungapped_.decode(
			read, qual, readi, readf, ref, refi, reff,
			snpPen, dseq, nedits, aedits, cedits, ccedits, rand);
	} else {
		score = gapped_.decode(
			read, qual, readi, readf, ref, refi, reff,
			snpPen, maxCost, readGaps, refGaps, edits,
			readOpenPen, readExtendPen, refOpenPen, refExtendPen,
			gapBarrier, dseq, nedits, aedits, cedits, ccedits, rand);
	}
	// Possibly trim the decoded nucleotides; calculate the decoded
	// qualities
	size_t cedIdx = 0;
	int lastQ = 0, curQ = 0;
	for(size_t i = 0; i < dseq.length(); i++) {
		// Color mismatches penalize quality
		bool mm = false;
		while(cedIdx < cedits.size() && cedits[cedIdx].pos == i) {
			assert_lt(i, dseq.length()-1); // 1 less color
			if(cedits[cedIdx++].isMismatch()) mm = true;
		}
		curQ = 0;
		if(i < qual.length() && read[i] != 4) {
			curQ = mmPenalty(maqRound, phredcToPhredq(qual[i]));
		}
		if(mm) curQ = -curQ;
		if(curQ + lastQ + (int)'!' > 126) {
			dqual.append(126);
		} else if (curQ + lastQ < 0) {
			dqual.append('!');
		} else {
			dqual.append(curQ + lastQ + (int)'!');
		}
		lastQ = curQ;
	}
	if(exEnds) {
		// Delete ends
		ASSERT_ONLY(size_t sz = dseq.length());
		dseq.trimEnd(1);
		dseq.remove(0);
		assert_eq(dseq.length()+2, sz);
		dqual.trimEnd(1);
		dqual.remove(0);
		for(int i = 0; i < (int)nedits.size(); i++) {
			assert(i == 0 || nedits[i].pos >= nedits[i-1].pos);
			// Account for the excluded end
			if(nedits[i].pos == 0) {
				// Chop
				nedits.remove(0);
				i--;
			} else if(nedits[i].pos > dseq.length()) {
				// Chop everyone else off
				nedits.resize(i);
				break;
			} else {
				nedits[i].pos--;
			}
		}
		for(int i = 0; i < (int)aedits.size(); i++) {
			assert(i == 0 || aedits[i].pos >= aedits[i-1].pos);
			// Account for the excluded end
			if(aedits[i].pos == 0) {
				// Chop
				aedits.remove(0);
				i--;
			} else if(aedits[i].pos > dseq.length()) {
				// Chop everyone else off
				aedits.resize(i);
				break;
			} else {
				aedits[i].pos--;
			}
		}
	}
#ifndef NDEBUG
	for(int i = 0; i < (int)nedits.size(); i++) {
		assert_lt(nedits[i].pos, dseq.length());
	}
	for(int i = 0; i < (int)aedits.size(); i++) {
		assert_lt(aedits[i].pos, dseq.length());
	}
	for(int i = 0; i < (int)cedits.size(); i++) {
		assert_lt(cedits[i].pos, readf-readi);
	}
	for(int i = 0; i < (int)ccedits.size(); i++) {
		assert_lt(ccedits[i].pos, readf-readi);
	}
#endif
	return score;
}

#ifdef MAIN_COLOR_DEC

#include <sstream>
#include <getopt.h>

static const char *short_opts = "s:m:r:d:i:";
static struct option long_opts[] = {
	{(char*)"snppen",  required_argument, 0, 's'},
	{(char*)"misspen", required_argument, 0, 'm'},
	{(char*)"seed",    required_argument, 0, 'r'},
	{(char*)"inserts", required_argument, 0, 'i'},
	{(char*)"deletes", required_argument, 0, 'd'}
};

static void printUsage(ostream& os) {
	os << "Usage: color_dec <read-color-seq> <ref-nuc-seq> [options]*" << endl;
	os << "Options:" << endl;
	os << "  -s/--snppen <int>   penalty incurred by SNP; used for decoding" << endl;
	os << "  -m/--misspen <int>  quality to use for read chars" << endl;
	os << "  -r/-seed <int>      seed for pseudo-random generator" << endl;
	os << "  -i/--inserts <int>  # inserts in the color-to-color alignment" << endl;
	os << "  -d/--deletes <int>  # deletes in the color-to-color alignment" << endl;
}

template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

int main(int argc, char **argv) {
	ColorspaceDecoder dec;
	int option_index = 0;
	int next_option;
	int snppen = 30;
	int misspen = 20;
	int refGaps = 0;
	int readGaps = 0;
	unsigned seed = 0;
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 's': snppen     = parse<int>(optarg); break;
			case 'm': misspen    = parse<int>(optarg); break;
			case 'r': seed       = parse<unsigned>(optarg); break;
			case 'i': readGaps   = parse<int>(optarg); break;
			case 'd': refGaps    = parse<int>(optarg); break;
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	srand(seed);
	if(argc - optind < 2) {
		cerr << "Not enough options" << endl;
		printUsage(cerr);
		exit(1);
	}
	BTDnaString read;
	BTString ref;
	if(isalpha(argv[optind][0])) {
		read.installChars(argv[optind]);
	} else {
		read.installColors(argv[optind]);
	}
	ref.install(argv[optind+1]);
	// Convert reference string to masks
	for(size_t i = 0; i < ref.length(); i++) {
		int num = 0;
		int alts[] = {4, 4, 4, 4};
		decodeNuc(toupper(ref[i]), num, alts);
		assert_leq(num, 4);
		assert_gt(num, 0);
		ref.set(0, i);
		for(int j = 0; j < num; j++) {
			ref.set(ref[i] | (1 << alts[j]), i);
		}
	}
	BTDnaString decoded;
	BTString quals;
	EList<Edit> nedits;
	EList<Edit> cedits;
	quals.resize(read.length(), misspen);
	int score = dec.decode(read, quals, 0, read.length(),
	                       ref, 0, ref.length(), snppen,
	                       readGaps, refGaps, decoded, nedits, cedits);

	cout << " Score: " << score << endl;
	cout << "   Read colors: " << endl;
	cout << "     ";
	for(size_t i = 0; i < read.length(); i++) {
		printColor((int)read[i]);
	}
	cout << endl;
	cout << "   Color alignment (decoded): " << endl;
	Edit::printQAlign(cout, "     ", read, cedits);
	cout << "   Nucleotide alignment (decoded): " << endl;
	Edit::printQAlign(cout, "     ", decoded, nedits);
}
#endif
