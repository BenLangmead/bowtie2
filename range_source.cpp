/*
 * range_source.cpp
 *
 *  Created on: Jan 25, 2010
 *      Author: Ben Langmead
 */

#include "range_source.h"

using namespace std;

// Requires that isJoinCell and depth are already set properly
#define VISITED(t,o) (isJoinCell && visit.isVisited(depth, qdep0, gapOff_+o, t))

static void printPadded(ostream& os, const string& s, int amt) {
	if(s.length() < 6) {
		for(size_t i = 0; i < amt - s.length(); i++) {
			os << "0";
		}
	}
	os << s;
}

/**
 * Initialize a new branch object with an empty path.
 */
bool Branch::init(
	AllocOnlyPool<Pos>& pool,
	Visit& visit,
	uint32_t id,      // id of this branch
	Branch *parent,   // pointer to parent branch (NULL -> this is root)
	uint32_t parentId,// id of parent branch
	Edit edit,        // edit at the root of this branch
	uint32_t nedits,  // # edits prior to this branch
	uint32_t qlen,    // length of read
	uint16_t *deps,   // edit budgets
	uint16_t rdepth,  // # read characters consumed so far
	uint16_t len,     // length of streak of matches so far
	uint16_t cost,    // cost so far
	uint16_t ham,     // penalty (non-stratum) part of cost so far
	uint32_t itop,    // top and bot of BW range at root of this branch
	uint32_t ibot,
	const EbwtParams& ep,
	const uint8_t* ebwt)
{
	assert(parent == NULL || edit.initialized());
	assert(parent != NULL || !edit.initialized());
	id_ = id;
	delayedCost_ = 0;
	memcpy(deps_, deps, 4 * sizeof(uint16_t));
	rdepth_ = rdepth;
	len_ = ilen_ = len;
	cost_ = cost;
	ham_ = ham;
	itop_ = top_ = itop;
	ibot_ = bot_ = ibot;
	edit_ = edit;
	assert(!edit.initialized() || bot_ > top_);
	numEdits_ = nedits;
	parent_ = parent;
	ASSERT_ONLY(parentId_ = parentId);
	// Prefetch the portion of the BW string that we'll soon be
	// searching
	prepped_ = false;
	assert(!active_);
	ASSERT_ONLY(active_ = true);
	prep(ep, ebwt);
	// Allocate RangeStates
	if(qlen - rdepth_ > 0) {
		ranges_ = pool.alloc(qlen - rdepth_); // allocated from the RangeStatePool
		if(ranges_ == NULL) {
			return false;
		}
		rangesSz_ = qlen - rdepth_;
	} else {
		ranges_ = NULL;
		rangesSz_ = 0;
	}
	curtailed_ = false;
	exhausted_ = false;
	delayedIncrease_ = false;
	// If we're starting with a non-zero length, that means we're
	// jumping over a bunch of unrevisitable positions.
	int i = (int)deps_[0];
	i = max(0, i - rdepth_);
	for(; i < rangesSz_; i++) {
		ranges_[i].eliminated_ = true;
		assert(eliminated(i));
	}
	gapOff_ = 0;
	if(parent != NULL) {
		gapOff_ = parent->gapOff_;
		if(edit_.isDelete()) {
			gapOff_--;
		} else if(edit_.isInsert()) {
			gapOff_++;
		}
	}
	assert(repOk(qlen));
	return true;
}

/**
 * Using randomness when picking from among multiple choices, pick
 * an edit to make.  Called from Range.splitBranch().  Assumes
 * elimination flags have already been set appropriately.
 */
Edit Branch::pickEdit(int i, RandomSource& rand, uint32_t& top, uint32_t& bot) {
	bool color = false;
	assert_lt(i, rangesSz_);
	Pos& r = ranges_[i];
	Edit ret;
	ret.type = EDIT_TYPE_MM;
	ret.pos = i + rdepth_;
	ret.chr = 0;
	ret.qchr = "ACGTN"[r.flags.qchr];
	assert(!eliminated(i));
	assert_neq(INVALID_COST, r.flags.costlo);
	assert(repOk());
	// Only need to pick randomly if there's a quality tie
	if(r.flags.costlo == r.flags.costlo2) {
		// Sum up range sizes and do a random weighted pick
		uint32_t tot = 0;
		bool candMm[] = {
			!r.flags.mmA && r.flags.qualA == r.flags.costlo,
			!r.flags.mmC && r.flags.qualC == r.flags.costlo,
			!r.flags.mmG && r.flags.qualG == r.flags.costlo,
			!r.flags.mmT && r.flags.qualT == r.flags.costlo };
		unsigned insCost = r.flags.ins_ext ? gInsExtend : gInsOpen;
		bool candIns[] = { false, false, false, false };
		if(insCost == r.flags.costlo) {
			candIns[0] = !r.flags.insA;
			candIns[1] = !r.flags.insC;
			candIns[2] = !r.flags.insG;
			candIns[3] = !r.flags.insT;
		}
		unsigned delCost = r.flags.del_ext ? gDelExtend : gDelOpen;
		bool candDel = (!r.flags.del && delCost == r.flags.costlo);
		for(int j = 0; j < 4; j++) {
			if(candMm[j]) {
				assert_gt(r.bots[j], r.tops[j]);
				tot += (r.bots[j] - r.tops[j]);
			}
			if(gGaps) {
				if(candIns[j]) {
					assert_gt(r.bots[j], r.tops[j]);
					tot += (r.bots[j] - r.tops[j]);
				}
			}
		}
		assert_gt(tot, 0);
		if(gGaps && candDel) {
			// Know what the next char is going to be?  Not
			// necessarily, because we don't know how long the
			// deletion will be.  Just let the chance of choosing a
			// deletion versus either a mismatch or an insertion be
			// 50/50.
			tot *= 2;
		}
		// Randomly choose among the substitutions and insertions,
		// weighted by BW range size
		uint32_t dart = rand.nextU32() % tot;
		for(int j = 0; j < 4; j++) {
			if(candMm[j]) {
				if(dart < (r.bots[j] - r.tops[j])) {
					top = r.tops[j]; bot = r.bots[j];
					r.elimMm(j);
					ret.chr = color ? "0123"[j] : "ACGT"[j];
					r.updateLo();
					return ret;
				}
				dart -= (r.bots[j] - r.tops[j]);
			}
			if(candIns[j]) {
				if(dart < (r.bots[j] - r.tops[j])) {
					top = r.tops[j]; bot = r.bots[j];
					r.elimIns(j);
					ret.type = EDIT_TYPE_INS;
					ret.chr = color ? "0123"[j] : "ACGT"[j];
					ret.qchr = '-';
					r.updateLo();
					return ret;
				}
				dart -= (r.bots[j] - r.tops[j]);
			}
		}
		assert_lt(dart, tot/2);
		r.flags.del = 1;
		ret.type = EDIT_TYPE_DEL;
		ret.chr = '-';
		r.updateLo();
		assert(repOk());
	} else {
		// Pick the first valid edit we see; it should be the last
		int chr = r.pickFirstEdit(ret);
		top = r.tops[chr]; bot = r.bots[chr];
		r.eliminated_ = (r.flags.costlo2 == INVALID_COST);
		if(!r.eliminated_) {
			r.updateLo();
		} else {
			r.flags.costlo = r.flags.costlo2 = INVALID_COST;
		}
		assert(r.repOk());
	}
	// Fixup for the case where a deletion was selected.  In this case,
	// top/bot are a function of information from the previous
	// position.
	if(ret.isDelete()) {
		if(i > ilen_) {
			int pc = ranges_[i-1].flags.qchr;
			assert_lt(pc, 4);
			top = ranges_[i-1].tops[pc];
			bot = ranges_[i-1].bots[pc];
			assert_gt(bot, top);
		} else {
			top = itop_;
			bot = ibot_;
			assert_gt(bot, top);
		}
	} else {
		assert_gt(bot, top);
	}
	// Can't do repOk() here; potentially have to wait for the caller
	// to re-set exhausted_ flag if we've just exhausted our options.
	return ret;
}

/**
 * For each as-yet-uneliminated outward edge, see if we can now
 * eliminate it because the path has been redundantly examined.
 */
bool Branch::eliminateVisited(size_t i, int qdepth, int qdep0, Visit& visit) {
	assert(ranges_ != NULL);
	assert_lt(i, rangesSz_);
	assert_geq(rdepth_ + i, deps_[0]);
	bool isJoinCell0 = visit.isJoinCell(qdepth, qdep0);
	bool isJoinCell1 = visit.isJoinCell(qdepth+1, qdep0);
	if(!isJoinCell0 && !isJoinCell1) return false;
	bool isJoinCell = isJoinCell1;
	Pos& r = ranges_[i];
	bool eliminated = false;
	//
	// Set which mismatches are eliminated
	//
	int depth = qdepth + 1;
	if(!r.flags.mmA && VISITED(r.tops[0], 0)) {
		r.flags.mmA = 1; eliminated = true;
	}
	if(!r.flags.mmC && VISITED(r.tops[1], 0)) {
		r.flags.mmC = 1; eliminated = true;
	}
	if(!r.flags.mmG && VISITED(r.tops[2], 0)) {
		r.flags.mmG = 1; eliminated = true;
	}
	if(!r.flags.mmT && VISITED(r.tops[3], 0)) {
		r.flags.mmT = 1; eliminated = true;
	}
	uint32_t top;
	if(i > 0) {
		int pc = ranges_[i-1].flags.qchr;
		assert_lt(pc, 4);
		top = ranges_[i-1].tops[pc];
	} else {
		top = itop_;
	}
	if(!r.flags.del && VISITED(top, -1)) {
		r.flags.del = 1; eliminated = true;
	}
	depth = qdepth;
	isJoinCell = isJoinCell0;
	if(!r.flags.insA && VISITED(r.tops[0], 1)) {
		r.flags.insA = 1; eliminated = true;
	}
	if(!r.flags.insC && VISITED(r.tops[1], 1)) {
		r.flags.insC = 1; eliminated = true;
	}
	if(!r.flags.insG && VISITED(r.tops[2], 1)) {
		r.flags.insG = 1; eliminated = true;
	}
	if(!r.flags.insT && VISITED(r.tops[3], 1)) {
		r.flags.insT = 1; eliminated = true;
	}
	if(eliminated) {
		r.updateLo();
		r.eliminated_ = (r.flags.costlo == INVALID_COST);
	}
	assert(r.repOk());
	return eliminated;
}

/**
 * Split off a new branch by selecting a good outgoing path and
 * creating a new Branch object for it and inserting that branch
 * into the priority queue.  Mark that outgoing path from the
 * parent branch as eliminated.  If the second-best outgoing path
 * costs more, add the difference to the cost of this branch.
 *
 * We're returning to a curtailed Branch after having explored other
 * paths; because some of the explored paths may have been redundant
 * with paths outgoing from this branch, we may have to eliminate some
 * more paths.
 */
Branch* Branch::splitBranch(AllocOnlyPool<Pos>& rpool,
                            AllocOnlyPool<Branch>& bpool,
                            Visit& visited,
                            const CoordMap& cmap,
                            RandomSource& rand,
                            uint32_t qlen, // full length of query
							int qdep0,
                            uint32_t qualLim,
                            int seedLen,
                            bool& triedGap,
                            const EbwtParams& ep,
                            const uint8_t* ebwt)
{
	assert(!exhausted_);
	assert(active_);
	assert(ranges_ != NULL);
	assert(curtailed_);
	assert(repOk(qlen));
	assert_neq(0xffff, nextIncrCost_);
	int tiedPositions[3];
	int numTiedPositions = 0;
	// Lowest marginal cost incurred by any outgoing edge
	uint16_t bestCost = 0xffff;
	// Next-lowest (can be tied)
	uint16_t nextCost = 0xffff;
	bool eliminatedVisited = false;
	int i = (int)deps_[0];
	i = max(0, i - rdepth_); // skip must-match positions
	// Iterate revisitable positions
	for(; i <= len_; i++) {
		// Still valid options for leaving out of this position?
		if(!eliminated(i)) {
			assert(ranges_[i].repOk());
			if(eliminateVisited(i, cmap.r2q[rdepth_ + i], qdep0, visited)) {
				// At least one outgoing path was eliminated, having
				// been explored in the meantime.
				eliminatedVisited = true;
				if(eliminated(i)) {
					assert(triedGap);
					continue;
				}
				ranges_[i].updateLo();
				assert_lt(ranges_[i].flags.costlo, INVALID_COST);
			}
			uint16_t stratum = (rdepth_ + i < seedLen) ? (1 << 14) : 0;
			uint16_t cost = stratum;
			cost |= ranges_[i].flags.costlo;
			if(cost < bestCost) {
				// Demote the old best to the next-best
				nextCost = bestCost;
				// Update the new best
				bestCost = cost;
				numTiedPositions = 1;
				tiedPositions[0] = i;
			} else if(cost == bestCost) {
				// As good as the best so far
				nextCost = cost;
				assert_gt(numTiedPositions, 0);
				if(numTiedPositions < 3) {
					tiedPositions[numTiedPositions++] = i;
				} else {
					tiedPositions[0] = tiedPositions[1];
					tiedPositions[1] = tiedPositions[2];
					tiedPositions[2] = i;
				}
			} else if(cost < nextCost) {
				// 'cost' isn't better than the best, but it is
				// better than the next-best
				nextCost = cost;
			}
			assert_geq(cost, bestCost);
			// Handle case where *this* position's next-best is
			// even better than the current next-best
			if(ranges_[i].flags.costlo2 < INVALID_COST) {
				nextCost = min<uint16_t>(
					nextCost, ranges_[i].flags.costlo2 | stratum);
			}
		}
	}
	if(bestCost != nextIncrCost_) {
		// This should only happen if our previous calculation of
		// nextIncrCost_ became invalid in the meantime due to the
		// best path having been investigated on another branch.
		assert(eliminatedVisited);
		assert_gt(bestCost, nextIncrCost_);
		if(bestCost == 0xffff) {
			// Branch is now exhausted
			exhausted_ = true;
		} else {
			cost_ = cost_ - nextIncrCost_ + bestCost;
			nextIncrCost_ = bestCost;
			// Caller should re-insert this branch into the queue using its
			// new cost
		}
		return NULL;
	}
	assert_eq(bestCost, nextIncrCost_);
	assert_gt(numTiedPositions, 0);
	int r = 0;
	if(numTiedPositions > 1) {
		r = rand.nextU32() % numTiedPositions;
	}
	int pos = tiedPositions[r];
	// Pick an edit from among the edits tied for lowest cost
	// (using randomness to break ties).  If the selected edit is
	// the last remaining one at this position, 'last' is set to
	// true.
	uint32_t top = 0, bot = 0;
	numChildren_++;
	Edit e = pickEdit(pos, rand, top, bot);
	// If this is a mismatch or an insertion, then top & bot will
	// be the BW range for the new, longer suffix.  If this is a
	// deletion, then top & bot will be the same as before.
	assert_gt(bot, top);
	// Create and initialize a new Branch
	assert_lt((bestCost >> 14), 4);
	uint32_t hamadd = (bestCost & ~0xc000);
	uint16_t depth = pos + rdepth_;
	assert_geq(depth, deps_[0]);
	uint16_t ndep[4] = {deps_[0], deps_[1], deps_[2], deps_[3]};
	for(int i = 0; i < 3; i++) {
		if(depth < deps_[i+1]) ndep[i] = deps_[i+1];
	}
	assert_eq((uint32_t)(cost_ & ~0xc000), (uint32_t)(ham_ + hamadd));
	Branch *newBranch = bpool.alloc();
	if(newBranch == NULL) {
		return NULL; // out of branch memory
	}
	ASSERT_ONLY(splits_++);
	// Create a new branch rooted at the selected outgoing path
	uint16_t newRdepth = rdepth_ + pos + (e.isInsert() ? 0 : 1);
	if(!newBranch->init(
		rpool, visited, bpool.lastId(), this, id_, e, numEdits_+1,
		qlen, ndep, newRdepth, 0, cost_, ham_ + hamadd, top, bot, ep,
		ebwt))
	{
		return NULL;
	}
	if(e.isDelete() && newBranch->rangesSz_ > 0) {
		assert(newBranch->ranges_ != NULL);
		// copy outgoing BW ranges from parent
		Pos *r = newBranch->backPos();
		for(int i = 0; i < 4; i++) {
			r->tops[i] = ranges_[pos].tops[i];
			r->bots[i] = ranges_[pos].bots[i];
		}
	}
	if(e.isGap()) triedGap = true;
	nextIncrCost_ = nextCost;
	if(bestCost != 0xffff && nextCost == 0xffff) {
		// No more valid outgoing paths after this one; remove
		// branch from the PathManager and mark it as exhausted.
		exhausted_ = true;
		// Free the ranges structure
		if(ranges_ != NULL) {
			assert_gt(rangesSz_, 0);
			if(rpool.free(ranges_, rangesSz_)) {
				ranges_ = NULL;
				rangesSz_ = 0;
			}
		}
		assert(repOk(qlen));
	}
	// Cost updates are delayed (until the next time this branch
	// comes to the top of the queue) so that queue re-insertion
	// can be lazy.
	else if(bestCost != nextCost) {
		// We exhausted the last outgoing edge at the current best
		// cost; update the best cost to be the next-best
		assert_gt(nextCost, bestCost);
		delayedCost_ = (cost_ - bestCost + nextCost);
		assert_gt(delayedCost_, cost_); // new cost must be higher
		delayedIncrease_ = true;
		assert(repOk(qlen));
	}
	return newBranch;
}

/**
 * Called when the most recent branch extension resulted in an
 * empty range or some other constraint violation (e.g., a
 * half-and-half constraint).
 */
void Branch::curtail(AllocOnlyPool<Pos>& rpool, int seedLen) {
	assert(active_);
	assert(!curtailed_);
	assert(!exhausted_);
	assert_eq(seedLen, deps_[3]);
	assert(cost_ < 0xc000 || deps_[0] >= seedLen);
	if(ranges_ == NULL) {
		exhausted_ = true;
		curtailed_ = true;
		return;
	}
	uint16_t lowestCost = 0xffff;
	// Iterate over positions in the path looking for the cost of
	// the lowest-cost non-eliminated position
	uint32_t eliminatedStretch = 0;
	int i = (int)deps_[0];
	i = max(0, i - rdepth_);
	for(; i <= len_; i++) {
		if(!eliminated(i)) {
			assert(ranges_[i].repOk());
			eliminatedStretch = 0;
			assert(rdepth_ + i >= seedLen || cost_ < 0xc000);
			uint16_t stratum = (rdepth_ + i < seedLen) ? (1 << 14) : 0;
			uint16_t cost = ranges_[i].flags.costlo | stratum;
			if(cost < lowestCost) lowestCost = cost;
		} else if(i < rangesSz_) {
			eliminatedStretch++;
		}
	}
	nextIncrCost_ = lowestCost;
	if(lowestCost > 0 && lowestCost != 0xffff) {
		// This branch's cost will change when curtailed; the
		// caller should re-insert it into the priority queue so
		// that the new cost takes effect.
		cost_ += lowestCost;
	} else if(lowestCost == 0xffff) {
		// This branch is totally exhausted; there are no more
		// valid outgoing paths from any positions within it.
		// Remove it from the PathManager and mark it as exhausted.
		// The caller should delete it.
		exhausted_ = true;
		if(ranges_ != NULL) {
			assert_gt(rangesSz_, 0);
			if(rpool.free(ranges_, rangesSz_)) {
				ranges_ = NULL;
				rangesSz_ = 0;
			}
		}
	} else {
		// Just mark it as curtailed and keep the same cost
	}
	if(ranges_ != NULL) {
		// Try to trim off no-longer-relevant elements of the
		// ranges_ array
		assert(!exhausted_);
		assert_gt(rangesSz_, 0);
		if(eliminatedStretch > 0) {
			// Decrement by one so that we can keep the last position
			// before.  We do this so that, if we attempt a deletion at
			// the following position, can still reach back and
			// determine the top/bot.
			eliminatedStretch--;
		}
		uint32_t trim = (rangesSz_ - len_ - 1) + eliminatedStretch;
		assert_leq(trim, rangesSz_);
		if(rpool.free(ranges_ + rangesSz_ - trim, trim)) {
			rangesSz_ -= trim;
			if(rangesSz_ == 0) {
				ranges_ = NULL;
			}
		}
	}
	curtailed_ = true;
}

/**
 * Set the elims to match the ranges in ranges_[len_], already
 * calculated by the caller.  We assume that at least one of the
 * ranges_ is non-empty.
 */
int Branch::installRanges(
	unsigned c,
	Visit& visit,
	int endoff,    // offset from one extreme end of the read (doesn't matter which)
	int fullqlen,  // full length of read (even if this is just a partial alignment)
	int qdepth,
	int qdep0,
	uint32_t qAllow,
	const uint8_t* qs,
	bool noInsert,
	bool overDelExtend)
{
	assert(active_);
	assert(!exhausted_);
	assert(ranges_ != NULL);
	assert_leq(c, 4);
	assert_geq(c, 0);
	Pos& r = ranges_[len_];
	r.flags.qchr = c;
	r.eliminated_ = true; // start with everything eliminated
	if((r.bots[0] == r.tops[0] &&
	    r.bots[1] == r.tops[1] &&
	    r.bots[2] == r.tops[2] &&
	    r.bots[3] == r.tops[3]) ||
	    rdepth_ + len_ < deps_[0])
	{
		// This position is eliminated, so we don't need to do any of
		// the bookkeeping below.  Note that we still have to set
		// r.flags.qchr (which we do above).
		return 0;
	}
	// Calculate our current offset from either side of the read; be
	// careful to accurately account for partial alignment state.
	int offFromNear = endoff;
	int offFromFar = fullqlen - endoff - 1;
	bool insExtend = (edit_.isInsert() && len_ == 0);
	bool delExtend = (edit_.isDelete() && len_ == 0) || overDelExtend;
	r.flags.ins_ext = insExtend;
	r.flags.del_ext = delExtend;
	int ret = 0;
	r.eliminate(); // set all elim flags to 1
	// Grab characters from last position for simple pruning of obvious
	// redundant paths.  We care only if gapped alignments are enabled
	// and if gaps are permitted at this query depth.  We disallow a
	// reference gap if both this query character and the previous are
	// the same.  We disallow a read gap if both this reference
	// character and the previous are the same.
	int prdc = -1; // read character from last position
	int prfc = -1; // ref character from last position
	int i = (int)deps_[0];
	i = max(0, i - rdepth_); // skip must-match positions
	if(gGaps && gAllowRedundant < 1 && visit.isJoinCell(qdepth, qdep0) && offFromNear > gGapBarrier+1) {
		if(len_ > i) {
			// Last pos was a match
			prdc = prfc = "ACGT"[ranges_[len_-1].flags.qchr];
		} else if(len_ == 0) {
			// Last pos was an edit
			prdc = edit_.qchr;
			prfc = edit_.chr;
		}
	}
	// if prdc and prfc are still -1, no outgoing paths will be pruned
	// for redundancy
	if(gVerbose) {
		cout << "installRanges(): " << len_ << ", [" << r.tops[0] << ", " << r.bots[0] << "], " << prdc << ", " << prfc << endl;
	}
	assert_lt(qs[0], INVALID_COST);
	assert_lt(qs[1], INVALID_COST);
	assert_lt(qs[2], INVALID_COST);
	assert_lt(qs[3], INVALID_COST);
	bool empty0 = !(r.bots[0] > r.tops[0]);
	bool empty1 = !(r.bots[1] > r.tops[1]);
	bool empty2 = !(r.bots[2] > r.tops[2]);
	bool empty3 = !(r.bots[3] > r.tops[3]);
	//
	// Set which mismatches are eliminated
	//
	// We can proceed on an A
	int depth = qdepth + 1;
	bool isJoinCell = visit.isJoinCell(depth, qdep0);
	if(c != 0 && !empty0 && qs[0] <= qAllow && !VISITED(r.tops[0], 0)) {
		r.flags.mmA = 0;
		ret++;
	}
	// We can proceed on a C
	if(c != 1 && !empty1 && qs[1] <= qAllow && !VISITED(r.tops[1], 0)) {
		r.flags.mmC = 0;
		ret++;
	}
	// We can proceed on a G
	if(c != 2 && !empty2 && qs[2] <= qAllow && !VISITED(r.tops[2], 0)) {
		r.flags.mmG = 0;
		ret++;
	}
	// We can proceed on a T
	if(c != 3 && !empty3 && qs[3] <= qAllow && !VISITED(r.tops[3], 0)) {
		r.flags.mmT = 0;
		ret++;
	}
	if(gGaps) {
		if(offFromNear > 0 &&  // no inserts before first char
		   offFromNear >= gGapBarrier && // respect near barrier
		   offFromFar  >= (gGapBarrier-1) && // respect far barrier
		   !noInsert) // override for the insert in middle of seed
		{
			//
			// Set which insertions are eliminated
			//
			depth = qdepth;
			isJoinCell = visit.isJoinCell(depth, qdep0);
			unsigned insCost = insExtend ? gInsExtend : gInsOpen;
			if(insCost <= qAllow) {
				if(!empty0 && (prdc == '-' || prfc != 'A') && !VISITED(r.tops[0], 1)) {
					r.flags.insA = 0;
					ret++;
				}
				if(!empty1 && (prdc == '-' || prfc != 'C') && !VISITED(r.tops[1], 1)) {
					r.flags.insC = 0;
					ret++;
				}
				if(!empty2 && (prdc == '-' || prfc != 'G') && !VISITED(r.tops[2], 1)) {
					r.flags.insG = 0;
					ret++;
				}
				if(!empty3 && (prdc == '-' || prfc != 'T') && !VISITED(r.tops[3], 1)) {
					r.flags.insT = 0;
					ret++;
				}
			}
		}
		bool nonRedundantDel = (prfc == '-' || "ACGTN"[c] != prdc);
		if(offFromNear >= gGapBarrier && // respect near barrier
		   offFromFar  >= gGapBarrier)   // respect far barrier
		{
			if(nonRedundantDel) {
				//
				// Set whether deletions are eliminated
				//
				depth = qdepth+1;
				isJoinCell = visit.isJoinCell(depth, qdep0);
				unsigned delCost = delExtend ? gDelExtend : gDelOpen;
				if(delCost <= qAllow) {
					uint32_t top;
					if(len_ > 0) {
						int pc = ranges_[len_-1].flags.qchr;
						assert_lt(pc, 4);
						top = ranges_[len_-1].tops[pc];
					} else {
						top = itop_;
					}
					if(!VISITED(top, -1)) {
						r.flags.del = 0;
						ret++;
					}
				}
			}
		}
	}
	if(ret > 0) r.eliminated_ = false;
	r.flags.costlo = r.flags.costlo2 = INVALID_COST;
	// Copy quals
	r.flags.qualA = qs[0];
	r.flags.qualC = qs[1];
	r.flags.qualG = qs[2];
	r.flags.qualT = qs[3];
	if(!r.eliminated_) {
		// Now that the quals are set and the elim flags are set,
		// determine best and second-best quals
		r.updateLo();
		assert_lt(r.flags.costlo, INVALID_COST);
		assert(ret == 1 || r.flags.costlo2 < INVALID_COST);
	}
	assert_leq(r.flags.costlo, r.flags.costlo2);
	assert(r.repOk());
	return ret;
}

/**
 * Pretty-print the Edits that lead to this branch.
 */
void Branch::printEdits(uint32_t qlen, bool ebwtFw) {
	if(!edit_.initialized()) return;
	cout << "Edit: ";
	Branch *cur = this;
	Edit e = cur->edit_;
	while(e.initialized()) {
		cout << (ebwtFw ? (qlen - e.pos - 1) : e.pos)
			 << ':' << (char)e.chr << ">" << (char)e.qchr;
		cur = cur->parent_;
		e = cur->edit_;
		if(e.initialized()) cout << " ";
	}
	cout << endl;
}

/**
 * Pretty-print the state of this branch.
 */
void Branch::print(const BTDnaString& qry,
                   uint16_t minCost,
                   ostream& out,
                   bool halfAndHalf,
                   int seedLen,
                   bool fw,
                   bool ebwtFw)
{
	size_t editidx = 0;
	const size_t qlen = qry.length();
	if(exhausted_)      out << "E ";
	else if(curtailed_) out << "C ";
	else                out << "  ";
	out << (ebwtFw ? "<" : ">");
	out << (fw ? "F " : "R ");
	stringstream ss;
	ss << cost_;
	printPadded(out, ss.str(), 6);
	out << ' ';
	stringstream ss2;
	ss2 << minCost;
	printPadded(out, ss2.str(), 6);
	if(halfAndHalf)      out << " h ";
	else if(seedLen > 0) out << " s ";
	else                 out << "   ";
	stringstream ss3;
	EList<Edit> edits;
	{
		Branch *cur = this;
		Edit e = cur->edit_;
		while(e.initialized()) {
			edits.push_back(e);
			cur = cur->parent_;
			e = cur->edit_;
		}
	}
	edits.reverse();
	if(edits.size() > 1) {
		assert_geq(edits[1].pos, edits[0].pos);
	}
	const size_t numEdits = edits.size();
	for(size_t i = 0; i < qlen; i++) {
		if(seedLen > 0 && (int)i == seedLen) {
			ss3 << "|";
		} else {
			ss3 << " ";
		}
		if(i > rdepth_ + len_) {
			// As-yet-unvisited
			ss3 << " =";
		} else if(i == rdepth_ + len_) {
			// May be inserts
			bool first = true;
			while(editidx < numEdits && edits[editidx].pos == i) {
				if(!first) ss3 << " ";
				assert(edits[editidx].isInsert());
				ss3 << "+" << (char)tolower(edits[editidx].chr);
				editidx++;
				first = false;
			}
		} else {
			bool first = true;
			bool ins = false, mm = false, del = false;
			// print edits
			while(editidx < numEdits && edits[editidx].pos == i) {
				if(!first) ss3 << " ";
				if(edits[editidx].isInsert()) ss3 << "+";
				else ss3 << " ";
				ss3 << (char)tolower(edits[editidx].chr);
				first = false;
				if(edits[editidx].isInsert()) ins = true;
				if(edits[editidx].isDelete()) del = true;
				if(edits[editidx].isMismatch()) {
					assert(!mm);
					mm = true;
				}
				editidx++;
			}
			assert(!ins || !del);
			assert(!mm || !del);
			// print reference char
			if(!mm && !del) {
				if(!first) ss3 << " ";
				ss3 << " " << (char)qry.toChar(qlen - i - 1);
			}
		}
	}
	assert_eq(editidx, edits.size());
	string s = ss3.str();
	if(ebwtFw) reverse(s.begin(), s.end());
	out << s << endl;
}
