/*
 * range_source.h
 */

#ifndef RANGE_SOURCE_H_
#define RANGE_SOURCE_H_

#include <stdint.h>
#include <vector>
#include <queue>
#include "ebwt.h"
#include "range.h"
#include "pool.h"
#include "edit.h"
#include "sstring.h"

enum AdvanceUntil {
	ADV_FOUND_RANGE = 1,
	ADV_COST_CHANGES,
	ADV_STEP
};

// 255 is max cost
#define INVALID_COST ((1<<8)-1)

/**
 * Helper class that relates the query and reference "depth" coordinate
 * systems.  This is helpful if the query has been transformed by some
 * edits that were determined in a partial alignment.  This is
 * currently filled in by EbwtRangeSource.setQuery(), though the
 * initialization should probably take place within this class.
 */
struct CoordMap {
	void clear() { r2q.clear(); q2r.clear(); }
	SStringExpandable<int> r2q; // ref depth to query depth
	SStringExpandable<int> q2r; // query depth to ref depth
};

/**
 * Class
 */
class Visit {
public:
	Visit() : visit_(32), loGapDepth_(0), hiGapDepth_(0), nGaps_(0), reset_(false) { }

	/**
	 * A "join cell" is a cell where it's possible for two separate
	 * alignments, where both alignments have the same characters on
	 * the reference side, to merge.  Such cells are only possible if
	 * gaps are permitted.  If the user has requested that redundant
	 * alignments be treated as distinct, valid alignments, we always
	 * return false.
	 */
	bool isJoinCell(int qdepth, int qdep0) const {
		assert(reset_);
		return qdepth >= qdep0 &&
		       nGaps_ > 0 &&
		       qdepth >= loGapDepth_ &&
		       qdepth <  hiGapDepth_;
	}

	/**
	 * Marks a <join-cell, top> pair as having been visited by
	 * inserting 'top' into a set.
	 */
	void setVisited(int qdepth, int qdep0, int gapOff, uint32_t top) {
		assert(reset_);
		assert(!isVisited(qdepth, qdep0, gapOff, top));
		setAndTestVisited(qdepth, qdep0, gapOff, top);
	}

	/**
	 * If the given <join-cell, top> pair has already been visited,
	 * this function returns false; otherwise, it marks the pair as
	 * visited and returns true.
	 */
	bool setAndTestVisited(int qdepth, int qdep0, int gapOff, uint32_t top) {
		assert(reset_);
		assert(isJoinCell(qdepth, qdep0));
		assert_leq(abs(gapOff), nGaps_);
		int qidx = qdepth - loGapDepth_;
		assert_geq(qidx, 0);
		assert_lt(qidx, (int)visit_.size());
		if(gapOff < 0) gapOff += visit_[qidx].size();
		assert_geq(gapOff, 0);
		assert_lt(gapOff, (int)visit_[qidx].size());
		return visit_[qidx][gapOff].insert(top);
	}

	/**
	 * Return true iff the given <join-cell, top> pair has already been
	 * visited.
	 */
	bool isVisited(int qdepth, int qdep0, int gapOff, uint32_t top) const {
		assert(reset_);
		assert(isJoinCell(qdepth, qdep0));
		assert_leq(abs(gapOff), nGaps_);
		int qidx = qdepth - loGapDepth_;
		assert_geq(qidx, 0);
		assert_lt(qidx, (int)visit_.size());
		if(gapOff < 0) gapOff += visit_[qidx].size();
		assert_geq(gapOff, 0);
		assert_lt(gapOff, (int)visit_[qidx].size());
		return visit_[qidx][gapOff].contains(top);
	}

	/**
	 * Reset the Visited data structure for a new search.
	 *
	 * loGapDep: alignment must be gapless as long as qdepth < loGapDep
	 * hiGapDep: alignment must be gapless as long as qdepth >= hiGapDep
	 * nGaps: largest possible number of total gaps in alignment,
	 * calculated by dividing the budget by the lowest of the 4 types
	 * of gap penalty (open/extend insert/delete)
	 */
	void reset(size_t loGapDep, size_t hiGapDep, int nGaps) {
		reset_ = true;
		loGapDepth_ = (uint32_t)loGapDep;
		hiGapDepth_ = (uint32_t)hiGapDep;
		nGaps_ = 0;
		if(!gGaps || gAllowRedundant >= 3) {
			// Leaving nGaps_ at 0 disables visit checking
			return;
		}
		nGaps_ = nGaps;
		if(nGaps > 0) {
			size_t mostref = nGaps*2 + 1;
			// Initialize the matrix of visited BW ranges before
			// resizing it and clear()ing all elements.
			visit_.resize(hiGapDepth_-loGapDepth_);
			for(int i = 0; i < (hiGapDepth_-loGapDepth_); i++) {
				visit_[i].resize(mostref);
				for(size_t j = 0; j < mostref; j++) {
					visit_[i][j].clear();
				}
			}
		} else {
			visit_.resize(0);
		}
	}

private:

	EList<EList<ESet<uint32_t> > > visit_;
	// lo/hiGapDepth_ are relatively tight bounds on the range of
	// BWT-SW *rows* where gaps (i.e. horizontal and vertical
	// movements) are permitted.
	int loGapDepth_;
	int hiGapDepth_;
	int nGaps_;
	bool reset_; // true -> reset() at least once already
};

/**
 * All per-position state, including the ranges calculated for each
 * character, the quality value at the position, and a set of flags
 * recording whether we've tried each way of proceeding from this
 * position.
 */
struct Pos {

	/**
	 * Pick an as-yet-uneliminated Edit by simply scanning each edit
	 * type and picking the first uneliminated one.  We'd generally
	 * rather do something random if there is more than one choice, but
	 * this is safe as long as there is exactly one.
	 */
	int pickFirstEdit(Edit& ret) {
		bool color = false;
		int chr = -1;
		ret.type = EDIT_TYPE_MM;
		unsigned insCost = flags.ins_ext ? gInsExtend : gInsOpen;
		unsigned delCost = flags.del_ext ? gDelExtend : gDelOpen;
		if       (flags.qualA == flags.costlo && flags.mmA == 0) {
			flags.mmA = 1;
			chr = 0;
		} else if(flags.qualC == flags.costlo && flags.mmC == 0) {
			flags.mmC = 1;
			chr = 1;
		} else if(flags.qualG == flags.costlo && flags.mmG == 0) {
			flags.mmG = 1;
			chr = 2;
		} else if(flags.qualT == flags.costlo && flags.mmT == 0) {
			flags.mmT = 1;
			chr = 3;
		} else if(insCost == flags.costlo && flags.insA == 0) {
			assert(gGaps);
			flags.insA = 1;
			ret.type = EDIT_TYPE_INS;
			chr = 0;
		} else if(insCost == flags.costlo && flags.insC == 0) {
			assert(gGaps);
			flags.insC = 1;
			ret.type = EDIT_TYPE_INS;
			chr = 1;
		} else if(insCost == flags.costlo && flags.insG == 0) {
			assert(gGaps);
			flags.insG = 1;
			ret.type = EDIT_TYPE_INS;
			chr = 2;
		} else if(insCost == flags.costlo && flags.insT == 0) {
			assert(gGaps);
			flags.insT = 1;
			ret.type = EDIT_TYPE_INS;
			chr = 3;
		} else {
			assert(gGaps);
			assert_eq(0, flags.del);
			assert_eq(delCost, flags.costlo);
			delCost += 0;
			flags.del = 1;
			ret.type = EDIT_TYPE_DEL;
		}
		if(chr == -1) {
			ret.chr = '-';
		} else {
			ret.chr = color ? "0123"[chr] : "ACGT"[chr];
		}
		if(ret.type == EDIT_TYPE_INS) {
			ret.qchr = '-';
		} else {
			ret.qchr = "ACGTN"[flags.qchr];
		}
		return chr; // chr is -1 iff deletion was chosen
	}

	/**
	 * Eliminate mismatch path for DNA char 'i'.
	 */
	void elimMm(int i) {
		switch(i) {
			case 0: flags.mmA = 1; break;
			case 1: flags.mmC = 1; break;
			case 2: flags.mmG = 1; break;
			case 3: flags.mmT = 1; break;
			default: throw 1; break;
		}
	}

	/**
	 * Eliminate insertion path for DNA char 'i'.
	 */
	void elimIns(int i) {
		switch(i) {
			case 0: flags.insA = 1; break;
			case 1: flags.insC = 1; break;
			case 2: flags.insG = 1; break;
			case 3: flags.insT = 1; break;
			default: throw 1; break;
		}
	}

	/**
	 * If elimflag == false, then update costlo and costlo2 w/r/t the
	 * cost c.
	 */
	void update(bool elimflag, unsigned c) {
		if(!elimflag) {
			if(c < flags.costlo) {
				flags.costlo2 = flags.costlo;
				flags.costlo = c;
			} else if(c == flags.costlo || c < flags.costlo2) {
				flags.costlo2 = c;
			}
		}
	}

	/**
	 * Assuming flags.qual* and all of the elimination flags are
	 * already set, set costlo and costlo2 to the additional cost
	 * incurred by the least and second-least costly as-yet-
	 * uneliminated paths.  If there is a tie, costlo will equal
	 * costlo2.
	 */
	void updateLo() {
		assert_lt(gInsExtend, INVALID_COST);
		assert_lt(gInsOpen,   INVALID_COST);
		assert_lt(gDelExtend, INVALID_COST);
		assert_lt(gDelOpen,   INVALID_COST);
		flags.costlo = flags.costlo2 = INVALID_COST;
		// Take mismatches into account
		update(flags.mmA, flags.qualA);
		update(flags.mmC, flags.qualC);
		update(flags.mmG, flags.qualG);
		update(flags.mmT, flags.qualT);
		if(gGaps) {
			// Take indels into account; if indels aren't enabled, then
			// flags.ins* and flags.del should be 1
			if(flags.ins_ext && gInsExtend < flags.costlo2) {
				// Gap extension; gaps in subject
				update(flags.insA, gInsExtend);
				update(flags.insC, gInsExtend);
				update(flags.insG, gInsExtend);
				update(flags.insT, gInsExtend);
			} else if(!flags.ins_ext && gInsOpen < flags.costlo2) {
				// Gap open; gaps in subject
				update(flags.insA, gInsOpen);
				update(flags.insC, gInsOpen);
				update(flags.insG, gInsOpen);
				update(flags.insT, gInsOpen);
			}
			if(flags.del_ext && gDelExtend < flags.costlo2) {
				update(flags.del, gDelExtend);
			} else if(!flags.del_ext && gDelOpen < flags.costlo2) {
				update(flags.del, gDelOpen);
			}
		} else {
			// gGaps is false, so these should all have been eliminated
			assert(flags.insA);
			assert(flags.insC);
			assert(flags.insG);
			assert(flags.insT);
			assert(flags.del);
		}
		assert(repOk());
	}

	/**
	 * Set all 9 elimination bits of the flags field to 1, indicating
	 * that all outgoing paths are eliminated.
	 */
	void eliminate() {
		flags.mmA = flags.mmC = flags.mmG = flags.mmT =
		flags.insA = flags.insC = flags.insG = flags.insT =
		flags.del = 1;
	}

	/**
	 * If elimflag == false, then update lo and lo2 w/r/t the cost c.
	 */
	void update2(bool elimflag, unsigned c, unsigned& lo, unsigned& lo2) const {
		if(!elimflag) {
			if(c < lo) { lo2 = lo; lo = c; }
			else if(c == lo || c < lo2) lo2 = c;
		}
	}

	/**
	 * Internal consistency check.  Basically just checks that lo and
	 * lo2 are set correctly.
	 */
	bool repOk() const {
		if(eliminated_) {
			return true;
		}
		// Uneliminated chars must have non-empty ranges
		if(!flags.mmA || !flags.insA) assert_gt(bots[0], tops[0]);
		if(!flags.mmC || !flags.insC) assert_gt(bots[1], tops[1]);
		if(!flags.mmG || !flags.insG) assert_gt(bots[2], tops[2]);
		if(!flags.mmT || !flags.insT) assert_gt(bots[3], tops[3]);
		if(!flags.del) {
			assert(bots[0] > tops[0] || bots[1] > tops[1] ||
			       bots[2] > tops[2] || bots[3] > tops[3]);
		}
		unsigned lo = INVALID_COST;
		unsigned lo2 = INVALID_COST;
		assert_lt(flags.qualA, INVALID_COST);
		assert_lt(flags.qualC, INVALID_COST);
		assert_lt(flags.qualG, INVALID_COST);
		assert_lt(flags.qualT, INVALID_COST);
		// Take mismatches into account
		update2(flags.mmA, flags.qualA, lo, lo2);
		update2(flags.mmC, flags.qualC, lo, lo2);
		update2(flags.mmG, flags.qualG, lo, lo2);
		update2(flags.mmT, flags.qualT, lo, lo2);
		// Take indels into account
		if(flags.ins_ext && gInsExtend <= flags.costlo2) {
			// Gap extension; gaps in subject
			update2(flags.insA, gInsExtend, lo, lo2);
			update2(flags.insC, gInsExtend, lo, lo2);
			update2(flags.insG, gInsExtend, lo, lo2);
			update2(flags.insT, gInsExtend, lo, lo2);
		} else if(!flags.ins_ext && gInsOpen <= flags.costlo2) {
			// Gap open; gaps in subject
			update2(flags.insA, gInsOpen, lo, lo2);
			update2(flags.insC, gInsOpen, lo, lo2);
			update2(flags.insG, gInsOpen, lo, lo2);
			update2(flags.insT, gInsOpen, lo, lo2);
		}
		if(flags.del_ext && gDelExtend <= flags.costlo2) {
			update2(flags.del, gDelExtend, lo, lo2);
		} else if(!flags.del_ext && gDelOpen <= flags.costlo2) {
			update2(flags.del, gDelOpen, lo, lo2);
		}
		assert_eq(lo, flags.costlo);
		assert_eq(lo2, flags.costlo2);
		return true;
	}

	// Outgoing ranges; if the position being described is not a
	// legitimate jumping-off point for a branch, tops[] and bots[]
	// will be filled with 0s and all possibilities in eq will be
	// eliminated
	uint32_t tops[4]; // A, C, G, T top offsets
	uint32_t bots[4]; // A, C, G, T bot offsets
	bool eliminated_;  // Whether all outgoing paths have been eliminated

	struct {
		uint64_t mmA      : 1; // A in ref aligns to non-A char in read
		uint64_t mmC      : 1; // C in ref aligns to non-C char in read
		uint64_t mmG      : 1; // G in ref aligns to non-G char in read
		uint64_t mmT      : 1; // T in ref aligns to non-T char in read
		uint64_t insA     : 1; // A insertion in reference w/r/t read
		uint64_t insC     : 1; // C insertion in reference w/r/t read
		uint64_t insG     : 1; // G insertion in reference w/r/t read
		uint64_t insT     : 1; // T insertion in reference w/r/t read
		uint64_t del      : 1; // deletion of read character
		uint64_t qualA    : 8; // quality penalty for picking A at this position
		uint64_t qualC    : 8; // quality penalty for picking C at this position
		uint64_t qualG    : 8; // quality penalty for picking G at this position
		uint64_t qualT    : 8; // quality penalty for picking T at this position
		uint64_t costlo   : 8; // lowest quality penalty at this position
		uint64_t costlo2  : 8; // 2nd-lowest quality penalty at this position
		uint64_t qchr     : 3; // query char at this position (invalid for del path)
		uint64_t ins_ext  : 1; // whether an insertion at this pos is an extension
		uint64_t del_ext  : 1; // whether a deletion at this pos is an extension
		uint64_t reserved : 2;
	} flags;
};

/**
 * Encapsulates a "branch" of the search space; i.e. all of the
 * information deduced by walking along a path with only matches, along
 * with information about the decisions that lead to the root of that
 * path.
 */
class Branch {
	typedef std::pair<uint32_t, uint32_t> U32Pair;
public:
	Branch() :
		delayedCost_(0), curtailed_(false), exhausted_(false),
		prepped_(false), delayedIncrease_(false)
	{
		ASSERT_ONLY(active_ = false);
		ASSERT_ONLY(incarnation_ = 0);
	}

	Branch* splitBranch(
		AllocOnlyPool<Pos>& rpool,
		AllocOnlyPool<Branch>& bpool,
		Visit& visit,
		const CoordMap& cmap,
		RandomSource& rand,
		uint32_t qlen,
		int qdep0,
		uint32_t qualLim,
		int seedLen,
		bool& triedGap,
		const EbwtParams& ep,
		const uint8_t* ebwt);

	void curtail(AllocOnlyPool<Pos>& rpool, int seedLen);

	int installRanges(
		unsigned c,
		Visit& visit,
		int endoff,
		int fullqlen,
		int qdepth,
		int qdep0,
		uint32_t qAllow,
		const uint8_t* qs,
		bool noInsert,
		bool overDelExtend);

	bool eliminateVisited(size_t i, int qdepth, int qdep0, Visit& visited);

	/**
	 * Initialize a new branch object with an empty path.
	 */
	bool init(AllocOnlyPool<Pos>& pool,
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
	          const uint8_t* ebwt);

	/**
	 * Depth of the deepest tip of the branch.
	 */
	int tipDepth() const { return (int)(rdepth_ + len_); }

	/**
	 * Return true iff all outgoing edges from position i have been
	 * eliminated.
	 */
	inline bool eliminated(int i) const {
		assert_geq(i, max(0, deps_[0] - rdepth_));
		if(i <= len_ && i < rangesSz_) {
			assert(ranges_ != NULL);
#ifndef NDEBUG
			if(!ranges_[i].eliminated_) {
				// Someone must be as-yet-uneliminated
				assert(!ranges_[i].flags.mmA ||
				       !ranges_[i].flags.mmC ||
				       !ranges_[i].flags.mmG ||
				       !ranges_[i].flags.mmT ||
				       !ranges_[i].flags.insA ||
				       !ranges_[i].flags.insC ||
				       !ranges_[i].flags.insG ||
				       !ranges_[i].flags.insT ||
				       !ranges_[i].flags.del);
				assert_lt(ranges_[i].flags.costlo, INVALID_COST);
			}
#endif
			return ranges_[i].eliminated_;
		}
		return true;
	}

	/**
	 * Return true iff we've just started to extend a Branch that split
	 * off from a deletion.
	 */
	bool justDeleted() const {
		return len_ == 0 && edit_.isDelete();
	}

	/**
	 * Using randomness when picking from among multiple choices, pick
	 * an edit to make.  Called from Branch.splitBranch().
	 */
	Edit pickEdit(int pos, RandomSource& rand, uint32_t& top, uint32_t& bot);

	/**
	 * Free a branch and all its contents.
	 */
	void free(AllocOnlyPool<Pos>& rpool, AllocOnlyPool<Branch>& bpool) {
		if(ranges_ != NULL) {
			assert_gt(rangesSz_, 0);
			rpool.free(ranges_, rangesSz_);
			ranges_ = NULL;
			rangesSz_ = 0;
		}
		if(parent_ != NULL) {
			assert_gt(parent_->numChildren_, 0);
			parent_->numChildren_--;
		}
		// TODO: This isn't great; we're simply not going to de-
		// allocate slots for Branches that still have children because
		// it's too hard to ensure the parent lives on the queue long
		// enough.  So we'll let the caller remove it from the queue
		// but we won't free its memory.
		if(numChildren_ == 0) {
			// Only place where bpool.free() is called.  Note that
			// bpool.free() zeroes out our memory.
			ASSERT_ONLY(uint32_t incarnation = incarnation_);
			bpool.free(this);
			ASSERT_ONLY(incarnation_ = incarnation+1);
			ASSERT_ONLY(active_ = false);
		} else {
			ASSERT_ONLY(active_ = false);
		}
		// Restore incarnation counter
	}

	/**
	 * Pretty-print the Edits that lead to this branch.
	 */
	void printEdits(uint32_t qlen, bool ebwtFw);

	/**
	 * Pretty-print the state of this branch.
	 */
	void print(const BTDnaString& qry,
	           uint16_t minCost,
	           std::ostream& out,
	           bool halfAndHalf,
	           int seedLen,
	           bool fw,
	           bool ebwtFw);

	/**
	 * Prep this branch for the next extension by calculating the
	 * SideLocus information and prefetching cache lines from the
	 * appropriate loci.
	 */
	void prep(const EbwtParams& ep, const uint8_t* ebwt) {
		assert(!prepped_);
		if(prepped_) return;
		if(bot_ > top_+1) {
			SideLocus::initFromTopBot(top_, bot_, ep, ebwt, ltop_, lbot_);
		} else if(bot_ > top_) {
			ltop_.initFromRow(top_, ep, ebwt);
			lbot_.invalidate();
		}
		prepped_ = true;
	}

	/**
	 * Get the furthest-out Pos.
	 */
	Pos* backPos() {
		assert(!exhausted_);
		assert(ranges_ != NULL);
		assert_lt(len_, rangesSz_);
		return &ranges_[len_];
	}

	/**
	 * Extend this branch by one position.
	 */
	void extend() {
		assert(active_);
		assert(!exhausted_);
		assert(!curtailed_);
		assert(ranges_ != NULL);
		assert(repOk());
#ifndef NDEBUG
		int chouts = 0;
		for(int i = 0; i < 4; i++) {
			if(ranges_[len_].bots[i] > ranges_[len_].tops[i]) chouts++;
		}
		assert_gt(chouts, 0);
#endif
		prepped_ = false;
		len_++;
	}

	/**
	 * Undo a recent extension.
	 */
	void retract() {
		assert_gt(len_ ,0);
		prepped_ = true;
		len_--;
	}

	/**
	 * Do an internal consistency check
	 */
	bool repOk(uint32_t qlen = 0) const{
		assert_leq(deps_[0], deps_[1]);
		assert_leq(deps_[1], deps_[2]);
		assert_leq(deps_[2], deps_[3]);
		// Make sure all the parent relationships look sane
		const Branch * cur = this;
		while(cur->parent_ != NULL) {
			Branch *p = cur->parent_;
			assert_eq(p->id_, cur->parentId_);
			if(p->edit_.initialized()) {
				assert_geq(cur->edit_.pos, p->edit_.pos);
				assert_gt(cur->numEdits_, 0);
			}
			if(!cur->edit_.initialized()) {
				assert_eq(cur->numEdits_, 0);
			}
			assert_eq(cur->numEdits_-1, p->numEdits_);
			cur = p;
		}
		// cur should now equal root
		assert(cur->parent_ == NULL); // root doesn't have parent
		assert_eq(0, cur->parentId_);
		assert(!cur->edit_.initialized()); // root doesn't have edit at its base
		if(qlen > 0) {
			assert_leq(rdepth_, qlen);
		}
		if(ranges_ != NULL) {
			int i = max((int)deps_[0] - rdepth_, 0);
			for(; i < min(len_, rangesSz_); i++) {
				assert(ranges_[i].repOk());
				if(i < len_) {
					assert_neq(4, ranges_[i].flags.qchr);
				}
			}
		}
		// Make sure individual eliminated() calls match the exhausted_
		// flag.
		if(exhausted_) {
			int i = (int)deps_[0];
			i = max(0, i - rdepth_); // skip must-match positions
			for(; i <= len_; i++) assert(eliminated(i));
		} else if(curtailed_) {
			int i = (int)deps_[0];
			i = max(0, i - rdepth_); // skip must-match positions
			bool elim = true;
			for(; i <= len_; i++) {
				if(!eliminated(i)) {
					elim = false;
					break;
				}
			}
			assert(!elim);
		}
		return true;
	}

	uint32_t id_;     // branch id; needed to make the ordering of
	                  // branches that are tied in the priority queue
	                  // totally unambiguous.  Otherwise, things start
	                  // getting non-deterministic.
	uint16_t deps_[4];// depths at which max # edits before are 0, 1, 2, 3
	uint16_t rdepth_; // offset in read space from root of search space
	int16_t  gapOff_; // additive offset due to gaps; = # inserts - # deletions
	uint16_t len_;    // length of the branch
	uint16_t ilen_;   // initial setting for len_
	uint16_t cost_;   // top 2 bits = stratum, bottom 14 = qual ham
	                  // it's up to Branch to keep this updated with the
	                  // cumulative cost of the best path leaving the
	                  // branch; if the branch hasn't been fully
	                  // extended yet, then that path will always be the
	                  // one that extends it by one more
	uint16_t ham_;    // quality-weighted hamming distance so far
	Pos *ranges_; // Allocated from the RangeStatePool
	uint16_t rangesSz_;// number of RangeStates allocated for ranges_
	uint32_t itop_;   // top of range at the root of the branch
	uint32_t ibot_;   // bot of range at the root of the branch
	uint32_t top_;    // last-calculated top offset (using ltop_)
	uint32_t bot_;    // last-calculated bot offset (using lbot_)
	SideLocus ltop_;  // SideLocus for calculating tops
	SideLocus lbot_;  // SideLocus for calculating bots
	Edit edit_;       // the edit at the root of this branch
	uint16_t numEdits_; // number of edits prior to this branch
	Branch *parent_;  // parent of this Branch
	uint32_t numChildren_; // number of Branches that split off from this one

	uint16_t delayedCost_;

	bool curtailed_;  // can't be extended anymore without using edits
	bool exhausted_;  // all outgoing edges exhausted, including all edits
	bool prepped_;    // whether SideLocus's are inited
	ASSERT_ONLY(uint32_t parentId_); // id of parent (for sanity-checking)
	ASSERT_ONLY(bool active_); // true -> Branch is initialized and in the mix
	ASSERT_ONLY(uint32_t incarnation_); // # times this bit of memory has been initialized
	ASSERT_ONLY(uint32_t splits_);      // # times a branch split off from this one
	uint32_t nextIncrCost_;
	bool delayedIncrease_;
};

/**
 * Order two Branches based on cost.
 */
class CostCompare {
public:
	/**
	 * true -> b before a
	 * false -> a before b
	 */
	bool operator()(const Branch* a, const Branch* b) const {
		bool aUnextendable = a->curtailed_ || a->exhausted_;
		bool bUnextendable = b->curtailed_ || b->exhausted_;
		// Branch with the best cost
		if(a->cost_ == b->cost_) {
			// If one or the other is curtailed, take the one that's
			// still getting extended
			if(bUnextendable && !aUnextendable) {
				// a still being extended, return false
				return false;
			}
			if(aUnextendable && !bUnextendable) {
				// b still being extended, return true
				return true;
			}
			// Either both are curtailed or both are still being
			// extended, pick based on which one is deeper
			if(a->tipDepth() != b->tipDepth()) {
				// Expression is true if b is deeper
				return a->tipDepth() < b->tipDepth();
			}
			// Keep things deterministic by providing an unambiguous
			// order using the id_ field
			assert_neq(b->id_, a->id_);
			return b->id_ < a->id_;
		} else {
			return b->cost_ < a->cost_;
		}
	}

	static bool equal(const Branch* a, const Branch* b) {
		return a->cost_ == b->cost_ && a->curtailed_ == b->curtailed_ && a->tipDepth() == b->tipDepth();
	}
};

/**
 * A priority queue for Branch objects; makes it easy to process
 * branches in a best-first manner by prioritizing branches with lower
 * cumulative costs over branches with higher cumulative costs.
 */
class BranchQueue {

	typedef std::pair<int, int> TIntPair;
	typedef std::priority_queue<Branch*, vector<Branch*>, CostCompare> TBranchQueue;

public:

	BranchQueue() : sz_(0), branchQ_(), patid_(0) { }

	/**
	 * Return the front (highest-priority) element of the queue.
	 */
	Branch *front(bool verbose) {
		Branch *b = branchQ_.top();
		if(gVerbose && verbose) {
			stringstream ss;
			ss << patid_ << ": Fronting " << b->id_ << ", " << b->cost_ << ", " << b->exhausted_ << ", " << b->curtailed_ << ", " << sz_ << "->" << sz_;
			glog.msg(ss.str());
		}
		return b;
	}

	/**
	 * Remove and return the front (highest-priority) element of the
	 * queue.
	 */
	Branch *dequeue() {
		Branch *b = branchQ_.top(); // get it
		branchQ_.pop(); // remove it
		if(gVerbose) {
			stringstream ss;
			ss << patid_ << ": DeQing " << b->id_ << ", " << b->cost_ << ", " << b->exhausted_ << ", " << b->curtailed_ << ", " << sz_ << "->" << (sz_-1);
			glog.msg(ss.str());
		}
		sz_--;
		return b;
	}

	/**
	 * Insert a new Branch into the sorted priority queue.
	 */
	void enqueue(Branch *b) {
#ifndef NDEBUG
		bool bIsBetter = empty() || !CostCompare()(b, branchQ_.top());
#endif
		if(gVerbose) {
			stringstream ss;
			ss << patid_ << ": EnQing " << b->id_ << ", " << b->cost_ << ", " << b->exhausted_ << ", " << b->curtailed_ << ", " << sz_ << "->" << (sz_+1);
			glog.msg(ss.str());
		}
		branchQ_.push(b);
#ifndef NDEBUG
		assert(bIsBetter  || branchQ_.top() != b || CostCompare::equal(branchQ_.top(), b));
		assert(!bIsBetter || branchQ_.top() == b || CostCompare::equal(branchQ_.top(), b));
#endif
		sz_++;
	}

	/**
	 * Empty the priority queue and reset the count.
	 */
	void reset(TReadId patid) {
		patid_ = patid;
		branchQ_ = TBranchQueue();
		sz_ = 0;
	}

	/**
	 * Return true iff the priority queue of branches is empty.
	 */
	bool empty() const {
		bool ret = branchQ_.empty();
		assert(ret || sz_ > 0);
		assert(!ret || sz_ == 0);
		return ret;
	}

	/**
	 * Return the number of Branches in the queue.
	 */
	uint32_t size() const {
		return sz_;
	}

#ifndef NDEBUG
	/**
	 * Consistency check.
	 */
	bool repOk(std::set<Branch*>& bset) {
		TIntPair pair = bestStratumAndHam(bset);
		Branch *b = branchQ_.top();
		assert_eq(pair.first, (b->cost_ >> 14));
		assert_eq(pair.second, (b->cost_ & ~0xc000));
		return true;
	}
#endif

protected:

#ifndef NDEBUG
	/**
	 * Return the stratum and quality-weight (sum of qualities of all
	 * edited positions) of the lowest-cost branch.
	 */
	TIntPair bestStratumAndHam(std::set<Branch*>& bset) const {
		TIntPair ret = make_pair(0xffff, 0xffff);
		std::set<Branch*>::iterator it;
		for(it = bset.begin(); it != bset.end(); it++) {
			Branch *b = *it;
			int stratum = b->cost_ >> 14;
			assert_lt(stratum, 4);
			int qual = b->cost_ & ~0xc000;
			if(stratum < ret.first ||
			   (stratum == ret.first && qual < ret.second))
			{
				ret.first = stratum;
				ret.second = qual;
			}
		}
		return ret;
	}
#endif

	uint32_t sz_;
	TBranchQueue branchQ_; // priority queue of branches
	TReadId patid_;
};

/**
 * A class that both contains Branches and determines how those
 * branches are extended to form longer paths.  The overall goal is to
 * find the best full alignment(s) as quickly as possible so that a
 * successful search can finish quickly.  A second goal is to ensure
 * that the most "promising" paths are tried first so that, if there is
 * a limit on the amount of effort spent searching before we give up,
 * we can be as sensitive as possible given that limit.
 *
 * The quality (or cost) of a given alignment path will ultimately be
 * configurable.  The default cost model is:
 *
 * 1. Mismatches incur one "edit" penalty and a "quality" penalty with
 *    magnitude equal to the Phred quality score of the read position
 *    involved in the edit (note that insertions into the read are a
 *    little trickier).
 * 2. Edit penalties are all more costly than any quality penalty; i.e.
 *    the policy sorts alignments first by edit penalty, then by
 *    quality penalty.
 * 3. For the Maq-like alignment policy, edit penalties saturate (don't
 *    get any greater) after leaving the seed region of the alignment.
 */
class PathManager {

public:

	PathManager(ChunkPool* cp, int *btCnt) :
		branchQ_(),
		cpool(cp),
		bpool(cpool, "branch"),
		rpool(cpool, "range state"),
		minCost(0),
		triedGap_(false),
		btCnt_(btCnt)
	{ }

	~PathManager() { }

	/**
	 * Return the "front" (highest-priority) branch in the collection.
	 */
	Branch* front(bool verbose) {
		assert(!empty());
		assert_gt(branchQ_.front(false)->deps_[3], 0);
		return branchQ_.front(verbose);
	}

	/**
	 * Pop the highest-priority (lowest cost) element from the
	 * priority queue.
	 */
	Branch* dequeue() {
		Branch* b = branchQ_.dequeue();
		assert_gt(b->deps_[3], 0);
#ifndef NDEBUG
		// Also remove it from the set
		assert(branchSet_.find(b) != branchSet_.end());
		ASSERT_ONLY(size_t setSz = branchSet_.size());
		branchSet_.erase(branchSet_.find(b));
		assert_eq(setSz-1, branchSet_.size());
		if(!branchQ_.empty()) {
			// Top shouldn't be b any more
			Branch *newtop = branchQ_.front(false);
			assert(b != newtop);
		}
#endif
		// Update this PathManager's cost
		minCost = branchQ_.front(false)->cost_;
		assert(repOk());
		return b;
	}

	/**
	 * Push a new element onto the priority queue.
	 */
	void enqueue(Branch *b) {
		assert(!b->exhausted_);
		assert_gt(b->deps_[3], 0);
		branchQ_.enqueue(b);
#ifndef NDEBUG
		// Also insert it into the set
		assert(branchSet_.find(b) == branchSet_.end());
		branchSet_.insert(b);
#endif
		// Update this PathManager's cost
		minCost = branchQ_.front(false)->cost_;
	}

	/**
	 * Return the number of active branches in the best-first
	 * BranchQueue.
	 */
	uint32_t size() {
		return branchQ_.size();
	}

	/**
	 * Reset the PathManager, clearing out the priority queue and
	 * resetting the RangeStatePool.
	 */
	void reset(TReadId patid) {
		branchQ_.reset(patid); // reset the priority queue
		assert(branchQ_.empty());
		bpool.reset();    // reset the Branch pool
		rpool.reset();    // reset the Pos pool
		assert(bpool.empty());
		assert(rpool.empty());
		ASSERT_ONLY(branchSet_.clear());
		assert_eq(0, branchSet_.size());
		assert_eq(0, branchQ_.size());
		minCost = 0;
	}

#ifndef NDEBUG
	/**
	 * Return true iff Branch b is in the priority queue;
	 */
	bool contains(Branch *b) const {
		bool ret = branchSet_.find(b) != branchSet_.end();
		assert(!ret || !b->exhausted_);
		return ret;
	}

	/**
	 * Do a consistenty-check on the collection of branches contained
	 * in this PathManager.
	 */
	bool repOk() {
		if(empty()) return true;
		assert(branchQ_.repOk(branchSet_));
		return true;
	}
#endif

	/**
	 * Return true iff the priority queue of branches is empty.
	 */
	bool empty() const {
		bool ret = branchQ_.empty();
		assert_eq(ret, branchSet_.empty());
		return ret;
	}

	/**
	 * Curtail the given branch, and possibly remove it from or
	 * re-insert it into the priority queue.
	 */
	void curtail(Branch *br, int seedLen) {
		assert(!br->exhausted_);
		assert(!br->curtailed_);
		uint16_t origCost = br->cost_;
		br->curtail(rpool, seedLen);
		assert(br->curtailed_);
		assert_geq(br->cost_, origCost);
		if(br->exhausted_) {
			assert(br == front(false));
			ASSERT_ONLY(Branch *popped =) dequeue();
			assert(popped == br);
			br->free(rpool, bpool);
		} else if(br->cost_ != origCost) {
			// Re-insert the newly-curtailed branch
			assert(br == front(false));
			Branch *popped = dequeue();
			assert(popped == br);
			enqueue(popped);
		} else {
			// Either there was no change in cost (e.g. a substitution
			// at a 0-quality position is possible) or there was a
			// delayed change in cost
		}
	}

	/**
	 * If the frontmost branch is a curtailed branch, split off an
	 * extendable branch and add it to the queue.
	 */
	bool splitAndPrep(Visit& visit,
	                  const CoordMap& cmap,
	                  RandomSource& rand,
	                  uint32_t qlen,
					  int qdep0,
	                  uint32_t qualLim,
	                  int seedLen,
	                  const EbwtParams& ep,
	                  const uint8_t* ebwt)
	{
		if(empty()) return true;
		// This counts as a backtrack
		if(btCnt_ != NULL && (*btCnt_ == 0)) {
			// Abruptly end search
			return false;
		}
		while(!empty()) {
			Branch *f = front(true);
			assert(!f->exhausted_);
			while(f->delayedIncrease_) {
				assert(!f->exhausted_);
				if(f->delayedIncrease_) {
					assert_neq(0, f->delayedCost_);
					ASSERT_ONLY(Branch *popped =) dequeue();
					assert(popped == f);
					f->cost_ = f->delayedCost_;
					f->delayedIncrease_ = false;
					f->delayedCost_ = 0;
					enqueue(f);
					assert(!empty());
				}
				f = front(true);
				assert(!f->exhausted_);
			}
			if(f->curtailed_) {
				uint16_t origCost = f->cost_;
				// This counts as a backtrack
				if(btCnt_ != NULL) {
					if(--(*btCnt_) == 0) {
						// Abruptly end search
						return false;
					}
				}
				Branch* newbr = f->splitBranch(
					rpool, bpool, visit, cmap, rand, qlen, qdep0,
				    qualLim, seedLen, triedGap_, ep, ebwt);
				assert(f->repOk());
				if(newbr == NULL && !f->exhausted_ && origCost != f->cost_) {
					// After eliminating visited outgoing paths, the
					// Branch's cost went up, so we have to re-insert
					// it and keep popping
					assert(!f->delayedIncrease_);
					ASSERT_ONLY(Branch *popped =) dequeue();
					assert(popped == f);
					enqueue(f);
					assert(!empty());
					continue;
				} else if(newbr == NULL && f->exhausted_) {
					// After eliminating visited outgoing paths, the
					// Branch was exhausted
					ASSERT_ONLY(Branch *popped =) dequeue();
					assert(popped == f);
					f->free(rpool, bpool);
					continue;
				} else if(newbr == NULL) {
					return false;
				}
				// If f is exhausted, get rid of it immediately
				if(f->exhausted_) {
					assert(!f->delayedIncrease_);
					ASSERT_ONLY(Branch *popped =) dequeue();
					assert(popped == f);
					f->free(rpool, bpool);
				}
				assert_eq(origCost, f->cost_);
				assert(newbr != NULL);
				enqueue(newbr);
				assert(newbr == front(false)); // already prepped in Branch.init(...)
			} else {
				prep(ep, ebwt);
			}
			break;
		}
		return true;
	}

	/**
	 * Return true iff the front element of the queue is prepped.
	 */
	bool prepped() {
		return front(false)->prepped_;
	}

	/**
	 * Return true iff we've tried a gap at some point so far.
	 */
	bool triedGap() const { return triedGap_; }

	/**
	 * Set flag indicating we've tried a gap.
	 */
	void setTriedGap() { triedGap_ = true; }

protected:

	/**
	 * Prep the next branch to be extended in advanceBranch().
	 */
	void prep(const EbwtParams& ep, const uint8_t* ebwt) {
		if(!branchQ_.empty()) {
			branchQ_.front(false)->prep(ep, ebwt);
		}
	}

	BranchQueue branchQ_; // priority queue for selecting lowest-cost Branch
	// set of branches in priority queue, for sanity checks
	ASSERT_ONLY(std::set<Branch*> branchSet_);

public:

	ChunkPool *cpool; // pool for generic chunks of memory
	AllocOnlyPool<Branch> bpool; // pool for allocating Branches
	AllocOnlyPool<Pos> rpool; // pool for allocating RangeStates
	/// The minimum possible cost for any alignments obtained by
	/// advancing further
	uint16_t minCost;

protected:
	/// Pointer to the aligner's per-read backtrack counter.  We
	/// increment it in splitBranch.
	bool triedGap_; // true iff a gap has been tried at some point in this search
	int *btCnt_;
};

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
class RangeSource {

public:
	RangeSource() :
		done(false), foundRange(false), curEbwt_(NULL) { }

	virtual ~RangeSource() { }

	/// Set query to find ranges for
	virtual void setQuery(Read& r, Range *partial) = 0;
	/// Set up the range search.
	virtual void initBranch(PathManager& pm) = 0;
	/// Advance the range search by one memory op
	virtual void advanceBranch(int until, uint16_t minCost, PathManager& pm) = 0;

	/// Return the last valid range found
	virtual Range& range() = 0;
	/// Return ptr to index this RangeSource is currently getting ranges from
	const Ebwt *curEbwt() const { return curEbwt_; }

	/// All searching w/r/t the current query is finished
	bool done;
	/// Set to true iff the last call to advance yielded a range
	bool foundRange;
protected:
	/// ptr to index this RangeSource is currently getting ranges from
	const Ebwt *curEbwt_;
};

/**
 * Abstract parent of RangeSourceDrivers.
 */
template<typename TRangeSource>
class RangeSourceDriver {

public:
	RangeSourceDriver(bool _done, uint32_t minCostAdjustment = 0) :
		foundRange(false), done(_done), minCostAdjustment_(minCostAdjustment)
	{
		minCost = minCostAdjustment_;
	}

	virtual ~RangeSourceDriver() { }

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc, Range *r) {
		// Clear our buffer of previously-dished-out top offsets
		setQueryImpl(patsrc, r);
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *r) = 0;

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advance(int until) {
		advanceImpl(until);
	}
	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advanceImpl(int until) = 0;
	/**
	 * Return the range found.
	 */
	virtual Range& range() = 0;

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	virtual bool mate1() const = 0;

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	virtual bool fw() const = 0;

	virtual void removeMate(int m) { }

	/// Set to true iff we just found a range.
	bool foundRange;

	/**
	 * Set to true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	bool done;

	/**
	 * The minimum "combined" stratum/qual cost that could possibly
	 * result from subsequent calls to advance() for this driver.
	 */
	uint16_t minCost;

	/**
	 * Adjustment to the minCost given by the caller that constructed
	 * the object.  This is useful if we know the lowest-cost branch is
	 * likely to cost less than the any of the alignments that could
	 * possibly result from advancing (e.g. when we're going to force a
	 * mismatch somewhere down the line).
	 */
	uint16_t minCostAdjustment_;
};

/**
 * A concrete driver wrapper for a single RangeSource.
 */
template<typename TRangeSource>
class SingleRangeSourceDriver : public RangeSourceDriver<TRangeSource> {

public:
	SingleRangeSourceDriver(
		EbwtSearchParams& params,
		TRangeSource* rs,
		bool fw,
		HitSink& sink,
		HitSinkPerThread* sinkPt,
		EList<SString<char> >& os,
		bool mate1,
		uint32_t minCostAdjustment,
		ChunkPool* pool,
		int *btCnt) :
		RangeSourceDriver<TRangeSource>(true, minCostAdjustment),
		len_(0), mate1_(mate1),
		sinkPt_(sinkPt),
		params_(params),
		fw_(fw), rs_(rs),
		ebwtFw_(rs_->curEbwt()->fw()),
		pm_(pool, btCnt)
	{
		assert(rs_ != NULL);
	}

	virtual ~SingleRangeSourceDriver() {
		delete rs_; rs_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *r) {
		this->done = false;
		pm_.reset(patsrc->patid());
		Read* buf = mate1_ ? &patsrc->bufa() : &patsrc->bufb();
		len_ = buf->length();
		assert_eq(buf->qual.length(), len_);
		rs_->setQuery(*buf, r);
		initRangeSource((fw_ == ebwtFw_) ? buf->qual : buf->qualRev,
		                buf->alts,
		                (fw_ == ebwtFw_) ? buf->altQual : buf->altQualRev);
		assert_gt(len_, 0);
		if(this->done) return;
		if(!rs_->done) {
			rs_->initBranch(pm_); // set up initial branch
		}
		uint16_t icost = (r != NULL) ? r->cost() : 0;
		this->minCost = max<uint16_t>(icost, this->minCostAdjustment_);
		this->done = rs_->done;
		this->foundRange = rs_->foundRange;
		if(!pm_.empty()) {
			assert(!pm_.front(false)->curtailed_);
			assert(!pm_.front(false)->exhausted_);
		}
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advanceImpl(int until) {
		if(this->done || pm_.empty()) {
			this->done = true;
			return;
		}
		assert(!pm_.empty());
		assert(!pm_.front(false)->curtailed_);
		assert(!pm_.front(false)->exhausted_);
		// Advance the RangeSource for the forward-oriented read
		ASSERT_ONLY(uint16_t oldMinCost = this->minCost);
		ASSERT_ONLY(uint16_t oldPmMinCost = pm_.minCost);
		rs_->advanceBranch(until, this->minCost, pm_);
		this->done = pm_.empty();
		if(pm_.minCost != 0) {
			this->minCost = max(pm_.minCost, this->minCostAdjustment_);
		} else {
			// pm_.minCost is 0 because we reset it due to exceptional
			// circumstances, like out-of-memory
		}
#ifndef NDEBUG
		{
			bool error = false;
			if(pm_.minCost != 0 && pm_.minCost < oldPmMinCost) {
				cerr << "PathManager's cost went down" << endl;
				error = true;
			}
			if(this->minCost < oldMinCost) {
				cerr << "this->minCost cost went down" << endl;
				error = true;
			}
			if(error) {
				cerr << "pm.minCost went from " << oldPmMinCost
				     << " to " << pm_.minCost << endl;
				cerr << "this->minCost went from " << oldMinCost
				     << " to " << this->minCost << endl;
				cerr << "this->minCostAdjustment_ == "
				     << this->minCostAdjustment_ << endl;
			}
			assert(!error);
		}
#endif
		this->foundRange = rs_->foundRange;
		if(!pm_.empty()) {
			assert(!pm_.front(false)->curtailed_);
			assert(!pm_.front(false)->exhausted_);
		}
	}

	/**
	 * Return the range found.
	 */
	virtual Range& range() {
		rs_->range().setMate1(mate1_);
		return rs_->range();
	}

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	bool mate1() const {
		return mate1_;
	}

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	bool fw() const {
		return fw_;
	}

	virtual void initRangeSource(const BTString& qual,
	                             int alts,
	                             const BTString* altQuals) = 0;

protected:

	// Progress state
	size_t len_;
	bool mate1_;

	// Temporary HitSink; to be deleted
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams& params_;
	bool              fw_;
	TRangeSource*     rs_; // delete this in destructor
	bool              ebwtFw_;
	PathManager       pm_;
};

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
template<typename TRangeSource>
class StubRangeSourceDriver : public RangeSourceDriver<TRangeSource> {

	typedef EList<RangeSourceDriver<TRangeSource>*> TRangeSrcDrPtrVec;

public:

	StubRangeSourceDriver() :
		RangeSourceDriver<TRangeSource>(false)
	{
		this->done = true;
		this->foundRange = false;
	}

	virtual ~StubRangeSourceDriver() { }

	/// Set query to find ranges for
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *r) { }

	/// Advance the range search by one memory op
	virtual void advanceImpl(int until) { }

	/// Return the last valid range found
	virtual Range& range() { throw 1; }

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	virtual bool mate1() const {
		return true;
	}

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	virtual bool fw() const {
		return true;
	}

};

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
template<typename TRangeSource>
class ListRangeSourceDriver : public RangeSourceDriver<TRangeSource> {

	typedef EList<RangeSourceDriver<TRangeSource>*> TRangeSrcDrPtrVec;

public:

	ListRangeSourceDriver(const TRangeSrcDrPtrVec& rss) :
		RangeSourceDriver<TRangeSource>(false),
		cur_(0), ham_(0), rss_(rss) /* copy */,
		patsrc_(NULL), seedRange_(NULL)
	{
		assert_gt(rss_.size(), 0);
		assert(!this->done);
	}

	virtual ~ListRangeSourceDriver() {
		for(size_t i = 0; i < rss_.size(); i++) {
			delete rss_[i];
		}
	}

	/// Set query to find ranges for
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *r) {
		cur_ = 0; // go back to first RangeSource in list
		this->done = false;
		rss_[cur_]->setQuery(patsrc, r);
		patsrc_ = patsrc; // so that we can call setQuery on the other elements later
		seedRange_ = r;
		this->done = (cur_ == rss_.size()-1) && rss_[cur_]->done;
		this->minCost = max(rss_[cur_]->minCost, this->minCostAdjustment_);
		this->foundRange = rss_[cur_]->foundRange;
	}

	/// Advance the range search by one memory op
	virtual void advanceImpl(int until) {
		assert(!this->done);
		assert_lt(cur_, rss_.size());
		if(rss_[cur_]->done) {
			// Move on to next RangeSourceDriver
			if(cur_ < rss_.size()-1) {
				rss_[++cur_]->setQuery(patsrc_, seedRange_);
				this->minCost = max(rss_[cur_]->minCost, this->minCostAdjustment_);
				this->foundRange = rss_[cur_]->foundRange;
			} else {
				// No RangeSources in list; done
				cur_ = 0xffffffff;
				this->done = true;
			}
		} else {
			// Advance current RangeSource
			rss_[cur_]->advance(until);
			this->done = (cur_ == rss_.size()-1 && rss_[cur_]->done);
			this->foundRange = rss_[cur_]->foundRange;
			this->minCost = max(rss_[cur_]->minCost, this->minCostAdjustment_);
		}
	}

	/// Return the last valid range found
	virtual Range& range() { return rss_[cur_]->range(); }

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	virtual bool mate1() const {
		return rss_[0]->mate1();
	}

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	virtual bool fw() const {
		return rss_[0]->fw();
	}

protected:

	uint32_t cur_;
	uint32_t ham_;
	TRangeSrcDrPtrVec rss_;
	PatternSourcePerThread* patsrc_;
	Range *seedRange_;
};

/**
 * A RangeSourceDriver that wraps a set of other RangeSourceDrivers and
 * chooses which one to advance at any given moment by picking one with
 * minimal "cumulative cost" so far.
 *
 * Note that costs have to be "adjusted" to account for the fact that
 * the alignment policy for the underlying RangeSource might force
 * mismatches.
 */
template<typename TRangeSource>
class CostAwareRangeSourceDriver : public RangeSourceDriver<TRangeSource> {

	typedef RangeSourceDriver<TRangeSource>* TRangeSrcDrPtr;
	typedef EList<TRangeSrcDrPtr> TRangeSrcDrPtrVec;

public:

	CostAwareRangeSourceDriver(
			const TRangeSrcDrPtrVec* rss,
			bool mixesReads) :
		RangeSourceDriver<TRangeSource>(false),
		rss_(), active_(),
		lastRange_(NULL), delayedRange_(NULL), patsrc_(NULL),
		mixesReads_(mixesReads)
	{
		if(rss != NULL) {
			rss_ = (*rss);
		}
		paired_ = false;
		this->foundRange = false;
		this->done = false;
		if(rss_.empty()) {
			return;
		}
		calcPaired();
		active_ = rss_;
		this->minCost = 0;
	}

	/// Destroy all underlying RangeSourceDrivers
	virtual ~CostAwareRangeSourceDriver() {
		const size_t rssSz = rss_.size();
		for(size_t i = 0; i < rssSz; i++) {
			delete rss_[i];
		}
		rss_.clear();
		active_.clear();
	}

	/// Set query to find ranges for
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *r) {
		this->done = false;
		this->foundRange = false;
		lastRange_ = NULL;
		delayedRange_ = NULL;
		ASSERT_ONLY(allTopsRc_.clear());
		patsrc_ = patsrc;
		rand_.init(patsrc->bufa().seed);
		const size_t rssSz = rss_.size();
		if(rssSz == 0) return;
		for(size_t i = 0; i < rssSz; i++) {
			// Assuming that all
			rss_[i]->setQuery(patsrc, r);
		}
		active_ = rss_;
		this->minCost = 0;
		sortActives();
	}

	/**
	 * Add a new RangeSource to the list and re-sort it.
	 */
	void addSource(TRangeSrcDrPtr p, Range *r) {
		assert(!this->foundRange);
		this->lastRange_ = NULL;
		this->delayedRange_ = NULL;
		this->done = false;
		if(patsrc_ != NULL) {
			p->setQuery(patsrc_, r);
		}
		rss_.push_back(p);
		active_.push_back(p);
		calcPaired();
		this->minCost = 0;
		sortActives();
	}

	/**
	 * Free and remove all contained RangeSources.
	 */
	void clearSources() {
		const size_t rssSz = rss_.size();
		for(size_t i = 0; i < rssSz; i++) {
			delete rss_[i];
		}
		rss_.clear();
		active_.clear();
		paired_ = false;
	}

	/**
	 * Return the number of RangeSources contained within.
	 */
	size_t size() const {
		return rss_.size();
	}

	/**
	 * Return true iff no RangeSources are contained within.
	 */
	bool empty() const {
		return rss_.empty();
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advance(int until) {
		ASSERT_ONLY(uint16_t precost = this->minCost);
		assert(!this->done);
		assert(!this->foundRange);
		until = max<int>(until, ADV_COST_CHANGES);
		advanceImpl(until);
		assert(!this->foundRange || lastRange_ != NULL);
		if(this->foundRange) {
			assert_eq(range().cost(), precost);
		}
	}

	/// Advance the range search
	virtual void advanceImpl(int until) {
		lastRange_ = NULL;
		ASSERT_ONLY(uint16_t iminCost = this->minCost);
		const size_t actSz = active_.size();
		assert(sortedActives());
		if(delayedRange_ != NULL) {
			assert_eq(iminCost, delayedRange_->cost());
			lastRange_ = delayedRange_;
			delayedRange_ = NULL;
			this->foundRange = true;
			assert_eq(range().cost(), iminCost);
			if(!active_.empty()) {
				assert_geq(active_[0]->minCost, this->minCost);
				this->minCost = max(active_[0]->minCost, this->minCost);
			} else {
				this->done = true;
			}
			return; // found a range
		}
		assert(delayedRange_ == NULL);
		if(mateEliminated() || actSz == 0) {
			// No more alternatoves; clear the active set and signal
			// we're done
			active_.clear();
			this->done = true;
			return;
		}
		// Advance lowest-cost RangeSourceDriver
		TRangeSrcDrPtr p = active_[0];
		uint16_t precost = p->minCost;
		assert(!p->done || p->foundRange);
		if(!p->foundRange) {
			p->advance(until);
		}
		bool needsSort = false;
		if(p->foundRange) {
			Range *r = &p->range();
			assert_eq(r->cost(), iminCost);
			needsSort = foundFirstRange(r); // may set delayedRange_; re-sorts active_
			assert_eq(lastRange_->cost(), iminCost);
			if(delayedRange_ != NULL) assert_eq(delayedRange_->cost(), iminCost);
			p->foundRange = false;
		}
		if(p->done || (precost != p->minCost) || needsSort) {
			sortActives();
			if(mateEliminated() || active_.empty()) {
				active_.clear();
				this->done = (delayedRange_ == NULL);
			}
		}
		assert(sortedActives());
		assert(lastRange_ == NULL || lastRange_->cost() == iminCost);
		assert(delayedRange_ == NULL || delayedRange_->cost() == iminCost);
	}

	/// Return the last valid range found
	virtual Range& range() {
		assert(lastRange_ != NULL);
		return *lastRange_;
	}

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	virtual bool mate1() const {
		return rss_[0]->mate1();
	}

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	virtual bool fw() const {
		return rss_[0]->fw();
	}

	virtual void removeMate(int m) {
		bool qmate1 = (m == 1);
		assert(paired_);
		for(size_t i = 0; i < active_.size(); i++) {
			if(active_[i]->mate1() == qmate1) {
				active_[i]->done = true;
			}
		}
		sortActives();
		assert(mateEliminated());
	}

protected:

	/**
	 * Set paired_ to true iff there are mate1 and mate2 range sources
	 * in the rss_ vector.
	 */
	void calcPaired() {
		const size_t rssSz = rss_.size();
		bool saw1 = false;
		bool saw2 = false;
		for(size_t i = 0; i < rssSz; i++) {
			if(rss_[i]->mate1()) saw1 = true;
			else saw2 = true;
		}
		assert(saw1 || saw2);
		paired_ = saw1 && saw2;
	}

	/**
	 * Return true iff one mate or the other has been eliminated.
	 */
	bool mateEliminated() {
		if(!paired_) return false;
		bool mate1sLeft = false;
		bool mate2sLeft = false;
		// If this RangeSourceDriver is done, shift everyone else
		// up and delete it
		const size_t rssSz = active_.size();
		for(size_t i = 0; i < rssSz; i++) {
			if(!active_[i]->done) {
				if(active_[i]->mate1()) mate1sLeft = true;
				else mate2sLeft = true;
			}
		}
		return !mate1sLeft || !mate2sLeft;
	}

	/**
	 * We found a range; check whether we should attempt to find a
	 * range of equal quality from the opposite strand so that we can
	 * resolve the strand bias.  Return true iff we modified the cost
	 * of one or more items after the first item.
	 */
	bool foundFirstRange(Range* r) {
		assert(r != NULL);
		this->foundRange = true;
		lastRange_ = r;
		if(gStrandFix) {
			// We found a range but there may be an equally good range
			// on the other strand; let's try to get it.
			const size_t sz = active_.size();
			for(size_t i = 1; i < sz; i++) {
				// Same mate, different orientation?
				if(rss_[i]->mate1() == r->mate1() && rss_[i]->fw() != r->fw()) {
					// Yes; see if it has the same cost
					TRangeSrcDrPtr p = active_[i];
					uint16_t minCost = max(this->minCost, p->minCost);
					if(minCost > r->cost()) break;
					// Yes, it has the same cost
					assert_eq(minCost, r->cost()); // can't be better
					// Advance it until it's done, we've found a range,
					// or its cost increases
					if(gVerbose) cout << " Looking for opposite range to avoid strand bias:" << endl;
					while(!p->done && !p->foundRange) {
						p->advance(ADV_COST_CHANGES);
						assert_geq(p->minCost, minCost);
						if(p->minCost > minCost) break;
					}
					if(p->foundRange) {
						// Found one!  Now we have to choose which one
						// to give out first; we choose randomly using
						// the size of the ranges as weights.
						delayedRange_ = &p->range();
						size_t tot = (delayedRange_->bot() - delayedRange_->top()) +
						             (lastRange_->bot()    - lastRange_->top());
						uint32_t rq = (uint32_t)(rand_.nextU32() % tot);
						// We picked this range, not the first one
						if(rq < (delayedRange_->bot() - delayedRange_->top())) {
							Range *tmp = lastRange_;
							lastRange_ = delayedRange_;
							delayedRange_ = tmp;
						}
						p->foundRange = false;
					}
					// Return true iff we need to force a re-sort
					return true;
				}
			}
			// OK, now we have a choice of two equally good ranges from
			// each strand.
		}
		return false;
	}

	/**
	 * Sort all of the RangeSourceDriver ptrs in the rss_ array so that
	 * the one with the lowest cumulative cost is at the top.  Break
	 * ties randomly.  Just do selection sort for now; we don't expect
	 * the list to be long.
	 */
	void sortActives() {
		TRangeSrcDrPtrVec& vec = active_;
		if(vec.empty()) return;
		size_t sz = vec.size();
		// Selection sort / removal outer loop
		for(size_t i = 0; i < sz;) {
			// Remove elements that we're done with
			if(vec[i]->done && !vec[i]->foundRange) {
				vec.remove(i);
				if(sz == 0) break;
				else sz--;
				continue;
			}
			uint16_t minCost = vec[i]->minCost;
			size_t minOff = i;
			// Selection sort inner loop
			for(size_t j = i+1; j < sz; j++) {
				if(vec[j]->done && !vec[j]->foundRange) {
					// We'll get rid of this guy later
					continue;
				}
				if(vec[j]->minCost < minCost) {
					minCost = vec[j]->minCost;
					minOff = j;
				} else if(vec[j]->minCost == minCost) {
					// Possibly randomly pick the other
					if(rand_.nextU32() & 0x1000) {
						minOff = j;
					}
				}
			}
			// Do the swap, if necessary
			if(i != minOff) {
				assert_leq(minCost, vec[i]->minCost);
				TRangeSrcDrPtr tmp = vec[i];
				vec[i] = vec[minOff];
				vec[minOff] = tmp;
			}
			i++;
		}
		if(delayedRange_ == NULL && !vec.empty()) {
			assert_geq(this->minCost, this->minCostAdjustment_);
			assert_geq(vec[0]->minCost, this->minCost);
			this->minCost = vec[0]->minCost;
		}
		assert(sortedActives());
	}

#ifndef NDEBUG
	/**
	 * Check that the rss_ array is sorted by minCost; assert if it's
	 * not.
	 */
	bool sortedActives() const {
		// Selection sort outer loop
		const TRangeSrcDrPtrVec& vec = active_;
		const size_t sz = vec.size();
		for(size_t i = 0; i < sz; i++) {
			assert(!vec[i]->done || vec[i]->foundRange);
			for(size_t j = i+1; j < sz; j++) {
				assert(!vec[j]->done || vec[j]->foundRange);
				assert_leq(vec[i]->minCost, vec[j]->minCost);
			}
		}
		if(delayedRange_ == NULL && sz > 0) {
			// Only assert this if there's no delayed range; if there's
			// a delayed range, the minCost is its cost, not the 0th
			// element's cost
			assert_leq(vec[0]->minCost, this->minCost);
		}
		return true;
	}
#endif

	/// List of all the drivers
	TRangeSrcDrPtrVec rss_;
	/// List of all the as-yet-uneliminated drivers
	TRangeSrcDrPtrVec active_;
	/// Whether the list of drivers contains drivers for both mates 1 and 2
	bool paired_;
	uint32_t randSeed_;
	/// The random seed from the Aligner, which we use to randomly break ties
	RandomSource rand_;
	Range *lastRange_;
	Range *delayedRange_;
	PatternSourcePerThread* patsrc_;
	bool mixesReads_;
	ASSERT_ONLY(std::set<int64_t> allTopsRc_);
};

#endif /* RANGE_SOURCE_H_ */
