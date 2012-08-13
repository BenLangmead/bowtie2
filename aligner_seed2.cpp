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

#include "aligner_seed2.h"
#include "assert_helpers.h"
#include "bt2_idx.h"

/**
 * Drive the process of descending from all search roots.
 */
void DescentDriver::go() {
}

/**
 * Recalculate our summary of the outgoing edges from this descent.  When
 * deciding what outgoing edges are legal, we abide by constraints.
 * Typically, they limit the total of the penalties accumulated so far, as
 * a function of distance from the search root.  E.g. a constraint might
 * disallow any gaps or mismatches within 20 ply of the search root, then
 * allow 1 mismatch within 30 ply, then allow up to 1 mismatch or 1 gap
 * within 40 ply, etc.
 *
 * Return the total number of valid outgoing edges found.
 *
 * TODO: Eliminate outgoing gap edges that are redundant with others.
 */
size_t Descent::recalcOutgoing(
	const DescentQuery& q,     // query string
	DescentConstraints& cons,  // constraints
	const Scoring& sc,         // scoring scheme
	EFactory<DescentPos>& pf)  // factory with DescentPoss
{
	assert(!empty());
	assert(out_.empty());
	assert(repOk(&q));
	// Get initial 5' and 3' offsets
	bool toward3p = (l2r_ == fw_);
	size_t off5p;
	assert_geq(al5pf_, al5pi_);
	if(toward3p) {
		off5p = al5pf_;
	} else {
		assert_gt(al5pi_, 0);
		off5p = al5pi_ - 1;
	}
	size_t off3p = q.len - off5p - 1;
	// By "depth" we essentially mean the number of characters already aligned
	size_t depth = al5pf_ - al5pi_;
	// Get gap penalties
	TScore pen_rdg_ex = sc.readGapExtend(), pen_rfg_ex = sc.refGapExtend();
	TScore pen_rdg_op = sc.readGapOpen(),   pen_rfg_op = sc.refGapOpen();
	// Top and bot in the direction of the descent
	TIndexOff top  = l2r_ ? topb_ : topf_;
	TIndexOff bot  = l2r_ ? botb_ : botf_;
	// Top and bot in the opposite direction
	TIndexOff topp = l2r_ ? topf_ : topb_;
	TIndexOff botp = l2r_ ? botf_ : botb_;
	assert_eq(botp - topp, bot - top);
	DescentEdge edge;
	size_t nout = 0;
	// Enumerate all outgoing edges, starting at the root and going out
	for(size_t i = 0; i < len_; i++) {
		assert_gt(bot, top);
		assert_eq(botp - topp, bot - top);
		size_t d = descid_ + i;
		TScore maxpen = cons[depth]; // maximum penalty total at this position
		assert_geq(maxpen, pen_);    // can't have already exceeded max penalty
		TScore diff = maxpen - pen_; // room we have left
		// Get pointer to SA ranges in the direction of descent
		const TIndexOff *t  = l2r_ ? pf[d].topb : pf[d].topf;
		const TIndexOff *b  = l2r_ ? pf[d].botb : pf[d].botf;
		const TIndexOff *tp = l2r_ ? pf[d].topf : pf[d].topb;
		const TIndexOff *bp = l2r_ ? pf[d].botf : pf[d].botb;
		// What are the read char / quality?
		std::pair<int, int> p = q.get(off5p, fw_);
		int c = p.first;
		assert_range(0, 4, c);
		// Only entertain edits if there is at least one type of edit left and
		// there is some penalty budget left
		if(!pf[d].flags.exhausted() && diff > 0) {
			// What would the penalty be if we mismatched at this position?
			// This includes the case where the mismatch is for an N in the
			// read.
			int q = p.second;
			assert_lt(q, 64);
			TScore pen_mm = sc.mm(c, q);
			if(pen_mm <= diff) {
				for(int j = 0; j < 4; j++) {
					if(j == c) continue; // Match, not mismatch
					if(b[j] <= t[j]) {
						continue; // No outgoing edge with this nucleotide
					}
					if(!pf[d].flags.mmExplore(j)) {
						continue; // Already been explored
					}
					TIndexOff width = b[j] - t[j];
					Edit edit((uint32_t)off5p, j, c, EDIT_TYPE_MM);
					DescentPriority pri(pen_ + pen_mm, depth, width, rootpri_);
					TIndexOff topf = pf[d].topf[j], botf = pf[d].botf[j];
					TIndexOff topb = pf[d].topb[j], botb = pf[d].botb[j];
					assert_eq(botb - topb, botf - topf);
					edge.init(edit, pri, topf, botf, topb, botb);
					out_.update(edge);
					nout++;
				}
			}
			bool gapsAllowed = (off5p >= sc.gapbar && off3p >= sc.gapbar);
			if(gapsAllowed) {
				bool rdex = false, rfex = false;
				if(i == 0) {
					// If we're at the root of the descent, and the descent
					// branched on a gap, then this could be scored as an
					// extension of that gap.
					if(pen_rdg_ex <= diff && edit_.isReadGap()) {
						// Extension of a read gap
						rdex = true;
						for(int j = 0; j < 4; j++) {
							if(b[j] <= t[j]) {
								continue; // No outgoing edge with this nucleotide
							}
							if(!pf[d].flags.rdgExplore(j)) {
								continue; // Already been explored
							}
							TIndexOff width = b[j] - t[j];
							Edit edit((uint32_t)off5p, j, '-', EDIT_TYPE_READ_GAP);
							DescentPriority pri(pen_ + pen_rdg_ex, depth, width, rootpri_);
							TIndexOff topf = pf[d].topf[j], botf = pf[d].botf[j];
							TIndexOff topb = pf[d].topb[j], botb = pf[d].botb[j];
							assert_eq(botb - topb, botf - topf);
							edge.init(edit, pri, topf, botf, topb, botb);
							out_.update(edge);
							nout++;
						}
					}
					if(pen_rfg_ex <= diff && edit_.isRefGap()) {
						// Extension of a reference gap
						rfex = true;
						if(!pf[d].flags.rfgExplore()) {
							continue; // Already been explored
						}
						TIndexOff width = bot - top;
						Edit edit((uint32_t)off5p, '-', c, EDIT_TYPE_REF_GAP);
						DescentPriority pri(pen_ + pen_rfg_ex, depth, width, rootpri_);
						TIndexOff topf = l2r_ ? topp : top;
						TIndexOff topb = l2r_ ? top : topp;
						TIndexOff botf = l2r_ ? botp : bot;
						TIndexOff botb = l2r_ ? bot : botp;
						edge.init(edit, pri, topf, botf, topb, botb);
						out_.update(edge);
						nout++;
					}
				}
				if(pen_rdg_op <= diff && !rdex) {
					// Opening a new read gap
					for(int j = 0; j < 4; j++) {
						if(b[j] <= t[j]) {
							continue; // No outgoing edge with this nucleotide
						}
						if(!pf[d].flags.rdgExplore(j)) {
							continue; // Already been explored
						}
						TIndexOff width = b[j] - t[j];
						Edit edit((uint32_t)off5p, j, '-', EDIT_TYPE_READ_GAP);
						DescentPriority pri(pen_ + pen_rdg_op, depth, width, rootpri_);
						TIndexOff topf = pf[d].topf[j], botf = pf[d].botf[j];
						TIndexOff topb = pf[d].topb[j], botb = pf[d].botb[j];
						assert_eq(botb - topb, botf - topf);
						edge.init(edit, pri, topf, botf, topb, botb);
						out_.update(edge);
						nout++;
					}
				}
				if(pen_rfg_op <= diff && !rfex) {
					// Opening a new reference gap
					TIndexOff width = bot - top;
					Edit edit((uint32_t)off5p, '-', c, EDIT_TYPE_REF_GAP);
					DescentPriority pri(pen_ + pen_rfg_op, depth, width, rootpri_);
					TIndexOff topf = l2r_ ? topp : top;
					TIndexOff topb = l2r_ ? top : topp;
					TIndexOff botf = l2r_ ? botp : bot;
					TIndexOff botb = l2r_ ? bot : botp;
					edge.init(edit, pri, topf, botf, topb, botb);
					out_.update(edge);
					nout++;
				}
			}
		}
		// Update off5p, off3p, depth
		if(toward3p) {
			assert_gt(off3p, 0);
			off5p++;
			off3p--;
		} else {
			assert_gt(off5p, 0);
			off3p++;
			off5p--;
		}
		depth++;
		// Update top and bot
		top = t[c], topp = tp[c];
		bot = b[c], botp = bp[c];
	}
	lastRecalc_ = (nout <= 5);
	return nout;
}

/**
 * Take the best outgoing edge and place it in the heap.  When deciding what outgoing
 * edges are possible, we have to abide by constraints.  Typically, these
 * constraints limit the total penalties accumulated so far versus distance
 * from the search root.  E.g. a constraint might disallow any gaps or
 * mismatches within 20 ply of the search root, then allow 1 mismatch within
 * 30 ply, then allow up to 1 mismatch or 1 gap within 40 ply, etc.
 */
void Descent::followBestOutgoing(
	const DescentQuery& q,     // query string
	const Ebwt& ebwtFw,        // forward index
	const Ebwt& ebwtBw,        // mirror index
	DescentConstraints& cons,  // constraints
	const Scoring& sc,         // scoring scheme
	EFactory<DescentPos>& pf,  // factory with DescentPoss
	EFactory<Descent>& df,     // factory with Descent
	EHeap<TDescentPair>& heap, // heap of descents
	DescentMetrics& met)       // metrics
{
	// We assume this descent has been popped off the heap.  We'll re-add it if
	// it hasn't been exhausted.
	assert(!empty());
	assert(!out_.empty());
	DescentPriority best = out_.bestPri();
	DescentEdge e = out_.rotate();
	if(out_.empty()) {
		if(!lastRecalc_) {
			recalcOutgoing(q, cons, sc, pf);
			assert(!empty());
		} else {
			assert(empty());
		}
	}
	if(!empty()) {
		// Re-insert this
		DescentPriority best2 = out_.bestPri();
		assert_leq(best, best2);
		heap.insert(make_pair(best2, descid_));
	}
	TDescentId id = df.alloc();
	// TODO: are we changing direction?
	DescentRoot r(
		e.e.pos,    // offset from 5' end
		l2r_,       // true -> left-to-right
		fw_,        // true -> aligning forward read, otherwise aligning revcomp
		q.len,      // query length
		rootpri_);  // root priority
	//df[id].init(
	//	q,          // query
	//	r,          // root
	//	cons,       // constraint scheme
	//	sc,         // scoring scheme
	//	e.topf,     // top of SA range in fw index
	//	e.botf,     // bottom of SA range in fw index
	//	e.topb,     // top of SA range in bw index
	//	e.botb,     // bottom of SA range in bw index
	//	e.e,        // edit
	//	descid_,    // parent ID
	//	id,         // self ID
	//	e.pri.pen,  // total penalties so far
	//	ebwtFw,     // forward index
	//	ebwtBw,     // mirror index
	//	pf,         // factory for descent poss
	//	heap,       // priority queue for descents
	//	met);       // metrics
}

/**
 * Given the forward and backward indexes, and given topf/botf/topb/botb, get
 * tloc, bloc ready for the next step.
 */
void Descent::nextLocsBi(
	const Ebwt& ebwtFw, // forward index
	const Ebwt& ebwtBw, // mirror index
	SideLocus& tloc,    // top locus
	SideLocus& bloc,    // bot locus
	TIndexOff topf,     // top in BWT
	TIndexOff botf,     // bot in BWT
	TIndexOff topb,     // top in BWT'
	TIndexOff botb)     // bot in BWT'
{
	assert_gt(botf, 0);
	// Which direction are we going in next?
	if(l2r_) {
		// Left to right; use BWT'
		if(botb - topb == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topb, ebwtBw.eh(), ebwtBw.ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topb, botb, ebwtBw.eh(), ebwtBw.ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	} else {
		// Right to left; use BWT
		if(botf - topf == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topf, ebwtFw.eh(), ebwtFw.ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topf, botf, ebwtFw.eh(), ebwtFw.ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	// Check if we should update the tracker with this refinement
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
}

/**
 * Advance this descent by following read matches as far as possible.
 *
 * This routine doesn't have to consider the whole gamut of constraints on
 * which outgoing edges can be followed.  If it is a root descent, it does have
 * to know how deep the no-edit constraint goes, though, so we can decide
 * whether using the ftab would potentially jump over relevant branch points.
 * Apart from that, though, we simply proceed as far as it can go by matching
 * characters in the query, irrespective of the constraints.
 * recalcOutgoing(...) and followBestOutgoing(...) do have to consider these
 * constraints, though.
 *
 * Conceptually, as we make descending steps, we have:
 * 1. Before each step, a single range indicating how we departed the previous
 *    step
 * 2. As part of each step, a quad of ranges indicating what range would result
 *    if we proceeded on an a, c, g ot t
 *
 * Return true iff it is possible to branch from this descent.  If we haven't
 * exceeded the no-branch depth.
 */
bool Descent::followMatches(
	const DescentQuery& q, // query string
	const Ebwt& ebwtFw,    // forward index
	const Ebwt& ebwtBw,    // mirror index
	EFactory<DescentPos>& pf,
	EHeap<TDescentPair>& heap,
	DescentMetrics& met)
{
	size_t nobranchDepth = 20;
	assert(repOk(&q));
	assert_eq(ebwtFw.eh().ftabChars(), ebwtBw.eh().ftabChars());
#ifndef NDEBUG
	for(int i = 0; i < 4; i++) {
		assert_eq(ebwtFw.fchr()[i], ebwtBw.fchr()[i]);
	}
#endif
	SideLocus tloc, bloc;
	TIndexOff topf = topf_, botf = botf_, topb = topb_, botb = botb_;
	bool toward3p = (l2r_ == fw_);
	size_t off5p;
	assert_geq(al5pf_, al5pi_);
	if(toward3p) {
		off5p = al5pf_;
	} else {
		assert_gt(al5pi_, 0);
		off5p = al5pi_ - 1;
	}
	size_t off3p = q.len - off5p - 1;
	assert_lt(off5p, q.len);
	assert_lt(off3p, q.len);
	bool firstPos = true;
	assert_eq(0, len_);
	if(root()) {
		// Check whether/how far we can jump using ftab
		int ftabLen = ebwtFw.eh().ftabChars();
		if(ftabLen > 1 && ftabLen <= nobranchDepth) {
			// Right-to-left
			size_t off_r2l = off5p;
			if(l2r_) {
				off_r2l += ftabLen;
				assert_lt(off_r2l, q.len);
			}
			ebwtFw.ftabLoHi(fw_ ? q.seq : q.seqrc, off_r2l, false, topf, botf);
			if(botf - topf == 0) {
				return false;
			}
			// Left-to-right
			size_t off_l2r = off5p;
			if(!l2r_) {
				assert_geq(off_l2r, ftabLen);
				off_l2r -= ftabLen;
			}
			ebwtBw.ftabLoHi(fw_ ? q.seq : q.seqrc, off_l2r, false, topb, botb);
			assert_eq(botf - topf, botb - topb);
			if(toward3p) {
				assert_geq(off3p, ftabLen);
				off5p += ftabLen; off3p -= ftabLen;
			} else {
				assert_geq(off5p, ftabLen);
				off5p -= ftabLen; off3p += ftabLen;
			}
			len_ += ftabLen;
			// Allocate DescentPos data structures and leave them empty.  We
			// jumped over them by doing our lookup in the ftab, so we have no
			// info about outgoing edges from them, besides the matching
			// outgoing edge from the last pos which is in topf/botf and
			// topb/botb.
			if(firstPos) {
				posid_ = pf.alloc();
				pf[posid_].reset();
				firstPos = false;
				for(int i = 1; i < ftabLen; i++) {
					size_t id = pf.alloc();
					pf[id].reset();
				}
			} else {
				for(int i = 0; i < ftabLen; i++) {
					size_t id = pf.alloc();
					pf[id].reset();
				}
			}
		} else {
			// Can't use ftab, use fchr instead
			int rdc = q.getc(off5p, fw_);
			assert_range(0, 3, rdc);
			topf = topb = ebwtFw.fchr()[rdc];
			botf = botb = ebwtFw.fchr()[rdc+1];
			if(botf - topf == 0) {
				return false;
			}
			if(toward3p) {
				assert_gt(off3p, 0);
				off5p++; off3p--;
			} else {
				assert_gt(off5p, 0);
				off5p--; off3p++;
			}
			len_++;
			// Allocate DescentPos data structure.  We could fill it with the
			// four ranges from fchr if we wanted to, but that will never be
			// relevant.
			if(firstPos) {
				posid_ = pf.alloc();
			} else {
				pf.alloc();
			}
		}
		assert_gt(botf, topf);
		assert_eq(botf - topf, botb - topb);
		// We just advanced either ftabLen characters, or 1 character,
		// depending on whether we used ftab or fchr.
		nextLocsBi(ebwtFw, ebwtBw, tloc, bloc, topf, botf, topb, botb);
		assert(tloc.valid());
	}
	assert(tloc.valid());
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
	TIndexOff t[4], b[4];   // dest BW ranges
	TIndexOff tp[4], bp[4]; // dest BW ranges for "prime" index
	ASSERT_ONLY(TIndexOff lasttot = botf - topf);
	bool fail = false, done = false;
	while(!fail && !done) {
		int rdc = q.getc(off5p, fw_);
		assert_range(0, 3, rdc); // TODO: handle Ns
		assert_gt(botf, topf);
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
		assert(tloc.valid());
		bool ltr = l2r_;
		const Ebwt& ebwt = ltr ? ebwtBw : ebwtFw;
		// Set up initial values for the primes
		if(ltr) {
			tp[0] = tp[1] = tp[2] = tp[3] = topf;
			bp[0] = bp[1] = bp[2] = bp[3] = botf;
		} else {
			tp[0] = tp[1] = tp[2] = tp[3] = topb;
			bp[0] = bp[1] = bp[2] = bp[3] = botb;
		}
		t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
		if(bloc.valid()) {
			// Range delimited by tloc/bloc has size >1.  If size == 1,
			// we use a simpler query (see if(!bloc.valid()) blocks below)
			met.bwops++;
			met.bwops_bi++;
			ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
			// t, b, tp and bp now filled
			ASSERT_ONLY(TIndexOff tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3]));
			ASSERT_ONLY(TIndexOff totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3]));
			assert_eq(tot, totp);
			assert_leq(tot, lasttot);
			ASSERT_ONLY(lasttot = tot);
			fail = (b[rdc] <= t[rdc]);
		} else {
			// Range delimited by tloc/bloc has size 1
			TIndexOff ntop = ltr ? topb : topf;
			met.bwops++;
			met.bwops_1++;
			int cc = ebwt.mapLF1(ntop, tloc);
			assert_range(-1, 3, cc);
			if(cc >= 0) {
				t[cc] = ntop; b[cc] = ntop+1;
			} else {
				fail = true;
			}
		}
		// failed?
		// Convert t, tp info tf, bf
		TIndexOff *tf = ltr ? tp : t, *bf = ltr ? bp : b;
		TIndexOff *tb = ltr ? t : tp, *bb = ltr ? b : bp;
		nextLocsBi(ebwtFw, ebwtBw, tloc, bloc, tf[rdc], bf[rdc], tb[rdc], bb[rdc]);
		// Update off5p
		if(toward3p) {
			assert_gt(off3p, 0);
			off5p++;
			off3p--;
			done = (off3p == 0);
		} else {
			assert_gt(off5p, 0);
			off5p--;
			off3p++;
			done = (off5p == 0);
		}
		// Allocate DescentPos data structure.
		if(firstPos) {
			posid_ = pf.alloc();
		} else {
			pf.alloc();
		}
		pf[posid_ + len_].reset();
		for(size_t i = 0; i < 4; i++) {
			pf[posid_ + len_].topf[i] = l2r_ ? tp[i] : t[i];
			pf[posid_ + len_].botf[i] = l2r_ ? bp[i] : b[i];
			pf[posid_ + len_].topb[i] = l2r_ ? t[i] : tp[i];
			pf[posid_ + len_].botb[i] = l2r_ ? b[i] : bp[i];
		}
		len_++;
	}
	if(done) {
		// Create new DescentRoot
	}
	assert(repOk(&q));
	return len_ > nobranchDepth;
}
