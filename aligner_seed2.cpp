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

#include <ctype.h>
#include "aligner_seed2.h"
#include "assert_helpers.h"
#include "bt2_idx.h"

/**
 * Drive the process of descending from all search roots.
 */
void DescentDriver::go(
    const Scoring& sc,    // scoring scheme
    const Ebwt& ebwtFw,   // forward index
    const Ebwt& ebwtBw,   // mirror index
    DescentMetrics& met)  // metrics
{
    // Convert DescentRoots to the initial Descents
    for(size_t i = 0; i < roots_.size(); i++) {
        size_t dfsz = df_.size();
        size_t pfsz = pf_.size();
        TDescentId id = df_.alloc();
        Edit e_null;
        assert(!e_null.inited());
        bool succ = df_[id].init(
            q_,        // query
            i,         // root and conf id
            sc,        // scoring scheme
            id,        // new Descent's id
            ebwtFw,    // forward index
            ebwtBw,    // mirror index
            df_,       // Descent factory
            pf_,       // DescentPos factory
            roots_,    // DescentRoots
            confs_,    // DescentConfs
            heap_,     // heap
            alsink_,   // alignment sink
            met);      // metrics
        if(!succ) {
            // Reclaim memory we had used for this descent and its DescentPos info
            df_.resize(dfsz);
            pf_.resize(pfsz);
        }
    }
    // Advance until some stopping condition
    bool stop = heap_.empty();
    while(!stop) {
        TDescentPair p = heap_.pop();
        df_[p.second].followBestOutgoing(
            q_,
            ebwtFw,
            ebwtBw,
            sc,
            pf_,
            df_,
            roots_,
            confs_,
            heap_,
            alsink_,
            met);
        stop = heap_.empty();
    }
}

/**
 * Initialize a new descent branching from the given descent via the given
 * edit.  Return false if the Descent has no outgoing edges (and can
 * therefore have its memory freed), true otherwise.
 */
bool Descent::init(
    const DescentQuery& q,          // query
    TRootId rid,                    // root id
    const Scoring& sc,              // scoring scheme
    TReadOff al5pi,                 // offset from 5' of 1st aligned char
    TReadOff al5pf,                 // offset from 5' of last aligned char
    TIndexOff topf,                 // SA range top in FW index
    TIndexOff botf,                 // SA range bottom in FW index
    TIndexOff topb,                 // SA range top in BW index
    TIndexOff botb,                 // SA range bottom in BW index
    bool l2r,                       // direction this descent will go in
    size_t descid,                  // my ID
    TDescentId parent,              // parent ID
    TScore pen,                     // total penalties so far
    const Edit& e,                  // edit for incoming edge
    const Ebwt& ebwtFw,             // forward index
    const Ebwt& ebwtBw,             // mirror index
    EFactory<Descent>& df,          // Descent factory
    EFactory<DescentPos>& pf,       // DescentPos factory
    const EList<DescentRoot>& rs,   // roots
    const EList<DescentConfig>& cs, // configs
    EHeap<TDescentPair>& heap,      // heap
    DescentAlignmentSink& alsink,   // alignment sink
    DescentMetrics& met)            // metrics
{
    rid_ = rid;
    al5pi_ = al5pi;
    al5pf_ = al5pf;
    l2r_ = l2r;
    topf_ = topf;
    botf_ = botf;
    topb_ = topb;
    botb_ = botb;
    descid_ = descid;
    parent_ = parent;
    pen_ = pen;
    posid_ = std::numeric_limits<size_t>::max();
    len_ = 0;
    out_.clear();
    edit_ = e;
    lastRecalc_ = true;
    bool branches = false, hitEnd = false, done = false;
    TIndexOff topf_new = 0, botf_new = 0, topb_new = 0, botb_new = 0;
    off5p_i_ = 0;
#ifndef NDEBUG
    size_t depth = al5pf_ - al5pi_ + 1;
    TScore maxpen = cs[rid_].cons[depth]; // maximum penalty total at this position
    assert_geq(maxpen, pen_);    // can't have already exceeded max penalty
#endif
    followMatches(
        q,
        ebwtFw,
        ebwtBw,
        df,
        pf,
        rs,
        cs,
        heap,
        alsink,
        met,
        branches,
        hitEnd,
        done,
        off5p_i_,
        topf_new,
        botf_new,
        topb_new,
        botb_new);
    bool bounceSucc = false;
    if(hitEnd && !done) {
        bounceSucc = bounce(
            q,
            topf_new,
            botf_new,
            topb_new,
            botb_new,
            ebwtFw,
            ebwtBw,
            sc,
            pf,
            df,
            rs,
            cs,
            heap,
            alsink,
            met);
    }
    // Calculate info about outgoing edges
    recalcOutgoing(q, sc, pf, rs, cs);
    if(!empty()) {
		heap.insert(make_pair(out_.bestPri(), descid)); // Add to heap
    }
#ifndef NDEBUG
    toStacked(cerr, q, ebwtFw, ebwtBw, pf, df, rs, cs);
#endif
    return !empty() || bounceSucc;
}

/**
 * Initialize a new descent beginning at the given root.  Return false if
 * the Descent has no outgoing edges (and can therefore have its memory
 * freed), true otherwise.
 */
bool Descent::init(
    const DescentQuery& q,          // query
    TRootId rid,                    // root id
    const Scoring& sc,              // scoring scheme
    size_t descid,                  // id of this Descent
    const Ebwt& ebwtFw,             // forward index
    const Ebwt& ebwtBw,             // mirror index
    EFactory<Descent>& df,          // Descent factory
    EFactory<DescentPos>& pf,       // DescentPos factory
    const EList<DescentRoot>& rs,   // roots
    const EList<DescentConfig>& cs, // configs
    EHeap<TDescentPair>& heap,      // heap
    DescentAlignmentSink& alsink,   // alignment sink
    DescentMetrics& met)            // metrics
{
    rid_ = rid;
    al5pi_ = rs[rid].off5p;
    al5pf_ = rs[rid].off5p;
    l2r_ = rs[rid].l2r;
    topf_ = botf_ = topb_ = botb_ = 0;
    descid_ = descid;
    parent_ = std::numeric_limits<size_t>::max();
    pen_ = 0;
    posid_ = std::numeric_limits<size_t>::max();
    len_ = 0;
    out_.clear();
    edit_.reset();
    lastRecalc_ = true;
    bool branches = false, hitEnd = false, done = false;
    TIndexOff topf_new = 0, botf_new = 0, topb_new = 0, botb_new = 0;
    off5p_i_ = 0;
    followMatches(
        q,
        ebwtFw,
        ebwtBw,
        df,
        pf,
        rs,
        cs,
        heap,
        alsink,
        met,
        branches,
        hitEnd,
        done,
        off5p_i_,
        topf_new,
        botf_new,
        topb_new,
        botb_new);
    bool bounceSucc = false;
    if(hitEnd && !done) {
        bounceSucc = bounce(
            q,
            topf_new,
            botf_new,
            topb_new,
            botb_new,
            ebwtFw,
            ebwtBw,
            sc,
            pf,
            df,
            rs,
            cs,
            heap,
            alsink,
            met);
    }
    // Calculate info about outgoing edges
    assert(empty());
    recalcOutgoing(q, sc, pf, rs, cs);
    if(!empty()) {
		heap.insert(make_pair(out_.bestPri(), descid)); // Add to heap
    }
#ifndef NDEBUG
    toStacked(cerr, q, ebwtFw, ebwtBw, pf, df, rs, cs);
#endif
    return !empty() || bounceSucc;
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
 * TODO: Eliminate outgoing gap edges that are redundant with others owing to
 *       the DNA sequence and the fact that we don't care to distinguish among
 *       "equivalent" homopolymer extensinos and retractions.
 */
size_t Descent::recalcOutgoing(
    const DescentQuery& q,           // query string
    const Scoring& sc,               // scoring scheme
    EFactory<DescentPos>& pf,        // factory with DescentPoss
    const EList<DescentRoot>& rs,    // roots
    const EList<DescentConfig>& cs)  // configs
{
    assert_eq(botf_ - topf_, botb_ - topb_);
	assert(out_.empty());
	assert(repOk(&q));
	// Get initial 5' and 3' offsets
    bool fw = rs[rid_].fw;
    float rootpri = rs[rid_].pri;
	bool toward3p = (l2r_ == fw);
	size_t off5p = off5p_i_;
	assert_geq(al5pf_, al5pi_);
    assert_geq(off5p, al5pi_);
	size_t off3p = q.length() - off5p - 1;
	// By "depth" we essentially mean the number of characters already aligned
	size_t depth;
    if(toward3p) {
        depth = off5p - al5pi_;
    } else {
        depth = al5pf_ - off5p;
    }
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
    size_t d = posid_;
	while(off5p >= al5pi_ && off5p <= al5pf_) {
        assert_lt(off5p, q.length());
        assert_lt(off3p, q.length());
		TScore maxpen = cs[rid_].cons[depth]; // maximum penalty total at this position
		assert_geq(maxpen, pen_);    // can't have already exceeded max penalty
		TScore diff = maxpen - pen_; // room we have left
		// Get pointer to SA ranges in the direction of descent
		const TIndexOff *t  = l2r_ ? pf[d].topb : pf[d].topf;
		const TIndexOff *b  = l2r_ ? pf[d].botb : pf[d].botf;
		const TIndexOff *tp = l2r_ ? pf[d].topf : pf[d].topb;
		const TIndexOff *bp = l2r_ ? pf[d].botf : pf[d].botb;
		// What are the read char / quality?
		std::pair<int, int> p = q.get(off5p, fw);
		int c = p.first;
		assert_range(0, 4, c);
		// Only entertain edits if there is at least one type of edit left and
		// there is some penalty budget left
		if(!pf[d].flags.exhausted() && diff > 0) {
			// What would the penalty be if we mismatched at this position?
			// This includes the case where the mismatch is for an N in the
			// read.
			int qq = p.second;
            assert_geq(qq, 0);
			TScore pen_mm = sc.mm(c, qq);
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
					Edit edit((uint32_t)off5p, (int)("ACGTN"[j]), (int)("ACGTN"[c]), EDIT_TYPE_MM);
					DescentPriority pri(pen_ + pen_mm, depth, width, rootpri);
					TIndexOff topf = pf[d].topf[j], botf = pf[d].botf[j];
					TIndexOff topb = pf[d].topb[j], botb = pf[d].botb[j];
                    assert(topf != 0 || botf != 0);
                    assert(topb != 0 || botb != 0);
					assert_eq(botb - topb, botf - topf);
					edge.init(edit, pri, d
#ifndef NDEBUG
                    , d, topf, botf, topb, botb
#endif
                    );
					out_.update(edge);
					nout++;
				}
			}
			bool gapsAllowed = (off5p >= sc.gapbar && off3p >= sc.gapbar);
			if(gapsAllowed) {
				bool rdex = false, rfex = false;
				if(off5p == off5p_i_) {
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
							Edit edit((uint32_t)off5p, (int)("ACGTN"[j]), '-', EDIT_TYPE_READ_GAP);
							DescentPriority pri(pen_ + pen_rdg_ex, depth, width, rootpri);
							TIndexOff topf = pf[d].topf[j], botf = pf[d].botf[j];
							TIndexOff topb = pf[d].topb[j], botb = pf[d].botb[j];
                            assert(topf != 0 || botf != 0);
                            assert(topb != 0 || botb != 0);
							assert_eq(botb - topb, botf - topf);
							edge.init(edit, pri, d
#ifndef NDEBUG
                            , d,
                            topf, botf, topb, botb
#endif
                            );
							out_.update(edge);
							nout++;
						}
					}
					if(pen_rfg_ex <= diff && edit_.isRefGap()) {
						// Extension of a reference gap
						rfex = true;
						if(pf[d].flags.rfgExplore()) {
                            TIndexOff width = bot - top;
                            Edit edit((uint32_t)off5p, '-', (int)("ACGTN"[c]), EDIT_TYPE_REF_GAP);
                            DescentPriority pri(pen_ + pen_rfg_ex, depth, width, rootpri);
                            TIndexOff topf = l2r_ ? topp : top;
                            TIndexOff topb = l2r_ ? top : topp;
                            TIndexOff botf = l2r_ ? botp : bot;
                            TIndexOff botb = l2r_ ? bot : botp;
                            assert(topf != 0 || botf != 0);
                            assert(topb != 0 || botb != 0);
                            edge.init(edit, pri, d
#ifndef NDEBUG
                            // It's a little unclear what the depth ought to be.
                            // Is it the depth we were at when we did the ref
                            // gap?  I.e. the depth of the flags where rfgExplore()
                            // returned true?  Or is it the depth where we can
                            // retrieve the appropriate top/bot?  We make it the
                            // latter, might wrap around, indicating we should get
                            // top/bot from the descent's topf_, ... fields.
                            , (d == posid_) ? std::numeric_limits<size_t>::max() : (d - 1),
                            topf, botf, topb, botb
#endif
                            );
                            out_.update(edge);
                            nout++;
						}
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
						Edit edit((uint32_t)off5p, (int)("ACGTN"[j]), '-', EDIT_TYPE_READ_GAP);
						DescentPriority pri(pen_ + pen_rdg_op, depth, width, rootpri);
						TIndexOff topf = pf[d].topf[j], botf = pf[d].botf[j];
						TIndexOff topb = pf[d].topb[j], botb = pf[d].botb[j];
                        assert(topf != 0 || botf != 0);
                        assert(topb != 0 || botb != 0);
						assert_eq(botb - topb, botf - topf);
						edge.init(edit, pri, d
#ifndef NDEBUG
                        , d, topf, botf, topb, botb
#endif
                        );
						out_.update(edge);
						nout++;
					}
				}
				if(pen_rfg_op <= diff && !rfex) {
					// Opening a new reference gap
                    if(pf[d].flags.rfgExplore()) {
                        TIndexOff width = bot - top;
                        Edit edit((uint32_t)off5p, '-', (int)("ACGTN"[c]), EDIT_TYPE_REF_GAP);
                        DescentPriority pri(pen_ + pen_rfg_op, depth, width, rootpri);
                        TIndexOff topf = l2r_ ? topp : top;
                        TIndexOff topb = l2r_ ? top : topp;
                        TIndexOff botf = l2r_ ? botp : bot;
                        TIndexOff botb = l2r_ ? bot : botp;
                        assert(topf != 0 || botf != 0);
                        assert(topb != 0 || botb != 0);
                        edge.init(edit, pri, d
#ifndef NDEBUG
                        // It's a little unclear what the depth ought to be.
                        // Is it the depth we were at when we did the ref
                        // gap?  I.e. the depth of the flags where rfgExplore()
                        // returned true?  Or is it the depth where we can
                        // retrieve the appropriate top/bot?  We make it the
                        // latter, might wrap around, indicating we should get
                        // top/bot from the descent's topf_, ... fields.
                        , (d == posid_) ? std::numeric_limits<size_t>::max() : (d - 1),
                        topf, botf, topb, botb
#endif
                        );
                        out_.update(edge);
                        nout++;
                    }
				}
			}
		}
		// Update off5p, off3p, depth
        d++;
		depth++;
        assert_leq(depth, al5pf_ - al5pi_ + 2);
        if(toward3p) {
            if(off3p == 0) {
                break;
            }
            off5p++;
            off3p--;
        } else {
            if(off5p == 0) {
                break;
            }
            off3p++;
            off5p--;
        }
		// Update top and bot
		top = t[c], topp = tp[c];
		bot = b[c], botp = bp[c];
	}
	lastRecalc_ = (nout <= 5);
    out_.best1.updateFlags(pf);
    out_.best2.updateFlags(pf);
    out_.best3.updateFlags(pf);
    out_.best4.updateFlags(pf);
    out_.best5.updateFlags(pf);
	return nout;
}

void Descent::print(
	std::ostream& os,
	const char *prefix,
	const DescentQuery& q,
	size_t trimLf,
	size_t trimRg,
	bool fw,
	const EList<Edit>& edits,
	BTDnaString& rf)
{
	const BTDnaString& read = fw ? *q.seq : *q.seqrc;
	size_t eidx = 0;
	os << prefix;
	// Print read
	for(size_t i = 0; i < read.length(); i++) {
		if(i < trimLf || i >= read.length() - trimRg) {
			os << (char)tolower(read.toChar(i));
			continue;
		}
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				ins = true;
				os << '-';
			} else if(edits[eidx].isRefGap()) {
				del = true;
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				os << read.toChar(i);
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				assert_eq((int)edits[eidx].qchr, read.toChar(i));
				os << (char)edits[eidx].qchr;
			}
			eidx++;
		}
		if(!del && !mm) {
			// Print read character
			os << read.toChar(i);
		}
	}
	os << endl;
	os << prefix;
	eidx = 0;
	// Print match bars
	for(size_t i = 0; i < read.length(); i++) {
		if(i < trimLf || i >= read.length() - trimRg) {
			os << ' ';
			continue;
		}
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				ins = true;
				os << ' ';
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << ' ';
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				os << ' ';
			}
			eidx++;
		}
		if(!del && !mm) os << '|';
	}
	os << endl;
	os << prefix;
	eidx = 0;
	// Print reference
	for(size_t i = 0; i < read.length(); i++) {
		if(i < trimLf || i >= read.length() - trimRg) {
			os << ' ';
			continue;
		}
		bool ins = false, del = false, mm = false;
		while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
				ins = true;
				rf.appendChar((char)edits[eidx].chr);
				os << (char)edits[eidx].chr;
			} else if(edits[eidx].isRefGap()) {
				del = true;
				os << '-';
			} else {
				mm = true;
				assert(edits[eidx].isMismatch());
				rf.appendChar((char)edits[eidx].chr);
				os << (char)edits[eidx].chr;
			}
			eidx++;
		}
		if(!del && !mm) {
			rf.append(read[i]);
			os << read.toChar(i);
		}
	}
	os << endl;
}

/**
 * Create a new Descent 
 */
bool Descent::bounce(
	const DescentQuery& q,          // query string
    TIndexOff topf,                 // SA range top in fw index
    TIndexOff botf,                 // SA range bottom in fw index
    TIndexOff topb,                 // SA range top in bw index
    TIndexOff botb,                 // SA range bottom in bw index
	const Ebwt& ebwtFw,             // forward index
	const Ebwt& ebwtBw,             // mirror index
	const Scoring& sc,              // scoring scheme
	EFactory<DescentPos>& pf,       // factory with DescentPoss
	EFactory<Descent>& df,          // factory with Descent
    const EList<DescentRoot>& rs,   // roots
    const EList<DescentConfig>& cs, // configs
	EHeap<TDescentPair>& heap,      // heap of descents
    DescentAlignmentSink& alsink,   // alignment sink
    DescentMetrics& met)            // metrics
{
    assert_gt(botf, topf);
    assert(al5pi_ == 0 || al5pf_ == q.length()-1);
    assert(!(al5pi_ == 0 && al5pf_ == q.length()-1));
    size_t dfsz = df.size();
    size_t pfsz = pf.size();
	TDescentId id = df.alloc();
    Edit e_null;
    assert(!e_null.inited());
	bool succ = df[id].init(
		q,         // query
        rid_,      // root id
		sc,        // scoring scheme
		al5pi_,    // new near-5' extreme
		al5pf_,    // new far-5' extreme
		topf,      // SA range top in FW index
		botf,      // SA range bottom in FW index
		topb,      // SA range top in BW index
		botb,      // SA range bottom in BW index
		!l2r_,     // direction this descent will go in; opposite from parent
		id,        // my ID
		descid_,   // parent ID
		pen_,      // total penalties so far - same as parent
		e_null,    // edit for incoming edge; uninitialized if bounced
		ebwtFw,    // forward index
		ebwtBw,    // mirror index
		df,        // Descent factory
		pf,        // DescentPos factory
        rs,        // DescentRoot list
        cs,        // DescentConfig list
		heap,      // heap
        alsink,    // alignment sink
		met);      // metrics
    if(!succ) {
        // Reclaim memory we had used for this descent and its DescentPos info
        df.resize(dfsz);
        pf.resize(pfsz);
    }
    return succ;
}

/**
 * Take the best outgoing edge and place it in the heap.  When deciding what
 * outgoing edges exist, abide by constraints in DescentConfig.  These
 * constraints limit total penalty accumulated so far versus distance from
 * search root.  E.g. a constraint might disallow any gaps or mismatches within
 * 20 ply of the root, then allow 1 mismatch within 30 ply, 1 mismatch or 1 gap
 * within 40 ply, etc.
 */
void Descent::followBestOutgoing(
	const DescentQuery& q,          // query string
	const Ebwt& ebwtFw,             // forward index
	const Ebwt& ebwtBw,             // mirror index
	const Scoring& sc,              // scoring scheme
	EFactory<DescentPos>& pf,       // factory with DescentPoss
	EFactory<Descent>& df,          // factory with Descent
    const EList<DescentRoot>& rs,   // roots
    const EList<DescentConfig>& cs, // configs
	EHeap<TDescentPair>& heap,      // heap of descents
    DescentAlignmentSink& alsink,   // alignment sink
	DescentMetrics& met)            // metrics
{
	// We assume this descent has been popped off the heap.  We'll re-add it if
	// it hasn't been exhausted.
	assert(!empty());
	assert(!out_.empty());
	DescentPriority best = out_.bestPri();
	DescentEdge e = out_.rotate();
	if(out_.empty()) {
		if(!lastRecalc_) {
			recalcOutgoing(q, sc, pf, rs, cs);
			assert(!empty());
		} else {
			assert(empty());
		}
	}
	if(!empty()) {
		// Re-insert this
        assert(!(out_.bestPri() < best));
		heap.insert(make_pair(out_.bestPri(), descid_));
	}
	// Allocate a new Descent object
    TReadOff edoff = e.e.pos; // 5' offset of edit
    assert_leq(edoff, al5pf_);
    assert_geq(edoff, al5pi_);
	TReadOff al5pi_new = al5pi_, al5pf_new = al5pf_;
    bool fw = rs[rid_].fw;
    bool toward3p = (l2r_ == fw);
    TReadOff doff; // edit's offset into this descent
    int chr = asc2dna[e.e.chr];
    // hitEnd is set to true iff this edit pushes us to the extreme 5' or 3'
    // end of the alignment
    bool hitEnd = false;
    // done is set to true iff this edit aligns the only remaining character of
    // the read
    bool done = false;
    if(toward3p) {
        // The 3' extreme of the new Descent is further in (away from the 3'
        // end) than the parent's.
        al5pf_new = doff = edoff;
        if(e.e.isReadGap()) {
            // We didn't actually consume the read character at 'edoff', so
            // retract al5pf_new by one position.  This doesn't effect the
            // "depth" (doff) of the SA range we took, though.
            assert_gt(al5pf_new, 0);
            al5pf_new--;
        }
        assert_lt(al5pf_new, q.length());
        hitEnd = (al5pf_new == q.length() - 1);
        done = (hitEnd && al5pi_new == 0);
        assert_geq(doff, off5p_i_);
        doff = doff - off5p_i_;
        assert_lt(doff, len_);
    } else {
        // The 5' extreme of the new Descent is further in (away from the 5'
        // end) than the parent's.
        al5pi_new = doff = edoff;
        if(e.e.isReadGap()) {
            // We didn't actually consume the read character at 'edoff', so
            // move al5pi_new closer to the 3' end by one position.  This
            // doesn't effect the "depth" (doff) of the SA range we took,
            // though.
            al5pi_new++;
        }
        hitEnd = (al5pi_new == 0);
        done = (hitEnd && al5pf_new == q.length());
        assert_geq(off5p_i_, doff);
        doff = off5p_i_ - doff;
        assert_lt(doff, len_);
    }
    bool l2r = l2r_; // gets overridden if we bounce
    if(done) {
        // Aligned the entire read end-to-end.  Presumably there's no need to
        // create a new Descent object.  We just report the alignment.
        alsink.reportAlignment(
            q,        // query
            descid_,  // Descent at the leaf
            e.e,      // extra edit, if necessary
            best.pen, // penalty
            df,       // factory with Descent
            pf,       // factory with DescentPoss
            rs,       // roots
            cs);      // configs
        return;
    } else if(hitEnd) {
        // Alignment finsihed extending in one direction
        l2r = !l2r;
    }
    size_t dfsz = df.size();
    size_t pfsz = pf.size();
	TDescentId id = df.alloc();
	TIndexOff topf, botf, topb, botb;
    size_t d = posid_ + doff;
    if(e.e.isRefGap()) {
        d--; // might underflow
        if(doff == 0) {
            topf = topf_;
            botf = botf_;
            topb = topb_;
            botb = botb_;
			d = std::numeric_limits<size_t>::max();
        } else {
            assert_gt(al5pf_new, 0);
            assert_gt(d, 0);
            chr = pf[d].c;
            topf = pf[d].topf[chr];
            botf = pf[d].botf[chr];
            topb = pf[d].topb[chr];
            botb = pf[d].botb[chr];
        }
    } else {
        // A read gap or a mismatch
        topf = pf[d].topf[chr];
        botf = pf[d].botf[chr];
        topb = pf[d].topb[chr];
        botb = pf[d].botb[chr];
    }
    assert_eq(d, e.d);
    assert_eq(topf, e.topf);
    assert_eq(botf, e.botf);
    assert_eq(topb, e.topb);
    assert_eq(botb, e.botb);
	bool succ = df[id].init(
		q,         // query
        rid_,      // root id
		sc,        // scoring scheme
		al5pi_new, // new near-5' extreme
		al5pf_new, // new far-5' extreme
		topf,      // SA range top in FW index
		botf,      // SA range bottom in FW index
		topb,      // SA range top in BW index
		botb,      // SA range bottom in BW index
		l2r,       // direction this descent will go in
		id,        // my ID
		descid_,   // parent ID
		best.pen,  // total penalties so far
		e.e,       // edit for incoming edge; uninitialized if bounced
		ebwtFw,    // forward index
		ebwtBw,    // mirror index
		df,        // Descent factory
		pf,        // DescentPos factory
        rs,        // DescentRoot list
        cs,        // DescentConfig list
		heap,      // heap
        alsink,    // alignment sink
		met);      // metrics
    if(!succ) {
        // Reclaim memory we had used for this descent and its DescentPos info
        df.resize(dfsz);
        pf.resize(pfsz);
    }
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
void Descent::followMatches(
	const DescentQuery& q,     // query string
	const Ebwt& ebwtFw,        // forward index
	const Ebwt& ebwtBw,        // mirror index
	EFactory<Descent>& df,     // Descent factory
	EFactory<DescentPos>& pf,  // DescentPos factory
    const EList<DescentRoot>& rs,   // roots
    const EList<DescentConfig>& cs, // configs
	EHeap<TDescentPair>& heap, // heap
    DescentAlignmentSink& alsink, // alignment sink
	DescentMetrics& met,       // metrics
    bool& branches,            // out: true -> there are > 0 ways to branch
    bool& hitEnd,              // out: true -> hit read end with non-empty range
    bool& done,                // out: true -> we made a full alignment
    TReadOff& off5p_i,         // out: initial 5' offset
    TIndexOff& topf_bounce,    // out: top of SA range for fw idx for bounce
    TIndexOff& botf_bounce,    // out: bot of SA range for fw idx for bounce
    TIndexOff& topb_bounce,    // out: top of SA range for bw idx for bounce
    TIndexOff& botb_bounce)    // out: bot of SA range for bw idx for bounce
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
    bool fw = rs[rid_].fw;
	bool toward3p = (l2r_ == fw);
	size_t off5p;
	assert_geq(al5pf_, al5pi_);
	if(toward3p) {
		if(al5pi_ == al5pf_) {
			off5p = off5p_i = 0;
		} else {
			off5p = off5p_i = (al5pf_ + 1);
		}
	} else {
		assert_gt(al5pi_, 0);
		if(al5pi_ == al5pf_) {
			off5p = off5p_i = (q.length() - 1);
		} else {
			off5p = off5p_i = (al5pi_ - 1);
		}
	}
	size_t off3p = q.length() - off5p - 1;
	assert_lt(off5p, q.length());
	assert_lt(off3p, q.length());
	bool firstPos = true;
	assert_eq(0, len_);
    
    // Set to true as soon as we encounter a branch point along this descent.
    branches = false;
    // hitEnd is set to true iff this edit pushes us to the extreme 5' or 3'
    // end of the alignment
    hitEnd = false;
    // done is set to true iff this edit aligns the only remaining character of
    // the read
    done = false;
	if(root()) {
        assert_eq(al5pi_, al5pf_);
		// Check whether/how far we can jump using ftab
		int ftabLen = ebwtFw.eh().ftabChars();
		if(ftabLen > 1 && ftabLen <= nobranchDepth) {
			// Right-to-left
			size_t off_r2l = off5p;
			if(l2r_) {
				off_r2l += ftabLen;
				assert_lt(off_r2l, q.length());
			}
			ebwtFw.ftabLoHi(fw ? *q.seq : *q.seqrc, off_r2l, false, topf, botf);
			if(botf - topf == 0) {
				return;
			}
			// Left-to-right
			size_t off_l2r = off5p;
			if(!l2r_) {
				assert_geq(off_l2r, ftabLen);
				off_l2r -= ftabLen;
			}
			ebwtBw.ftabLoHi(fw ? *q.seq : *q.seqrc, off_l2r, false, topb, botb);
			assert_eq(botf - topf, botb - topb);
			if(toward3p) {
				assert_geq(off3p, ftabLen);
				off5p += ftabLen; off3p -= ftabLen;
			} else {
				assert_geq(off5p, ftabLen);
				off5p -= ftabLen; off3p += ftabLen;
			}
			len_ += ftabLen;
            if(toward3p) {
                al5pf_ += (ftabLen - 1); // -1 accounts for inclusive al5pf_
                if(al5pf_ == q.length() - 1) {
                    hitEnd = true;
                    done = (al5pi_ == 0);
                }
            } else {
                assert_geq(al5pi_, ftabLen);
                al5pi_ -= ftabLen;
                if(al5pi_ == 0) {
                    hitEnd = true;
                    done = (al5pf_ == q.length()-1);
                }
            }
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
			int rdc = q.getc(off5p, fw);
			assert_range(0, 3, rdc);
			topf = topb = ebwtFw.fchr()[rdc];
			botf = botb = ebwtFw.fchr()[rdc+1];
			if(botf - topf == 0) {
				return;
			}
			if(toward3p) {
				assert_gt(off3p, 0);
				off5p++; off3p--;
			} else {
				assert_gt(off5p, 0);
				off5p--; off3p++;
			}
			len_++;
            if(toward3p) {
                if(al5pf_ == q.length()-1) {
                    hitEnd = true;
                    done = (al5pi_ == 0);
                }
            } else {
                if(al5pi_ == 0) {
                    hitEnd = true;
                    done = (al5pf_ == q.length()-1);
                }
            }
			// Allocate DescentPos data structure.  We could fill it with the
			// four ranges from fchr if we wanted to, but that will never be
			// relevant.
			if(firstPos) {
				posid_ = pf.alloc();
                firstPos = false;
			} else {
				pf.alloc();
			}
		}
		assert_gt(botf, topf);
		assert_eq(botf - topf, botb - topb);
	}
    if(done) {
        Edit eempty;
        alsink.reportAlignment(
            q,        // query
            descid_,  // Descent at the leaf
            eempty,   // extra edit, if necessary
            pen_,     // penalty
            df,       // factory with Descent
            pf,       // factory with DescentPoss
            rs,       // roots
            cs);      // configs
        return;
    } else if(hitEnd) {
        assert_gt(botf, topf);
        topf_bounce = topf;
        botf_bounce = botf;
        topb_bounce = topb;
        botb_bounce = botb;
        return; // Bounced
    }
    // We just advanced either ftabLen characters, or 1 character,
    // depending on whether we used ftab or fchr.
    nextLocsBi(ebwtFw, ebwtBw, tloc, bloc, topf, botf, topb, botb);
    assert(tloc.valid());
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
	TIndexOff t[4], b[4];   // dest BW ranges
	TIndexOff tp[4], bp[4]; // dest BW ranges for "prime" index
	ASSERT_ONLY(TIndexOff lasttot = botf - topf);
	bool fail = false;
	while(!fail && !hitEnd) {
        assert(!done);
		int rdc = q.getc(off5p, fw);
		assert_range(0, 3, rdc); // TODO: handle Ns
		assert_gt(botf, topf);
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
		assert(tloc.valid());
        TIndexOff width = botf - topf;
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
            if(b[rdc] - t[rdc] < width) {
                branches = true;
            }
		} else {
			// Range delimited by tloc/bloc has size 1
			TIndexOff ntop = ltr ? topb : topf;
			met.bwops++;
			met.bwops_1++;
			int cc = ebwt.mapLF1(ntop, tloc);
			assert_range(-1, 3, cc);
            fail = (cc != rdc);
            if(fail) {
                branches = true;
            }
			if(cc >= 0) {
				t[cc] = ntop; b[cc] = ntop+1;
			}
		}
        topf = ltr ? tp[rdc] : t[rdc];
        botf = ltr ? bp[rdc] : b[rdc];
        topb = ltr ? t[rdc] : tp[rdc];
        botb = ltr ? b[rdc] : bp[rdc];
        // If 'fail' is true, we failed to align this read character.  We still
        // install the SA ranges into the DescentPos and increment len_ in this
        // case.
        
		// Convert t, tp, b, bp info tf, bf, tb, bb
		TIndexOff *tf = ltr ? tp : t, *bf = ltr ? bp : b;
		TIndexOff *tb = ltr ? t : tp, *bb = ltr ? b : bp;
		// Allocate DescentPos data structure.
		if(firstPos) {
			posid_ = pf.alloc();
            firstPos = false;
		} else {
			pf.alloc();
		}
		pf[posid_ + len_].reset();
        pf[posid_ + len_].c = rdc;
		for(size_t i = 0; i < 4; i++) {
			pf[posid_ + len_].topf[i] = l2r_ ? tp[i] : t[i];
			pf[posid_ + len_].botf[i] = l2r_ ? bp[i] : b[i];
			pf[posid_ + len_].topb[i] = l2r_ ? t[i] : tp[i];
			pf[posid_ + len_].botb[i] = l2r_ ? b[i] : bp[i];
		}
		if(!fail) {
			len_++;
			if(toward3p) {
				al5pf_++;
				off5p++;
				off3p--;
				if(al5pf_ == q.length() - 1) {
					hitEnd = true;
					done = (al5pi_ == 0);
				}
			} else {
				assert_gt(al5pi_, 0);
				al5pi_--;
				off5p--;
				off3p++;
				if(al5pi_ == 0) {
					hitEnd = true;
					done = (al5pf_ == q.length() - 1);
				}
			}
		}
        if(!fail && !hitEnd) {
            nextLocsBi(ebwtFw, ebwtBw, tloc, bloc, tf[rdc], bf[rdc], tb[rdc], bb[rdc]);
        }
	}
    if(done) {
        Edit eempty;
        alsink.reportAlignment(
            q,        // query
            descid_,  // Descent at the leaf
            eempty,   // extra edit, if necessary
            pen_,     // total penalty
            df,       // factory with Descent
            pf,       // factory with DescentPoss
            rs,       // roots
            cs);      // configs
        return;
    } else if(hitEnd) {
        assert_gt(botf, topf);
        topf_bounce = topf;
        botf_bounce = botf;
        topb_bounce = topb;
        botb_bounce = botb;
        return; // Bounced
    }
    assert(repOk(&q));
}

#ifdef ALIGNER_SEED2_MAIN

#include <string>
#include "sstring.h"

using namespace std;

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {

    EList<string> strs;
    //                            GCTATATAGCGCGCTCGCATCATTTTGTGT
    strs.push_back(string("CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"
                          "NNNNNNNNNN"
                          "CATGTCAGCTATATAGCGCGCTCGCATCATTTTGTGTGTAAACCA"));
    //                            GCTATATAGCGCGCTTGCATCATTTTGTGT
    //                                           ^
    bool packed = false;
    int color = 0;
	pair<Ebwt*, Ebwt*> ebwts = Ebwt::fromStrings<SString<char> >(
		strs,
		packed,
		color,
		REF_READ_REVERSE,
		Ebwt::default_bigEndian,
		Ebwt::default_lineRate,
		Ebwt::default_offRate,
		Ebwt::default_ftabChars,
		".aligner_seed2.cpp.tmp",
		Ebwt::default_useBlockwise,
		Ebwt::default_bmax,
		Ebwt::default_bmaxMultSqrt,
		Ebwt::default_bmaxDivN,
		Ebwt::default_dcv,
		Ebwt::default_seed,
		false,  // verbose
		false,  // autoMem
		false); // sanity
    
    ebwts.first->loadIntoMemory (color, -1, true, true, true, true, false);
    ebwts.second->loadIntoMemory(color,  1, true, true, true, true, false);
    
    {
        DescentMetrics mets;
        DescentDriver dr;
        
        // Set up the read
        BTDnaString seq ("GCTATATAGCGCGCTCGCATCATTTTGTGT", true);
        BTString    qual("ABCDEFGHIabcdefghiABCDEFGHIabc");
        BTDnaString seqrc = seq;
        BTString    qualrc = qual;
        seqrc.reverseComp();
        qualrc.reverse();
        dr.initRead(
            seq,
            qual,
            seqrc,
            qualrc);
        
        // Set up the DescentConfig
        DescentConfig conf;
        conf.cons.init(SIMPLE_FUNC_LINEAR, 0.0, 1.0);
        conf.expol = DESC_EX_NONE;
        
        // Set up the search roots
        dr.addRoot(
            conf,   // DescentConfig
            0,      // 5' offset into read of root
            true,   // left-to-right?
            true,   // forward?
            0.0);   // root priority
        
        // Do the search
        Scoring sc = Scoring::base1();
        dr.go(sc, *ebwts.first, *ebwts.second, mets);
    }
    
    delete ebwts.first;
    delete ebwts.second;
}
#endif

