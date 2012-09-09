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

#ifndef ALIGNER_SEED2_H_
#define ALIGNER_SEED2_H_

/**
 * The user of the DescentDriver class specifies a collection of search roots.
 * Logic for picking these search roots is located elsewhere, not in this
 * module.  The search roots are annotated with a priority score, which 
 *
 * The heap is a min-heap over pairs, where the first element of each pair is
 * the score associated with a descent and the second element of each pair is
 * the descent ID.
 */

#include <stdint.h>
#include <utility>
#include "bt2_idx.h"
#include "simple_func.h"
#include "scoring.h"
#include "edit.h"
#include "ds.h"

typedef uint32_t TIndexOff;
typedef size_t   TReadOff;
typedef int64_t  TScore;
typedef float    TRootPri;
typedef size_t   TDescentId;
typedef size_t   TRootId;

/**
 * enum encapsulating a few different policies for how we might extend descents
 * in the direction opposite from their primary direction.  
 */
enum {
	// Never extened in the direction opposite from the primary.  Just go in
	// the primary direction until the bounce.
	DESC_EX_NONE = 1,
	
	// When we're finished extending out the matches for a descent, try to
	// extend in the opposite direction in a way that extends all branches
	// simultaneously.  The Descent.nex_ field contains the number of positions
	// we were able to extend through in this way.
	DESC_EX_FROM_1ST_BRANCH = 2,
	
	// Each time we add an edge to the summary, extend it in the opposite
	// direction.  The DescentEdge.nex field contains the number of positions
	// we were able to extend through, and this in turn gets propagated to
	// Descent.nex_ if and when we branch from the DescentEdge.
	DESC_EX_EACH_EDGE = 3
};

/**
 * Counters to keep track of how much work is being done.
 */
struct DescentMetrics {

	DescentMetrics() { reset(); }

	void reset() {
		bwops = bwops_1 = bwops_bi = recalc = branch = branch_mm =
		branch_del = branch_ins = heap_max = descent_max = descentpos_max = 
		nex = 0;
	}

	uint64_t bwops;          // # FM Index opbs
	uint64_t bwops_1;        // # LF1 FM Index opbs
	uint64_t bwops_bi;       // # BiEx FM Index opbs
	uint64_t recalc;         // # times outgoing edge summary was recalculated
	uint64_t branch;         // # times we descended from another descent
	uint64_t branch_mm;      // # times branch was on a mismatch
	uint64_t branch_del;     // # times branch was on a deletion
	uint64_t branch_ins;     // # times branch was on a insertion
	uint64_t heap_max;       // maximum size of Descent heap
	uint64_t descent_max;    // maximum size of Descent factory
	uint64_t descentpos_max; // maximum size of DescentPos factory
	uint64_t nex;            // # extensions
};

/**
 * Priority used to rank which descent we should branch from next.  Right now,
 * priority is governed by a 4-tuple.  From higher to lower priority:
 *
 *  1. Penalty accumulated so far
 *  2. Depth into the search space, including extensions
 *  3. Width of the SA range (i.e. uniqueness)
 *  4. Root priority
 */
struct DescentPriority {

	DescentPriority() { reset(); }

	DescentPriority(
		TScore pen_,
		size_t depth_,
		TIndexOff width_,
		float rootpri_)
	{
		pen = pen_;
		depth = depth_;
		width = width_;
		rootpri = rootpri_;
	}
	
	/**
	 * Initialize new DescentPriority.
	 */
	void init(TScore pen_, size_t depth_, TIndexOff width_, float rootpri_) {
		pen = pen_;
		depth = depth_;
		width = width_;
		rootpri = rootpri_;
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		width = 0;
	}
	
	/**
	 * Return true iff DescentPriority is initialized.
	 */
	bool inited() const {
		return width > 0;
	}
	
	/**
	 * Return true iff this priority is prior to given priority.
	 */
	bool operator<(const DescentPriority& o) const {
		assert(inited());
		assert(o.inited());
		// 1st priority: penalty accumulated so far
		if(pen < o.pen) return true;
		if(pen > o.pen) return false;
		// 2nd priority: depth into the search space, including extensions
		if(depth > o.depth) return true;
		if(depth < o.depth) return false;
		// 3rd priority: width of the SA range (i.e. uniqueness)
		if(width < o.width) return true;
		if(width > o.width) return false;
		// 4th priority: root priority
		if(rootpri > o.rootpri) return true;
		return false;
	}

	/**
	 * Return true iff this priority is prior to or equal to given priority.
	 */
	bool operator<=(const DescentPriority& o) const {
		assert(inited());
		assert(o.inited());
		// 1st priority: penalty accumulated so far
		if(pen < o.pen) return true;
		if(pen > o.pen) return false;
		// 2nd priority: depth into the search space, including extensions
		if(depth > o.depth) return true;
		if(depth < o.depth) return false;
		// 3rd priority: width of the SA range (i.e. uniqueness)
		if(width < o.depth) return true;
		if(width > o.width) return false;
		// 4th priority: root priority
		if(rootpri > o.rootpri) return true;
		return true;
	}

	/**
	 * Return true iff this priority is prior to or equal to given priority.
	 */
	bool operator==(const DescentPriority& o) const {
		assert(inited());
		assert(o.inited());
		return pen == o.pen && depth == o.depth && width == o.width && rootpri == o.rootpri;
	}

	TScore pen;      // total penalty accumulated so far
	size_t depth;    // depth from root of descent
	TIndexOff width; // width of the SA range
	float  rootpri;  // priority of the root
};

static inline std::ostream& operator<<(
	std::ostream& os,
	const DescentPriority& o)
{
	os << "[" << o.pen << ", " << o.depth << ", " << o.width << ", " << o.rootpri << "]";
	return os;
}

static inline std::ostream& operator<<(
	std::ostream& os,
	const std::pair<DescentPriority, TDescentId>& o)
{
	os << "{[" << o.first.pen << ", " << o.first.depth << ", "
	   << o.first.width << ", " << o.first.rootpri << "], " << o.second << "}";
	return os;
}

typedef std::pair<DescentPriority, TDescentId> TDescentPair;

/**
 * Encapsulates the constraints limiting which outgoing edges are permitted.
 */
struct DescentConstraints {

	DescentConstraints() { reset(); }
	
	/**
	 * Initialize with new constraint function.
	 */
	DescentConstraints(int type, double C, double L) {
        init(type, C, L);
	}
    
    /**
     * Initialize with given function.
     */
    void init(int type, double C, double L) {
		double mn = 0.0;
		double mx = std::numeric_limits<double>::max();
		f.init(type, mn, mx, C, L);
		scs.resize(1024);
		for(size_t i = 0; i < 1024; i++) {
			scs[i] = f.f<TScore>(i);
		}
    }
    
    /**
     * Reset to uninitialized state.
     */
    void reset() {
        scs.clear();
    }
    
    /**
     * Return true iff the DescentConstraints has been initialized.
     */
    bool inited() const {
        return !scs.empty();
    }
	
	/**
	 * Get the maximum penalty total for depth 'off'.
	 */
	const TScore operator[](TReadOff off) const {
		if(off > scs.size()) {
			size_t oldsz = scs.size();
            EList<TScore>& scs_rw = const_cast<EList<TScore>&>(scs);
			scs_rw.resize((size_t)(off * 1.5 + 0.5));
			for(size_t i = oldsz; i < scs.size(); i++) {
				scs_rw[i] = f.f<TScore>(i);
			}
		}
		return scs[off];
	}

	EList<TScore> scs; // precalculated scores
	SimpleFunc    f;   // maximum penalties w/r/t depth in the search tree
};

/**
 * Encapsulates settings governing how we descent.
 */
struct DescentConfig {

    DescentConfig() { reset(); }
    
    /**
     * Reset the DescentConfig to an uninitialized state.
     */
    void reset() { expol = 0; }
    
    /**
     * Return true iff this DescentConfig is initialized.
     */
    bool inited() const { return expol != 0; }

    DescentConstraints cons; // constraints
    int expol; // extend policy
};

/**
 * A search root.  Consists of an offset from the 5' end read and flags
 * indicating (a) whether we're initially heading left-to-right or
 * right-to-left, and (b) whether we're examining the read or its reverse
 * complement.
 *
 * A root also comes with a priority ("pri") score indicating how promising it
 * is as a root.  Promising roots have long stretches of high-quality,
 * non-repetitive nucleotides in the first several ply of the search tree.
 * Also, roots beginning at the 5' end of the read may receive a higher
 * priority.
 */
struct DescentRoot {

	DescentRoot() { reset(); }

	DescentRoot(size_t off5p_, bool l2r_, bool fw_, size_t len, float pri_) {
		init(off5p_, l2r_, fw_, len, pri_);
	}
	
	/**
	 * Reset this DescentRoot to uninitialized state.
	 */
	void reset() {
		off5p = std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Return true iff this DescentRoot is uninitialized.
	 */
	bool inited() const {
		return off5p == std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Initialize a new descent root.
	 */
	void init(size_t off5p_, bool l2r_, bool fw_, size_t len, float pri_) {
		off5p = off5p_;
		l2r = l2r_;
		fw = fw_;
		pri = pri_;
		assert_lt(off5p, len);
	}

	TReadOff off5p;   // root origin offset, expressed as offset from 5' end
	bool     l2r;     // true -> move in left-to-right direction
	bool     fw;      // true -> work with forward read, false -> revcomp
	float    pri;     // priority of seed
};

/**
 * Set of flags indicating outgoing edges we've tried from a DescentPos.
 */
struct DescentPosFlags {

	DescentPosFlags() { reset(); }
	
	/**
	 * Set all flags to 1, indicating all outgoing edges are yet to be
	 * explored.
	 */
	void reset() {
		mm_a = mm_c = mm_g = mm_t = rdg_a = rdg_c = rdg_g = rdg_t = rfg = 1;
		reserved = 0;
	}
	
	/**
	 * Return true iff all outgoing edges have already been explored.
	 */
	bool exhausted() const {
		return ((uint16_t*)this)[0] == 0;
	}
	
	/**
	 * Return false iff the specified mismatch has already been explored.
	 */
	bool mmExplore(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			return mm_a;
		} else if(c == 1) {
			return mm_c;
		} else if(c == 2) {
			return mm_g;
		} else {
			return mm_t;
		}
	}

	/**
	 * Try to explore a mismatch.  Return false iff it has already been
	 * explored.
	 */
	bool mmSet(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			bool ret = mm_a; mm_a = 0; return ret;
		} else if(c == 1) {
			bool ret = mm_c; mm_c = 0; return ret;
		} else if(c == 2) {
			bool ret = mm_g; mm_g = 0; return ret;
		} else {
			bool ret = mm_t; mm_t = 0; return ret;
		}
	}

	/**
	 * Return false iff specified read gap has already been explored.
	 */
	bool rdgExplore(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			return rdg_a;
		} else if(c == 1) {
			return rdg_c;
		} else if(c == 2) {
			return rdg_g;
		} else {
			return rdg_t;
		}
	}

	/**
	 * Try to explore a read gap.  Return false iff it has already been
	 * explored.
	 */
	bool rdgSet(int c) {
		assert_range(0, 3, c);
		if(c == 0) {
			bool ret = rdg_a; rdg_a = 0; return ret;
		} else if(c == 1) {
			bool ret = rdg_c; rdg_c = 0; return ret;
		} else if(c == 2) {
			bool ret = rdg_g; rdg_g = 0; return ret;
		} else {
			bool ret = rdg_t; rdg_t = 0; return ret;
		}
	}

	/**
	 * Return false iff the reference gap has already been explored.
	 */
	bool rfgExplore() {
		return rfg;
	}

	/**
	 * Try to explore a reference gap.  Return false iff it has already been
	 * explored.
	 */
	bool rfgSet() {
		bool ret = rfg; rfg = 0; return ret;
	}

	uint16_t mm_a     : 1;
	uint16_t mm_c     : 1;
	uint16_t mm_g     : 1;
	uint16_t mm_t     : 1;

	uint16_t rdg_a    : 1;
	uint16_t rdg_c    : 1;
	uint16_t rdg_g    : 1;
	uint16_t rdg_t    : 1;

	uint16_t rfg      : 1;
	
	uint16_t reserved : 7;
};

/**
 * FM Index state associated with a single position in a descent.  For both the
 * forward and backward indexes, it stores the four SA ranges corresponding to
 * the four nucleotides.
 */
struct DescentPos {

	/**
	 * Reset all tops and bots to 0.
	 */
	void reset() {
		topf[0] = topf[1] = topf[2] = topf[3] = 0;
		botf[0] = botf[1] = botf[2] = botf[3] = 0;
		topb[0] = topb[1] = topb[2] = topb[3] = 0;
		botb[0] = botb[1] = botb[2] = botb[3] = 0;
        c = -1;
		flags.reset();
	}
    
    /**
     * Return true iff DescentPos has been initialized.
     */
    bool inited() const {
        return c >= 0;
    }
	
	TIndexOff       topf[4]; // SA range top indexes in fw index
	TIndexOff       botf[4]; // SA range bottom indexes (exclusive) in fw index
	TIndexOff       topb[4]; // SA range top indexes in bw index
	TIndexOff       botb[4]; // SA range bottom indexes (exclusive) in bw index
    char            c;       // read char that would yield match
	DescentPosFlags flags;   // flags 
};

/**
 * Encapsulates an edge outgoing from a descent.
 */
struct DescentEdge {

	DescentEdge() { reset(); }

	DescentEdge(
		Edit e_,
		DescentPriority pri_,
        size_t posFlag_,
		TReadOff nex_
#ifndef NDEBUG
        ,
        size_t d_,
		TIndexOff topf_,
		TIndexOff botf_,
		TIndexOff topb_,
		TIndexOff botb_
#endif
        )
	{
		init(e_, pri_, posFlag_
#ifndef NDEBUG
        , d_, topf_, botf_, topb_, botb_
#endif
        );
	}

	/**
	 * Return true iff edge is initialized.
	 */
	bool inited() const { return e.inited(); }

	/**
	 * Reset to uninitialized state.
	 */
	void reset() { e.reset(); }
	
	/**
	 * Initialize DescentEdge given 5' offset, nucleotide, and priority.
	 */
	void init(
		Edit e_,
		DescentPriority pri_,
        size_t posFlag_
#ifndef NDEBUG
        ,
        size_t d_,
		TIndexOff topf_,
		TIndexOff botf_,
		TIndexOff topb_,
		TIndexOff botb_
#endif
        )
	{
		e = e_;
		pri = pri_;
        posFlag = posFlag_;
#ifndef NDEBUG
        d = d_;
		topf = topf_;
		botf = botf_;
		topb = topb_;
		botb = botb_;
#endif
	}

    /**
     * Update flags to show this edge as visited.
     */
    void updateFlags(EFactory<DescentPos>& pf) {
        if(inited()) {
            if(e.isReadGap()) {
                assert_neq('-', e.chr);
                pf[posFlag].flags.rdgSet(asc2dna[e.chr]);
            } else if(e.isRefGap()) {
                pf[posFlag].flags.rfgSet();
            } else {
                assert_neq('-', e.chr);
                pf[posFlag].flags.mmSet(asc2dna[e.chr]);
            }
        }
    }
	
	/**
	 * Return true iff this edge has higher priority than the given edge.
	 */
	bool operator<(const DescentEdge& o) const {
        if(inited() && !o.inited()) {
            return true;
        } else if(!inited()) {
            return false;
        }
		return pri < o.pri;
	}

	DescentPriority pri; // priority of the edge
	TReadOff nex;        // # extends possible from this edge
    size_t posFlag;      // depth of DescentPos where flag should be set


#ifndef NDEBUG
    // This can be recreated by looking at the edit, the paren't descent's
    // len_, al5pi_, al5pf_.  I have it here so we can sanity check.
    size_t d;
	TIndexOff topf, botf, topb, botb;
#endif

	Edit e;
};

/**
 * Encapsulates an incomplete summary of the outgoing edges from a descent.  We
 * don't try to store information about all outgoing edges, because doing so
 * will generally be wasteful.  We'll typically only try a handful of them per
 * descent.
 */
class DescentOutgoing {

public:

	/**
	 * Return the best edge and rotate in preparation for next call.
	 */
	DescentEdge rotate() {
		DescentEdge tmp = best1;
        assert(!(best2 < tmp));
		best1 = best2;
        assert(!(best3 < best2));
		best2 = best3;
        assert(!(best4 < best3));
		best3 = best4;
        assert(!(best5 < best4));
		best4 = best5;
		best5.reset();
		return tmp;
	}
	
	/**
	 * Given a potental outgoing edge, place it where it belongs in the running
	 * list of best 5 outgoing edges from this descent.
	 */
	void update(DescentEdge e) {
		if(!best1.inited()) {
			best1 = e;
		} else if(e < best1) {
			best5 = best4;
			best4 = best3;
			best3 = best2;
			best2 = best1;
			best1 = e;
		} else if(!best2.inited()) {
			best2 = e;
		} else if(e < best2) {
			best5 = best4;
			best4 = best3;
			best3 = best2;
			best2 = e;
		} else if(!best3.inited()) {
			best3 = e;
		} else if(e < best3) {
			best5 = best4;
			best4 = best3;
			best3 = e;
		} else if(!best4.inited()) {
			best4 = e;
		} else if(e < best4) {
			best5 = best4;
			best4 = e;
		}  else if(!best5.inited() || e < best5) {
			best5 = e;
		}
	}
	
	/**
	 * Clear all the outgoing edges stored here.
	 */
	void clear() {
		best1.reset();
		best2.reset();
		best3.reset();
		best4.reset();
		best5.reset();
	}
	
	/**
	 * Return true iff there are no outgoing edges currently represented in
	 * this summary.  There may still be outgoing edges, they just haven't
	 * been added to the summary.
	 */
	bool empty() const {
		return !best1.inited();
	}
	
	/**
	 * Return the DescentPriority of the best outgoing edge.
	 */
	DescentPriority bestPri() const {
		assert(!empty());
		return best1.pri;
	}

	DescentEdge best1; // best
	DescentEdge best2; // 2nd-best
	DescentEdge best3; // 3rd-best
	DescentEdge best4; // 4th-best
	DescentEdge best5; // 5th-best
};

/**
 * Encapsulates the string we're matching during our descent search.
 */
struct DescentQuery {

	DescentQuery() { reset(); }

	/**
	 * Get the nucleotide and quality value at the given offset from 5' end.
	 * If 'fw' is false, get the reverse complement.
	 */
	std::pair<int, int> get(TReadOff off5p, bool fw) const {
		assert_lt(off5p, length());
		int c = (int)(*seq)[off5p];
        int q = (*qual)[off5p];
        assert_geq(q, 33);
		return make_pair((!fw && c < 4) ? (c ^ 3) : c, q - 33);
	}
	
	/**
	 * Get the nucleotide at the given offset from 5' end.
	 * If 'fw' is false, get the reverse complement.
	 */
	int getc(TReadOff off5p, bool fw) const {
		assert_lt(off5p, length());
		int c = (int)(*seq)[off5p];
		return (!fw && c < 4) ? (c ^ 3) : c;
	}
	
	/**
	 * Get the quality value at the given offset from 5' end.
	 */
	int getq(TReadOff off5p) const {
		assert_lt(off5p, length());
        int q = (*qual)[off5p];
        assert_geq(q, 33);
		return q-33;
	}
	
	/**
	 * Initialize.
	 */
	void init(
		const BTDnaString& seq_,
		const BTString&    qual_,
		const BTDnaString& seqrc_,
		const BTString&    qualrc_)
	{
		seq = &seq_;
		qual = &qual_;
		seqrc = &seqrc_;
		qualrc = &qualrc_;
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		seq = NULL;
	}
	
	/**
	 * Return true iff DescentQuery is initialized.
	 */
	bool inited() const {
		return seq != NULL;
	}
    
    size_t length() const {
        assert(inited());
        return seq->length();
    }

    const BTDnaString* seq;
    const BTString* qual;

    const BTDnaString* seqrc;
    const BTString* qualrc;
};

class DescentAlignmentSink;

/**
 * Encapsulates a descent through a search tree, along a path of matches.
 * Descents that are part of the same alignment form a chain.  Two aligments
 * adjacent in the chain are connected either by an edit, or by a switch in
 * direction.  Because a descent might have a different direction from the
 * DescentRoot it ultimately came from, it has its own 'l2r' field, which might
 * differ from the root's.
 */
class Descent {

public:

	Descent() { reset(); }

	/**
	 * Initialize a new descent branching from the given descent via the given
	 * edit.  Return false if the Descent has no outgoing edges (and can
     * therefore have its memory freed), true otherwise.
	 */
	bool init(
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
		DescentMetrics& met);           // metrics

	/**
	 * Initialize a new descent beginning at the given root.  Return false if
     * the Descent has no outgoing edges (and can therefore have its memory
     * freed), true otherwise.
	 */
	bool init(
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
        DescentMetrics& met);           // metrics
	
	/**
	 * Return true iff this Descent has been initialized.
	 */
	bool inited() const {
		return descid_ != std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
        lastRecalc_ = true;
		descid_ = std::numeric_limits<size_t>::max();
	}
	
	/**
	 * Return true iff this Descent is a search root.
	 */
	bool root() const {
		return parent_ == std::numeric_limits<TDescentId>::max();
	}
	
	/**
	 * Take the best outgoing edge and follow it.
	 */
	void followBestOutgoing(
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
        DescentMetrics& met);           // metrics
	
	/**
	 * Return true iff no outgoing edges from this descent remain unexplored.
	 */
	bool empty() const { return lastRecalc_ && out_.empty(); }
	
	/**
	 * Turn the current descent and its chain of parents into a stacked
	 * alignment string.
	 */
#ifndef NDEBUG
	void toStacked(
        std::ostream& os,
        const DescentQuery& q,
        const Ebwt& ebwtFw,             // forward index
        const Ebwt& ebwtBw,             // mirror index
        EFactory<DescentPos>& pf,       // factory with DescentPoss
        EFactory<Descent>& df,          // factory with Descent
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs) // configs
    {
        // Take just the portion of the read that has aligned up until this
        // point
		TDescentId cur = descid_;
        edittmp_.clear();
        size_t nuninited = 0;
		while(cur != std::numeric_limits<TDescentId>::max()) {
            if(!df[cur].edit_.inited()) {
                nuninited++;
                assert_leq(nuninited, 2);
            } else {
                edittmp_.push_back(df[cur].edit_);
            }
			cur = df[cur].parent_;
		}
        edittmp_.sort();
		const bool fw = rs[rid_].fw;
        // Read on top, ref on bottom.  Printing this partial alignment is a
		// little complicated since it might not involve all the characters
		// from the read.  What's needed here is a function that prints all the
		// read characters with the aligned characters in uppercase and the
		// unaligned in lowercase.  And it just prints the reference characters
		// involved in the alignment so far.
		size_t len = q.length();
		assert_lt(al5pf_, len);
		size_t trimLf = fw ? al5pi_ : (len - al5pf_ - 1);
		size_t trimRg = fw ? (len - al5pf_ - 1) : al5pi_;
		BTDnaString& rf = tmprfdnastr_;
		rf.clear();
		print(os, "", q, trimLf, trimRg, fw, edittmp_, rf);
        os << std::endl;
		ASSERT_ONLY(uint32_t toptmp = 0);
		ASSERT_ONLY(uint32_t bottmp = 0);
        // Check that the edited string occurs in the reference
		if(!ebwtFw.contains(rf, &toptmp, &bottmp)) {
			os << rf << std::endl;
			assert(false);
		}
	}
#endif
	
	/**
	 * Check whether this Descent is redundant with the given Descent.
	 */
	bool redundantTo(const Descent& o) const {
		// Return true iff this Descent is redundant with the given Descent and
		// the other descent should be preferred.
		//
		// Problem: we don't have a priority calculated in this class
		return false;
	}
	
	/**
	 * Return true iff the Descent is internally consistent.
	 */
	bool repOk(const DescentQuery *q) const {
		assert( root() ||  edit_.inited());
		assert(!root() || !edit_.inited());
		assert_eq(botf_ - topf_, botb_ - topb_);
		if(q != NULL) {
			assert_leq(len_, q->length());
		}
		return true;
	}

protected:

#ifndef NDEBUG
	void print(
		std::ostream& os,
		const char *prefix,
        const DescentQuery& q,
		size_t trimLf,
		size_t trimRg,
		bool fw,
		const EList<Edit>& edits,
		BTDnaString& rf);
#endif

    bool bounce(
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
        DescentMetrics& met);           // metrics

	/**
	 * Given the forward and backward indexes, and given topf/botf/topb/botb,
	 * get tloc, bloc ready for the next step.
	 */
	void nextLocsBi(
		const Ebwt& ebwtFw, // forward index
		const Ebwt& ebwtBw, // mirror index
		SideLocus& tloc,    // top locus
		SideLocus& bloc,    // bot locus
		uint32_t topf,      // top in BWT
		uint32_t botf,      // bot in BWT
		uint32_t topb,      // top in BWT'
		uint32_t botb);     // bot in BWT'

	/**
	 * Advance this descent by following read matches as far as possible.
	 */
    void followMatches(
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
        TIndexOff& botb_bounce);   // out: bot of SA range for bw idx for bounce

	/**
	 * Recalculate our summary of the outgoing edges from this descent.  When
	 * deciding what outgoing edges are legal, we abide by constraints.
	 * Typically, they limit the total of the penalties accumulated so far, as
	 * a function of distance from the search root.  E.g. a constraint might
	 * disallow any gaps or mismatches within 20 ply of the search root, then
	 * allow 1 mismatch within 30 ply, then allow up to 1 mismatch or 1 gap
	 * within 40 ply, etc.
	 */
	size_t recalcOutgoing(
		const DescentQuery& q,           // query string
		const Scoring& sc,               // scoring scheme
		EFactory<DescentPos>& pf,        // factory with DescentPoss
        const EList<DescentRoot>& rs,    // roots
        const EList<DescentConfig>& cs); // configs

    TRootId         rid_;         // root id

	TReadOff        al5pi_;       // lo offset from 5' end of aligned read char
	TReadOff        al5pf_;       // hi offset from 5' end of aligned read char
	bool            l2r_;         // left-to-right?
    TReadOff        off5p_i_;     // offset we started out at for this descent

	TIndexOff       topf_, botf_; // incoming SA range w/r/t forward index
	TIndexOff       topb_, botb_; // incoming SA range w/r/t forward index

	size_t          descid_;      // ID of this descent
	TDescentId      parent_;      // ID of parent descent
	TScore          pen_;         // total penalties accumulated so far
	size_t          posid_;       // ID of 1st elt of the DescentPos factory w/
	                              // descent pos info for this descent
	size_t          len_;         // length of stretch of matches
	DescentOutgoing out_;         // summary of outgoing edges
	Edit            edit_;        // edit joining this descent with parent
	bool            lastRecalc_;  // set by recalcOutgoing if out edges empty

#ifndef NDEBUG
    EList<Edit>     edittmp_;
	BTDnaString     tmprfdnastr_;
#endif
};


/**
 * Class that accepts alignments found during descent.
 */
class DescentAlignmentSink {

public:

    /**
     * If this is the final descent in a complete end-to-end alignment, report
     * the alignment.
     */
    void reportAlignment(
        const DescentQuery& q,          // query string
        TDescentId id,
        const Edit& e,
        TScore pen,
        EFactory<Descent>& df,          // factory with Descent
        EFactory<DescentPos>& pf,       // factory with DescentPoss
        const EList<DescentRoot>& rs,   // roots
        const EList<DescentConfig>& cs) // configs
    {
        std::cerr << "Found an alignment with total penalty = " << pen << std::endl;
    }
    
    /**
     * Reset to uninitialized state.
     */
    void reset() {
    }

protected:

};

/**
 * Class responsible for advancing all the descents.  The initial descents may
 * emanate from several different locations in the read.  Note that descents
 * may become redundant with each other, and should then be eliminated.
 */
class DescentDriver {
public:

	DescentDriver() { reset(); }
	
	/**
	 * Initialize driver with respect to a new read.
	 */
	void initRead(
		const BTDnaString& seq,
		const BTString& qual,
		const BTDnaString& seqrc,
		const BTString& qualrc)
	{
		reset();
		q_.init(seq, qual, seqrc, qualrc);
	}
	
	/**
	 * Add a new search root, which might (a) prefer to move in a left-to-right
	 * direction, and might (b) be with respect to the read or its reverse
	 * complement.
	 */
	void addRoot(
        const DescentConfig& conf,
        TReadOff off,
        bool l2r,
        bool fw,
        float pri)
    {
        confs_.push_back(conf);
		roots_.push_back(DescentRoot(off, l2r, fw, q_.length(), pri));
	}
	
	/**
	 * Clear the Descent driver so that we're ready to re-start seed alignment
	 * for the current read.
	 */
	void reset() {
		df_.clear();     // clear Descents
		pf_.clear();     // clear DescentPoss
		heap_.clear();   // clear Heap
		roots_.clear();  // clear roots
        alsink_.reset(); // clear alignment sink
	}

	/**
	 * Perform seed alignment.
	 */
	void go(
        const Scoring& sc,    // scoring scheme
		const Ebwt& ebwtFw,   // forward index
		const Ebwt& ebwtBw,   // mirror index
        DescentMetrics& met); // metrics

	/**
	 * Return true iff this DescentDriver is well formed.  Throw an assertion
	 * otherwise.
	 */
	bool repOk() {
		return true;
	}

protected:

	DescentQuery         q_;      // query nucleotide and quality strings
	EFactory<Descent>    df_;     // factory holding all the Descents, which
	                              // must be referred to by ID
	EFactory<DescentPos> pf_;     // factory holding all the DescentPoss, which
	                              // must be referred to by ID
	EList<DescentRoot>   roots_;  // search roots
    EList<DescentConfig> confs_;  // configuration params for each root
	EHeap<TDescentPair>  heap_;   // priority queue of Descents
    DescentAlignmentSink alsink_; // alignment sink
};

#endif
