#ifndef HIT_H_
#define HIT_H_

#include <stdint.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <utility>
#include "alphabet.h"
#include "assert_helpers.h"
#include "spinlock.h"
#include "threading.h"
#include "bitset.h"
#include "tokenize.h"
#include "pat.h"
#include "formats.h"
#include "filebuf.h"
#include "edit.h"
#include "refmap.h"
#include "annot.h"
#include "sstring.h"
#include "ds.h"
#include "ref_coord.h"
#include "aligner_result.h"
#include "unique.h"
#include "read_sink.h"

/**
 * Classes for dealing with reporting alignments.
 */

using namespace std;

extern void printAlSumm(
	uint64_t al,
	uint64_t un,
	uint64_t mx,
	uint64_t rep,
	uint64_t repp,
	bool sample,
	bool hadoopOut);

/// Constants for the various output modes
enum output_types {
	OUTPUT_FULL = 1,
	OUTPUT_CONCISE,
	OUTPUT_BINARY,
	OUTPUT_CHAIN,
	OUTPUT_SAM,
	OUTPUT_NONE
};

/// Names of the various output modes
static const std::string output_type_names[] = {
	"Invalid!",
	"Full",
	"Concise",
	"Binary",
	"None"
};

typedef pair<uint32_t,uint32_t> U32Pair;

/**
 * Print the given string.  If chopws = true, print only up to and not
 * including the first space or tab.  Useful for printing reference
 * names.
 */
static inline void printUptoWs(
	std::ostream& os,
	const std::string& str,
	bool chopws)
{
	if(!chopws) {
		os << str;
	} else {
		size_t pos = str.find_first_of(" \t");
		if(pos != string::npos) {
			os << str.substr(0, pos);
		} else {
			os << str;
		}
	}
}

/**
 * Print the given string.  If chopws = true, print only up to and not
 * including the first space or tab.  Useful for printing reference
 * names.
 */
static inline void printUptoWs(
	OutFileBuf& o,
	const std::string& str,
	bool chopws)
{
	if(!chopws) {
		o.writeString(str);
	} else {
		size_t pos = str.find_first_of(" \t");
		if(pos != string::npos) {
			o.writeString(str.substr(0, pos));
		} else {
			o.writeString(str);
		}
	}
}

/**
 * Print the given string.  If ws = true, print only up to and not
 * including the first space or tab.  Useful for printing reference
 * names.
 */
template<typename T>
static inline void printUptoWs(
	OutFileBuf& o,
	const T& str,
	bool chopws)
{
	if(!chopws) {
		o.writeString(str);
	} else {
		size_t len = str.length();
		for(size_t i = 0; i < len; i++) {
			if(str[i] != ' ' && str[i] != '\t') {
				o.write(str[i]);
			} else {
				break;
			}
		}
	}
}

/**
 * Encapsulates a hit, including a text-id/text-offset pair, a pattern
 * id, and a boolean indicating whether it matched as its forward or
 * reverse-complement version.
 */
class Hit {
public:
	Hit() : edits(16), aedits(16), cedits(16), ccedits(16), stratum(-1) { }

	Hit(const Hit& o) {
		*this = o;
	}
	
	/**
	 * Initialize this Ent to represent the given hit.
	 */
	void toHitSetEnt(HitSetEnt& hse, AnnotationMap *amap) const {
		hse.h = h;
		hse.fw = fw ? 1 : 0;
		hse.oms = oms;
		hse.stratum = stratum;
		hse.cost = cost;
		hse.edits = edits;
	}

	/**
	 * Convert a list of Hit objects to a single HitSet object suitable
	 * for chaining.
	 */
	static void toHitSet(const EList<Hit>& hits, HitSet& hs,
	                     AnnotationMap *amap)
	{
		if(hits.empty()) return;
		// Initialize HitSet
		hs.name = hits.front().name;
		hs.seq  = hits.front().seq;
		hs.qual = hits.front().quals;
		hs.color = hits.front().color;
		if(!hits.front().fw) {
			// Re-reverse
			hs.seq.reverseComp(hs.color);
			hs.qual.reverse();
		}
		// Convert hits to entries
		hs.ents.resize(hits.size());
		for(size_t i = 0; i < hs.ents.size(); i++) {
			hits[i].toHitSetEnt(hs.ents[i], amap);
		}
	}

	/**
	 * Populate a vector of Hits with hits from the given HitSet.
	 */
	static void fromHitSet(EList<Hit>& hits, const HitSet& hs) {
		assert(hs.sorted());
		hits.resize(hs.size());
		for(size_t i = 0; i < hs.size(); i++) {
			hits[i].h.first = hs[i].h.first;
			hits[i].h.second = hs[i].h.second;
			hits[i].oms = hs[i].oms;
			hits[i].fw = hs[i].fw;
			hits[i].stratum = hs[i].stratum;
			hits[i].cost = hs[i].cost;
			hits[i].mate = 0; //hs[i].mate;
			hits[i].name = hs.name;
			hits[i].seq = hs.seq;
			hits[i].quals = hs.qual;
			hits[i].color = hs.color;
			hits[i].edits = hs[i].edits;
			if(!hs[i].fw) {
				hits[i].seq.reverseComp(hs.color);
				hits[i].quals.reverse();
			}
		}
	}

	U32Pair      h;       /// reference index & offset
	TReadId      patid;   /// read index
	BTString     name;    /// read name
	BTDnaString  seq;     /// read sequence
	BTDnaString  cseq;    /// original color sequence, not decoded
	BTString     quals;   /// read qualities
	BTString     cquals;  /// original color qualities, not decoded
	EList<Edit>  edits;   /// decoded nucleotide edits
	EList<Edit>  aedits;  /// resolutions for ambiguous reference characters
	EList<Edit>  cedits;  /// decoded color edits (just miscalls)
	EList<Edit>  ccedits; /// color-to-color edits
	uint32_t     oms;     /// # of other possible mappings; 0 -> this is unique
	bool         fw;      /// orientation of read in alignment
	int8_t       stratum; /// stratum of hit (= mismatches in seed)
	uint16_t     cost;    /// total cost, factoring in stratum and quality penalty
	int8_t       cstratum;/// stratum of colorspace-to-colospace alignment
	uint16_t     ccost;   /// total cost of colorspace-to-colorspace alignment
	uint8_t      mate;    /// matedness; 0 = not a mate
	                      ///            1 = upstream mate
	                      ///            2 = downstream mate
	Hit*         pmate;   /// pointer to mate Hit
	bool         color;   /// read is in colorspace?
	int          readGaps;/// how many gaps in the read (i.e. inserts)?
	int          refGaps; /// how many gaps in the reference (i.e. deletes)?
	uint32_t     seed;    /// pseudo-random seed for aligned read
	
	/**
	 * Check that Hit is internally consistent.
	 */
	bool repOk() const {
#ifndef NDEBUG
		if(pmate == NULL) {
			assert_eq(0, mate);
		} else {
			assert(mate == 1 || mate == 2);
			assert(pmate->mate == 1 || pmate->mate == 2);
			assert_neq(pmate->mate, mate);
		}
		if(stratum > 0) {
			if(color) {
//				assert_geq((int)ccedits.size(), stratum);
			} else {
				assert_geq((int)edits.size(), stratum);
			}
		}
		for(size_t i = 0; i < edits.size(); i++) assert(edits[i].repOk());
		for(size_t i = 0; i < cedits.size(); i++) assert(cedits[i].repOk());
		for(size_t i = 0; i < ccedits.size(); i++) assert(ccedits[i].repOk());
#endif
		return true;
	}

	size_t length() const { return seq.length(); }

	bool operator<(const Hit& b) const {
		if(cost < b.cost) return true;
		if(cost > b.cost) return false;
		if(h < b.h) return true;
		if(h > b.h) return false;
		if(fw < b.fw) return true;
		if(fw > b.fw) return false;
		return false;
	}

	Hit& operator = (const Hit &other) {
		this->h       = other.h;
		this->patid   = other.patid;
		this->name    = other.name;
		this->seq     = other.seq;
		this->cseq    = other.cseq;
		this->quals   = other.quals;
		this->cquals  = other.cquals;
		this->edits   = other.edits;
		this->aedits  = other.aedits;
		this->cedits  = other.cedits;
		this->ccedits = other.ccedits;
		this->oms     = other.oms;
		this->fw      = other.fw;
		this->stratum = other.stratum;
		this->cost    = other.cost;
		this->cstratum= other.cstratum;
		this->ccost   = other.ccost;
		this->mate    = other.mate;
		this->pmate   = other.pmate;
		this->color   = other.color;
		this->readGaps= other.readGaps;
		this->refGaps = other.refGaps;
		this->seed    = other.seed;
		return *this;
	}
};

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b);

/**
 * Needed for asserts involving pairs 
 */
template<typename T1, typename T2>
ostream& operator<< (ostream& out, const std::pair<T1, T2>& c) {
	out << c.first << "," << c.second;
	return out;
}

/**
 * Encapsulates an object that accepts hits, optionally retains them in
 * a vector, and does something else with them according to
 * descendent's implementation of pure virtual member reportHitImpl().
 */
class HitSink {

public:

	explicit HitSink(
		OutFileBuf* out,
		ReadSink* readSink,
		const EList<string>* refnames) :
		outs_(),
		deleteOuts_(false),
		refnames_(refnames),
		numWrappers_(0),
		locks_(),
		first_(true),
		numAligned_(0llu),
		numUnaligned_(0llu),
		numMaxed_(0llu),
		numReported_(0llu),
		numReportedPaired_(0llu),
		quiet_(false),
		ssmode_(ios_base::out),
		readSink_(readSink)
	{
		outs_.push_back(out);
		locks_.resize(1);
		MUTEX_INIT(locks_[0]);
		MUTEX_INIT(mainlock_);
	}

	/**
	 * Destroy HitSinkobject;
	 */
	virtual ~HitSink() {
		closeOuts();
		if(deleteOuts_) {
			// Delete all non-NULL output streams
			for(size_t i = 0; i < outs_.size(); i++) {
				if(outs_[i] != NULL) {
					delete outs_[i];
					outs_[i] = NULL;
				}
			}
		}
	}

	/**
	 * Call this whenever this HitSink is wrapped by a new
	 * HitSinkPerThread.  This helps us keep track of whether the main
	 * lock or any of the per-stream locks will be contended.
	 */
	void addWrapper() {
		numWrappers_++;
	}

	/**
	 * Called by concrete subclasses to figure out which elements of
	 * the outs_/locks_ array to use when outputting the alignment.
	 */
	size_t refIdxToStreamIdx(size_t refIdx) {
		if(refIdx >= outs_.size()) return 0;
		return refIdx;
	}

	/**
	 * Append a single hit to the given output stream.
	 */
	virtual void append(ostream& o, const Hit& h) = 0;

	/**
	 * Report a batch of hits; all in the given vector.
	 */
	virtual void reportHits(EList<Hit>& hs) {
		reportHits(hs, 0, hs.size());
	}

	/**
	 * Report a batch of hits from a vector, perhaps subsetting it.
	 */
	virtual void reportHits(EList<Hit>& hs, size_t start, size_t end) {
		assert_geq(end, start);
		if(end-start == 0) return;
		bool paired = hs[start].mate > 0;
		char buf[4096];
		lock(0);
		for(size_t i = start; i < end; i++) {
			const Hit& h = hs[i];
			ostringstream ss(ssmode_);
			ss.rdbuf()->pubsetbuf(buf, 4096);
			append(ss, h);
			out(h.h.first).writeChars(buf, (size_t)ss.tellp());
		}
		unlock(0);
		mainlock();
		commitHits(hs);
		first_ = false;
		numAligned_++;
		if(paired) numReportedPaired_ += (end-start);
		else       numReported_ += (end-start);
		mainunlock();
	}

	void commitHit(const Hit& hit) { }

	void commitHits(const EList<Hit>& hits) { }

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary.
	 */
	void finish(bool hadoopOut, bool sampleMax) {
		// Close output streams
		closeOuts();
		if(!quiet_) {
			printAlSumm(
				numAligned_,
				numUnaligned_,
				numMaxed_,
				numReported_,
				numReportedPaired_,
				sampleMax,
				hadoopOut);
		}
	}

	/// Returns the alignment output stream; if the stream needs to be
	/// created, create it
	OutFileBuf& out(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		if(outs_[strIdx] == NULL) {
			assert(deleteOuts_);
			ostringstream oss;
			oss << "ref";
			if     (strIdx < 10)    oss << "0000";
			else if(strIdx < 100)   oss << "000";
			else if(strIdx < 1000)  oss << "00";
			else if(strIdx < 10000) oss << "0";
			oss << strIdx << ".map";
			outs_[strIdx] = new OutFileBuf(oss.str().c_str(), ssmode_ == ios_base::binary);
		}
		assert(outs_[strIdx] != NULL);
		return *(outs_[strIdx]);
	}

	/**
	 * Lock the monolithic lock for this HitSink.  This is useful when,
	 * for example, outputting a read to an unaligned-read file.
	 */
	void mainlock() { MUTEX_LOCK(mainlock_); }

	/**
	 * Unlock the monolithic lock for this HitSink.  This is useful
	 * when, for example, outputting a read to an unaligned-read file.
	 */
	void mainunlock() { MUTEX_UNLOCK(mainlock_); }

	/**
	 * Report a maxed-out read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportMaxed(EList<Hit>& hs, PatternSourcePerThread& p) {
		mainlock();
		numMaxed_++;
		mainunlock();
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(PatternSourcePerThread& p) {
		mainlock();
		numUnaligned_++;
		mainunlock();
	}

	void dumpMaxed(const Read* m1, const Read* m2, TReadId rdid) {
		if(readSink_ != NULL) readSink_->dumpMaxed(m1, m2, rdid);
	}
	
	void dumpUnal(const Read* m1, const Read* m2, TReadId rdid) {
		if(readSink_ != NULL) readSink_->dumpUnal(m1, m2, rdid);
	}
	
	void dumpAlign(const Read* m1, const Read* m2, TReadId rdid) {
		if(readSink_ != NULL) readSink_->dumpAlign(m1, m2, rdid);
	}

protected:

	/// Implementation of hit-report
	virtual void reportHit(const Hit& h) {
		mainlock();
		commitHit(h);
		first_ = false;
		if(h.mate > 0) numReportedPaired_++;
		else           numReported_++;
		numAligned_++;
		mainunlock();
	}

	/**
	 * Close (and flush) all OutFileBufs.
	 */
	void closeOuts() {
		// Flush and close all non-NULL output streams
		for(size_t i = 0; i < outs_.size(); i++) {
			if(outs_[i] != NULL && !outs_[i]->closed()) {
				outs_[i]->close();
			}
		}
	}

	/**
	 * Lock the output buffer for the output stream for reference with
	 * index 'refIdx'.  By default, hits for all references are
	 * directed to the same output stream, but if --refout is
	 * specified, each reference has its own reference stream.
	 */
	void lock(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		MUTEX_LOCK(locks_[strIdx]);
	}

	/**
	 * Lock the output buffer for the output stream for reference with
	 * index 'refIdx'.  By default, hits for all references are
	 * directed to the same output stream, but if --refout is
	 * specified, each reference has its own reference stream.
	 */
	void unlock(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		MUTEX_UNLOCK(locks_[strIdx]);
	}

	EList<OutFileBuf*>   outs_;        /// the alignment output stream(s)
	bool                 deleteOuts_;  /// Whether to delete elements of outs_ upon exit
	const EList<string>* refnames_;    /// map from reference indexes to names
	int                  numWrappers_; /// # threads owning a wrapper for this HitSink
	EList<MUTEX_T>       locks_;       /// pthreads mutexes for per-file critical sections
	MUTEX_T              mainlock_;    /// pthreads mutexes for fields of this object

	volatile bool     first_;       /// true -> first hit hasn't yet been reported
	volatile uint64_t numAligned_;  /// # reads with >= 1 alignment
	volatile uint64_t numUnaligned_;/// # reads with no alignments
	volatile uint64_t numMaxed_;    /// # reads with # alignments exceeding -m ceiling
	volatile uint64_t numReported_; /// # single-ended alignments reported
	volatile uint64_t numReportedPaired_; /// # paired-end alignments reported
	bool quiet_;  /// true -> don't print alignment stats at the end
	ios_base::openmode ssmode_;     /// output mode for stringstreams
	
	ReadSink* readSink_;
};

/**
 * A per-thread wrapper for a HitSink.  Incorporates state that a
 * single search thread cares about.
 */
class HitSinkPerThread {
public:
	HitSinkPerThread(HitSink& sink, uint32_t max, uint32_t n) :
		_sink(sink),
		_bestRemainingStratum(0),
		_hits(),
		hits_(),
		hitsForThisRead_(),
		_max(max),
		_n(n)
	{
		_sink.addWrapper();
		assert_gt(_n, 0);
	}

	virtual ~HitSinkPerThread() { }

	/// Return the vector of retained hits
	EList<Hit>& retainedHits()   { return _hits; }
	
	/**
	 * Clear out buffered information about hits and hit endpoints.
	 */
	void clearHits() {
		hits_.clear();
		if(gGaps) {
			lfws_.clear();
			rfws_.clear();
			lrcs_.clear();
			rrcs_.clear();
			lls_.clear();
			lrs_.clear();
			rls_.clear();
			rrs_.clear();
		}
	}

	/// Finalize current read
	virtual uint32_t finishRead(PatternSourcePerThread& p, bool report, bool dump) {
		uint32_t ret = finishReadImpl();
		_bestRemainingStratum = 0;
		if(!report) {
			clearHits();
			return 0;
		}
		bool maxed = (ret > _max);
		bool unal = (ret == 0);
		if(dump && (unal || maxed)) {
			// Either no reportable hits were found or the number of
			// reportable hits exceeded the -m limit specified by the
			// user
			assert(ret == 0 || ret > _max);
			if(maxed) _sink.dumpMaxed(&p.bufa(), &p.bufb(), p.patid());
			else      _sink.dumpUnal(&p.bufa(), &p.bufb(), p.patid());
		}
		ret = 0;
		if(maxed) {
			// Report that the read maxed-out; useful for chaining output
			if(dump) _sink.reportMaxed(hits_, p);
			clearHits();
		} else if(unal) {
			// Report that the read failed to align; useful for chaining output
			if(dump) _sink.reportUnaligned(p);
		} else {
			// Flush buffered hits
			assert_gt(hits_.size(), 0);
			if(hits_.size() > _n) {
				hits_.resize(_n);
			}
			_sink.reportHits(hits_);
			_sink.dumpAlign(&p.bufa(), &p.bufb(), p.patid());
			ret = (uint32_t)hits_.size();
			clearHits();
		}
		assert_eq(0, hits_.size());
		assert(lfws_.empty());
		assert(rfws_.empty());
		assert(lrcs_.empty());
		assert(rrcs_.empty());
		assert(lls_.empty());
		assert(lrs_.empty());
		assert(rls_.empty());
		assert(rrs_.empty());
		return ret;
	}

	virtual uint32_t finishReadImpl() = 0;

	/**
	 * Implementation for hit reporting; update per-thread _hits and
	 * _numReportableHits variables and call the master HitSink to do the actual
	 * reporting
	 */
	virtual void bufferHit(const Hit& h, int stratum) {
		assert(h.repOk());
		assert_eq(0, h.mate);
#ifndef NDEBUG
		// Ensure all buffered hits have the same patid
		for(size_t i = 1; i < hits_.size(); i++) {
			assert_eq(hits_[0].patid, hits_[i].patid);
		}
#endif
		hits_.push_back(h);
		if(gGaps && gAllowRedundant < 2) {
			// A subtlety is that this can spuriously blow away
			// mates from paired-end alignments.  A paired-end
			// alignment isn't redundant unless the entire pair is
			// redundant
			Coord left(h.h.first, h.h.second, h.fw);
			Coord right(h.h.first, h.h.second + h.length() + h.readGaps - h.refGaps, h.fw);
			if(h.fw) {
				assert(!lfws_.contains(left));
				assert(!rfws_.contains(right));
				lfws_.insert(left);
				rfws_.insert(right);
			} else {
				assert(!lrcs_.contains(left));
				assert(!rrcs_.contains(right));
				lrcs_.insert(left);
				rrcs_.insert(right);
			}
		}
	}

	/**
	 * Implementation for hit reporting; update per-thread _hits and
	 * _numReportableHits variables and call the master HitSink to do the actual
	 * reporting
	 */
	virtual void bufferHitPair(const Hit& hL, int stratumL, const Hit& hR, int stratumR) {
		assert(hL.repOk());
		assert(hR.repOk());
		const Hit& h = hL;
		assert(hL.mate == 1 || hL.mate == 2);
		assert(hR.mate == 1 || hR.mate == 2);
		assert_neq(hL.mate, hR.mate);
#ifndef NDEBUG
		// Ensure all buffered hits have the same patid
		for(size_t i = 1; i < hits_.size(); i++) {
			assert_eq(hits_[0].patid, hits_[i].patid);
		}
#endif
		if(gGaps && gAllowRedundant < 2) {
			const Hit& h1 = (h.mate == 1 ? h : *h.pmate);
			const Hit& h2 = (h.mate == 2 ? h : *h.pmate);
			uint64_t left1  = h1.h.second;
			uint64_t right1 = h1.h.second + h1.length() + h1.readGaps - h1.refGaps;
			uint64_t left2  = h2.h.second;
			uint64_t right2 = h2.h.second + h2.length() + h2.readGaps - h2.refGaps;
			Interval ll(h1.h.first, h2.h.first, left1,  left2);
			Interval lr(h1.h.first, h2.h.first, left1,  right2);
			Interval rl(h1.h.first, h2.h.first, right1, left2);
			Interval rr(h1.h.first, h2.h.first, right1, right2);
			assert(!lls_.contains(ll));
			assert(!lrs_.contains(lr));
			assert(!rls_.contains(rl));
			assert(!rrs_.contains(rr));
			lls_.insert(ll);
			lrs_.insert(lr);
			rls_.insert(rl);
			rrs_.insert(rr);
		}
		hits_.push_back(hL);
		hits_.push_back(hR);
		hits_[hits_.size()-2].pmate = &hits_[hits_.size()-1];
		hits_[hits_.size()-1].pmate = &hits_[hits_.size()-2];
	}
	
	/**
	 * Concrete subclasses override this to (possibly) report a hit and
	 * return true iff the caller should continue to report more hits.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		assert(h.repOk());
		assert_eq(0, h.mate);
		if(!gGaps) return true;
		else {
			Coord left(h.h.first, h.h.second, h.fw);
			Coord right(h.h.first, h.h.second + h.length() + h.readGaps - h.refGaps, h.fw);
			if(h.fw) {
				if(lfws_.contains(left) || rfws_.contains(right)) {
					return false;
				}
			} else {
				if(lrcs_.contains(left) || rrcs_.contains(right)) {
					return false;
				}
			}
			return true;
		}
	}
	
	/**
	 * Concrete subclasses override this to (possibly) report a hit and
	 * return true iff the caller should continue to report more hits.
	 */
	virtual bool reportHitPair(const Hit& h, const Hit& hR, int stratumL, int stratumR) {
		assert(h.repOk());
		assert_range(1, 2, (int)h.mate);
		if(!gGaps) return true;
		else {			
			const Hit& h1 = (h.mate == 1 ? h : *h.pmate);
			const Hit& h2 = (h.mate == 2 ? h : *h.pmate);
			uint64_t left1  = h1.h.second;
			uint64_t right1 = h1.h.second + h1.length() + h1.readGaps - h1.refGaps;
			uint64_t left2  = h2.h.second;
			uint64_t right2 = h2.h.second + h2.length() + h2.readGaps - h2.refGaps;
			Interval ll(h1.h.first, h2.h.first, left1,  left2);
			Interval lr(h1.h.first, h2.h.first, left1,  right2);
			Interval rl(h1.h.first, h2.h.first, right1, left2);
			Interval rr(h1.h.first, h2.h.first, right1, right2);
			if(lls_.contains(ll) || lrs_.contains(lr) ||
			   rls_.contains(rl) || rrs_.contains(rr))
			{
				return false;
			}
			return true;
		}
	}

	/**
	 * Return true if there are no more reportable hits.
	 */
	bool finishedWithStratum(int stratum) {
		bool ret = finishedWithStratumImpl(stratum);
		_bestRemainingStratum = stratum+1;
		return ret;
	}

	/**
	 * Use the given set of hits as a starting point.  By default, we don't
	 */
	virtual bool setHits(HitSet& hs) {
		if(!hs.empty()) {
			cerr << "Error: default setHits() called with non-empty HitSet" << endl;
			throw 1;
		}
		return false;
	}

	/**
	 * Return true if there are no reportable hits with the given cost
	 * (or worse).
	 */
	virtual bool irrelevantCost(uint16_t cost) {
		return false;
	}

	/**
	 * Concrete subclasses override this to determine whether the
	 * search routine should keep searching after having finished
	 * reporting all alignments at the given stratum.
	 */
	virtual bool finishedWithStratumImpl(int stratum) = 0;

	/// The mhits maximum
	uint32_t overThresh() { return _max; }

	/// Whether this thread, for this read, knows that we have already
	/// exceeded the mhits maximum
	bool exceededOverThresh() { return hitsForThisRead_ > _max; }

	/// Return whether we span strata
	virtual bool spanStrata() = 0;

	/// Return whether we report only the best possible hits
	virtual bool best() = 0;

	/**
	 * Return true iff there are currently no buffered hits.
	 */
	bool empty() const {
		return hits_.empty();
	}

	/**
	 * Return the number of currently buffered hits.
	 */
	bool size() const {
		return hits_.size();
	}

	/**
	 * Return max # hits to report (*2 in paired-end mode because mates
	 * count separately)
	 */
	virtual uint32_t maxHits() {
		return _n;
	}

protected:
	HitSink&    _sink; /// Ultimate destination of reported hits
	/// Least # mismatches in alignments that will be reported in the
	/// future.  Updated by the search routine.
	int         _bestRemainingStratum;
	/// # hits reported to this HitSink so far (not all of which were
	/// necesssary reported to _sink)
	EList<Hit> _hits; /// Repository for retained hits
	/// Buffered hits, to be reported and flushed at end of read-phase
	EList<Hit> hits_;
	// For detecting redundant unpaired alignments, these holds sets of
	// previously-seen offsets for extreme-left and extreme-right ends
	// of reported alignments
	ESet<Coord> lfws_;
	ESet<Coord> rfws_;
	ESet<Coord> lrcs_;
	ESet<Coord> rrcs_;
	// For detecting redundant paired alignments, these holds sets of
	// previously-seen L1/L2, L1/R2, L2/R1, R1/R2 offset combos.
	ESet<Interval> lls_;
	ESet<Interval> lrs_;
	ESet<Interval> rls_;
	ESet<Interval> rrs_;
	// Following variables are declared in the parent but maintained in
	// the concrete subcalsses
	uint32_t hitsForThisRead_; /// # hits for this read so far
	uint32_t _max; /// don't report any hits if there were > _max
	uint32_t _n;   /// report at most _n hits
};

/**
 * Abstract parent factory for HitSinkPerThreads.
 */
class HitSinkPerThreadFactory {
public:
	virtual ~HitSinkPerThreadFactory() { }
	virtual HitSinkPerThread* create() const = 0;
	virtual HitSinkPerThread* createMult(uint32_t m) const = 0;

	/// Free memory associated with a per-thread hit sink
	virtual void destroy(HitSinkPerThread* sink) const {
		assert(sink != NULL);
		// Free the HitSinkPerThread
		delete sink;
	}
};

/**
 * Report first N good alignments encountered; trust search routine
 * to try alignments in something approximating a best-first order.
 * Best used in combination with a stringent alignment policy.
 */
class NGoodHitSinkPerThread : public HitSinkPerThread {

public:
	NGoodHitSinkPerThread(
			HitSink& sink,
			uint32_t n,
			uint32_t max) :
				HitSinkPerThread(sink, max, n)
	{ }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return false; // we settle for "good" hits
	}

	/// Finalize current read
	virtual uint32_t finishReadImpl() {
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		return ret;
	}

	/**
	 * Report and then return true if we've already reported N good
	 * hits.  Ignore the stratum - it's not relevant for finding "good"
	 * hits.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		assert(h.repOk());
		if(!HitSinkPerThread::reportHit(h, stratum)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		hitsForThisRead_++;
		if(hitsForThisRead_ > _max) {
			return true; // done - report nothing
		}
		//if(hitsForThisRead_ <= _n) {
			// Only report hit if we haven't
			bufferHit(h, stratum);
		//}
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits and max isn't set; stop!
		}
		return false; // not at N or max yet; keep going
	}

	/**
	 * Report and then return true if we've already reported N good
	 * hits.  Ignore the stratum - it's not relevant for finding "good"
	 * hits.
	 */
	virtual bool reportHitPair(const Hit& hL, const Hit& hR, int stratumL, int stratumR) {
		if(!HitSinkPerThread::reportHitPair(hL, hR, stratumL, stratumR)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		hitsForThisRead_++;
		if(hitsForThisRead_ > _max) {
			return true; // done - report nothing
		}
		hitsForThisRead_++;
		// Only report hit if we haven't
		bufferHitPair(hL, stratumL, hR, stratumR);
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits and max isn't set; stop!
		}
		return false; // not at N or max yet; keep going
	}
	
	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return false; }
};

/**
 * Concrete factory for FirstNGoodHitSinkPerThreads.
 */
class NGoodHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NGoodHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max) :
			sink_(sink),
			n_(n),
			max_(max)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NGoodHitSinkPerThread(sink_, n_, max_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		uint32_t n = n_ * (n_ == 0xffffffff ? 1 : m);
		return new NGoodHitSinkPerThread(sink_, n, max);
	}

private:
	HitSink& sink_;
	uint32_t n_;
	uint32_t max_;
};

/**
 * Report the first N best alignments encountered in a single
 * alignment stratum assuming that we're receiving the alignments in
 * best-first order.
 */
class NBestFirstStratHitSinkPerThread : public HitSinkPerThread {

public:
	NBestFirstStratHitSinkPerThread(
			HitSink& sink,
			uint32_t n,
			uint32_t max,
			uint32_t mult) :
				HitSinkPerThread(sink, max, n),
				bestStratum_(999), mult_(mult)
	{ }

	/**
	 * false -> we do not allow strata to be spanned
	 */
	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	/**
	 * true -> we report best hits
	 */
	virtual bool best() {
		return true;
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		if(!HitSinkPerThread::reportHit(h, stratum)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		// This hit is within th best possible remaining stratum,
		// so it should definitely count
		hitsForThisRead_++;
		// It doesn't exceed the limit, so buffer it
		if(stratum < bestStratum_) {
			bestStratum_ = stratum;
		}
		if(hitsForThisRead_ > _max) {
			return true; // done - report nothing
		}
		//if(hitsForThisRead_ <= _n) {
			bufferHit(h, stratum);
		//}
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits; stop!
		}
		return false; // not at N yet; keep going
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHitPair(const Hit& hL, const Hit& hR, int stratumL, int stratumR) {
		if(!HitSinkPerThread::reportHitPair(hL, hR, stratumL, stratumR)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		// This hit is within th best possible remaining stratum,
		// so it should definitely count
		hitsForThisRead_++;
		// It doesn't exceed the limit, so buffer it
		if(stratumL < bestStratum_) {
			bestStratum_ = stratumL;
		}
		if(stratumR < bestStratum_) {
			bestStratum_ = stratumR;
		}
		if(hitsForThisRead_ > _max) {
			return true; // done - report nothing
		}
		hitsForThisRead_++;
		bufferHitPair(hL, stratumL, hR, stratumR);
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits; stop!
		}
		return false; // not at N yet; keep going
	}
	
	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual uint32_t finishReadImpl() {
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		bestStratum_ = 999;
		const size_t sz = hits_.size();
		for(size_t i = 0; i < sz; i++) {
			// Set 'oms' according to the number of other alignments
			// at this stratum
			hits_[i].oms = (uint32_t)(sz / mult_) - 1;
		}
		return ret;
	}

	/**
	 * If we had any alignments at all and we're now moving on to a new
	 * stratum, then we're done.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
		return hitsForThisRead_ > 0;
	}

	/**
	 * If there have been any hits reported so far, classify any
	 * subsequent alignments with higher strata as irrelevant.
	 */
	virtual bool irrelevantCost(uint16_t cost) {
		if(hitsForThisRead_) {
			// irrelevant iff at worse stratum
			return ((int)cost >> 14) > bestStratum_;
		}
		return false;
	}

private:

	int bestStratum_; /// best stratum observed so far
	uint32_t mult_; /// number of batched-up alignments
};

/**
 * Concrete factory for NBestStratHitSinkPerThread.
 */
class NBestFirstStratHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NBestFirstStratHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max) :
			sink_(sink),
			n_(n),
			max_(max)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NBestFirstStratHitSinkPerThread(sink_, n_, max_, 1);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		uint32_t n = n_ * (n_ == 0xffffffff ? 1 : m);
		return new NBestFirstStratHitSinkPerThread(sink_, n, max, m);
	}

private:
	HitSink& sink_;
	uint32_t n_;
	uint32_t max_;
};

/**
 *
 */
class ChainingHitSinkPerThread : public HitSinkPerThread {
public:

	ChainingHitSinkPerThread(HitSink& sink,
	                         uint32_t n,
	                         uint32_t max,
	                         bool strata,
	                         uint32_t mult) :
	HitSinkPerThread(sink, max, n),
	mult_(mult), strata_(strata), cutoff_(0xffff)
	{
		hs_ = NULL;
		hsISz_ = 0;
	}

	/**
	 * Return true iff we're allowed to span strata.
	 */
	virtual bool spanStrata() { return !strata_; }

	/**
	 * true -> we report best hits
	 */
	virtual bool best() { return true; }

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		if(!HitSinkPerThread::reportHit(h, stratum)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		assert(hs_->sorted());
		assert(hs_ != NULL);
		assert_eq(h.stratum, stratum);
		assert_eq(1, mult_);
		assert(consistentStrata());
		assert(!irrelevantCost(h.cost));

		if(!hs_->empty() && strata_ && stratum < hs_->front().stratum) {
			hs_->clear();
			hits_.clear();
			hitsForThisRead_ = 0;
		}
		assert(consistentStrata());

		size_t replPos = 0;
		if(!hs_->empty() && hs_->tryReplacing(h.h, h.fw, h.cost, replPos)) {
			if(replPos != 0xffffffff) {
				// Replaced an existing hit
				assert_lt(replPos, hits_.size());
				hits_[replPos] = h;
				hs_->sort();
			}
			// Size didn't change, so no need to check against _max and _n
			assert(hs_->sorted());
			assert(consistentStrata());
		} else {
			// Added a new hit
			hs_->expand();
			hs_->back().h = h.h;
			hs_->back().fw = h.fw;
			hs_->back().stratum = h.stratum;
			hs_->back().cost = h.cost;
			hitsForThisRead_++;
			if(hs_->size() > _max) {
				assert_eq(hs_->size(), hits_.size());
				return true; // done - report nothing
			}
			hits_.push_back(h);
			if(hsISz_ == 0 &&
			   hs_->size() == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				assert_eq(hs_->size(), hits_.size());
				return true; // already reported N good hits; stop!
			}
			hs_->sort();
			assert(consistentStrata());
		}
		assert_eq(hs_->size(), hits_.size());
		updateCutoff();
		return false; // not at N yet; keep going
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHitPair(const Hit& hL, const Hit& hR, int stratumL, int stratumR) {
		if(!HitSinkPerThread::reportHitPair(hL, hR, stratumL, stratumR)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		throw 1;
		return true;
	}
	
	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual uint32_t finishReadImpl() {
		assert_eq(1, mult_);
		assert(hs_ != NULL);
		assert(consistentStrata());
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		if(!hs_->empty() && hs_->size() < _n) {
			const size_t sz = hits_.size();
			for(size_t i = 0; i < sz; i++) {
				// Set 'oms' according to the number of other alignments
				// at this stratum
				hits_[i].oms = (uint32_t)(sz / mult_) - 1;
			}
		}
		//hits_.sort();
		if(hs_->size() > _n) {
			hits_.resize(_n);
		}
		assert(consistentStrata());
		// Make sure that a chained, maxed-out read gets treated as
		// maxed-out and not as unaligned
		if(hs_->empty() && hs_->maxedStratum != -1) {
			assert(hits_.empty());
			// Boy, this is stupid.  Need to switch to just using HitSet
			// internally all the time.
			hits_.resize(_max+1);
			for(size_t i = 0; i < _max+1; i++) {
				hits_[i].stratum = hs_->maxedStratum;
			}
			ret = _max+1;
		} else if(!hs_->empty() && hs_->maxedStratum != -1) {
			assert_lt(hs_->front().stratum, hs_->maxedStratum);
		}
		return ret;
	}

	/**
	 * Set the initial set of Hits.
	 */
	virtual bool setHits(HitSet& hs) {
		hs_ = &hs;
		assert(hs_ != NULL);
		hsISz_ = hs.size();
		cutoff_ = 0xffff;
		hitsForThisRead_ = (uint32_t)hs.size();
		assert(hits_.empty());
		assert_geq(hs.maxedStratum, -1);
		assert_lt(hs.maxedStratum, 4);
		if(!hs.empty()) {
			assert_eq(-1, hs.maxedStratum);
			hs.sort();
			Hit::fromHitSet(hits_, hs);
			assert(!hits_.empty());
			assert_leq(hits_.size(), _max);
			assert(consistentStrata());
		} else if(hs.maxedStratum != -1) {
			if(hs.maxedStratum == 0) {
				cutoff_ = 0;
				// Already done
				return true;
			}
			cutoff_ = (hs.maxedStratum << 14);
		}
		assert_eq(hs_->size(), hits_.size());
		updateCutoff();
		return false;
	}

	/**
	 * If we had any alignments at all and we're now moving on to a new
	 * stratum, then we're done.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
		assert(false);
		return false;
	}

	/**
	 * If there have been any hits reported so far, classify any
	 * subsequent alignments with higher strata as irrelevant.
	 */
	virtual bool irrelevantCost(uint16_t cost) {
		if(cutoff_ == 0) return false;
		return cost > cutoff_;
	}

protected:

#ifndef NDEBUG
	/**
	 * Sanity check that, if we're in strata_ mode, all 'stratum's
	 * should be the same among hits.
	 */
	bool consistentStrata() {
		if(hs_->empty() || !strata_) return true;
		int stratum = hs_->front().stratum;
		for(size_t i = 1; i < hs_->size(); i++) {
			assert_eq(stratum, (*hs_)[i].stratum);
		}
		return true;
	}
#endif

	/**
	 * Update the lowest relevant cost.
	 */
	void updateCutoff() {
		ASSERT_ONLY(uint16_t origCutoff = cutoff_);
		assert(hs_->sorted());
		assert(hs_->empty() || hs_->back().cost >= hs_->front().cost);
		bool atCapacity = (hs_->size() >= _n);
		if(atCapacity && (_max == 0xffffffff || _max < _n)) {
			cutoff_ = min(hs_->back().cost, cutoff_);
		}
		if(strata_ && !hs_->empty()) {
			uint16_t sc = hs_->back().cost;
			sc = ((sc >> 14) + 1) << 14;
			cutoff_ = min(cutoff_, sc);
			assert_leq(cutoff_, origCutoff);
		}
	}

	HitSet *hs_;
	size_t hsISz_;
	uint32_t mult_;
	bool strata_; /// true -> reporting is stratified
	uint16_t cutoff_; /// the smallest irrelevant cost
};

/**
 * Concrete factory for ChainingHitSinkPerThread.
 */
class ChainingHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	ChainingHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max,
			bool strata) :
			sink_(sink),
			n_(n),
			max_(max),
			strata_(strata)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new ChainingHitSinkPerThread(sink_, n_, max_, strata_, 1);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		uint32_t n = n_ * (n_ == 0xffffffff ? 1 : m);
		return new ChainingHitSinkPerThread(sink_, n, max, strata_, m);
	}

private:
	HitSink& sink_;
	uint32_t n_;
	uint32_t max_;
	bool strata_;
};

/**
 * Report all valid alignments.
 */
class AllHitSinkPerThread : public HitSinkPerThread {

public:
	AllHitSinkPerThread(
			HitSink& sink,
	        uint32_t max) :
		    HitSinkPerThread(sink, max, 0xffffffff) { }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and always return true; we're finiding all hits so that
	 * search routine should always continue.
	 */
	virtual bool reportHitPair(const Hit& hL, const Hit& hR, int stratumL, int stratumR) {
		if(!HitSinkPerThread::reportHitPair(hL, hR, stratumL, stratumR)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		hitsForThisRead_++;
		if(hitsForThisRead_ > _max) {
			return true; // done - report nothing
		}
		hitsForThisRead_++;
		bufferHitPair(hL, stratumL, hR, stratumR);
		return false; // reporting all; always keep going
	}
	
	/**
	 * Report and always return true; we're finiding all hits so that
	 * search routine should always continue.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		if(!HitSinkPerThread::reportHit(h, stratum)) {
			// Hit rejected because one of its ending anchors overlaps
			// with a previously reported alignment
			return false;
		}
		hitsForThisRead_++;
		if(hitsForThisRead_ > _max) {
			return true; // done - report nothing
		}
		bufferHit(h, stratum);
		return false; // reporting all; always keep going
	}

	/**
	 * Finalize; do nothing because we haven't buffered anything
	 */
	virtual uint32_t finishReadImpl() {
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		return ret;
	}

	/**
	 * Always return false; search routine should not stop.
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return false; }
};

/**
 * Concrete factory for AllHitSinkPerThread.
 */
class AllHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	AllHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t max) :
			sink_(sink),
			max_(max)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new AllHitSinkPerThread(sink_, max_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		return new AllHitSinkPerThread(sink_, max);
	}

private:
	HitSink& sink_;
	uint32_t max_;
};

/**
 * Sink that prints lines like this:
 * pat-name \t [-|+] \t ref-name \t ref-off \t pat \t qual \t #-alt-hits \t mm-list
 */
class VerboseHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream VerboseHitSink (default)
	 */
	VerboseHitSink(
		OutFileBuf* out,
		ReadSink* readSink,
		const EList<std::string>* refnames,
		int offBase,
		bool colorSeq,
		bool colorQual,
		bool printCost,
		const EList<bool>& suppressOuts,
		ReferenceMap *rmap,
		AnnotationMap *amap,
		bool fullRef,
		bool sampleMax,
		int partition = 0) :
		HitSink(out, readSink, refnames),
		partition_(partition),
		offBase_(offBase),
		colorSeq_(colorSeq),
		colorQual_(colorQual),
		cost_(printCost),
		suppress_(suppressOuts),
		fullRef_(fullRef),
		sampleMax_(sampleMax),
		rmap_(rmap),
		amap_(amap)
	{ }

	// In hit.cpp
	static void append(
		ostream& ss,
		const Hit& h,
		const EList<string>* refnames,
		ReferenceMap *rmap,
		AnnotationMap *amap,
		bool fullRef,
		int partition,
		int offBase,
		bool colorSeq,
		bool colorQual,
		bool cost,
		const EList<bool>& suppress);

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(ostream& ss, const Hit& h) {
		VerboseHitSink::append(
			ss,
			h,
			refnames_,
			rmap_,
			amap_,
			fullRef_,
			partition_,
			offBase_,
			colorSeq_,
			colorQual_,
			cost_,
			suppress_);
	}

	/**
	 * See hit.cpp
	 */
	virtual void reportMaxed(EList<Hit>& hs, PatternSourcePerThread& p);

protected:

	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream.
	 */
	virtual void reportHit(const Hit& h) {
		reportHit(h, true);
	}

	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream.
	 */
	virtual void reportHit(const Hit& h, bool count) {
		if(count) HitSink::reportHit(h);
		ostringstream ss;
		append(ss, h);
		// Make sure to grab lock before writing to output stream
		lock(h.h.first);
		out(h.h.first).writeString(ss.str());
		unlock(h.h.first);
	}

private:
	int      partition_;   /// partition size, or 0 if partitioning is disabled
	int      offBase_;     /// Add this to reference offsets before outputting.
	                       /// (An easy way to make things 1-based instead of
	                       /// 0-based)
	bool     colorSeq_;    /// true -> print colorspace alignment sequence in colors
	bool     colorQual_;   /// true -> print colorspace quals as originals, not decoded
	bool     cost_;        /// true -> print statum and cost
	EList<bool> suppress_; /// output fields to suppress
	bool     fullRef_;     /// print full reference name
	bool     sampleMax_;   /// true -> user specified -M
	ReferenceMap *rmap_;   /// mapping to reference coordinate system.
	AnnotationMap *amap_;  ///
};

/**
 * Sink for outputting alignments in a binary format.
 */
class ChainingHitSink : public HitSink {
public:

	/**
	 * Construct a single-stream BinaryHitSink (default)
	 */
	ChainingHitSink(
		OutFileBuf* out,
		ReadSink *readSink,
		const EList<std::string>* refnames,
		bool strata,
		AnnotationMap *amap) :
		HitSink(out, readSink, refnames),
		amap_(amap),
		strata_(strata)
	{
		ssmode_ |= ios_base::binary;
	}

	/**
	 * Report a batch of hits.
	 */
	virtual void reportHits(EList<Hit>& hs);

	/**
	 * Append a binary alignment to the output stream corresponding to
	 * the reference sequence involved.
	 */
	virtual void append(ostream& o, const Hit& h) {
		cerr << "Error: ChainingHitSink::append() not implemented" << endl;
		throw 1;
	}

	/**
	 * See hit.cpp
	 */
	virtual void reportMaxed(EList<Hit>& hs, PatternSourcePerThread& p);

	/**
	 * See hit.cpp
	 */
	virtual void reportUnaligned(PatternSourcePerThread& p);

protected:
	AnnotationMap *amap_;
	bool strata_;
};

/**
 * Sink that does nothing.
 */
class StubHitSink : public HitSink {
public:
	StubHitSink() : HitSink(new OutFileBuf(".tmp"), NULL, NULL) { }
	virtual void append(ostream& o, const Hit& h) { }
};

#endif /*HIT_H_*/
