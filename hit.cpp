#include "hit.h"
#include "hit_set.h"
#include "search_globals.h"
#include "sstring.h"
#include "util.h"
#include "unique.h"
#include "aligner_seed.h"

using namespace std;

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b) {
    return a.h < b.h;
}

/**
 * Print a friendly summary of:
 *
 *  1. How many reads were aligned and had one or more alignments
 *     reported
 *  2. How many reads exceeded the -m or -M ceiling and therefore had
 *     their alignments suppressed or sampled
 *  3. How many reads failed to align entirely
 *
 * Optionally print a series of Hadoop streaming-style counter updates
 * with similar information.
 */
void printAlSumm(
	uint64_t al,
	uint64_t un,
	uint64_t mx,
	uint64_t rep,
	uint64_t repp,
	bool sample,
	bool hadoopOut)
{
	// Print information about how many unpaired and/or paired
	// reads were aligned.
	uint64_t tot = al + un + mx;
	double alPct = 0.0, unalPct = 0.0, maxPct = 0.0;
	if(tot > 0) {
		alPct   = 100.0 * (double)al / (double)tot;
		unalPct = 100.0 * (double)un / (double)tot;
		maxPct  = 100.0 * (double)mx / (double)tot;
	}
	cerr << "# reads processed: " << tot << endl;
	cerr << "# reads with at least one reported alignment: "
		 << al << " (" << fixed << setprecision(2)
		 << alPct << "%)" << endl;
	cerr << "# reads that failed to align: "
		 << un << " (" << fixed << setprecision(2)
		 << unalPct << "%)" << endl;
	if(mx > 0) {
		if(sample) {
			cerr << "# reads with alignments sampled due to -M: "
				 << mx << " (" << fixed << setprecision(2)
				 << maxPct << "%)" << endl;
		} else {
			cerr << "# reads with alignments suppressed due to -m: "
				 << mx << " (" << fixed << setprecision(2)
				 << maxPct << "%)" << endl;
		}
	}
	if((rep + repp) == 0) {
		assert_eq(0llu, rep);
		cerr << "No alignments" << endl;
	}
	else if(repp > 0 && rep == 0) {
		cerr << "Reported " << (repp >> 1)
			 << " paired-end alignments" << endl;
	}
	else if(rep > 0 && repp == 0) {
		cerr << "Reported " << rep << " alignments" << endl;
	}
	else {
		assert_gt(rep, 0);
		assert_gt(repp, 0);
		cerr << "Reported " << (repp >> 1)
			 << " paired-end alignments and " << rep
			 << " unpaired alignments" << endl;
	}
	if(hadoopOut) {
		cerr << "reporter:counter:Bowtie,Reads with reported alignments," << al   << endl;
		cerr << "reporter:counter:Bowtie,Reads with no alignments,"       << un   << endl;
		cerr << "reporter:counter:Bowtie,Reads exceeding -m limit,"       << mx   << endl;
		cerr << "reporter:counter:Bowtie,Unpaired alignments reported,"   << rep  << endl;
		cerr << "reporter:counter:Bowtie,Paired alignments reported,"     << repp << endl;
	}
}

/**
 * Report a batch of hits to a chaining file.
 */
void ChainingHitSink::reportHits(EList<Hit>& hs) {
	size_t hssz = hs.size();
	assert_gt(hssz, 0);
	assert_eq(0, hs[0].mate);
	// Convert EList<Hit> into HitSet
	{
		// Critical section for output stream 0
		HitSet s;
		Hit::toHitSet(hs, s, amap_);
		lock(0);
		s.serialize(out(0));
		unlock(0);
	}
	{
		// Global critical section
		mainlock();
		commitHits(hs); // Commit to recalibration table
		first_ = false;
		numReported_ += hssz;
		numAligned_++;
		mainunlock();
	}
}

/**
 * Report a maxed-out read.  Typically we do nothing, but we might
 * want to print a placeholder when output is chained.
 */
void ChainingHitSink::reportMaxed(EList<Hit>& hs, PatternSourcePerThread& p) {
	HitSink::reportMaxed(hs, p);
	assert(!hs.empty());
	int8_t loStrat = (strata_ ? hs.front().stratum : 0);
	HitSet s;
	s.fromRead(p.bufa());
	//p.bufa().toHitSet(s);
	s.maxedStratum = loStrat;
	lock(0);
	s.serialize(out(0));
	unlock(0);
}

/**
 * Report an unaligned read.  Typically we do nothing, but we might
 * want to print a placeholder when output is chained.
 */
void ChainingHitSink::reportUnaligned(PatternSourcePerThread& p) {
	HitSink::reportUnaligned(p);
	// Read is unaligned; just report a huge starting stratum
	HitSet s;
	s.fromRead(p.bufa());
	//p.bufa().toHitSet(s);
	lock(0);
	s.serialize(out(0));
	unlock(0);
}

/**
 * Report a maxed-out read.
 */
void VerboseHitSink::reportMaxed(EList<Hit>& hs, PatternSourcePerThread& p) {
	HitSink::reportMaxed(hs, p);
	if(sampleMax_) {
		RandomSource rand;
		rand.init(p.bufa().seed);
		assert_gt(hs.size(), 0);
		bool paired = hs.front().mate > 0;
		size_t num = 1;
		if(paired) {
			num = 0;
			int bestStratum = 999;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat < bestStratum) {
					bestStratum = strat;
					num = 1;
				} else if(strat == bestStratum) {
					num++;
				}
			}
			assert_leq(num, hs.size());
			uint32_t r = rand.nextU32() % num;
			num = 0;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat == bestStratum) {
					if(num == r) {
						hs[i].oms = hs[i+1].oms = hs.size()/2;
						reportHits(hs, i, i+2);
						break;
					}
					num++;
				}
			}
			assert_eq(num, r);
		} else {
			for(size_t i = 1; i < hs.size(); i++) {
				assert_geq(hs[i].stratum, hs[i-1].stratum);
				if(hs[i].stratum == hs[i-1].stratum) num++;
				else break;
			}
			assert_leq(num, hs.size());
			uint32_t r = rand.nextU32() % num;
			Hit& h = hs[r];
			h.oms = hs.size();
			reportHit(h, false);
		}
	}
}

/**
 * Append a verbose, readable hit to the given output stream.
 */
void VerboseHitSink::append(
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
	const Bitset& suppress)
{
	assert(h.repOk());
	bool spill = false;
	int spillAmt = 0;
	uint32_t pdiv = 0xffffffff;
	uint32_t pmod = 0xffffffff;
	do {
		bool dospill = false;
		if(spill) {
			// The read spilled over a partition boundary and so
			// needs to be printed more than once
			spill = false;
			dospill = true;
			spillAmt++;
		}
		assert(!spill);
		size_t field = 0;
		bool firstfield = true;
		if(partition != 0) {
			int pospart = abs(partition);
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// Output a partitioning key
				// First component of the key is the reference index
				if(refnames != NULL && rmap != NULL) {
					printUptoWs(ss, rmap->getName(h.h.first), !fullRef);
				} else if(refnames != NULL && h.h.first < refnames->size()) {
					printUptoWs(ss, (*refnames)[h.h.first], !fullRef);
				} else {
					ss << h.h.first;
				}
			}
			ostringstream ss2, ss3;
			// Next component of the key is the partition id
			if(!dospill) {
				pdiv = (h.h.second + offBase) / pospart;
				pmod = (h.h.second + offBase) % pospart;
			}
			assert_neq(0xffffffff, pdiv);
			assert_neq(0xffffffff, pmod);
			if(dospill) assert_gt(spillAmt, 0);
			ss2 << (pdiv + (dospill ? spillAmt : 0));
			if(partition > 0 &&
			   (pmod + h.length()) >= ((uint32_t)pospart * (spillAmt + 1))) {
				// Spills into the next partition so we need to
				// output another alignment for that partition
				spill = true;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// Print partition id with leading 0s so that Hadoop
				// can do lexicographical sort (modern Hadoop versions
				// seen to support numeric)
				string s2 = ss2.str();
				size_t partDigits = 1;
				if(pospart >= 10) partDigits++;
				if(pospart >= 100) partDigits++;
				if(pospart >= 1000) partDigits++;
				if(pospart >= 10000) partDigits++;
				if(pospart >= 100000) partDigits++;
				for(size_t i = s2.length(); i < (10-partDigits); i++) {
					ss << "0";
				}
				ss << s2.c_str();
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// Print offset with leading 0s
				ss3 << (h.h.second + offBase);
				string s3 = ss3.str();
				for(size_t i = s3.length(); i < 9; i++) {
					ss << "0";
				}
				ss << s3;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.fw? "+":"-");
			}
			// end if(partition != 0)
		} else {
			assert(!dospill);
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.name;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.fw? '+' : '-');
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// .first is text id, .second is offset
				if(refnames != NULL && rmap != NULL) {
					printUptoWs(ss, rmap->getName(h.h.first), !fullRef);
				} else if(refnames != NULL && h.h.first < refnames->size()) {
					printUptoWs(ss, (*refnames)[h.h.first], !fullRef);
				} else {
					ss << h.h.first;
				}
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.h.second + offBase);
			}
			// end else clause of if(partition != 0)
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			const BTDnaString* pat = &h.seq;
			if(h.color && colorSeq) pat = &h.cseq;
			ss << *pat;
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			const BTString* qual = &h.quals;
			if(h.color && colorQual) qual = &h.cquals;
			ss << *qual;
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			ss << h.oms;
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			// Look for SNP annotations falling within the alignment
			map<int, char> snpAnnots;
			const size_t len = h.seq.length();
			if(amap != NULL) {
				AnnotationMap::Iter ai = amap->lower_bound(h.h);
				for(; ai != amap->end(); ai++) {
					assert_geq(ai->first.first, h.h.first);
					if(ai->first.first != h.h.first) {
						// Different chromosome
						break;
					}
					if(ai->first.second >= h.h.second + len) {
						// Doesn't fall into alignment
						break;
					}
					if(ai->second.first != 'S') {
						// Not a SNP annotation
						continue;
					}
					size_t off = ai->first.second - h.h.second;
					if(!h.fw) off = len - off - 1;
					snpAnnots[off] = ai->second.second;
				}
			}
			// Output edit column
			const EList<Edit> *es = &h.edits;
			if(gColorEdit) {
				es = &h.ccedits;
			}
			for(size_t i = 0; i < es->size(); i++) {
				const Edit& e = (*es)[i];
				assert(i == es->size()-1 || e.pos <= (*es)[i+1].pos);
				assert_neq(e.chr, e.qchr);
				if(i > 0) ss << ',';
				ss << e.pos << ':' << (char)e.chr;
				while((*es)[i].isInsert() &&
				      i+1 < es->size() &&
				      (*es)[i+1].isInsert() &&
				      (*es)[i+1].pos == (*es)[i].pos)
				{
					i++;
					ss << (char)(*es)[i].chr;
					assert_eq('-', (char)(*es)[i].qchr);
				}
				ss << '>' << (char)e.qchr;
//				if(snpAnnots.find(i) != snpAnnots.end()) {
//					if (!firstmm) ss << ",";
//					ss << i; // position
//					char qryChar = (h.fw ? h.seq[i] : h.seq[h.seq.length()-i-1]);
//					ss << "S:" << snpAnnots[i] << ">" << qryChar;
//					firstmm = false;
//				}
			}
			if(es->empty()) ss << '-';
		}
		if(partition != 0) {
			// Fields addded as of Crossbow 0.1.4
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.mate;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.name;
			}
		}
		if(cost) {
			// Stratum
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.stratum;
			}
			// Cost
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)(h.cost & ~0xc000);
			}
		}
		if(gShowSeed) {
			// Seed
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.seed;
			}
		}
		ss << endl;
	} while(spill);
}

/**
 * Print a list of edits to an OutFileBuf.
 */
static void printEdits(
	const EList<Edit>& es,
	size_t len,
	bool exEnds,
	OutFileBuf& o)
{
	char buf[1024];
	size_t elen = es.size();
	bool first = true;
	int posAdj = exEnds ? -1 : 0;
	for(size_t i = 0; i < elen; i++) {
		const Edit& e = es[i];
		assert(i == elen-1 || e.pos <= es[i+1].pos);
		assert_neq(e.chr, e.qchr);
		assert_lt(e.pos, len);
		if(exEnds && (e.pos == 0 || e.pos == len-1)) {
			continue;
		}
		if(!first) o.write(',');
		first = false;
		itoa10(e.pos + posAdj, buf);
		o.writeChars(buf);
		o.write(':');
		o.write(e.chr);
		while(
			es[i].isInsert() &&
			i+1 < elen &&
			es[i+1].isInsert() &&
			es[i+1].pos == es[i].pos)
		{
			i++;
			o.write(es[i].chr);
			assert_eq('-', (char)es[i].qchr);
		}
		o.write('>');
		o.write(e.qchr);
	}
	if(es.empty()) o.write('-');
}

/**
 * Initialize the wrapper with a new read pair and return an
 * integer >= -1 indicating which stage the aligner should start
 * at.  If -1 is returned, the aligner can skip the read entirely.
 * at.  If .  Checks if the new read pair is identical to the
 * previous pair.  If it is, then we return the id of the first
 * stage to run.
 */
int MSHitSinkWrap::nextRead(
	// One of the other of rd1, rd2 will = NULL if read is unpaired
	const Read* rd1,      // new mate #1
	const Read* rd2,      // new mate #2
	TReadId rdid,         // read ID for new pair
	bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
	assert(!init_);
	assert(rd1 != NULL || rd2 != NULL);
	init_ = true;
	bool same = false;
	if(rd1_ != NULL || rd2_ != NULL) {
		// This is not the first time the sink was initialized with
		// a read.  Check if new read/pair is identical to previous
		// read/pair
		if((rd1_ == NULL) == (rd1 == NULL) &&
		   (rd2_ == NULL) == (rd2 == NULL))
		{
			bool m1same = (rd1 == NULL && rd1_ == NULL);
			if(!m1same) {
				assert(rd1 != NULL);
				assert(rd1_ != NULL);
				m1same = Read::same(
					rd1->patFw,  // new seq
					rd1->qual,   // new quals
					rd1_->patFw, // old seq
					rd1_->qual,  // old quals
					qualitiesMatter);
			}
			if(m1same) {
				bool m2same = (rd2 == NULL && rd2_ == NULL);
				if(!m2same) {
					m2same = Read::same(
						rd2->patFw,  // new seq
						rd2->qual,   // new quals
						rd2_->patFw, // old seq
						rd2_->qual,  // old quals
						qualitiesMatter);
				}
				same = m2same;
			}
		}
	}
	// Keep copy of new read, so that we can compare it with the
	// next one
	if(rd1 != NULL) {
		rd1buf_ = *rd1;
		rd1_ = &rd1buf_;
	} else rd1_ = NULL;
	if(rd2 != NULL) {
		rd2buf_ = *rd2;
		rd2_ = &rd2buf_;
	} else rd2_ = NULL;
	rdid_ = rdid;
	if(same) {
		if(short_) {
			// We were short circuited, so we start from the
			// beginning of the short-circuited stage
			if(lastStage_ > -1) return lastStage_ + 1;
		} else {
			// We were not short circuited, so we can skip the read
			// entirely
			return -1;
		}
	}
	// Caller must now align the read
	maxed1_ = false;
	maxed2_ = false;
	maxedPair_ = false;
	short_ = false;
	lastStage_ = -1;
	best_ = std::numeric_limits<THitInt>::min();
	rs1_.clear();
	rs2_.clear();
	rs1u_.clear();
	rs2u_.clear();
	assert(empty());
	assert(!maxed());
	// Start from the first stage
	return 0;
}

/**
 * Finish reporting for the read in rd1_ or in rd2_, depending on
 * how 'one' is set.
 *
 * Called by finishRead.
 */
void MSHitSinkWrap::finishGroup(
	bool paired,        // true iff alns being reported are paired
	bool condord,       // true iff paired-end alns are concordant
	bool one,           // true iff unpaired alns are from mate 1
	RandomSource& rnd,  // pseudo-random generator
	uint64_t& al,       // counter to inc for reads that align
	uint64_t& mx,       // counter to inc for reads that align repetitively
	uint64_t& un) const // counter to inc for reads that don't align
{
	const Read *rda, *rdb;
	const EList<AlnRes> *rsa, *rsb;
	bool maxed;
	if(paired) {
		assert(rd1_ != NULL);
		assert(rd2_ != NULL);
		rda = rd1_;
		rdb = rd2_;
		assert(rs1_.size() == rs2_.size());
		rsa = &rs1_;
		rsb = &rs2_;
		maxed = maxedPair_;
	} else {
		rda = one ? rd1_ : rd2_;
		assert(rda != NULL);
		rdb = NULL;
		rsa = one ? &rs1u_ : &rs2u_;
		rsb = NULL;
		maxed = one ? maxed1_ : maxed2_;
	}
	assert(rda != NULL);
	assert(rsa != NULL);
	if(!rsa->empty()) {
		AlnSetSumm summ(rda, rdb, rsa, rsb);
		if(maxed && !rp_.msample) {
			g_.reportMaxed(rda, rdb, rdid_, rsa, rsb, summ);
			mx++;
			al++;
		} else {
			// Make any random selections that may be necessary
			selectAlnsToReport(*rsa, const_cast<EList<bool>& >(select_), rnd);
			if(maxed) {
				// Read exceded -M ceiling
				assert(rp_.msample);
				mx++;
			}
			al++;
			g_.reportHits(rda, rdb, rdid_, select_, rsa, rsb, maxed, summ);
		}
	} else {
		g_.reportUnaligned(rda, rdb, rdid_, true);
		un++;
	}
}

/**
 * Inform global, shared MSHitSink object that we're finished with
 * this read.  The global MSHitSink is responsible for updating
 * counters, creating the output record, and delivering the record
 * to the appropriate output stream.
 */
void MSHitSinkWrap::finishRead(
	const SeedResults *sr1,
	const SeedResults *sr2,
	RandomSource&      rnd,
	ReportingMetrics&  met,
	bool suppressSeedSummary, // = true
	bool suppressAlignments)  // = false
{
	if(!suppressSeedSummary) {
		if(sr1 != NULL) {
			assert(rd1_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(*rd1_, rdid_, *sr1, true);
		} else if(rd1_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(*rd1_, rdid_, true);
		}
		if(sr2 != NULL) {
			assert(rd2_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(*rd2_, rdid_, *sr2, true);
		} else if(rd2_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(*rd2_, rdid_, true);
		}
	}
	if(!suppressAlignments) {
		// First decide whether we're reporting from among the
		// paired-end or the unpaired alignments
		bool paired = false;  // true -> there are 1 or more
							  // paired-end alignments to report
		bool concord = false; // true -> there are 1 or more
							  // condordant pairs to report
		// Was the input paired?
		if(rd1_ != NULL && rd2_ != NULL) {
			// Do we have any concordant paired-end alignments?
			if(!rs1_.empty()) {
				assert(!rs2_.empty());
				// Yes
				paired = true;
				concord = true;
			} else {
				assert(rs2_.empty());
				// Do we have one or more discordant paired-end
				// alignments?  Currently, we just check to see if
				// we have unique unpaired alignments for both
				// mates.
				if(prepareDiscordants()) {
					// Yes, we have a discordant paired-end alignment
					assert(!rs1_.empty());
					assert_eq(rs1_.size(), rs2_.size());
					paired = true;
					concord = false;
				}
			}
		}
		if(paired) {
			// There are paired-end alignments - report those
			finishGroup(paired, concord,  true, rnd, met.al, met.max, met.unal);
		} else {
			// No paired-end alignments - report unpaired
			finishGroup(paired, concord,  true, rnd, met.al, met.max, met.unal);
			finishGroup(paired, concord, false, rnd, met.al, met.max, met.unal);
		}
	}
	init_ = false;
}

/**
 * Called by the aligner when a new unpaired or paired alignment is
 * discovered in the given stage.  This function checks whether the
 * addition of this alignment causes the reporting policy to be
 * violated (by meeting or exceeding the limits set by -k, -m, -M),
 * in which case true is returned immediately and the aligner is
 * short circuited.  Otherwise, the alignment is tallied and false
 * is returned.
 */
bool MSHitSinkWrap::report(
	int stage,
	const AlnRes* rs1,
	const AlnRes* rs2)
{
	assert(rs1 != NULL || rs2 != NULL);
	assert(rs1 == NULL || !rs1->empty());
	assert(rs2 == NULL || !rs2->empty());
	assert(rs1 == NULL || rs1->repOk());
	assert(rs2 == NULL || rs2->repOk());
	bool paired = (rs1 != NULL && rs2 != NULL);
	bool one = (rs1 != NULL);
	bool maxed = (paired ? maxedPair_ : (one ? maxed1_ : maxed2_));
	const AlnRes* rsa = one ? rs1 : rs2;
	const AlnRes* rsb = one ? rs2 : rs1;
	// We shouldn't have *already* exceeded the mhits ceiling
	// before we got here
	assert(!maxed);
	assert_leq(rs1_.size(),  (size_t)rp_.mhits);
	assert_leq(rs2_.size(),  (size_t)rp_.mhits);
	assert_leq(rs1u_.size(), (size_t)rp_.mhits);
	assert_leq(rs2u_.size(), (size_t)rp_.mhits);
	assert_gt(stage, lastStage_);
	// Copy the alignment result into rs1_/rs2_/rs1u_/rs2u_
	size_t newsz = 0;
	if(paired) {
		rs1_.push_back(*rs1);
		rs2_.push_back(*rs2);
		newsz = rs1_.size();
	} else {
		if(one) {
			rs1u_.push_back(*rs1);
			newsz = rs1u_.size();
		} else {
			rs2u_.push_back(*rs2);
			newsz = rs2u_.size();
		}
	}
	assert_leq(newsz-1, (size_t)rp_.mhits);
	assert(rp_.mhitsSet() || newsz <= (size_t)rp_.khits);
	// Were -m/-M or -k ceilings exceeded?
	if(rp_.mhitsSet() && newsz > (size_t)rp_.mhits) {
		// Set the appropriate maxed field
		if(paired) {
			maxed = maxedPair_ = true;
		} else {
			if(one) {
				maxed = maxed1_ = true;
			} else {
				maxed = maxed2_ = true;
			}
		}
	}
	// Did we max out a category that allows us to short-circuit
	// the alignment process?
	bool inputPaired = (rd1_ != NULL && rd2_ != NULL);
	bool maxShortCircuit =
		(inputPaired && maxedPair_) ||
		(!inputPaired && maxed);
		
	bool kShortCircuit = !rp_.mhitsSet();
	if(kShortCircuit) {
		assert_leq(newsz, rp_.khits);
		if(inputPaired) {
			kShortCircuit = (paired && newsz == rp_.khits);
		} else {
			kShortCircuit = (newsz == rp_.khits);
		}
	}
	bool done = false;
	if(maxShortCircuit || kShortCircuit) {
		// Yes, a ceiling (either -k or -m/-M) was exceeded so we
		// short circuit the aligner here.  We expect that the
		// aligner will do no additional work.  If the next read is
		// identical to this one, we still have to redo stage
		// 'stage' because here we only looked at a subset of
		// alignments before short-circuiting.
		short_ = true;
		done = true;
	}
	// Tally overall alignment score
	THitInt score = rsa->score().score();
	if(rsb != NULL) score += rsb->score().score();
	// Update best score so far
	if(score < best_) {
		best_ = score;
	}
	assert(done || !this->maxed());
	return done;
}

#define NOT_SUPPRESSED !suppress_.test(field++)
#define BEGIN_FIELD { \
	if(firstfield) firstfield = false; \
	else o.write('\t'); \
}

/**
 *
 */
void MSHitSink::reportSeedSummary(
	const Read&        rd,
	TReadId            rdid,
	const SeedResults& rs,
	bool               getLock)
{
	ThreadSafe ts(&locks_[0], getLock);
	appendSeedSummary(
		*(outs_[0]),           // output stream to write to
		rd,                    // read
		rdid,                  // read id
		rs.numOffs()*2,        // # seeds tried
		rs.nonzeroOffsets(),   // # seeds with non-empty results
		rs.numRanges(),        // # ranges for all seed hits
		rs.numElts(),          // # elements for all seed hits
		rs.numOffs(),          // # seeds tried from fw read
		rs.nonzeroOffsetsFw(), // # seeds with non-empty results from fw read
		rs.numRangesFw(),      // # ranges for seed hits from fw read
		rs.numEltsFw(),        // # elements for seed hits from fw read
		rs.numOffs(),          // # seeds tried from rc read
		rs.nonzeroOffsetsRc(), // # seeds with non-empty results from fw read
		rs.numRangesRc(),      // # ranges for seed hits from fw read
		rs.numEltsRc());       // # elements for seed hits from fw read
}

/**
 *
 */
void MSHitSink::reportEmptySeedSummary(
	const Read&        rd,
	TReadId            rdid,
	bool               getLock)
{
	ThreadSafe ts(&locks_[0], getLock);
	appendSeedSummary(
		*(outs_[0]),           // output stream to write to
		rd,                    // read
		rdid,                  // read id
		0,                     // # seeds tried
		0,                     // # seeds with non-empty results
		0,                     // # ranges for all seed hits
		0,                     // # elements for all seed hits
		0,                     // # seeds tried from fw read
		0,                     // # seeds with non-empty results from fw read
		0,                     // # ranges for seed hits from fw read
		0,                     // # elements for seed hits from fw read
		0,                     // # seeds tried from rc read
		0,                     // # seeds with non-empty results from fw read
		0,                     // # ranges for seed hits from fw read
		0);                    // # elements for seed hits from fw read
}

/**
 * Append a batch of unresolved seed alignment summary results (i.e.
 * seed alignments where all we know is the reference sequence aligned
 * to and its SA range, not where it falls in the reference
 * sequence) to the given output stream in Bowtie's seed-sumamry
 * verbose-mode format.
 *
 * The seed summary format is:
 *
 *  - One line per read
 *  - A typical line consists of a set of tab-delimited fields:
 *
 *    1. Read name
 *    2. Total number of seeds extracted from the read
 *    3. Total number of seeds that aligned to the reference at
 *       least once (always <= field 2)
 *    4. Total number of distinct BW ranges found in all seed hits
 *       (always >= field 3)
 *    5. Total number of distinct BW elements found in all seed
 *       hits (always >= field 4)
 *    6-9.:   Like 2-5. but just for seeds extracted from the
 *            forward representation of the read
 *    10-13.: Like 2-5. but just for seeds extracted from the
 *            reverse-complement representation of the read
 *
 *    Note that fields 6 and 10 should add to field 2, 7 and 11
 *    should add to 3, etc.
 *
 *  - Lines for reads that are filtered out for any reason (e.g. too
 *    many Ns) have columns 2 through 13 set to 0.
 */
void MSHitSink::appendSeedSummary(
	OutFileBuf&   o,
	const Read&   rd,
	const TReadId rdid,
	size_t        seedsTried,
	size_t        nonzero,
	size_t        ranges,
	size_t        elts,
	size_t        seedsTriedFw,
	size_t        nonzeroFw,
	size_t        rangesFw,
	size_t        eltsFw,
	size_t        seedsTriedRc,
	size_t        nonzeroRc,
	size_t        rangesRc,
	size_t        eltsRc)
{
	char buf[1024];
	size_t field = 0;
	bool firstfield = true;
	//
	// Read name
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		printUptoWs(o, rd.name, true);
	}
	
	//
	// Total number of seeds tried
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(seedsTried, buf);
		o.writeChars(buf);
	}
	//
	// Total number of seeds tried where at least one range was found.
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(nonzero, buf);
		o.writeChars(buf);
	}
	//
	// Total number of ranges found
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(ranges, buf);
		o.writeChars(buf);
	}
	//
	// Total number of elements found
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(elts, buf);
		o.writeChars(buf);
	}
	
	//
	// The same four numbers, but only for seeds extracted from the
	// forward read representation.
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(seedsTriedFw, buf);
		o.writeChars(buf);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(nonzeroFw, buf);
		o.writeChars(buf);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(rangesFw, buf);
		o.writeChars(buf);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(eltsFw, buf);
		o.writeChars(buf);
	}

	//
	// The same four numbers, but only for seeds extracted from the
	// reverse complement read representation.
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(seedsTriedRc, buf);
		o.writeChars(buf);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(nonzeroRc, buf);
		o.writeChars(buf);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(rangesRc, buf);
		o.writeChars(buf);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		itoa10(eltsRc, buf);
		o.writeChars(buf);
	}
	o.write('\n');
}

/**
 * Append a single hit to the given output stream in Bowtie's
 * verbose-mode format.
 */
void MSVerboseHitSink::appendMate(
	OutFileBuf&   o,
	const Read&   rd,
	const Read*   rdo,
	const TReadId rdid,
	const AlnRes& rs,
	const AlnRes* rso,
	const AlnSetSumm& summ)
{
	bool spill = false;
	int spillAmt = 0;
	int offAdj = ((rd.color && exEnds_) ? 1 : 0);
	size_t rdlen = rd.length();
	TRefOff pdiv = std::numeric_limits<TRefOff>::max();
	uint32_t pmod = 0xffffffff;
	char buf[1024];
	do {
		bool dospill = false;
		if(spill) {
			// The read spilled over a partition boundary in a
			// previous iteration and so needs to be printed again
			// in this iteration
			spill = false;
			dospill = true;
			spillAmt++;
		}
		assert(!spill);
		size_t field = 0;
		bool firstfield = true;
		if(partition_ != 0) {
			int pospart = abs(partition_);
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				// Output a partitioning key
				// First component of the key is the reference index
				if(refnames_ != NULL && rmap_ != NULL) {
					printUptoWs(o, rmap_->getName(rs.refid()), !fullRef_);
				} else if(refnames_ != NULL && rs.refid() < refnames_->size()) {
					printUptoWs(o, (*refnames_)[rs.refid()], !fullRef_);
				} else {
					itoa10<TRefId>(rs.refid(), buf);
					o.writeChars(buf);
				}
			}
			// Next component of the key is the partition id
			if(!dospill) {
				pdiv = (rs.refoff() + offAdj + offBase_) / pospart;
				pmod = (rs.refoff() + offAdj + offBase_) % pospart;
			}
			assert_neq(std::numeric_limits<TRefOff>::max(), pdiv);
			assert_neq(0xffffffff, pmod);
			if(dospill) assert_gt(spillAmt, 0);
			char * buf2 = itoa10<TRefOff>(pdiv + (dospill ? spillAmt : 0), buf);
			size_t pdivLen = buf2 - buf;
			assert_gt(pdivLen, 0);
			// TODO: is colorspace being properly accounted for in the
			// following?
			if(partition_ > 0 &&
			   (pmod + rdlen) >= ((uint32_t)pospart * (spillAmt + 1))) {
				// Spills into the next partition so we need to
				// output another alignment for that partition
				spill = true;
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				// Print partition id with leading 0s so that Hadoop
				// can do lexicographical sort
				size_t partDigits = 1;
				if(pospart >= 10) partDigits++;
				if(pospart >= 100) partDigits++;
				if(pospart >= 1000) partDigits++;
				if(pospart >= 10000) partDigits++;
				if(pospart >= 100000) partDigits++;
				for(size_t i = pdivLen; i < (10-partDigits); i++) {
					o.write('0');
				}
				o.writeChars(buf);
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				// Print offset with leading 0s
				buf2 = itoa10<TRefOff>(rs.refoff() + offAdj + offBase_, buf);
				size_t offLen = buf2 - buf;
				assert_gt(offLen, 0);
				for(size_t i = offLen; i < 9; i++) o.write('0');
				o.writeChars(buf);
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				o.write(rs.refcoord().fw() ? '+' : '-');
			}
			// end if(partition != 0)
		} else {
			assert(!dospill);
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				printUptoWs(o, rd.name, true);
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				o.write(rs.refcoord().fw() ? '+' : '-');
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				// .first is text id, .second is offset
				if(refnames_ != NULL && rmap_ != NULL) {
					printUptoWs(o, rmap_->getName(rs.refid()), !fullRef_);
				} else if(refnames_ != NULL && rs.refid() < refnames_->size()) {
					printUptoWs(o, (*refnames_)[rs.refid()], !fullRef_);
				} else {
					itoa10<TRefId>(rs.refid(), buf);
					o.writeChars(buf);
				}
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				itoa10<TRefOff>(rs.refoff() + offAdj + offBase_, buf);
				o.writeChars(buf);
			}
			// end else clause of if(partition != 0)
		}
		// Set to false until we decode the colorspace alignment, at which point 
		bool decoded = false;
		if(!suppress_.test(field++)) {
			if(firstfield) firstfield = false;
			else o.write('\t');
			bool printColors = rd.color && colorSeq_;
			if(rd.color && !colorSeq_) {
				// decode colorspace alignment
				rs.decodedNucsAndQuals(rd, dseq_, dqual_);
				decoded = true;
			}
			rs.printSeq(
				rd,
				&dseq_,
				printColors,
				exEnds_,
				o);
		}
		if(!suppress_.test(field++)) {
			if(firstfield) firstfield = false;
			else o.write('\t');
			bool printColors = rd.color && colorQual_;
			if(rd.color && !decoded && !colorQual_) {
				// decode colorspace alignment
				rs.decodedNucsAndQuals(rd, dseq_, dqual_);
				decoded = true;
			}
			rs.printQuals(
				rd,
				&dqual_,
				printColors,
				exEnds_,
				o);
		}
		if(!suppress_.test(field++)) {
			if(firstfield) firstfield = false;
			else o.write('\t');
			itoa10<TMapq>(mapq_.mapq(summ), buf);
			o.writeChars(buf);
		}
		if(!suppress_.test(field++)) {
			if(firstfield) firstfield = false;
			else o.write('\t');
			// If ends are being excluded, we need to subtract 1 from
			// .pos's of ned and aed, and exclude elements at the
			// extreme ends.
			printEdits(
				rs.ned(),                   // edits to print
				rdlen + (rd.color ? 1 : 0), // length of read string that edits refer to
				rd.color && exEnds_,        // true -> exclude edits at ends and adjust poss
				o);                         // output stream
		}
		if(partition_ != 0) {
			// Fields addded as of Crossbow 0.1.4
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				itoa10(rd.mate, buf);
				o.writeChars(buf);
			}
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				printUptoWs<BTString>(o, rd.name, true);
			}
		}
		if(printCost_) {
			// Cost
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				itoa10(rs.score().penalty(), buf);
				o.writeChars(buf);
			}
		}
		if(printParams_) {
			// Cost
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				itoa10(rs.seedmms(), buf);
				o.writeChars(buf);
				o.write(',');
				itoa10(rs.seedlen(), buf);
				o.writeChars(buf);
				o.write(',');
				itoa10(rs.seedival(), buf);
				o.writeChars(buf);
				o.write(',');
				itoa10(rs.penceil(), buf);
				o.writeChars(buf);
			}
		}
		if(gShowSeed) {
			// Seed
			if(!suppress_.test(field++)) {
				if(firstfield) firstfield = false;
				else o.write('\t');
				itoa10(rd.seed, buf);
				o.writeChars(buf);
			}
		}
		o.write('\n');
	} while(spill);
}
