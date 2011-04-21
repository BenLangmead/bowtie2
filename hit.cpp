#include <iostream>
#include <iomanip>
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
			uint32_t r = (uint32_t)(rand.nextU32() % num);
			num = 0;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat == bestStratum) {
					if(num == r) {
						hs[i].oms = hs[i+1].oms = (uint32_t)(hs.size()/2);
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
			uint32_t r = (uint32_t)(rand.nextU32() % num);
			Hit& h = hs[r];
			h.oms = (uint32_t)hs.size();
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
	const EList<bool>& suppress)
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
		uint32_t field = 0;
		bool firstfield = true;
		if(partition != 0) {
			int pospart = abs(partition);
			if(!suppress[field++]) {
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
			if(!suppress[field++]) {
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
			if(!suppress[field++]) {
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
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.fw? "+":"-");
			}
			// end if(partition != 0)
		} else {
			assert(!dospill);
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.name;
			}
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.fw? '+' : '-');
			}
			if(!suppress[field++]) {
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
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.h.second + offBase);
			}
			// end else clause of if(partition != 0)
		}
		if(!suppress[field++]) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			const BTDnaString* pat = &h.seq;
			if(h.color && colorSeq) pat = &h.cseq;
			ss << *pat;
		}
		if(!suppress[field++]) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			const BTString* qual = &h.quals;
			if(h.color && colorQual) qual = &h.cquals;
			ss << *qual;
		}
		if(!suppress[field++]) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			ss << h.oms;
		}
		if(!suppress[field++]) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			// Look for SNP annotations falling within the alignment
			map<size_t, char> snpAnnots;
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
				while((*es)[i].isReadGap() &&
				      i+1 < es->size() &&
				      (*es)[i+1].isReadGap() &&
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
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.mate;
			}
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.name;
			}
		}
		if(cost) {
			// Stratum
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.stratum;
			}
			// Cost
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)(h.cost & ~0xc000);
			}
		}
		if(gShowSeed) {
			// Seed
			if(!suppress[field++]) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.seed;
			}
		}
		ss << endl;
	} while(spill);
}

