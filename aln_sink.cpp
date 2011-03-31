//
//  aln_sink.cpp
//

#include <iomanip>
#include "aln_sink.h"
#include "aligner_seed.h"
#include "util.h"

using namespace std;


/**
 * Initialize state machine with a new read.  The state we start in depends
 * on whether it's paired-end or unpaired.
 */
void ReportingState::nextRead(bool paired) {
	paired_ = paired;
	if(paired) {
		state_ = CONCORDANT_PAIRS;
		doneConcord_ = false;
		doneDiscord_ = p_.discord ? false : true;
		doneUnpair1_ = p_.mixed   ? false : true;
		doneUnpair2_ = p_.mixed   ? false : true;
		exitConcord_ = ReportingState::EXIT_DID_NOT_EXIT;
		exitDiscord_ = p_.discord ?
			ReportingState::EXIT_DID_NOT_EXIT :
			ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair1_ = p_.mixed   ?
			ReportingState::EXIT_DID_NOT_EXIT :
			ReportingState::EXIT_DID_NOT_ENTER;
		exitUnpair2_ = p_.mixed   ?
			ReportingState::EXIT_DID_NOT_EXIT :
			ReportingState::EXIT_DID_NOT_ENTER;
	} else {
		// Unpaired
		state_ = UNPAIRED;
		doneConcord_ = true;
		doneDiscord_ = true;
		doneUnpair1_ = false;
		doneUnpair2_ = true;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER; // not relevant
		exitDiscord_ = ReportingState::EXIT_DID_NOT_ENTER; // not relevant
		exitUnpair1_ = ReportingState::EXIT_DID_NOT_EXIT;
		exitUnpair2_ = ReportingState::EXIT_DID_NOT_ENTER; // not relevant
	}
	doneUnpair_ = doneUnpair1_ && doneUnpair2_;
	done_ = false;
	nconcord_ = ndiscord_ = nunpair1_ = nunpair2_ = 0;
}

/**
 * Caller uses this member function to indicate that one additional
 * concordant alignment has been found.
 */
bool ReportingState::foundConcordant() {
	assert(paired_);
	assert_geq(state_, ReportingState::CONCORDANT_PAIRS);
	assert(!doneConcord_);
	nconcord_++;
	areDone(nconcord_, doneConcord_, exitConcord_);
	// No need to search for discordant alignments if there are one or more
	// concordant alignments.
	doneDiscord_ = true;
	exitDiscord_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
	if(doneConcord_) {
		// TODO: If we're finished looking for concordant alignments, do we
		// have to continue on to search for unpaired alignments?  Only if
		// our exit from the concordant stage is EXIT_SHORT_CIRCUIT_m.  If
		// it's EXIT_SHORT_CIRCUIT_k or EXIT_SHORT_CIRCUIT_M or
		// EXIT_WITH_ALIGNMENTS, we can skip unpaired.
		assert_neq(ReportingState::EXIT_NO_ALIGNMENTS, exitConcord_);
		if(exitConcord_ != ReportingState::EXIT_SHORT_CIRCUIT_m) {
			if(!doneUnpair1_) {
				doneUnpair1_ = true;
				exitUnpair1_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
			}
			if(!doneUnpair2_) {
				doneUnpair2_ = true;
				exitUnpair2_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
			}
		}
	}
	updateDone();
	return done();
}

/**
 * Caller uses this member function to indicate that one additional
 * discordant alignment has been found.
 */
bool ReportingState::foundDiscordant() {
	assert(paired_);
	assert_gt(state_, ReportingState::NO_READ);
	ndiscord_++;
	// There can only be one discordant alignment per paired-end read, so
	// there's no need to search for any more.
	assert(!doneDiscord_);
	doneDiscord_ = true;
	exitDiscord_ = ReportingState::EXIT_WITH_ALIGNMENTS;
	// If there are any discordant alignments found, there can't be any
	// unpaired alignments reported.
	if(!doneUnpair1_) {
		doneUnpair1_ = true;
		exitUnpair1_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
	}
	if(!doneUnpair2_) {
		doneUnpair2_ = true;
		exitUnpair2_ = ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED;
	}
	updateDone();
	return done();
}

/**
 * Caller uses this member function to indicate that one additional
 * discordant alignment has been found.
 */
bool ReportingState::foundUnpaired(bool mate1) {
	assert_gt(state_, ReportingState::NO_READ);
	// Note: it's not right to assert !doneUnpair1_/!doneUnpair2_ here.
	// Even if we're done with finding 
	if(mate1) {
		nunpair1_++;
		// Did we just finish with this mate?
		if(!doneUnpair1_) {
			areDone(nunpair1_, doneUnpair1_, exitUnpair1_);
			if(nunpair1_ > 1 || nunpair2_ > 1) {
				doneDiscord_ = true;
				exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
			}
			if(doneUnpair1_) {
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		} else {
			// If we've already ruled out reporting unpaired alignments for
			// this mate, then this discovery must be in aid of finding
			// paired-end alignments.
			assert(!doneConcord_ || !doneDiscord_);
		}
	} else {
		nunpair2_++;
		// Did we just finish with this mate?
		if(!doneUnpair2_) {
			areDone(nunpair2_, doneUnpair2_, exitUnpair2_);
			if(nunpair1_ > 1 || nunpair2_ > 1) {
				doneDiscord_ = true;
				exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
			}
			if(doneUnpair2_) {
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		} else {
			// If we've already ruled out reporting unpaired alignments for
			// this mate, then this discovery must be in aid of finding
			// paired-end alignments.
			assert(!doneConcord_ || !doneDiscord_);
		}
	}
	return done();
}
	
/**
 * Called to indicate that the aligner has finished searching for
 * alignments.  This gives us a chance to finalize our state.
 *
 * TODO: Keep track of short-circuiting information.
 */
void ReportingState::finish() {
	if(!doneConcord_) {
		doneConcord_ = true;
		exitConcord_ =
			((nconcord_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitConcord_, EXIT_DID_NOT_EXIT);
	if(!doneUnpair1_) {
		doneUnpair1_ = true;
		exitUnpair1_ =
			((nunpair1_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitUnpair1_, EXIT_DID_NOT_EXIT);
	if(!doneUnpair2_) {
		doneUnpair2_ = true;
		exitUnpair2_ =
			((nunpair2_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert_gt(exitUnpair2_, EXIT_DID_NOT_EXIT);
	if(!doneDiscord_) {
		// Check if the unpaired alignments should be converted to a single
		// discordant paired-end alignment.
		assert_eq(0, ndiscord_);
		if(nconcord_ == 0 && nunpair1_ == 1 && nunpair2_ == 1) {
			convertUnpairedToDiscordant();
		}
		doneDiscord_ = true;
		exitDiscord_ =
			((ndiscord_ > 0) ?
				ReportingState::EXIT_WITH_ALIGNMENTS :
				ReportingState::EXIT_NO_ALIGNMENTS);
	}
	assert(!paired_ || exitDiscord_ > ReportingState::EXIT_DID_NOT_EXIT);
	doneUnpair_ = done_ = true;
	assert(done());
}
	
/**
 * Populate given counters with the number of various kinds of alignments
 * to report for this read.  Concordant alignments are preferable to (and
 * mutually exclusive with) discordant alignments, and paired-end
 * alignments are preferable to unpaired alignments.
 *
 * The caller also needs some additional information for the case where a
 * pair or unpaired read aligns repetitively.  If the read is paired-end
 * and the paired-end has repetitive concordant alignments, that should be
 * reported, and 'pairMax' is set to true to indicate this.  If the read is
 * paired-end, does not have any conordant alignments, but does have
 * repetitive alignments for one or both mates, then that should be
 * reported, and 'unpair1Max' and 'unpair2Max' are set accordingly.
 *
 * Note that it's possible in the case of a paired-end read for the read to
 * have repetitive concordant alignments, but for one mate to have a unique
 * unpaired alignment.
 */
void ReportingState::getReport(
	uint64_t& nconcordAln, // # concordant alignments to report
	uint64_t& ndiscordAln, // # discordant alignments to report
	uint64_t& nunpair1Aln, // # unpaired alignments for mate #1 to report
	uint64_t& nunpair2Aln, // # unpaired alignments for mate #2 to report
	bool& pairMax,         // repetitive concordant alignments
	bool& unpair1Max,      // repetitive alignments for mate #1
	bool& unpair2Max)      // repetitive alignments for mate #2
	const
{
	nconcordAln = ndiscordAln = nunpair1Aln = nunpair2Aln = 0;
	pairMax = unpair1Max = unpair2Max = false;
	assert_gt(p_.khits, 0);
	assert_gt(p_.mhits, 0);
	if(paired_) {
		// Do we have 1 or more concordant alignments to report?
		if(exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
			// k at random
			assert_geq(nconcord_, (uint64_t)p_.khits);
			nconcordAln = p_.khits;
			return;
		} else if(exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_M) {
			assert(p_.msample);
			assert_gt(nconcord_, 0);
			pairMax = true;  // repetitive concordant alignments
			nconcordAln = 1; // 1 at random
			return;
		} else if(exitConcord_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
			assert_gt(nconcord_, 0);
			// <= k at random
			nconcordAln = min<uint64_t>(nconcord_, p_.khits);
			return;
		}
		if(exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_m) {
			assert(!p_.msample);
			pairMax = true;  // repetitive concordant alignments
		} else {
			assert(!p_.mhitsSet() || nconcord_ <= (uint64_t)p_.mhits+1);
		}
		
		// Do we have a discordant alignment to report?
		if(exitDiscord_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
			// Report discordant
			assert(p_.discord);
			ndiscordAln = 1;
			return;
		}
	}
	
	assert_neq(ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED, exitUnpair1_);
	assert_neq(ReportingState::EXIT_SHORT_CIRCUIT_TRUMPED, exitUnpair2_);

	if((paired_ && !p_.mixed) || nunpair1_ + nunpair2_ == 0) {
		return;
	}

	// Do we have 1 or more alignments for mate #1 to report?
	if(exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
		// k at random
		assert_geq(nunpair1_, (uint64_t)p_.khits);
		nunpair1Aln = p_.khits;
	} else if(exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_M) {
		assert(p_.msample);
		assert_gt(nunpair1_, 0);
		unpair1Max = true;  // repetitive alignments for mate #1
		nunpair1Aln = 1; // 1 at random
	} else if(exitUnpair1_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
		assert_gt(nunpair1_, 0);
		// <= k at random
		nunpair1Aln = min<uint64_t>(nunpair1_, (uint64_t)p_.khits);
	}
	if(exitUnpair1_ == ReportingState::EXIT_SHORT_CIRCUIT_m) {
		assert(!p_.msample);
		unpair1Max = true;  // repetitive alignments for mate #1
	} else {
		assert(!p_.mhitsSet() || paired_ || nunpair1_ <= (uint64_t)p_.mhits+1);
	}

	// Do we have 2 or more alignments for mate #2 to report?
	if(exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
		// k at random
		nunpair2Aln = p_.khits;
	} else if(exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_M) {
		assert(p_.msample);
		assert_gt(nunpair2_, 0);
		unpair2Max = true;  // repetitive alignments for mate #1
		nunpair2Aln = 1; // 1 at random
	} else if(exitUnpair2_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
		assert_gt(nunpair2_, 0);
		// <= k at random
		nunpair2Aln = min<uint64_t>(nunpair2_, (uint64_t)p_.khits);
	}
	if(exitUnpair2_ == ReportingState::EXIT_SHORT_CIRCUIT_m) {
		assert(!p_.msample);
		unpair2Max = true;  // repetitive alignments for mate #2
	} else {
		assert(!p_.mhitsSet() || paired_ || nunpair2_ <= (uint64_t)p_.mhits+1);
	}
}

/**
 * Given the number of alignments in a category, check whether we
 * short-circuited out of the category.  Set the done and exit arguments to
 * indicate whether and how we short-circuited.
 */
inline void ReportingState::areDone(
	uint64_t cnt,    // # alignments in category
	bool& done,      // out: whether we short-circuited out of category
	int& exit) const // out: if done, how we short-circuited (-k? -m? etc)
{
	assert(!done);
	// Have we exceeded the -k limit?
	assert_gt(p_.khits, 0);
	assert_gt(p_.mhits, 0);
	if(cnt >= (uint64_t)p_.khits && !p_.mhitsSet()) {
		done = true;
		exit = ReportingState::EXIT_SHORT_CIRCUIT_k;
	}
	// Have we exceeded the -m or -M limit?
	else if(p_.mhitsSet() && cnt > (uint64_t)p_.mhits) {
		done = true;
		exit = p_.msample ?
			ReportingState::EXIT_SHORT_CIRCUIT_M :
			ReportingState::EXIT_SHORT_CIRCUIT_m;
	}
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
void AlnSink::printAlSumm(
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
 * Return true iff the read in rd1/rd2 matches the last read handled, which
 * should still be in rd1_/rd2_.
 */
bool AlnSinkWrap::sameRead(
	// One of the other of rd1, rd2 will = NULL if read is unpaired
	const Read* rd1,      // new mate #1
	const Read* rd2,      // new mate #2
	bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
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
	return same;
}

/**
 * Initialize the wrapper with a new read pair and return an integer >= -1
 * indicating which stage the aligner should start at.  If -1 is returned, the
 * aligner can skip the read entirely.  Checks if the new read pair is
 * identical to the previous pair.  If it is, then we return the id of the
 * first stage to run.
 */
int AlnSinkWrap::nextRead(
	// One of the other of rd1, rd2 will = NULL if read is unpaired
	const Read* rd1,      // new mate #1
	const Read* rd2,      // new mate #2
	TReadId rdid,         // read ID for new pair
	bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
	assert(!init_);
	assert(rd1 != NULL || rd2 != NULL);
	init_ = true;
	bool same = sameRead(rd1, rd2, qualitiesMatter);
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
#if 0
		if(short_) {
			// We were short circuited, so we start from the
			// beginning of the short-circuited stage
			if(lastStage_ > -1) return lastStage_ + 1;
		} else {
			// We were not short circuited, so we can skip the read
			// entirely
			return -1;
		}
#endif
		// don't skip anything.  Need to think harder about how to do this.
	}
	// Caller must now align the read
	maxed1_ = false;
	maxed2_ = false;
	maxedOverall_ = false;
	best_ = std::numeric_limits<THitInt>::min();
	rs1_.clear();     // clear out paired-end alignments
	rs2_.clear();     // clear out paired-end alignments
	rs1u_.clear();    // clear out unpaired alignments for mate #1
	rs2u_.clear();    // clear out unpaired alignments for mate #2
	st_.nextRead(readIsPair()); // reset state
	assert(empty());
	assert(!maxed());
	// Start from the first stage
	return 0;
}

/**
 * Inform global, shared AlnSink object that we're finished with this read.
 * The global AlnSink is responsible for updating counters, creating the output
 * record, and delivering the record to the appropriate output stream.
 */
void AlnSinkWrap::finishRead(
	const SeedResults *sr1,
	const SeedResults *sr2,
	RandomSource&      rnd,
	ReportingMetrics&  met,
	bool suppressSeedSummary, // = true
	bool suppressAlignments)  // = false
{
	assert(init_);
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
		// Ask the ReportingState what to report
		st_.finish();
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		st_.getReport(
			nconcord,
			ndiscord,
			nunpair1,
			nunpair2,
			pairMax,
			unpair1Max,
			unpair2Max);
		assert_leq(nconcord, rs1_.size());
		assert_leq(nunpair1, rs1u_.size());
		assert_leq(nunpair2, rs2u_.size());
		assert_leq(ndiscord, 1);
		assert_gt(rp_.khits, 0);
		assert_gt(rp_.mhits, 0);
		assert(!pairMax    || rs1_.size()  >= (uint64_t)rp_.mhits);
		assert(!unpair1Max || rs1u_.size() >= (uint64_t)rp_.mhits);
		assert(!unpair2Max || rs2u_.size() >= (uint64_t)rp_.mhits);
		// Report concordant paired-end alignments if possible
		if(nconcord > 0) {
			AlnSetSumm concordSumm(rd1_, rd2_, &rs1_, &rs2_);
			AlnFlags flags(ALN_FLAG_PAIR_CONCORD, pairMax, pairMax);
			// Possibly select a random subset
			selectAlnsToReport(rs1_, nconcord, select_, rnd);
			g_.reportHits(
				rd1_,
				rd2_,
				rdid_,
				select_,
				&rs1_,
				&rs2_,
				pairMax,
				concordSumm,
				flags);
		}
		// Report concordant paired-end alignments if possible
		else if(ndiscord > 0) {
			ASSERT_ONLY(bool ret =) prepareDiscordants();
			assert(ret);
			assert_eq(1, rs1_.size());
			assert_eq(1, rs2_.size());
			AlnSetSumm discordSumm(rd1_, rd2_, &rs1_, &rs2_);
			AlnFlags flags(ALN_FLAG_PAIR_DISCORD, false, pairMax);
			selectAlnsToReport(rs1_, ndiscord, select_, rnd);
			g_.reportHits(
				rd1_,
				rd2_,
				rdid_,
				select_,
				&rs1_,
				&rs2_,
				pairMax,
				discordSumm,
				flags);
		}
		// Report unpaired alignments if possbile
		else if(nunpair1 > 0 || nunpair2 > 0) {
			if(nunpair1 > 0) {
				AlnSetSumm unpair1Summ(rd1_, NULL, &rs1u_, NULL);
				AlnFlags flags(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR :
						ALN_FLAG_PAIR_UNPAIRED,
					unpair1Max,
					pairMax);
				// Possibly select a random subset
				selectAlnsToReport(rs1u_, nunpair1, select_, rnd);
				g_.reportHits(
					rd1_,
					NULL,
					rdid_,
					select_,
					&rs1u_,
					NULL,
					unpair1Max,
					unpair1Summ,
					flags);
			}
			if(nunpair2 > 0) {
				AlnSetSumm unpair2Summ(NULL, rd2_, NULL, &rs2u_);
				AlnFlags flags(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR :
						ALN_FLAG_PAIR_UNPAIRED,
					unpair2Max,
					pairMax);
				// Possibly select a random subset
				selectAlnsToReport(rs2u_, nunpair2, select_, rnd);
				g_.reportHits(
					NULL,
					rd2_,
					rdid_,
					select_,
					NULL,
					&rs2u_,
					unpair2Max,
					unpair2Summ,
					flags);
			}
		}
		// No alignments to report
		else if(pairMax || unpair1Max || unpair2Max) {
			// Pair aligned repetitively & concordantly?
			if(pairMax) {
				AlnSetSumm concordSumm(rd1_, rd2_, &rs1_, &rs2_);
				AlnFlags flags(ALN_FLAG_PAIR_CONCORD, pairMax, pairMax);
				g_.reportMaxed(
					rd1_,
					rd2_,
					rdid_,
					&rs1_,
					&rs2_,
					concordSumm,
					flags);
			}
			// Either mate aligned repetitively?
			else {
				if(unpair1Max) {
					AlnSetSumm unpair1Summ(rd1_, NULL, &rs1u_, NULL);
					AlnFlags flags(
						readIsPair() ?
							ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR :
							ALN_FLAG_PAIR_UNPAIRED,
						unpair1Max,
						pairMax);
					g_.reportMaxed(
						rd1_,
						NULL,
						rdid_,
						&rs1u_,
						NULL,
						unpair1Summ,
						flags);
				}
				if(unpair2Max) {
					AlnSetSumm unpair2Summ(NULL, rd2_, NULL, &rs2u_);
					AlnFlags flags(
						readIsPair() ?
							ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR :
							ALN_FLAG_PAIR_UNPAIRED,
						unpair2Max,
						pairMax);
					g_.reportMaxed(
						NULL,
						rd2_,
						rdid_,
						NULL,
						&rs2u_,
						unpair2Summ,
						flags);
				}
			}
		}
		// No alignments at all
		else {
			if(rd1_ != NULL) {
				AlnSetSumm summ(rd1_, NULL, NULL, NULL);
				AlnFlags flags(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR :
						ALN_FLAG_PAIR_UNPAIRED,
					false,
					false);
				g_.reportUnaligned(rd1_, rd2_, rdid_, summ, flags, true);
			}
			if(rd2_ != NULL) {
				AlnSetSumm summ(NULL, rd2_, NULL, NULL);
				AlnFlags flags(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR :
						ALN_FLAG_PAIR_UNPAIRED,
					false,
					false);
				g_.reportUnaligned(rd1_, rd2_, rdid_, summ, flags, true);
			}
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
bool AlnSinkWrap::report(
	int stage,
	const AlnRes* rs1,
	const AlnRes* rs2)
{
	assert(init_);
	assert(rs1 != NULL || rs2 != NULL);
	assert(rs1 == NULL || !rs1->empty());
	assert(rs2 == NULL || !rs2->empty());
	assert(rs1 == NULL || rs1->repOk());
	assert(rs2 == NULL || rs2->repOk());
	bool paired = (rs1 != NULL && rs2 != NULL);
	bool one = (rs1 != NULL);
	const AlnRes* rsa = one ? rs1 : rs2;
	const AlnRes* rsb = one ? rs2 : rs1;
	assert(!st_.done());
	if(paired) {
		assert(readIsPair());
		st_.foundConcordant();
		rs1_.push_back(*rs1);
		rs2_.push_back(*rs2);
	} else {
		st_.foundUnpaired(one);
		if(one) {
			rs1u_.push_back(*rs1);
		} else {
			rs2u_.push_back(*rs2);
		}
	}
	// Tally overall alignment score
	THitInt score = rsa->score().score();
	if(rsb != NULL) score += rsb->score().score();
	// Update best score so far
	if(score < best_) {
		best_ = score;
	}
	return st_.done();
}

/**
 * If there is a configuration of unpaired alignments that fits our
 * criteria for there being one or more discordant alignments, then
 * shift the discordant alignments over to the rs1_/rs2_ lists, clear the
 * rs1u_/rs2u_ lists and return true.  Otherwise, return false.
 */
bool AlnSinkWrap::prepareDiscordants() {
	if(rs1u_.size() == 1 && rs2u_.size() == 1) {
		assert(rs1_.empty());
		assert(rs2_.empty());
		rs1_.push_back(rs1u_[0]);
		rs2_.push_back(rs2u_[0]);
		rs1u_.clear();
		rs2u_.clear();
		return true;
	}
	return false;
}

/**
 * Given that rs is already populated with alignments, consider the
 * alignment policy and make random selections where necessary.  E.g. if we
 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
 * random.  We "select" an alignment by setting the parallel entry in the
 * 'select' list to true.
 */
void AlnSinkWrap::selectAlnsToReport(
	const EList<AlnRes>& rs,     // alignments to select from
	uint64_t             num,    // number of alignments to select
	EList<bool>&         select, // list to put results in
	RandomSource&        rnd)
	const
{
	assert(init_);
	assert(repOk());
	uint64_t sz = rs.size();
	if(sz < 1) return;
	select.resize(sz);
	if(num < sz) {
		// Select a random offset into the list of alignments
		uint32_t off = rnd.nextU32() % (uint32_t)sz;
		size_t take = (size_t)num;
		// Now take rp_.khits elements starting at that offset,
		// wrapping back to 0 if necessary, and leave the rest.
		for(size_t i = 0; i < sz; i++) {
			off++;
			if(off == sz) off = 0;
			select[off] = i < take;
		}
	} else {
		// Select them all!  No randomness needed.
		for(size_t i = 0; i < sz; i++) {
			select[i] = true;
		}
	}
}

#define NOT_SUPPRESSED !suppress_[field++]
#define BEGIN_FIELD { \
	if(firstfield) firstfield = false; \
	else o.write('\t'); \
}
#define WRITE_TAB { \
	if(firstfield) firstfield = false; \
	else o.write('\t'); \
}
#define WRITE_NUM(o, x) { \
	itoa10(x, buf); \
	o.writeChars(buf); \
}

/**
 * Print a seed summary to the first output stream in the outs_ list.
 */
void AlnSink::reportSeedSummary(
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
 * Print an empty seed summary to the first output stream in the outs_ list.
 */
void AlnSink::reportEmptySeedSummary(
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
void AlnSink::appendSeedSummary(
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
	uint32_t field = 0;
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
		WRITE_NUM(o, seedsTried);
	}
	//
	// Total number of seeds tried where at least one range was found.
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, nonzero);
	}
	//
	// Total number of ranges found
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, ranges);
	}
	//
	// Total number of elements found
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, elts);
	}
	
	//
	// The same four numbers, but only for seeds extracted from the
	// forward read representation.
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, seedsTriedFw);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, nonzeroFw);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, rangesFw);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, eltsFw);
	}

	//
	// The same four numbers, but only for seeds extracted from the
	// reverse complement read representation.
	//
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, seedsTriedRc);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, nonzeroRc);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, rangesRc);
	}
	if(NOT_SUPPRESSED) {
		BEGIN_FIELD;
		WRITE_NUM(o, eltsRc);
	}
	o.write('\n');
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
 * Append a single hit to the given output stream in Bowtie's
 * verbose-mode format.
 */
void AlnSinkVerbose::appendMate(
	OutFileBuf&   o,
	const Read&   rd,
	const Read*   rdo,
	const TReadId rdid,
	const AlnRes* rs,
	const AlnRes* rso,
	const AlnSetSumm& summ,
	const AlnFlags& flags)
{
	if(rs == NULL && !printPlaceholders_) return;
	bool spill = false;
	int spillAmt = 0;
	int offAdj = ((rd.color && exEnds_) ? 1 : 0);
	if(rs == NULL) offAdj = 0;
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
		uint32_t field = 0;
		bool firstfield = true;
		if(partition_ != 0) {
			int pospart = abs(partition_);
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				if(rs != NULL) {
					// Output a partitioning key
					// First component of the key is the reference index
					if(refnames_ != NULL && rmap_ != NULL) {
						printUptoWs(o, rmap_->getName(rs->refid()), !fullRef_);
					} else if(refnames_ != NULL && rs->refid() < refnames_->size()) {
						printUptoWs(o, (*refnames_)[rs->refid()], !fullRef_);
					} else {
						itoa10<TRefId>(rs->refid(), buf);
						o.writeChars(buf);
					}
				} else {
					o.write('*');
				}
			}
			TRefOff off = rs != NULL ? rs->refoff() : 0;
			pdiv = pmod = 0;
			size_t pdivLen = 0;
			char * buf2;
			if(rs != NULL) {
				// Next component of the key is the partition id
				if(!dospill) {
					pdiv = (uint32_t)((rs->refoff() + offAdj + offBase_) / pospart);
					pmod = (uint32_t)((rs->refoff() + offAdj + offBase_) % pospart);
				}
				assert_neq(std::numeric_limits<TRefOff>::max(), pdiv);
				assert_neq(0xffffffff, pmod);
				if(dospill) assert_gt(spillAmt, 0);
				// TODO: is colorspace being properly accounted for in the
				// following?
				if(partition_ > 0 &&
				   (pmod + rdlen) >= ((uint32_t)pospart * (spillAmt + 1))) {
					// Spills into the next partition so we need to
					// output another alignment for that partition
					spill = true;
				}
			}
			buf2 = itoa10<TRefOff>(pdiv + (dospill ? spillAmt : 0), buf);
			pdivLen = buf2 - buf;
			assert_gt(pdivLen, 0);
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
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
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				// Print offset with leading 0s
				int base = (rs == NULL ? 0 : offBase_);
				buf2 = itoa10<TRefOff>(off + offAdj + base, buf);
				size_t offLen = buf2 - buf;
				assert_gt(offLen, 0);
				for(size_t i = offLen; i < 9; i++) o.write('0');
				o.writeChars(buf);
			}
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				if(rs == NULL) {
					o.write('*');
				} else {
					o.write(rs->refcoord().fw() ? '+' : '-');
				}
			}
			// end if(partition != 0)
		} else {
			assert(!dospill);
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				printUptoWs(o, rd.name, true);
			}
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				if(rs == NULL) {
					o.write('*');
				} else {
					o.write(rs->refcoord().fw() ? '+' : '-');
				}
			}
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				// .first is text id, .second is offset
				if(rs != NULL) {
					if(refnames_ != NULL && rmap_ != NULL) {
						printUptoWs(o, rmap_->getName(rs->refid()), !fullRef_);
					} else if(refnames_ != NULL && rs->refid() < refnames_->size()) {
						printUptoWs(o, (*refnames_)[rs->refid()], !fullRef_);
					} else {
						itoa10<TRefId>(rs->refid(), buf);
						o.writeChars(buf);
					}
				} else {
					o.write('*');
				}
			}
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				if(rs != NULL) {
					TRefOff off = rs->refoff();
					itoa10<TRefOff>(off + offAdj + offBase_, buf);
					o.writeChars(buf);
				} else {
					o.write('*');
				}
			}
			// end else clause of if(partition != 0)
		}
		// Set to false until we decode the colorspace alignment, at which point 
		bool decoded = false;
		if(NOT_SUPPRESSED) {
			WRITE_TAB;
			bool printColors = rd.color && colorSeq_;
			if(rs != NULL && rd.color && !colorSeq_) {
				// decode colorspace alignment
				rs->decodedNucsAndQuals(rd, dseq_, dqual_);
				decoded = true;
			}
			if(rs != NULL) {
				rs->printSeq(
					rd,
					&dseq_,
					printColors,
					exEnds_,
					o);
			} else {
				// Print the read
				o.writeChars(rd.patFw.toZBuf());
			}
		}
		if(NOT_SUPPRESSED) {
			WRITE_TAB;
			bool printColors = rd.color && colorQual_;
			if(rs != NULL && rd.color && !decoded && !colorQual_) {
				// decode colorspace alignment
				rs->decodedNucsAndQuals(rd, dseq_, dqual_);
				decoded = true;
			}
			if(rs != NULL) {
				rs->printQuals(
					rd,
					&dqual_,
					printColors,
					exEnds_,
					o);
			} else {
				// Print the quals
				o.writeChars(rd.qual.toZBuf());
			}
		}
		if(NOT_SUPPRESSED) {
			WRITE_TAB;
			if(rs != NULL) {
				itoa10<TMapq>(mapq_.mapq(summ), buf);
				o.writeChars(buf);
			} else o.write('*');
		}
		if(NOT_SUPPRESSED) {
			WRITE_TAB;
			// If ends are being excluded, we need to subtract 1 from
			// .pos's of ned and aed, and exclude elements at the
			// extreme ends.
			if(rs != NULL) {
				printEdits(
					rs->ned(),                  // edits to print
					rdlen + (rd.color ? 1 : 0), // length of read string that edits refer to
					rd.color && exEnds_,        // true -> exclude edits at ends and adjust poss
					o);                         // output stream
			} else o.write('*');
		}
		if(partition_ != 0) {
			// Fields addded as of Crossbow 0.1.4
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				WRITE_NUM(o, rd.mate);
			}
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				printUptoWs<BTString>(o, rd.name, true);
			}
		}
		if(printFlags_) {
			// Print alignment flags, including:
			//
			// a. Whether this is a (i) half a concordant paired-end alignment,
			//    (ii) half a discordant paired-end alignment, (iii) an
			//    unpaired alignment
			// b. Whether the alignment is (i) itself repetitive, or (ii) is
			//    associated with a paired-end read that has repetitive
			//    concordant alignments
			// c. Whether alignment was found using BW-DP or Mate-DP
			//
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				bool first = true;
				if(flags.maxed()) {
					first = false;
					o.writeChars("XM:1");
				}
				if(flags.maxedPair()) {
					if(!first) o.write(',');
					first = false;
					o.writeChars("XP:1");
				}
				if(!first) o.write(',');
				o.writeChars("XT:");
				if(flags.pairing() == ALN_FLAG_PAIR_CONCORD) {
					o.writeChars("CP");
				} else if(flags.pairing() == ALN_FLAG_PAIR_DISCORD) {
					o.writeChars("DP");
				} else if(flags.pairing() == ALN_FLAG_PAIR_UNPAIRED_FROM_PAIR) {
					o.writeChars("UP");
				} else if(flags.pairing() == ALN_FLAG_PAIR_UNPAIRED) {
					o.writeChars("UU");
				}
			}
		}
		if(printCost_) {
			// Cost
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				if(rs != NULL) {
					WRITE_NUM(o, rs->score().penalty());
				} else o.write('*');
			}
		}
		if(printParams_) {
			if(NOT_SUPPRESSED) {
				WRITE_TAB;
				if(rs != NULL) {
					WRITE_NUM(o, rs->seedmms());  o.write(',');
					WRITE_NUM(o, rs->seedlen());  o.write(',');
					WRITE_NUM(o, rs->seedival()); o.write(',');
					WRITE_NUM(o, rs->penceil());
				} else o.write('*');
			}
		}
		o.write('\n');
	} while(spill);
}

#ifdef ALN_SINK_MAIN

#include <iostream>

bool testDones(
	const ReportingState& st,
	bool done1,
	bool done2,
	bool done3,
	bool done4,
	bool done5,
	bool done6)
{
	assert(st.doneConcordant()    == done1);
	assert(st.doneDiscordant()    == done2);
	assert(st.doneUnpaired(true)  == done3);
	assert(st.doneUnpaired(false) == done4);
	assert(st.doneUnpaired()      == done5);
	assert(st.done()              == done6);
	assert(st.repOk());
	return true;
}

int main(void) {
	cerr << "Case 1 (simple unpaired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			0,      // mhits
			0,      // pengap
			false,  // msample
			false,  // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(false); // unpaired read
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, true, true, true, true));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(2, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (simple unpaired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			false,  // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(false); // unpaired read
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, false, true, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 3 (simple paired 1) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			false,  // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(!unpair1Max); // because !mixed
		assert(!unpair2Max); // because !mixed
	}
	cerr << "PASSED" << endl;

	cerr << "Case 4 (simple paired 2) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, false, false, false));
		st.foundUnpaired(false);
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, false, true, true, true, true, false));
		st.foundConcordant();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(4, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(4, st.numUnpaired1());
		assert_eq(4, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(unpair1Max);
		assert(unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 5 (potential discordant after concordant) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		st.foundUnpaired(false);
		st.foundConcordant();
		assert(testDones(st, false, true, false, false, false, false));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(1, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(1, st.numUnpaired1());
		assert_eq(1, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(1, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;

	cerr << "Case 6 (true discordant) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			2,      // khits
			3,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		st.foundUnpaired(false);
		assert(testDones(st, false, false, false, false, false, false));
		st.finish();
		assert(testDones(st, true, true, true, true, true, true));
		assert_eq(0, st.numConcordant());
		assert_eq(1, st.numDiscordant());
		assert_eq(0, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		assert(st.repOk());
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(1, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max);
		assert(!unpair2Max);
	}
	cerr << "PASSED" << endl;
}

#endif /*def ALN_SINK_MAIN*/
