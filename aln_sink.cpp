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

#include <iomanip>
#include <limits>
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
		// If we're finished looking for concordant alignments, do we have to
		// continue on to search for unpaired alignments?  Only if our exit
		// from the concordant stage is EXIT_SHORT_CIRCUIT_M.  If it's
		// EXIT_SHORT_CIRCUIT_k or EXIT_WITH_ALIGNMENTS, we can skip unpaired.
		assert_neq(ReportingState::EXIT_NO_ALIGNMENTS, exitConcord_);
		if(exitConcord_ != ReportingState::EXIT_SHORT_CIRCUIT_M) {
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
 * Caller uses this member function to indicate that one additional unpaired
 * mate alignment has been found for the specified mate.
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
			if(doneUnpair1_) {
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		}
		if(nunpair1_ > 1) {
			doneDiscord_ = true;
			exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
		}
	} else {
		nunpair2_++;
		// Did we just finish with this mate?
		if(!doneUnpair2_) {
			areDone(nunpair2_, doneUnpair2_, exitUnpair2_);
			if(doneUnpair2_) {
				doneUnpair_ = doneUnpair1_ && doneUnpair2_;
				updateDone();
			}
		}
		if(nunpair2_ > 1) {
			doneDiscord_ = true;
			exitDiscord_ = ReportingState::EXIT_NO_ALIGNMENTS;
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
			if(p_.mixed) {
				unpair1Max = nunpair1_ > (uint64_t)p_.mhits;
				unpair2Max = nunpair2_ > (uint64_t)p_.mhits;
			}
			// Not sure if this is OK
			nconcordAln = 1; // 1 at random
			return;
		} else if(exitConcord_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
			assert_gt(nconcord_, 0);
			// <= k at random
			nconcordAln = min<uint64_t>(nconcord_, p_.khits);
			return;
		}
		assert(!p_.mhitsSet() || nconcord_ <= (uint64_t)p_.mhits+1);
		
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
		// Unpaired alignments either not reportable or non-existant
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
	assert(!p_.mhitsSet() || paired_ || nunpair1_ <= (uint64_t)p_.mhits+1);

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
	assert(!p_.mhitsSet() || paired_ || nunpair2_ <= (uint64_t)p_.mhits+1);
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
		assert(p_.msample);
		exit = ReportingState::EXIT_SHORT_CIRCUIT_M;
	}
}

static std::ostream& printPct(
	std::ostream& os,
	uint64_t num,
	uint64_t denom)
{
	double pct = 0.0f;
	if(denom != 0) { pct = 100.0 * (double)num / (double)denom; }
	os << fixed << setprecision(2) << pct << '%';
	return os;
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
	const ReportingMetrics& met,
	size_t repThresh,   // threshold for uniqueness, or max if no thresh
	bool discord,       // looked for discordant alignments
	bool mixed,         // looked for unpaired alignments where paired failed?
	bool hadoopOut)     // output Hadoop counters?
{
	// NOTE: there's a filtering step at the very beginning, so everything
	// being reported here is post filtering

	bool canRep = repThresh != MAX_SIZE_T;
	if(hadoopOut) {
		cerr << "reporter:counter:Bowtie,Reads processed," << met.nread << endl;
	}
	uint64_t totread = met.nread;
	if(totread > 0) {
		cerr << "" << met.nread << " reads; of these:" << endl;
	} else {
		assert_eq(0, met.npaired);
		assert_eq(0, met.nunpaired);
		cerr << "" << totread << " reads" << endl;
	}
	uint64_t totpair = met.npaired;
	if(totpair > 0) {
		// Paired output
		cerr << "  " << totpair << " (";
		printPct(cerr, totpair, totread);
		cerr << ") were paired; of these:" << endl;

		// Concordants
		cerr << "    " << met.nconcord_0 << " (";
		printPct(cerr, met.nconcord_0, met.npaired);
		cerr << ") aligned concordantly 0 times" << endl;
		if(canRep) {
			// Print the number that aligned concordantly exactly once
			assert_eq(met.nconcord_uni, met.nconcord_uni1+met.nconcord_uni2);
			cerr << "    " << met.nconcord_uni1 << " (";
			printPct(cerr, met.nconcord_uni1, met.npaired);
			cerr << ") aligned concordantly exactly 1 time" << endl;
			
			// Print the number that aligned concordantly more than once but
			// fewer times than the limit
			
			cerr << "    " << met.nconcord_uni2+met.nconcord_rep << " (";
			printPct(cerr, met.nconcord_uni2+met.nconcord_rep, met.npaired);
			cerr << ") aligned concordantly >1 times" << endl;
		} else {
			// Print the number that aligned concordantly exactly once
			assert_eq(met.nconcord_uni, met.nconcord_uni1+met.nconcord_uni2);
			cerr << "    " << met.nconcord_uni1 << " (";
			printPct(cerr, met.nconcord_uni1, met.npaired);
			cerr << ") aligned concordantly exactly 1 time" << endl;

			// Print the number that aligned concordantly more than once
			cerr << "    " << met.nconcord_uni2 << " (";
			printPct(cerr, met.nconcord_uni2, met.npaired);
			cerr << ") aligned concordantly >1 times" << endl;
		}
		if(discord) {
			// TODO: what about discoardant and on separate chromosomes?
		
			// Bring out the unaligned pair total so we can subtract discordants
			cerr << "    ----" << endl;
			cerr << "    " << met.nconcord_0
			     << " pairs aligned concordantly 0 times; of these:" << endl;
			// Discordants
			cerr << "      " << met.ndiscord << " (";
			printPct(cerr, met.ndiscord, met.nconcord_0);
			cerr << ") aligned discordantly 1 time" << endl;
		}
		uint64_t ncondiscord_0 = met.nconcord_0 - met.ndiscord;
		if(mixed) {
			// Bring out the unaligned pair total so we can subtract discordants
			cerr << "    ----" << endl;
			cerr << "    " << ncondiscord_0
			     << " pairs aligned 0 times concordantly or discordantly; of these:" << endl;
			cerr << "      " << (ncondiscord_0 * 2) << " mates make up the pairs; of these:" << endl;
			cerr << "        " << met.nunp_0_0 << " " << "(";
			printPct(cerr, met.nunp_0_0, ncondiscord_0 * 2);
			cerr << ") aligned 0 times" << endl;
			if(canRep) {
				// Print the number that aligned exactly once
				assert_eq(met.nunp_0_uni, met.nunp_0_uni1+met.nunp_0_uni2);
				cerr << "        " << met.nunp_0_uni1 << " (";
				printPct(cerr, met.nunp_0_uni1, ncondiscord_0 * 2);
				cerr << ") aligned exactly 1 time" << endl;

				// Print the number that aligned more than once but fewer times
				// than the limit
				cerr << "        " << met.nunp_0_uni2+met.nunp_0_rep << " (";
				printPct(cerr, met.nunp_0_uni2+met.nunp_0_rep, ncondiscord_0 * 2);
				cerr << ") aligned >1 times" << endl;
			} else {
				// Print the number that aligned exactly once
				assert_eq(met.nunp_0_uni, met.nunp_0_uni1+met.nunp_0_uni2);
				cerr << "        " << met.nunp_0_uni1 << " (";
				printPct(cerr, met.nunp_0_uni1, ncondiscord_0 * 2);
				cerr << ") aligned exactly 1 time" << endl;

				// Print the number that aligned more than once but fewer times
				// than the limit
				cerr << "        " << met.nunp_0_uni2 << " (";
				printPct(cerr, met.nunp_0_uni2, ncondiscord_0 * 2);
				cerr << ") aligned >1 times" << endl;
			}
			
			//if(canRep) {
			//	// Bring out the repetitively aligned pair total so we can subtract discordants
			//	cerr << "    ----" << endl;
			//	cerr << "    " << met.nconcord_rep
			//		 << " pairs aligned concordantly >" << repThresh
			//		 << " times; of these:" << endl;
			//	cerr << "      " << (met.nconcord_rep * 2) << " mates make up the pairs; of these:" << endl;
			//	
			//	cerr << "        " << met.nunp_rep_0 << " (";
			//	printPct(cerr, met.nunp_rep_0, met.nconcord_rep * 2);
			//	cerr << ") aligned 0 times" << endl;
			//	
			//	cerr << "        " << met.nunp_rep_uni << " (";
			//	printPct(cerr, met.nunp_rep_uni, met.nconcord_rep * 2);
			//	cerr << ") aligned >0 and <=" << repThresh << " times" << endl;
			//	
			//	cerr << "        " << met.nunp_rep_rep << " (";
			//	printPct(cerr, met.nunp_rep_rep, met.nconcord_rep * 2);
			//	cerr << ") aligned >" << repThresh << " times" << endl;
			//}
		}
	}
	uint64_t totunpair = met.nunpaired;
	if(totunpair > 0) {
		// Unpaired output
		cerr << "  " << totunpair << " (";
		printPct(cerr, totunpair, totread);
		cerr << ") were unpaired; of these:" << endl;
		
		cerr << "    " << met.nunp_0 << " (";
		printPct(cerr, met.nunp_0, met.nunpaired);
		cerr << ") aligned 0 times" << endl;
		if(hadoopOut) {
			cerr << "reporter:counter:Bowtie 2,Unpaired reads with 0 alignments,"
			     << met.nunpaired << endl;
		}
		
		if(canRep) {
			// Print the number that aligned exactly once
			assert_eq(met.nunp_uni, met.nunp_uni1+met.nunp_uni2);
			cerr << "    " << met.nunp_uni1 << " (";
			printPct(cerr, met.nunp_uni1, met.nunpaired);
			cerr << ") aligned exactly 1 time" << endl;

			// Print the number that aligned more than once but fewer times
			// than the limit
			cerr << "    " << met.nunp_uni2+met.nunp_rep << " (";
			printPct(cerr, met.nunp_uni2+met.nunp_rep, met.nunpaired);
			cerr << ") aligned >1 times" << endl;
		} else {
			// Print the number that aligned exactly once
			assert_eq(met.nunp_uni, met.nunp_uni1+met.nunp_uni2);
			cerr << "    " << met.nunp_uni1 << " (";
			printPct(cerr, met.nunp_uni1, met.nunpaired);
			cerr << ") aligned exactly 1 time" << endl;

			// Print the number that aligned more than once
			cerr << "    " << met.nunp_uni2 << " (";
			printPct(cerr, met.nunp_uni2, met.nunpaired);
			cerr << ") aligned >1 times" << endl;
		}
	}
	uint64_t tot_al_cand = totunpair + totpair*2;
	uint64_t tot_al =
		(met.nconcord_uni + met.nconcord_rep)*2 +
		(met.ndiscord)*2 +
		met.nunp_0_uni +
		met.nunp_0_rep + 
		met.nunp_uni +
		met.nunp_rep;
	assert_leq(tot_al, tot_al_cand);
	printPct(cerr, tot_al, tot_al_cand);
	cerr << " overall alignment rate" << endl;
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
		//rd1buf_ = *rd1;
		//rd1_ = &rd1buf_;
		rd1_ = rd1;
	} else rd1_ = NULL;
	if(rd2 != NULL) {
		//rd2buf_ = *rd2;
		//rd2_ = &rd2buf_;
		rd2_ = rd2;
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
	bestPair_ = best2Pair_ =
	bestUnp1_ = best2Unp1_ =
	bestUnp2_ = best2Unp2_ = std::numeric_limits<THitInt>::min();
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
 *
 * What gets reported for a paired-end alignment?
 *
 * 1. If there are reportable concordant alignments, report those and stop
 * 2. If there are reportable discordant alignments, report those and stop
 * 3. If unpaired alignments can be reported:
 *    3a. Report 
 #
 * Update metrics.  Only ambiguity is: what if a pair aligns repetitively and
 * one of its mates aligns uniquely?
 *
 * 	uint64_t al;   // # mates w/ >= 1 reported alignment
 *  uint64_t unal; // # mates w/ 0 alignments
 *  uint64_t max;  // # mates withheld for exceeding -M/-m ceiling
 *  uint64_t al_concord;  // # pairs w/ >= 1 concordant alignment
 *  uint64_t al_discord;  // # pairs w/ >= 1 discordant alignment
 *  uint64_t max_concord; // # pairs maxed out
 *  uint64_t unal_pair;   // # pairs where neither mate aligned
 */
void AlnSinkWrap::finishRead(
	const SeedResults *sr1,         // seed alignment results for mate 1
	const SeedResults *sr2,         // seed alignment results for mate 2
	bool               exhaust1,    // mate 1 exhausted?
	bool               exhaust2,    // mate 2 exhausted?
	bool               nfilt1,      // mate 1 N-filtered?
	bool               nfilt2,      // mate 2 N-filtered?
	bool               scfilt1,     // mate 1 score-filtered?
	bool               scfilt2,     // mate 2 score-filtered?
	bool               lenfilt1,    // mate 1 length-filtered?
	bool               lenfilt2,    // mate 2 length-filtered?
	bool               qcfilt1,     // mate 1 qc-filtered?
	bool               qcfilt2,     // mate 2 qc-filtered?
	bool               sortByScore, // prioritize alignments by score
	RandomSource&      rnd,         // pseudo-random generator
	ReportingMetrics&  met,         // reporting metrics
	const PerReadMetrics& prm,      // per-read metrics
	bool suppressSeedSummary,       // = true
	bool suppressAlignments)        // = false
{
	obuf_.clear();
	OutputQueueMark qqm(g_.outq(), obuf_, rdid_, threadid_);
	assert(init_);
	if(!suppressSeedSummary) {
		if(sr1 != NULL) {
			assert(rd1_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(obuf_, *rd1_, rdid_, threadid_, *sr1, true);
		} else if(rd1_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(obuf_, *rd1_, rdid_, true);
		}
		if(sr2 != NULL) {
			assert(rd2_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(obuf_, *rd2_, rdid_, threadid_, *sr2, true);
		} else if(rd2_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(obuf_, *rd2_, rdid_, true);
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
		met.nread++;
		if(readIsPair()) {
			met.npaired++;
		} else {
			met.nunpaired++;
		}
		// Report concordant paired-end alignments if possible
		if(nconcord > 0) {
			AlnSetSumm concordSumm(
				rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
				exhaust1, exhaust2, -1, -1);
			// Possibly select a random subset
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs1_, &rs2_, nconcord, select1_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs1_, nconcord, select1_, rnd);
			}
			assert_lt(off, rs1_.size());
			const AlnRes *rs1 = &rs1_[off];
			const AlnRes *rs2 = &rs2_[off];
			AlnFlags flags1(
				ALN_FLAG_PAIR_CONCORD_MATE1,
				st_.params().mhitsSet(),
				unpair1Max,
				pairMax,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,       // primary
				true,       // opp aligned
				rs2->fw()); // opp fw
			AlnFlags flags2(
				ALN_FLAG_PAIR_CONCORD_MATE2,
				st_.params().mhitsSet(),
				unpair2Max,
				pairMax,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				false,      // primary
				true,       // opp aligned
				rs1->fw()); // opp fw
			// Issue: we only set the flags once, but some of the flags might
			// vary from pair to pair among the pairs we're reporting.  For
			// instance, whether the a given mate aligns to the forward strand.
			SeedAlSumm ssm1, ssm2;
			sr1->toSeedAlSumm(ssm1);
			sr2->toSeedAlSumm(ssm2);
			for(size_t i = 0; i < rs1_.size(); i++) {
				rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
				rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
				assert_eq(abs(rs1_[i].fragmentLength()), abs(rs2_[i].fragmentLength()));
			}
			assert(!select1_.empty());
			g_.reportHits(
				obuf_,
				threadid_,
				rd1_,
				rd2_,
				rdid_,
				select1_,
				NULL,
				&rs1_,
				&rs2_,
				pairMax,
				concordSumm,
				ssm1,
				ssm2,
				&flags1,
				&flags2,
				prm,
				mapq_);
			if(pairMax) {
				met.nconcord_rep++;
			} else {
				met.nconcord_uni++;
				assert(!rs1_.empty());
				if(rs1_.size() == 1) {
					met.nconcord_uni1++;
				} else {
					met.nconcord_uni2++;
				}
			}
			init_ = false;
			//g_.outq().finishRead(obuf_, rdid_, threadid_);
			return;
		}
		// Report concordant paired-end alignments if possible
		else if(ndiscord > 0) {
			ASSERT_ONLY(bool ret =) prepareDiscordants();
			assert(ret);
			assert_eq(1, rs1_.size());
			assert_eq(1, rs2_.size());
			AlnSetSumm discordSumm(
				rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
				exhaust1, exhaust2, -1, -1);
			const AlnRes *rs1 = &rs1_[0];
			const AlnRes *rs2 = &rs2_[0];
			AlnFlags flags1(
				ALN_FLAG_PAIR_DISCORD_MATE1,
				st_.params().mhitsSet(),
				false,
				pairMax,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,       // primary
				true,       // opp aligned
				rs2->fw()); // opp fw
			AlnFlags flags2(
				ALN_FLAG_PAIR_DISCORD_MATE2,
				st_.params().mhitsSet(),
				false,
				pairMax,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				false,      // primary
				true,       // opp aligned
				rs1->fw()); // opp fw
			SeedAlSumm ssm1, ssm2;
			sr1->toSeedAlSumm(ssm1);
			sr2->toSeedAlSumm(ssm2);
			for(size_t i = 0; i < rs1_.size(); i++) {
				rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
				rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
				assert(rs1_[i].isFraglenSet() == rs2_[i].isFraglenSet());
				assert(!rs1_[i].isFraglenSet() || abs(rs1_[i].fragmentLength()) == abs(rs2_[i].fragmentLength()));
			}
			ASSERT_ONLY(size_t off);
			if(sortByScore) {
				// Sort by score then pick from low to high
				ASSERT_ONLY(off =) selectByScore(&rs1_, &rs2_, ndiscord, select1_, rnd);
			} else {
				// Select subset randomly
				ASSERT_ONLY(off =) selectAlnsToReport(rs1_, ndiscord, select1_, rnd);
			}
			assert_eq(0, off);
			assert(!select1_.empty());
			g_.reportHits(
				obuf_,
				threadid_,
				rd1_,
				rd2_,
				rdid_,
				select1_,
				NULL,
				&rs1_,
				&rs2_,
				pairMax,
				discordSumm,
				ssm1,
				ssm2,
				&flags1,
				&flags2,
				prm,
				mapq_);
			met.nconcord_0++;
			met.ndiscord++;
			init_ = false;
			//g_.outq().finishRead(obuf_, rdid_, threadid_);
			return;
		}
		// If we're at this point, at least one mate failed to align.
		// BTL: That's not true.  It could be that there are no concordant
		// alignments but both mates have unpaired alignments, with one of
		// the mates having more than one.
		//assert(nunpair1 == 0 || nunpair2 == 0);
		assert(!pairMax);

		// Update counters given that one mate didn't align
		if(readIsPair()) {
			met.nconcord_0++;
		}
		if(rd1_ != NULL) {
			if(nunpair1 > 0) {
				// Update counters
				if(readIsPair()) {
					if(unpair1Max) met.nunp_0_rep++;
					else {
						met.nunp_0_uni++;
						assert(!rs1u_.empty());
						if(rs1u_.size() == 1) {
							met.nunp_0_uni1++;
						} else {
							met.nunp_0_uni2++;
						}
					}
				} else {
					if(unpair1Max) met.nunp_rep++;
					else {
						met.nunp_uni++;
						assert(!rs1u_.empty());
						if(rs1u_.size() == 1) {
							met.nunp_uni1++;
						} else {
							met.nunp_uni2++;
						}
					}
				}
			} else if(unpair1Max) {
				// Update counters
				if(readIsPair())   met.nunp_0_rep++;
				else               met.nunp_rep++;
			} else {
				// Update counters
				if(readIsPair())   met.nunp_0_0++;
				else               met.nunp_0++;
			}
		}
		if(rd2_ != NULL) {
			if(nunpair2 > 0) {
				// Update counters
				if(readIsPair()) {
					if(unpair2Max) met.nunp_0_rep++;
					else {
						assert(!rs2u_.empty());
						met.nunp_0_uni++;
						if(rs2u_.size() == 1) {
							met.nunp_0_uni1++;
						} else {
							met.nunp_0_uni2++;
						}
					}
				} else {
					if(unpair2Max) met.nunp_rep++;
					else {
						assert(!rs2u_.empty());
						met.nunp_uni++;
						if(rs2u_.size() == 1) {
							met.nunp_uni1++;
						} else {
							met.nunp_uni2++;
						}
					}
				}
			} else if(unpair2Max) {
				// Update counters
				if(readIsPair())   met.nunp_0_rep++;
				else               met.nunp_rep++;
			} else {
				// Update counters
				if(readIsPair())   met.nunp_0_0++;
				else               met.nunp_0++;
			}
		}
		
		const AlnRes *repRs1 = NULL, *repRs2 = NULL;
		AlnSetSumm summ1, summ2;
		AlnFlags flags1, flags2;
		TRefId refid = -1; TRefOff refoff = -1;
		bool rep1 = rd1_ != NULL && nunpair1 > 0;
		bool rep2 = rd2_ != NULL && nunpair2 > 0;

		// This is the preliminary if statement for mate 1 - here we're
		// gathering some preliminary information, making it possible to call
		// g_.reportHits(...) with information about both mates potentially
		if(rep1) {
			// Mate 1 aligned at least once
			summ1.init(
				rd1_, NULL, NULL, NULL, &rs1u_, NULL,
				exhaust1, exhaust2, -1, -1);
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs1u_, NULL, nunpair1, select1_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs1u_, nunpair1, select1_, rnd);
			}
			repRs1 = &rs1u_[off];
		} else if(rd1_ != NULL) {
			// Mate 1 failed to align - don't do anything yet.  First we want
			// to collect information on mate 2 in case that factors into the
			// summary
			assert(!unpair1Max);
		}
		
		if(rep2) {
			summ2.init(
				NULL, rd2_, NULL, NULL, NULL, &rs2u_,
				exhaust1, exhaust2, -1, -1);
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs2u_, NULL, nunpair2, select2_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs2u_, nunpair2, select2_, rnd);
			}
			repRs2 = &rs2u_[off];
		} else if(rd2_ != NULL) {
			// Mate 2 failed to align - don't do anything yet.  First we want
			// to collect information on mate 1 in case that factors into the
			// summary
			assert(!unpair2Max);
		}

		// Now set up flags
		if(rep1) {
			// Initialize flags.  Note: We want to have information about how
			// the other mate aligned (if it did) at this point
			flags1.init(
				readIsPair() ?
					ALN_FLAG_PAIR_UNPAIRED_MATE1 :
					ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				unpair1Max,
				pairMax,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,   // primary
				repRs2 != NULL,                    // opp aligned
				repRs2 == NULL || repRs2->fw());   // opp fw
			for(size_t i = 0; i < rs1u_.size(); i++) {
				rs1u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE1, NULL, flags1);
			}
		}
		if(rep2) {
			// Initialize flags.  Note: We want to have information about how
			// the other mate aligned (if it did) at this point
			flags2.init(
				readIsPair() ?
					ALN_FLAG_PAIR_UNPAIRED_MATE2 :
					ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				unpair2Max,
				pairMax,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				true,   // primary
				repRs1 != NULL,                  // opp aligned
				repRs1 == NULL || repRs1->fw()); // opp fw
			for(size_t i = 0; i < rs2u_.size(); i++) {
				rs2u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE2, NULL, flags2);
			}
		}
		
		// Now report mate 1
		if(rep1) {
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			assert(!select1_.empty());
			g_.reportHits(
				obuf_,
				threadid_,
				rd1_,
				repRs2 != NULL ? rd2_ : NULL,
				rdid_,
				select1_,
				repRs2 != NULL ? &select2_ : NULL,
				&rs1u_,
				repRs2 != NULL ? &rs2u_ : NULL,
				unpair1Max,
				summ1,
				ssm1,
				ssm2,
				&flags1,
				repRs2 != NULL ? &flags2 : NULL,
				prm,
				mapq_);
			assert_lt(select1_[0], rs1u_.size());
			refid = rs1u_[select1_[0]].refid();
			refoff = rs1u_[select1_[0]].refoff();
		}
		
		// Now report mate 2
		if(rep2 && !rep1) {
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			assert(!select2_.empty());
			g_.reportHits(
				obuf_,
				threadid_,
				rd2_,
				repRs1 != NULL ? rd1_ : NULL,
				rdid_,
				select2_,
				repRs1 != NULL ? &select1_ : NULL,
				&rs2u_,
				repRs1 != NULL ? &rs1u_ : NULL,
				unpair2Max,
				summ2,
				ssm1,
				ssm2,
				&flags2,
				repRs1 != NULL ? &flags1 : NULL,
				prm,
				mapq_);
			assert_lt(select2_[0], rs2u_.size());
			refid = rs2u_[select2_[0]].refid();
			refoff = rs2u_[select2_[0]].refoff();
		}
		
		if(rd1_ != NULL && nunpair1 == 0) {
			if(nunpair2 > 0) {
				assert_neq(-1, refid);
				summ1.init(
					rd1_, NULL, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, refid, refoff);
			} else {
				summ1.init(
					rd1_, NULL, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, -1, -1);
			}
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			flags1.init(
				readIsPair() ?
					ALN_FLAG_PAIR_UNPAIRED_MATE1 :
					ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				false,
				false,
				nfilt1,
				scfilt1,
				lenfilt1,
				qcfilt1,
				st_.params().mixed,
				true,           // primary
				repRs2 != NULL, // opp aligned
				(repRs2 != NULL) ? repRs2->fw() : false); // opp fw
			g_.reportUnaligned(
				obuf_,      // string to write output to
				threadid_,
				rd1_,    // read 1
				NULL,    // read 2
				rdid_,   // read id
				summ1,   // summ
				ssm1,    // 
				ssm2,
				&flags1, // flags 1
				NULL,    // flags 2
				prm,
				mapq_,   // MAPQ calculator
				true);   // get lock?
		}
		if(rd2_ != NULL && nunpair2 == 0) {
			if(nunpair1 > 0) {
				assert_neq(-1, refid);
				summ2.init(
					NULL, rd2_, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, refid, refoff);
			} else {
				summ2.init(
					NULL, rd2_, NULL, NULL, NULL, NULL,
					exhaust1, exhaust2, -1, -1);
			}
			SeedAlSumm ssm1, ssm2;
			if(sr1 != NULL) sr1->toSeedAlSumm(ssm1);
			if(sr2 != NULL) sr2->toSeedAlSumm(ssm2);
			flags2.init(
				readIsPair() ?
					ALN_FLAG_PAIR_UNPAIRED_MATE2 :
					ALN_FLAG_PAIR_UNPAIRED,
				st_.params().mhitsSet(),
				false,
				false,
				nfilt2,
				scfilt2,
				lenfilt2,
				qcfilt2,
				st_.params().mixed,
				true,           // primary
				repRs1 != NULL, // opp aligned
				(repRs1 != NULL) ? repRs1->fw() : false); // opp fw
			g_.reportUnaligned(
				obuf_,      // string to write output to
				threadid_,
				rd2_,    // read 1
				NULL,    // read 2
				rdid_,   // read id
				summ2,   // summ
				ssm1,
				ssm2,
				&flags2, // flags 1
				NULL,    // flags 2
				prm,
				mapq_,   // MAPQ calculator
				true);   // get lock?
		}
	} // if(suppress alignments)
	init_ = false;
	//g_.outq().finishRead(obuf_, rdid_, threadid_);
	return;
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
	//assert(!st_.done());
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
	TAlScore score = rsa->score().score();
	if(rsb != NULL) score += rsb->score().score();
	// Update best score so far
	if(paired) {
		if(score > bestPair_) {
			best2Pair_ = bestPair_;
			bestPair_ = score;
		} else if(score > best2Pair_) {
			best2Pair_ = score;
		}
	} else {
		if(one) {
			if(score > bestUnp1_) {
				best2Unp1_ = bestUnp1_;
				bestUnp1_ = score;
			} else if(score > best2Unp1_) {
				best2Unp1_ = score;
			}
		} else {
			if(score > bestUnp2_) {
				best2Unp2_ = bestUnp2_;
				bestUnp2_ = score;
			} else if(score > best2Unp2_) {
				best2Unp2_ = score;
			}
		}
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
		//rs1u_.clear();
		//rs2u_.clear();
		return true;
	}
	return false;
}

/**
 * rs1 (possibly together with rs2 if reads are paired) are populated with
 * alignments.  Here we prioritize them according to alignment score, and
 * some randomness to break ties.  Priorities are returned in the 'select'
 * list.
 */
size_t AlnSinkWrap::selectByScore(
	const EList<AlnRes>* rs1,    // alignments to select from (mate 1)
	const EList<AlnRes>* rs2,    // alignments to select from (mate 2, or NULL)
	uint64_t             num,    // number of alignments to select
	EList<size_t>&       select, // prioritized list to put results in
	RandomSource&        rnd)
	const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	assert(rs1 != NULL);
	size_t sz = rs1->size(); // sz = # alignments found
	assert_leq(num, sz);
	if(sz < num) {
		num = sz;
	}
	// num = # to select
	if(sz < 1) {
		return 0;
	}
	select.resize((size_t)num);
	// Use 'selectBuf_' as a temporary list for sorting purposes
	EList<std::pair<TAlScore, size_t> >& buf =
		const_cast<EList<std::pair<TAlScore, size_t> >& >(selectBuf_);
	buf.resize(sz);
	// Sort by score.  If reads are pairs, sort by sum of mate scores.
	for(size_t i = 0; i < sz; i++) {
		buf[i].first = (*rs1)[i].score().score();
		if(rs2 != NULL) {
			buf[i].first += (*rs2)[i].score().score();
		}
		buf[i].second = i; // original offset
	}
	buf.sort(); buf.reverse();
	for(size_t i = 0; i < num; i++) { select[i] = selectBuf_[i].second; }
	// Returns index of the representative alignment, but in 'select' also
	// returns the indexes of the next best selected alignments in order by
	// score.
	return selectBuf_[0].second;
}

/**
 * Given that rs is already populated with alignments, consider the
 * alignment policy and make random selections where necessary.  E.g. if we
 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
 * random.  We "select" an alignment by setting the parallel entry in the
 * 'select' list to true.
 *
 * Return the "representative" alignment.  This is simply the first one
 * selected.  That will also be what SAM calls the "primary" alignment.
 */
size_t AlnSinkWrap::selectAlnsToReport(
	const EList<AlnRes>& rs,     // alignments to select from
	uint64_t             num,    // number of alignments to select
	EList<size_t>&       select, // list to put results in
	RandomSource&        rnd)
	const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	size_t sz = rs.size();
	if(sz < num) {
		num = sz;
	}
	if(sz < 1) {
		return 0;
	}
	select.resize((size_t)num);
	if(sz == 1) {
		assert_eq(1, num);
		select[0] = 0;
		return 0;
	}
	// Select a random offset into the list of alignments
	uint32_t off = rnd.nextU32() % (uint32_t)sz;
	uint32_t offOrig = off;
	// Now take elements starting at that offset, wrapping around to 0 if
	// necessary.  Leave the rest.
	for(size_t i = 0; i < num; i++) {
		select[i] = off;
		off++;
		if(off == sz) {
			off = 0;
		}
	}
	return offOrig;
}

#define NOT_SUPPRESSED !suppress_[field++]
#define BEGIN_FIELD { \
	if(firstfield) firstfield = false; \
	else o.append('\t'); \
}
#define WRITE_TAB { \
	if(firstfield) firstfield = false; \
	else o.append('\t'); \
}
#define WRITE_NUM(o, x) { \
	itoa10(x, buf); \
	o.append(buf); \
}

/**
 * Print a seed summary to the first output stream in the outs_ list.
 */
void AlnSink::reportSeedSummary(
	BTString&          o,
	const Read&        rd,
	TReadId            rdid,
	size_t             threadId,
	const SeedResults& rs,
	bool               getLock)
{
	appendSeedSummary(
		o,                     // string to write to
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
	BTString&          o,
	const Read&        rd,
	TReadId            rdid,
	size_t             threadId,
	bool               getLock)
{
	appendSeedSummary(
		o,                     // string to append to
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
	BTString& s,
	const T& str,
	bool chopws)
{
	size_t len = str.length();
	for(size_t i = 0; i < len; i++) {
		if(!chopws || (str[i] != ' ' && str[i] != '\t')) {
			s.append(str[i]);
		} else {
			break;
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
	BTString&     o,
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
	bool firstfield = true;
	//
	// Read name
	//
	BEGIN_FIELD;
	printUptoWs(o, rd.name, true);
	
	//
	// Total number of seeds tried
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTried);

	//
	// Total number of seeds tried where at least one range was found.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, nonzero);

	//
	// Total number of ranges found
	//
	BEGIN_FIELD;
	WRITE_NUM(o, ranges);

	//
	// Total number of elements found
	//
	BEGIN_FIELD;
	WRITE_NUM(o, elts);
	
	//
	// The same four numbers, but only for seeds extracted from the
	// forward read representation.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedFw);

	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroFw);

	BEGIN_FIELD;
	WRITE_NUM(o, rangesFw);

	BEGIN_FIELD;
	WRITE_NUM(o, eltsFw);

	//
	// The same four numbers, but only for seeds extracted from the
	// reverse complement read representation.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedRc);

	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroRc);

	BEGIN_FIELD;
	WRITE_NUM(o, rangesRc);

	BEGIN_FIELD;
	WRITE_NUM(o, eltsRc);

	o.append('\n');
}

/**
 * Append a single hit to the given output stream in Bowtie's
 * verbose-mode format.
 */
void AlnSinkSam::appendMate(
	BTString&     o,           // append to this string
	const Read&   rd,
	const Read*   rdo,
	const TReadId rdid,
	AlnRes* rs,
	AlnRes* rso,
	const AlnSetSumm& summ,
	const SeedAlSumm& ssm,
	const SeedAlSumm& ssmo,
	const AlnFlags& flags,
	const PerReadMetrics& prm,
	const Mapq& mapqCalc)
{
	if(rs == NULL && samc_.omitUnalignedReads()) {
		return;
	}
	char buf[1024];
	StackedAln staln;
	if(rs != NULL) {
		rs->initStacked(rd, staln);
		staln.leftAlign(false /* not past MMs */);
	}
	int offAdj = 0;
	// QNAME
	samc_.printReadName(o, rd.name, flags.partOfPair());
	o.append('\t');
	// FLAG
	int fl = 0;
	if(flags.partOfPair()) {
		fl |= SAM_FLAG_PAIRED;
		if(flags.alignedConcordant()) {
			fl |= SAM_FLAG_MAPPED_PAIRED;
 		}
		if(!flags.mateAligned()) {
			// Other fragment is unmapped
			fl |= SAM_FLAG_MATE_UNMAPPED;
		}
		fl |= (flags.readMate1() ?
			SAM_FLAG_FIRST_IN_PAIR : SAM_FLAG_SECOND_IN_PAIR);
		if(flags.mateAligned() && rso != NULL) {
			if(!rso->fw()) {
				fl |= SAM_FLAG_MATE_STRAND;
			}
		}
	}
	if(!flags.isPrimary()) {
		fl |= SAM_FLAG_NOT_PRIMARY;
	}
	if(rs != NULL && !rs->fw()) {
		fl |= SAM_FLAG_QUERY_STRAND;
	}
	if(rs == NULL) {
		// Failed to align
		fl |= SAM_FLAG_UNMAPPED;
	}
	itoa10<int>(fl, buf);
	o.append(buf);
	o.append('\t');
	// RNAME
	if(rs != NULL) {
		samc_.printRefNameFromIndex(o, (size_t)rs->refid());
		o.append('\t');
	} else {
		if(summ.orefid() != -1) {
			// Opposite mate aligned but this one didn't - print the opposite
			// mate's RNAME and POS as is customary
			assert(flags.partOfPair());
			samc_.printRefNameFromIndex(o, (size_t)summ.orefid());
		} else {		
			// No alignment
			o.append('*');
		}
		o.append('\t');
	}
	// POS
	// Note: POS is *after* soft clipping.  I.e. POS points to the
	// upstream-most character *involved in the clipped alignment*.
	if(rs != NULL) {
		itoa10<int64_t>(rs->refoff()+1+offAdj, buf);
		o.append(buf);
		o.append('\t');
	} else {
		if(summ.orefid() != -1) {
			// Opposite mate aligned but this one didn't - print the opposite
			// mate's RNAME and POS as is customary
			assert(flags.partOfPair());
			itoa10<int64_t>(summ.orefoff()+1+offAdj, buf);
			o.append(buf);
		} else {
			// No alignment
			o.append('0');
		}
		o.append('\t');
	}
	// MAPQ
	mapqInps_[0] = '\0';
	if(rs != NULL) {
		itoa10<TMapq>(mapqCalc.mapq(
			summ, flags, rd.mate < 2, rd.length(),
			rdo == NULL ? 0 : rdo->length(), mapqInps_), buf);
		o.append(buf);
		o.append('\t');
	} else {
		// No alignment
		o.append("0\t");
	}
	// CIGAR
	if(rs != NULL) {
		staln.buildCigar(false);
		staln.writeCigar(&o, NULL);
		o.append('\t');
	} else {
		// No alignment
		o.append("*\t");
	}
	// RNEXT
	if(rs != NULL && flags.partOfPair()) {
		if(rso != NULL && rs->refid() != rso->refid()) {
			samc_.printRefNameFromIndex(o, (size_t)rso->refid());
			o.append('\t');
		} else {
			o.append("=\t");
		}
	} else if(summ.orefid() != -1) {
		// The convention if this mate fails to align but the other doesn't is
		// to copy the mate's details into here
		o.append("=\t");
	} else {
		o.append("*\t");
	}
	// PNEXT
	if(rs != NULL && flags.partOfPair()) {
		if(rso != NULL) {
			itoa10<int64_t>(rso->refoff()+1, buf);
			o.append(buf);
			o.append('\t');
		} else {
			// The convenstion is that if this mate aligns but the opposite
			// doesn't, we print this mate's offset here
			itoa10<int64_t>(rs->refoff()+1, buf);
			o.append(buf);
			o.append('\t');
		}
	} else if(summ.orefid() != -1) {
		// The convention if this mate fails to align but the other doesn't is
		// to copy the mate's details into here
		itoa10<int64_t>(summ.orefoff()+1, buf);
		o.append(buf);
		o.append('\t');
	} else {
		o.append("0\t");
	}
	// ISIZE
	if(rs != NULL && rs->isFraglenSet()) {
		itoa10<int64_t>(rs->fragmentLength(), buf);
		o.append(buf);
		o.append('\t');
	} else {
		// No fragment
		o.append("0\t");
	}
	// SEQ
	if(!flags.isPrimary() && samc_.omitSecondarySeqQual()) {
		o.append('*');
	} else {
		// Print the read
		if(rd.patFw.length() == 0) {
			o.append('*');
		} else {
			if(rs == NULL || rs->fw()) {
				o.append(rd.patFw.toZBuf());
			} else {
				o.append(rd.patRc.toZBuf());
			}
		}
	}
	o.append('\t');
	// QUAL
	if(!flags.isPrimary() && samc_.omitSecondarySeqQual()) {
		o.append('*');
	} else {
		// Print the quals
		if(rd.qual.length() == 0) {
			o.append('*');
		} else {
			if(rs == NULL || rs->fw()) {
				o.append(rd.qual.toZBuf());
			} else {
				o.append(rd.qualRev.toZBuf());
			}
		}
	}
	o.append('\t');
	//
	// Optional fields
	//
	if(rs != NULL) {
		samc_.printAlignedOptFlags(
			o,           // output buffer
			true,        // first opt flag printed is first overall?
			rd,          // read
			*rs,         // individual alignment result
			staln,       // stacked alignment
			flags,       // alignment flags
			summ,        // summary of alignments for this read
			ssm,         // seed alignment summary
			prm,         // per-read metrics
			mapqInps_);  // inputs to MAPQ calculation
	} else {
		samc_.printEmptyOptFlags(
			o,           // output buffer
			true,        // first opt flag printed is first overall?
			rd,          // read
			flags,       // alignment flags
			summ,        // summary of alignments for this read
			ssm,         // seed alignment summary
			prm);        // per-read metrics
	}
	o.append('\n');
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

	cerr << "Case 7 (unaligned pair & uniquely aligned mate, mixed-mode) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,      // khits
			1,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		// assert(st.doneConcordant()    == done1);
		// assert(st.doneDiscordant()    == done2);
		// assert(st.doneUnpaired(true)  == done3);
		// assert(st.doneUnpaired(false) == done4);
		// assert(st.doneUnpaired()      == done5);
		// assert(st.done()              == done6);
		st.foundUnpaired(true);
		assert(testDones(st, false, false, false, false, false, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, false, false, false));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
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

	cerr << "Case 8 (unaligned pair & uniquely aligned mate, NOT mixed-mode) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,      // khits
			1,      // mhits
			0,      // pengap
			false,  // msample
			true,   // discord
			false); // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		// assert(st.doneConcordant()    == done1);
		// assert(st.doneDiscordant()    == done2);
		// assert(st.doneUnpaired(true)  == done3);
		// assert(st.doneUnpaired(false) == done4);
		// assert(st.doneUnpaired()      == done5);
		// assert(st.done()              == done6);
		st.foundUnpaired(true);
		assert(testDones(st, false, false, true, true, true, false));
		st.foundUnpaired(true);
		assert(testDones(st, false, true, true, true, true, false));
		assert_eq(0, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(0, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(0, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(!pairMax);
		assert(!unpair1Max); // not really relevant
		assert(!unpair2Max); // not really relevant
	}
	cerr << "PASSED" << endl;

	cerr << "Case 9 (repetitive pair, only one mate repetitive) ... ";
	{
		uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
		bool pairMax = false, unpair1Max = false, unpair2Max = false;
		ReportingParams rp(
			1,      // khits
			1,      // mhits
			0,      // pengap
			true,   // msample
			true,   // discord
			true);  // mixed
		ReportingState st(rp);
		st.nextRead(true); // unpaired read
		// assert(st.doneConcordant()    == done1);
		// assert(st.doneDiscordant()    == done2);
		// assert(st.doneUnpaired(true)  == done3);
		// assert(st.doneUnpaired(false) == done4);
		// assert(st.doneUnpaired()      == done5);
		// assert(st.done()              == done6);
		st.foundConcordant();
		assert(st.repOk());
		st.foundUnpaired(true);
		assert(st.repOk());
		st.foundUnpaired(false);
		assert(st.repOk());
		assert(testDones(st, false, true, false, false, false, false));
		assert(st.repOk());
		st.foundConcordant();
		assert(st.repOk());
		st.foundUnpaired(true);
		assert(st.repOk());
		assert(testDones(st, true, true, true, false, false, false));
		assert_eq(2, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(1, st.numUnpaired2());
		st.foundUnpaired(false);
		assert(st.repOk());
		assert(testDones(st, true, true, true, true, true, true));		
		assert_eq(2, st.numConcordant());
		assert_eq(0, st.numDiscordant());
		assert_eq(2, st.numUnpaired1());
		assert_eq(2, st.numUnpaired2());
		st.finish();
		st.getReport(nconcord, ndiscord, nunpair1, nunpair2,
		             pairMax, unpair1Max, unpair2Max);
		assert_eq(1, nconcord);
		assert_eq(0, ndiscord);
		assert_eq(0, nunpair1);
		assert_eq(0, nunpair2);
		assert(pairMax);
		assert(unpair1Max); // not really relevant
		assert(unpair2Max); // not really relevant
	}
	cerr << "PASSED" << endl;
}

#endif /*def ALN_SINK_MAIN*/
