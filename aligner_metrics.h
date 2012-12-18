/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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

#ifndef ALIGNER_METRICS_H_
#define ALIGNER_METRICS_H_

#include <math.h>
#include <iostream>
#include "alphabet.h"
#include "timer.h"
#include "sstring.h"

using namespace std;

/**
 * Borrowed from http://www.johndcook.com/standard_deviation.html,
 * which in turn is borrowed from Knuth.
 */
class RunningStat {
public:
	RunningStat() : m_n(0), m_tot(0.0) { }

	void clear() {
		m_n = 0;
		m_tot = 0.0;
	}

	void push(float x) {
		m_n++;
		m_tot += x;
		// See Knuth TAOCP vol 2, 3rd edition, page 232
		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		} else {
			m_newM = m_oldM + (x - m_oldM)/m_n;
			m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
			// set up for next iteration
			m_oldM = m_newM;
			m_oldS = m_newS;
		}
	}

	int num() const {
		return m_n;
	}

	double tot() const {
		return m_tot;
	}

	double mean() const {
		return (m_n > 0) ? m_newM : 0.0;
	}

	double variance() const {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	double stddev() const {
		return sqrt(variance());
	}

private:
	int m_n;
	double m_tot;
	double m_oldM, m_newM, m_oldS, m_newS;
};

/**
 * Encapsulates a set of metrics that we would like an aligner to keep
 * track of, so that we can possibly use it to diagnose performance
 * issues.
 */
class AlignerMetrics {

public:

	AlignerMetrics() :
		curBacktracks_(0),
		curBwtOps_(0),
		first_(true),
		curIsLowEntropy_(false),
		curIsHomoPoly_(false),
		curHadRanges_(false),
		curNumNs_(0),
		reads_(0),
		homoReads_(0),
		lowEntReads_(0),
		hiEntReads_(0),
		alignedReads_(0),
		unalignedReads_(0),
		threeOrMoreNReads_(0),
		lessThanThreeNRreads_(0),
		bwtOpsPerRead_(),
		backtracksPerRead_(),
		bwtOpsPerHomoRead_(),
		backtracksPerHomoRead_(),
		bwtOpsPerLoEntRead_(),
		backtracksPerLoEntRead_(),
		bwtOpsPerHiEntRead_(),
		backtracksPerHiEntRead_(),
		bwtOpsPerAlignedRead_(),
		backtracksPerAlignedRead_(),
		bwtOpsPerUnalignedRead_(),
		backtracksPerUnalignedRead_(),
		bwtOpsPer0nRead_(),
		backtracksPer0nRead_(),
		bwtOpsPer1nRead_(),
		backtracksPer1nRead_(),
		bwtOpsPer2nRead_(),
		backtracksPer2nRead_(),
		bwtOpsPer3orMoreNRead_(),
		backtracksPer3orMoreNRead_(),
		timer_(cout, "", false)
		{ }

	void printSummary() {
		if(!first_) {
			finishRead();
		}
		cout << "AlignerMetrics:" << endl;
		cout << "  # Reads:             " << reads_ << endl;
		float hopct = (reads_ > 0) ? (((float)homoReads_)/((float)reads_)) : (0.0f);
		hopct *= 100.0f;
		cout << "  % homo-polymeric:    " << (hopct) << endl;
		float lopct = (reads_ > 0) ? ((float)lowEntReads_/(float)(reads_)) : (0.0f);
		lopct *= 100.0f;
		cout << "  % low-entropy:       " << (lopct) << endl;
		float unpct = (reads_ > 0) ? ((float)unalignedReads_/(float)(reads_)) : (0.0f);
		unpct *= 100.0f;
		cout << "  % unaligned:         " << (unpct) << endl;
		float npct = (reads_ > 0) ? ((float)threeOrMoreNReads_/(float)(reads_)) : (0.0f);
		npct *= 100.0f;
		cout << "  % with 3 or more Ns: " << (npct) << endl;
		cout << endl;
		cout << "  Total BWT ops:    avg: " << bwtOpsPerRead_.mean() << ", stddev: " << bwtOpsPerRead_.stddev() << endl;
		cout << "  Total Backtracks: avg: " << backtracksPerRead_.mean() << ", stddev: " << backtracksPerRead_.stddev() << endl;
		time_t elapsed = timer_.elapsed();
		cout << "  BWT ops per second:    " << (bwtOpsPerRead_.tot()/elapsed) << endl;
		cout << "  Backtracks per second: " << (backtracksPerRead_.tot()/elapsed) << endl;
		cout << endl;
		cout << "  Homo-poly:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerHomoRead_.mean() << ", stddev: " << bwtOpsPerHomoRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerHomoRead_.mean() << ", stddev: " << backtracksPerHomoRead_.stddev() << endl;
		cout << "  Low-entropy:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerLoEntRead_.mean() << ", stddev: " << bwtOpsPerLoEntRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerLoEntRead_.mean() << ", stddev: " << backtracksPerLoEntRead_.stddev() << endl;
		cout << "  High-entropy:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerHiEntRead_.mean() << ", stddev: " << bwtOpsPerHiEntRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerHiEntRead_.mean() << ", stddev: " << backtracksPerHiEntRead_.stddev() << endl;
		cout << endl;
		cout << "  Unaligned:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerUnalignedRead_.mean() << ", stddev: " << bwtOpsPerUnalignedRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerUnalignedRead_.mean() << ", stddev: " << backtracksPerUnalignedRead_.stddev() << endl;
		cout << "  Aligned:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerAlignedRead_.mean() << ", stddev: " << bwtOpsPerAlignedRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerAlignedRead_.mean() << ", stddev: " << backtracksPerAlignedRead_.stddev() << endl;
		cout << endl;
		cout << "  0 Ns:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPer0nRead_.mean() << ", stddev: " << bwtOpsPer0nRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPer0nRead_.mean() << ", stddev: " << backtracksPer0nRead_.stddev() << endl;
		cout << "  1 N:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPer1nRead_.mean() << ", stddev: " << bwtOpsPer1nRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPer1nRead_.mean() << ", stddev: " << backtracksPer1nRead_.stddev() << endl;
		cout << "  2 Ns:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPer2nRead_.mean() << ", stddev: " << bwtOpsPer2nRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPer2nRead_.mean() << ", stddev: " << backtracksPer2nRead_.stddev() << endl;
		cout << "  >2 Ns:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPer3orMoreNRead_.mean() << ", stddev: " << bwtOpsPer3orMoreNRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPer3orMoreNRead_.mean() << ", stddev: " << backtracksPer3orMoreNRead_.stddev() << endl;
		cout << endl;
	}

	/**
	 *
	 */
	void nextRead(const BTDnaString& read) {
		if(!first_) {
			finishRead();
		}
		first_ = false;
		//float ent = entropyDna5(read);
		float ent = 0.0f;
		curIsLowEntropy_ = (ent < 0.75f);
		curIsHomoPoly_ = (ent < 0.001f);
		curHadRanges_ = false;
		curBwtOps_ = 0;
		curBacktracks_ = 0;
		// Count Ns
		curNumNs_ = 0;
		const size_t len = read.length();
		for(size_t i = 0; i < len; i++) {
			if((int)read[i] == 4) curNumNs_++;
		}
	}

	/**
	 *
	 */
	void setReadHasRange() {
		curHadRanges_ = true;
	}

	/**
	 * Commit the running statistics for this read to
	 */
	void finishRead() {
		reads_++;
		if(curIsHomoPoly_) homoReads_++;
		else if(curIsLowEntropy_) lowEntReads_++;
		else hiEntReads_++;
		if(curHadRanges_) alignedReads_++;
		else unalignedReads_++;
		bwtOpsPerRead_.push((float)curBwtOps_);
		backtracksPerRead_.push((float)curBacktracks_);
		// Drill down by entropy
		if(curIsHomoPoly_) {
			bwtOpsPerHomoRead_.push((float)curBwtOps_);
			backtracksPerHomoRead_.push((float)curBacktracks_);
		} else if(curIsLowEntropy_) {
			bwtOpsPerLoEntRead_.push((float)curBwtOps_);
			backtracksPerLoEntRead_.push((float)curBacktracks_);
		} else {
			bwtOpsPerHiEntRead_.push((float)curBwtOps_);
			backtracksPerHiEntRead_.push((float)curBacktracks_);
		}
		// Drill down by whether it aligned
		if(curHadRanges_) {
			bwtOpsPerAlignedRead_.push((float)curBwtOps_);
			backtracksPerAlignedRead_.push((float)curBacktracks_);
		} else {
			bwtOpsPerUnalignedRead_.push((float)curBwtOps_);
			backtracksPerUnalignedRead_.push((float)curBacktracks_);
		}
		if(curNumNs_ == 0) {
			lessThanThreeNRreads_++;
			bwtOpsPer0nRead_.push((float)curBwtOps_);
			backtracksPer0nRead_.push((float)curBacktracks_);
		} else if(curNumNs_ == 1) {
			lessThanThreeNRreads_++;
			bwtOpsPer1nRead_.push((float)curBwtOps_);
			backtracksPer1nRead_.push((float)curBacktracks_);
		} else if(curNumNs_ == 2) {
			lessThanThreeNRreads_++;
			bwtOpsPer2nRead_.push((float)curBwtOps_);
			backtracksPer2nRead_.push((float)curBacktracks_);
		} else {
			threeOrMoreNReads_++;
			bwtOpsPer3orMoreNRead_.push((float)curBwtOps_);
			backtracksPer3orMoreNRead_.push((float)curBacktracks_);
		}
	}

	// Running-total of the number of backtracks and BWT ops for the
	// current read
	uint32_t curBacktracks_;
	uint32_t curBwtOps_;

protected:

	bool first_;

	// true iff the current read is low entropy
	bool curIsLowEntropy_;
	// true if current read is all 1 char (or very close)
	bool curIsHomoPoly_;
	// true iff the current read has had one or more ranges reported
	bool curHadRanges_;
	// number of Ns in current read
	int curNumNs_;

	// # reads
	uint32_t reads_;
	// # homo-poly reads
	uint32_t homoReads_;
	// # low-entropy reads
	uint32_t lowEntReads_;
	// # high-entropy reads
	uint32_t hiEntReads_;
	// # reads with alignments
	uint32_t alignedReads_;
	// # reads without alignments
	uint32_t unalignedReads_;
	// # reads with 3 or more Ns
	uint32_t threeOrMoreNReads_;
	// # reads with < 3 Ns
	uint32_t lessThanThreeNRreads_;

	// Distribution of BWT operations per read
	RunningStat bwtOpsPerRead_;
	RunningStat backtracksPerRead_;

	// Distribution of BWT operations per homo-poly read
	RunningStat bwtOpsPerHomoRead_;
	RunningStat backtracksPerHomoRead_;

	// Distribution of BWT operations per low-entropy read
	RunningStat bwtOpsPerLoEntRead_;
	RunningStat backtracksPerLoEntRead_;

	// Distribution of BWT operations per high-entropy read
	RunningStat bwtOpsPerHiEntRead_;
	RunningStat backtracksPerHiEntRead_;

	// Distribution of BWT operations per read that "aligned" (for
	// which a range was arrived at - range may not have necessarily
	// lead to an alignment)
	RunningStat bwtOpsPerAlignedRead_;
	RunningStat backtracksPerAlignedRead_;

	// Distribution of BWT operations per read that didn't align
	RunningStat bwtOpsPerUnalignedRead_;
	RunningStat backtracksPerUnalignedRead_;

	// Distribution of BWT operations/backtracks per read with no Ns
	RunningStat bwtOpsPer0nRead_;
	RunningStat backtracksPer0nRead_;

	// Distribution of BWT operations/backtracks per read with one N
	RunningStat bwtOpsPer1nRead_;
	RunningStat backtracksPer1nRead_;

	// Distribution of BWT operations/backtracks per read with two Ns
	RunningStat bwtOpsPer2nRead_;
	RunningStat backtracksPer2nRead_;

	// Distribution of BWT operations/backtracks per read with three or
	// more Ns
	RunningStat bwtOpsPer3orMoreNRead_;
	RunningStat backtracksPer3orMoreNRead_;

	Timer timer_;
};

#endif /* ALIGNER_METRICS_H_ */
