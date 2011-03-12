/*
 *  aligner_counters.h
 *  bowtie-beta1
 *
 *  Created by Benjamin Langmead on 9/18/10.
 *  Copyright 2010 Johns Hopkins University. All rights reserved.
 *
 */

#ifndef ALIGNER_COUNTERS_H_
#define ALIGNER_COUNTERS_H_

#include "aligner_seed.h"
#include "aligner_sw.h"

/**
 * A set of counters for characterizing the work done by the seed
 * aligner.
 */
struct ReadCounters {
	SACounters   seed;
	SwCounters   sw;
};

/**
 * Abstract parent for a class with a method that gets passed every
 * set of counters for every join attempt.
 */
class ReadCounterSink {
public:
	ReadCounterSink() { MUTEX_INIT(lock_); }
	virtual ~ReadCounterSink() { }
	/**
	 * Grab the lock and call abstract member reportCountersImpl()
	 */
	virtual void reportCounters(const ReadCounters& c) {
		ThreadSafe(&this->lock_);
		reportCountersImpl(c);
	}
protected:
	virtual void reportCountersImpl(const ReadCounters& c) = 0;
	MUTEX_T lock_;
};

/**
 * Write each per-read set of counters to an output stream using a
 * simple record-per-line tab-delimited format.
 */
class StreamTabReadCounterSink : public ReadCounterSink {
public:
	StreamTabReadCounterSink(std::ostream& os) : ReadCounterSink(), os_(os) { }
	virtual ~StreamTabReadCounterSink() { }
protected:
	virtual void reportCountersImpl(
		const ReadCounters& c)
	{
		os_ << "\n"; // avoid 'endl' b/c flush is unnecessary
	}
	std::ostream& os_;
};

#endif /*ALIGNER_COUNTERS_H_*/
