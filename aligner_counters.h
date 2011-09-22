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
