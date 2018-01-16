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

#ifndef OUTQ_H_
#define OUTQ_H_

#include "assert_helpers.h"
#include "ds.h"
#include "sstring.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"
#include <vector>

static std::string compose_multisam_name(const std::string& fn, int idx) {
	assert_geq(idx, 0);
	const size_t len = fn.length();
	std::ostringstream os;
	if(fn.substr(len-4, 4) == ".sam") {
		os << fn.substr(0, len-3) << "part" << (idx+1) << ".sam";
		return os.str();
	}
	os << fn << ".part" << (idx+1);
	return os.str();
}

/**
 * Encapsulates a list of lines of output.  If the earliest as-yet-unreported
 * read has id N and Bowtie 2 wants to write a record for read with id N+1, we
 * resize the lines_ and committed_ lists to have at least 2 elements (1 for N,
 * 1 for N+1) and return the BTString * associated with the 2nd element.  When
 * the user calls commit() for the read with id N, 
 */
class OutputQueue {

	static const size_t NFLUSH_THRESH = 8;

public:

	OutputQueue(
		const std::string& ofn, // empty -> stdin
		size_t output_buffer_size,
		bool reorder,
		size_t nthreads,
		bool threadSafe,
		int perThreadBufSize,
		int nmulti_output,
		TReadId rdid = 0) :
		ofhs_(),
		obufs_(),
		cur_(rdid),
		nfinished_(0),
		nflushed_(0),
		lines_(RES_CAT),
		started_(RES_CAT),
		finished_(RES_CAT),
		reorder_(reorder),
		threadSafe_(threadSafe),
		mutex_global_(),
		mutexes_(),
		nmulti_output_(nmulti_output),
		nthreads_(nthreads),
		perThreadBufSize_(perThreadBufSize)
	{
		nstarted_=0;
		assert(nthreads_ <= 2 || threadSafe);
		if(!reorder)
		{
			perThreadBuf.resize(nthreads_, std::vector<BTString>(perThreadBufSize_));
			perThreadCounter = new int[nthreads_];
			size_t i = 0;
			for(i=0;i<nthreads_;i++)
			{
				perThreadCounter[i] = 0;
			}
		}
		assert_gt(nmulti_output_, 0);
		if(ofn.empty() || ofn == "/dev/null") {
			nmulti_output_ = 1;
		}
		std::string ofn_l = ofn;
		for(int i = 0; i < nmulti_output_; i++) {
			if(nmulti_output_ > 1) {
				ofn_l = compose_multisam_name(ofn, i);
			}
			FILE *ofh = fopen(ofn_l.c_str(), "w");
			if(ofh == NULL) {
				std::cerr << "Error: Could not open alignment output file "
				<< ofn_l << std::endl;
				throw 1;
			}
			char *obuf = new char[output_buffer_size];
			int ret = setvbuf(ofh, obuf, _IOFBF, output_buffer_size);
			if(ret != 0) {
				std::cerr << "Warning: Could not allocate the proper "
				<< "buffer size for output file stream. "
				<< "Return value = " << ret << std::endl;
			}
			obufs_.push_back(obuf);
			ofhs_.push_back(ofh);
			mutexes_.push_back(new MUTEX_T());
		}
	}

	~OutputQueue() {
		for(int i = 0; i < (int)ofhs_.size(); i++) {
			if(ofhs_[i] != NULL) {
				delete[] obufs_[i];
				delete mutexes_[i];
				fclose(ofhs_[i]);
				ofhs_[i] = NULL;
			}
		}
	}

	/**
	 * Caller is telling us that they're about to write output record(s) for
	 * the read with the given id.
	 */
	void beginRead(TReadId rdid, size_t threadId);
	
	/**
	 * Writer is finished writing to 
	 */
	void finishRead(const BTString& rec, TReadId rdid, size_t threadId);
	
	/**
	 * Return the number of records currently being buffered.
	 */
	size_t size() const {
		return lines_.size();
	}
	
	/**
	 * Return the number of records that have been flushed so far.
	 */
	TReadId numFlushed() const {
		return nflushed_;
	}

	/**
	 * Return the number of records that have been started so far.
	 */
	TReadId numStarted() const {
		return nstarted_;
	}

	/**
	 * Return the number of records that have been finished so far.
	 */
	TReadId numFinished() const {
		return nfinished_;
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	void writeString(const BTString& s, int outidx);
	
	/**
	 * Write already-committed lines starting from cur_.
	 */
	void flush(bool force = false, bool getLock = true);

protected:

	std::vector<FILE *> ofhs_;
	std::vector<char *> obufs_;
	TReadId         cur_;
#ifdef WITH_TBB
	tbb::atomic<TReadId> nstarted_;
#else
	TReadId         nstarted_;
#endif
	TReadId         nfinished_;
	TReadId         nflushed_;
	EList<BTString> lines_;
	EList<bool>     started_;
	EList<bool>     finished_;
	bool            reorder_;
	bool            threadSafe_;

	MUTEX_T         mutex_global_;  // for reorder bookkeeping
	std::vector<MUTEX_T*> mutexes_;
	int             nmulti_output_;
	
	// used for output read buffer	
	size_t nthreads_;
	std::vector<vector<BTString> > perThreadBuf;
	int* 		perThreadCounter;
	int perThreadBufSize_;

private:

	void flushImpl(bool force);
	void beginReadImpl(TReadId rdid, size_t threadId);
	void finishReadImpl(const BTString& rec, TReadId rdid, size_t threadId);
};

class OutputQueueMark {
public:
	OutputQueueMark(
		OutputQueue& q,
		const BTString& rec,
		TReadId rdid,
		size_t threadId) :
		q_(q),
		rec_(rec),
		rdid_(rdid),
		threadId_(threadId)
	{
		q_.beginRead(rdid, threadId);
	}
	
	~OutputQueueMark() {
		q_.finishRead(rec_, rdid_, threadId_);
	}
	
protected:
	OutputQueue& q_;
	const BTString& rec_;
	TReadId rdid_;
	size_t threadId_;
};

#endif
