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

#include "outq.h"

/**
 * Caller is telling us that they're about to write output record(s) for
 * the read with the given id.
 */
void OutputQueue::beginReadImpl(TReadId rdid, size_t threadId) {
#ifdef WITH_TBB
	nstarted_.fetch_and_add(1);
#else
	nstarted_++;
#endif
	if(reorder_) {
		assert_geq(rdid, cur_);
		assert_eq(lines_.size(), finished_.size());
		assert_eq(lines_.size(), started_.size());
		if(rdid - cur_ >= lines_.size()) {
			// Make sure there's enough room in lines_, started_ and finished_
			size_t oldsz = lines_.size();
			lines_.resize(rdid - cur_ + 1);
			started_.resize(rdid - cur_ + 1);
			finished_.resize(rdid - cur_ + 1);
			for(size_t i = oldsz; i < lines_.size(); i++) {
				started_[i] = finished_[i] = false;
			}
		}
		started_[rdid - cur_] = true;
		finished_[rdid - cur_] = false;
	}
}

void OutputQueue::beginRead(TReadId rdid, size_t threadId) {
	if(reorder_ && threadSafe_) {
		ThreadSafe ts(mutex_m);
		beginReadImpl(rdid, threadId);
	} else {
		beginReadImpl(rdid, threadId);
	}
}

/**
 * Writer is finished writing to 
 */
void OutputQueue::finishReadImpl(const BTString& rec, TReadId rdid, size_t threadId) {
	if(reorder_) {
		assert_geq(rdid, cur_);
		assert_eq(lines_.size(), finished_.size());
		assert_eq(lines_.size(), started_.size());
		assert_lt(rdid - cur_, lines_.size());
		assert(started_[rdid - cur_]);
		assert(!finished_[rdid - cur_]);
		lines_[rdid - cur_] = rec;
		nfinished_++;
		finished_[rdid - cur_] = true;
		flush(false, false); // don't force; already have lock
	} else {
		// obuf_ is the OutFileBuf for the output file
		int i = 0;
		for(i=0; i < perThreadBufSize_; i++)
		{
			obuf_.writeString(perThreadBuf[threadId][i]);
			//TODO: turn these into atomics
			nfinished_++;
			nflushed_++;
		}
		perThreadCounter[threadId] = 0;
	}
}

void OutputQueue::finishRead(const BTString& rec, TReadId rdid, size_t threadId) {
	if(reorder_ || perThreadCounter[threadId] >= perThreadBufSize_) {
		if(threadSafe_) {
			ThreadSafe ts(mutex_m);
			finishReadImpl(rec, rdid, threadId);
		} else {
			finishReadImpl(rec, rdid, threadId);
		}
	}
	if(!reorder_) {
		perThreadBuf[threadId][perThreadCounter[threadId]++] = rec;
	}
}

/**
 * Write already-finished lines starting from cur_.
 */
void OutputQueue::flushImpl(bool force) {
	if(!reorder_) {
		size_t i = 0;
		int j = 0;
		for(i=0;i<nthreads_;i++)
		{
			for(j=0;j<perThreadCounter[i];j++)
			{
				obuf_.writeString(perThreadBuf[i][j]);
				nfinished_++;
				nflushed_++;
			}
			perThreadCounter[i]=0;
		}
		return;
	}
	size_t nflush = 0;
	while(nflush < finished_.size() && finished_[nflush]) {
		assert(started_[nflush]);
		nflush++;
	}
	// Waiting until we have several in a row to flush cuts down on copies
	// (but requires more buffering)
	if(force || nflush >= NFLUSH_THRESH) {
		for(size_t i = 0; i < nflush; i++) {
			assert(started_[i]);
			assert(finished_[i]);
			obuf_.writeString(lines_[i]);
		}
		lines_.erase(0, nflush);
		started_.erase(0, nflush);
		finished_.erase(0, nflush);
		cur_ += nflush;
		nflushed_ += nflush;
	}
}

/**
 * Write already-finished lines starting from cur_.
 */
void OutputQueue::flush(bool force, bool getLock) {
	if(getLock && threadSafe_) {
		ThreadSafe ts(mutex_m);
		flushImpl(force);
	} else {
		flushImpl(force);
	}
}

#ifdef OUTQ_MAIN

#include <iostream>

using namespace std;

int main(void) {
	cerr << "Case 1 (one thread) ... ";
	{
		OutFileBuf ofb;
		OutputQueue oq(ofb, false);
		assert_eq(0, oq.numFlushed());
		assert_eq(0, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(1);
		assert_eq(0, oq.numFlushed());
		assert_eq(1, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(3);
		assert_eq(0, oq.numFlushed());
		assert_eq(2, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(2);
		assert_eq(0, oq.numFlushed());
		assert_eq(3, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(3, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.beginRead(0);
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(0, oq.numFinished());
		oq.finishRead(0);
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.flush();
		assert_eq(0, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.flush(true);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(1, oq.numFinished());
		oq.finishRead(2);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(2, oq.numFinished());
		oq.flush(true);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(2, oq.numFinished());
		oq.finishRead(1);
		assert_eq(1, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(3, oq.numFinished());
		oq.flush(true);
		assert_eq(3, oq.numFlushed());
		assert_eq(4, oq.numStarted());
		assert_eq(3, oq.numFinished());
	}
	cerr << "PASSED" << endl;

	cerr << "Case 2 (one thread) ... ";
	{
		OutFileBuf ofb;
		OutputQueue oq(ofb, false);
		BTString& buf1 = oq.beginRead(0);
		BTString& buf2 = oq.beginRead(1);
		BTString& buf3 = oq.beginRead(2);
		BTString& buf4 = oq.beginRead(3);
		BTString& buf5 = oq.beginRead(4);
		assert_eq(5, oq.numStarted());
		assert_eq(0, oq.numFinished());
		buf1.install("A\n");
		buf2.install("B\n");
		buf3.install("C\n");
		buf4.install("D\n");
		buf5.install("E\n");
		oq.finishRead(4);
		oq.finishRead(1);
		oq.finishRead(0);
		oq.finishRead(2);
		oq.finishRead(3);
		oq.flush(true);
		assert_eq(5, oq.numFlushed());
		assert_eq(5, oq.numStarted());
		assert_eq(5, oq.numFinished());
		ofb.flush();
	}
	cerr << "PASSED" << endl;
	return 0;
}

#endif /*def ALN_SINK_MAIN*/
