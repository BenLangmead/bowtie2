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

#ifndef INDEX_LOAD_H_
#define INDEX_LOAD_H_

/**
 * index_load.h:
 *
 * Reads the large arrays of an index (forward and mirror BWT, SA sample,
 * bit-packed reference) on helper threads while the main thread goes on
 * reading the rest of the index files in the usual order.  Reading the
 * arrays one after another is mostly the kernel copying from the page
 * cache, one file at a time; on several threads the copies overlap.
 *
 * Each deferred array is read in pieces of CHUNK bytes, queued in file
 * order.  The helpers take pieces from the queue, and so does the main
 * thread once it reaches finish(), so the threads stay busy until the
 * last piece whatever the arrays' sizes.  A piece is read through its own
 * stream, opened when the piece is taken and closed when it is read, so
 * only a few streams are open at once; the first piece uses the array's
 * own stream, opened by tryDefer().
 *
 * The caller allocates each array exactly as a serial load does, then
 * calls tryDefer() where the serial load would have read it.  tryDefer()
 * hands the read over only when a stream on the file can be opened and
 * positioned and the file holds the whole array; otherwise it returns
 * false and the caller reads the array itself, so a short or unreadable
 * file fails with the serial load's message at the serial load's point.
 * A piece whose own stream cannot be opened later is read in finish()
 * through the array's stream.  finish() waits for every read; if one
 * failed, it prints the message the serial load would have printed for
 * the first failed array (its byte counts rebuilt from the first failed
 * piece) and throws 1.  Nothing may read a deferred array before finish().
 */

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "assert_helpers.h"

class IndexLoader {
public:

	/**
	 * How the serial load reads an array, which decides the message
	 * printed when the read fails.
	 */
	enum Kind {
		LOOP_EBWT = 0, // Ebwt::readIntoMemory's loop over _ebwt[]
		LOOP_OFFS,     // Ebwt::readIntoMemory's loop over _offs[]
		ONE_REF        // BitPairReference's single fread of the .4 file
	};

	/**
	 * Bytes per piece.  64 MB splits an hg38 index into 56 pieces, so the
	 * threads finish within about one piece's read of each other.
	 */
	static const uint64_t CHUNK = (uint64_t)64 << 20;

	/**
	 * maxThreads: the most helper threads at once (the main thread reads
	 * too, in finish()); 0 reads nothing in the background (tryDefer()
	 * always returns false).
	 */
	explicit IndexLoader(int maxThreads) :
		maxThreads_(maxThreads < 0 ? 0 : maxThreads),
		live_(0),
		next_(0)
	{ }

	/**
	 * Waits for every read still in flight, so that no helper writes
	 * into an array after an exception has freed it.  Never throws.
	 */
	~IndexLoader() {
		joinAll();
		for(size_t i = 0; i < arrays_.size(); i++) {
			if(arrays_[i]->file != NULL) fclose(arrays_[i]->file);
			delete arrays_[i];
		}
		for(size_t i = 0; i < pieces_.size(); i++) {
			delete pieces_[i];
		}
	}

	bool enabled() const { return maxThreads_ > 0; }

	/**
	 * Reads len bytes of the file at path, starting at byte off, into
	 * dest on the helper threads.  Returns false, and changes nothing the
	 * caller can see (errno included), when the read cannot be handed
	 * off: the loader is disabled, a stream cannot be opened or
	 * positioned, the file is shorter than off + len, or no helper
	 * thread can be started.
	 */
	bool tryDefer(Kind kind, const std::string& path, uint64_t off, void *dest, uint64_t len) {
		if(!enabled()) return false;
		const int savedErrno = errno;
		FILE *f = fopen(path.c_str(), "rb");
		bool ok = f != NULL;
		if(ok) ok = fseeko(f, 0, SEEK_END) == 0;
		if(ok) {
			const off_t end = ftello(f);
			ok = end >= 0 && (uint64_t)end >= off && (uint64_t)end - off >= len;
		}
		if(ok) ok = fseeko(f, (off_t)off, SEEK_SET) == 0;
		if(!ok) {
			if(f != NULL) fclose(f);
			errno = savedErrno;
			return false;
		}
		const uint64_t chunk = CHUNK;
		// One piece for an empty array, so that it is still checked
		const uint64_t npieces = (len == 0) ? 1 : (len + chunk - 1) / chunk;
		std::unique_lock<std::mutex> lk(mu_);
		const size_t firstPiece = pieces_.size();
		Array *a = NULL;
		try {
			a = new Array();
			a->kind = kind;
			a->path = path;
			a->off = off;
			a->len = len;
			a->first = firstPiece;
			a->n = (size_t)npieces;
			arrays_.reserve(arrays_.size() + 1);
			pieces_.reserve(pieces_.size() + (size_t)npieces);
			for(uint64_t k = 0; k < npieces; k++) {
				Piece *p = new Piece();
				p->array = arrays_.size();
				p->off = k * chunk;
				p->dest = (char *)dest + p->off;
				p->len = (len - p->off < chunk) ? (len - p->off) : chunk;
				pieces_.push_back(p); // room reserved above
			}
			a->file = f; // from here on the array closes it
			arrays_.push_back(a); // room reserved above
		} catch(...) {
			dropFrom(firstPiece);
			delete a;
			fclose(f);
			errno = savedErrno;
			return false;
		}
		// Start helpers for the new pieces, up to maxThreads_ in all
		while(live_ < maxThreads_ && (size_t)live_ < pieces_.size() - next_) {
			bool started = false;
			try {
				// Room first, so that the started thread is never
				// dropped by a failed push_back
				threads_.reserve(threads_.size() + 1);
				threads_.push_back(std::thread(&IndexLoader::worker, this));
				started = true;
			} catch(...) { }
			if(!started) break;
			live_++; // the new helper waits for mu_ before looking
		}
		if(live_ == 0) {
			// No helper to take the pieces: the caller reads the array
			dropFrom(firstPiece);
			arrays_.pop_back();
			delete a; // closes nothing: the file is closed here
			fclose(f);
			errno = savedErrno;
			return false;
		}
		lk.unlock();
		errno = savedErrno;
		return true;
	}

	/**
	 * Reads pieces on the main thread until none is left, waits for the
	 * helpers, then reads any piece whose own stream could not be opened
	 * through its array's stream.  If any read failed, prints the serial
	 * load's message for the first failed array (in the order they were
	 * deferred, which is the serial load's order) and throws 1.
	 */
	void finish() {
		const int savedErrno = errno;
		take(false);
		joinAll();
		for(size_t i = 0; i < pieces_.size(); i++) {
			Piece& p = *pieces_[i];
			if(!p.unopened) continue;
			Array& a = *arrays_[p.array];
			if(fseeko(a.file, (off_t)(a.off + p.off), SEEK_SET) != 0) {
				p.thrown = true;
				p.failed = true;
				continue;
			}
			run(p, a.kind, a.file, i + 1 == a.first + a.n);
		}
		for(size_t i = 0; i < arrays_.size(); i++) {
			fclose(arrays_[i]->file);
			arrays_[i]->file = NULL;
		}
		errno = savedErrno;
		for(size_t i = 0; i < arrays_.size(); i++) {
			const Array& a = *arrays_[i];
			for(size_t k = a.first; k < a.first + a.n; k++) {
				if(!pieces_[k]->failed) continue;
				std::cerr << message(a, *pieces_[k]) << std::endl;
				throw 1;
			}
		}
	}

private:

	/** One deferred array. */
	struct Array {
		Array() : kind(LOOP_EBWT), file(NULL), off(0), len(0), first(0), n(0) { }
		Kind kind;
		std::string path;
		FILE *file;   // opened by tryDefer(), positioned at off
		uint64_t off; // the array's offset in the file
		uint64_t len;
		size_t first; // its first piece in pieces_
		size_t n;     // number of pieces
	};

	/** One piece of a deferred array. */
	struct Piece {
		Piece() : array(0), dest(NULL), off(0), len(0), unopened(false),
			failed(false), ioErr(false), thrown(false), before(0),
			lastRead(0) { }
		size_t array;
		char *dest;
		uint64_t off;      // offset of the piece in its array
		uint64_t len;
		bool unopened;     // its own stream could not be opened
		bool failed;
		bool ioErr;        // failed on ferror() (else at end of file)
		bool thrown;       // failed on an exception
		uint64_t before;   // bytes read by calls before the failing one
		uint64_t lastRead; // what the failing fread() returned
	};

	/**
	 * The message the serial load prints for array a when piece p is the
	 * first piece of a that failed.  Every earlier piece read in full, so
	 * a serial read gets to p's offset and then sees what p saw: its
	 * first fread() returns everything up to the error or the end of the
	 * file, with the whole array outstanding.
	 */
	static std::string message(const Array& a, const Piece& p) {
		if(p.thrown) {
			return "Error reading index file " + a.path;
		}
		std::ostringstream os;
		if(a.kind == ONE_REF) {
			os << "Only read " << (p.off + p.before + p.lastRead)
			   << " bytes (out of " << a.len
			   << ") from reference index file " << a.path.c_str();
			return os.str();
		}
		const char *what = (a.kind == LOOP_EBWT) ?
			"Error reading _ebwt[] array: " :
			"Error reading block of _offs[] array: ";
		if(!p.ioErr) {
			return std::string(what) + "no more data";
		}
		// is_fread_err() sets gLastIOErrMsg to this text
		if(p.before == 0) {
			os << what << (p.off + p.lastRead) << ", " << a.len << "Error Reading File!";
		} else {
			os << what << p.lastRead << ", " << (a.len - p.off - p.before) << "Error Reading File!";
		}
		return os.str();
	}

	/**
	 * Reads one piece from file, positioned at the piece, with the serial
	 * loop's failure tests, keeping what message() needs instead of
	 * printing.  last: the piece ends its array.
	 */
	static void run(Piece& p, Kind kind, FILE *file, bool last) {
		uint64_t bytesLeft = p.len;
		char *dest = p.dest;
		while(bytesLeft > 0) {
			size_t r = fread((void *)dest, 1, bytesLeft, file);
			if(ferror(file)) {
				p.lastRead = r;
				p.ioErr = true;
				p.failed = true;
				return;
			} else if(r == 0) {
				p.failed = true;
				return;
			}
			dest += r;
			bytesLeft -= r;
			p.before += r;
		}
#ifndef NDEBUG
		if(last && kind == ONE_REF) {
			// Make sure there's no more
			char c;
			size_t ret = fread(&c, 1, 1, file);
			assert_eq(0, ret); // should have failed
		}
#endif
	}

	/**
	 * Reads piece p of array a: the first piece through the array's
	 * stream, any other through a stream of its own.  A piece whose stream
	 * cannot be opened or positioned is left for finish().
	 */
	static void runPiece(Piece& p, const Array& a, bool first, bool last) {
		if(first) {
			run(p, a.kind, a.file, last);
			return;
		}
		FILE *f = fopen(a.path.c_str(), "rb");
		if(f == NULL || fseeko(f, (off_t)(a.off + p.off), SEEK_SET) != 0) {
			if(f != NULL) fclose(f);
			p.unopened = true;
			return;
		}
		run(p, a.kind, f, last);
		fclose(f);
	}

	/**
	 * Takes pieces from the queue and reads them until none is left; a
	 * helper then counts itself out.
	 */
	void take(bool helper) {
		std::unique_lock<std::mutex> lk(mu_);
		while(next_ < pieces_.size()) {
			const size_t i = next_++;
			Piece *p = pieces_[i];
			const Array *a = arrays_[p->array];
			const bool first = (i == a->first);
			const bool last = (i + 1 == a->first + a->n);
			lk.unlock();
			try {
				runPiece(*p, *a, first, last);
			} catch(...) {
				p->thrown = true;
				p->failed = true;
			}
			lk.lock();
		}
		if(helper) live_--;
	}

	void worker() { take(true); }

	/** Removes the pieces from index first on, none of them taken. */
	void dropFrom(size_t first) {
		while(pieces_.size() > first) {
			delete pieces_.back();
			pieces_.pop_back();
		}
	}

	void joinAll() {
		for(size_t i = 0; i < threads_.size(); i++) {
			if(threads_[i].joinable()) threads_[i].join();
		}
	}

	int maxThreads_;
	int live_;      // helpers that will still look at the queue
	size_t next_;   // first piece not yet taken
	std::mutex mu_;
	std::vector<Array*> arrays_;
	std::vector<Piece*> pieces_;
	std::vector<std::thread> threads_;
};

#endif /* INDEX_LOAD_H_ */
