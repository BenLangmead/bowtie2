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

#ifndef FILEBUF_H_
#define FILEBUF_H_

#include <iostream>
#include <fstream>
#include <string>
#include <condition_variable>
#include <thread>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdexcept>
#include "assert_helpers.h"
#include <errno.h>
#include <stdlib.h>
#include <zlib.h>
#ifdef WITH_ZSTD
#include "zstd_decompress.h"
#endif

/**
 * Simple, fast helper for determining if a character is a newline.
 */
static inline bool isnewline(int c) {
	return c == '\r' || c == '\n';
}

/**
 * Simple, fast helper for determining if a character is a non-newline
 * whitespace character.
 */
static inline bool isspace_notnl(int c) {
	return isspace(c) && !isnewline(c);
}

/**
 * Simple wrapper for a FILE*, istream or ifstream that reads it in chunks
 * using fread and keeps those chunks in a buffer.  It also services calls to
 * get(), peek() and gets() from the buffer, reading in additional chunks when
 * necessary.
 *
 * Helper functions do things like parse strings, numbers, and FASTA records.
 *
 *
 */
class FileBuf {
public:
	FileBuf() {
		init();
	}

	FileBuf(FILE *in) {
		init();
		_in = in;
		assert(_in != NULL);
	}

	FileBuf(gzFile in) {
		init();
		_zIn = in;
		assert(_zIn != NULL);
	}
#ifdef WITH_ZSTD
	FileBuf(zstdStrm *zstdIn) {
		init();
		_zstdIn = zstdIn;
		assert(_zstdIn != NULL);
	}
#endif

	FileBuf(std::ifstream *inf) {
		init();
		_inf = inf;
		assert(_inf != NULL);
	}

	FileBuf(std::istream *ins) {
		init();
		_ins = ins;
		assert(_ins != NULL);
	}


	~FileBuf() {
		close();
	}

	/**
	 * Return true iff there is a stream ready to read.
	 */
	bool isOpen() {
		return _in != NULL || _inf != NULL || _ins != NULL;
	}

	/**
	 * Close the input stream (if that's possible)
	 */
	void close() {
		if(_in != NULL && _in != stdin) {
			fclose(_in);
		} else if(_inf != NULL) {
			_inf->close();
		} else if(_zIn != NULL) {
			gzclose(_zIn);
#ifdef WITH_ZSTD
		} else if(_zstdIn != NULL) {
			zstdClose(_zstdIn);
#endif
		} else {
			// can't close _ins
		}
	}

	/**
	 * Get the next character of input and advance.
	 */
	int get() {
		assert(_in != NULL || _zIn != NULL || _inf != NULL || _ins != NULL);
#ifdef WITH_ZSTD
		assert(_zstdIn != NULL);
#endif
		int c = peek();
		if(c != -1) {
			_cur++;
			if(_lastn_cur < LASTN_BUF_SZ) _lastn_buf[_lastn_cur++] = c;
		}
		return c;
	}

	/**
	 * Return true iff all input is exhausted.
	 */
	bool eof() {
		return (_cur == _buf_sz) && _done;
	}

	/**
	 * Initialize the buffer with a new C-style file.
	 */
	void newFile(FILE *in) {
		_in = in;
		_zIn = NULL;
		_inf = NULL;
		_ins = NULL;
#ifdef WITH_ZSTD
		_zstdIn = NULL;
#endif
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Initialize the buffer with a new gz file.
	 */
	void newFile(gzFile in) {
		_in = NULL;
		_zIn = in;
		_inf = NULL;
		_ins = NULL;
#ifdef WITH_ZSTD
		_zstdIn = NULL;
#endif
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

#ifdef WITH_ZSTD
	/**
	 * Initialize the buffer with a new ZSTD file.
	 */
	void newFile(zstdStrm *s) {
		_in = NULL;
		_zIn = NULL;
		_inf = NULL;
		_ins = NULL;
		_zstdIn = s;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}
#endif
	/**
	 * Initialize the buffer with a new ifstream.
	 */
	void newFile(std::ifstream *__inf) {
		_in = NULL;
		_zIn = NULL;
		_inf = __inf;
		_ins = NULL;
#ifdef WITH_ZSTD
		_zstdIn = NULL;
#endif
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Initialize the buffer with a new istream.
	 */
	void newFile(std::istream *__ins) {
		_in = NULL;
		_zIn = NULL;
		_inf = NULL;
		_ins = __ins;
#ifdef WITH_ZSTD
		_zstdIn = NULL;
#endif
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Restore state as though we just started reading the input
	 * stream.
	 */
	void reset() {
		if(_inf != NULL) {
			_inf->clear();
			_inf->seekg(0, std::ios::beg);
		} else if(_ins != NULL) {
			_ins->clear();
			_ins->seekg(0, std::ios::beg);
		} else if (_zIn != NULL) {
			gzrewind(_zIn);
#ifdef WITH_ZSTD
		} else if (_zstdIn) {
			zstdRewind(_zstdIn);
#endif
		} else {
			rewind(_in);
		}
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Peek at the next character of the input stream without
	 * advancing.  Typically we can simple read it from the buffer.
	 * Occasionally we'll need to read in a new buffer's worth of data.
	 */
	int peek() {
		assert(_in != NULL || _zIn != NULL || _inf != NULL || _ins != NULL);
#ifdef WITH_ZSTD
		assert(_zstdIn != NULL);
#endif
		assert_leq(_cur, _buf_sz);
		if(_cur == _buf_sz) {
			if(_done) {
				// We already exhausted the input stream
				return -1;
			}
			// Read a new buffer's worth of data
			else {
				// Get the next chunk
				if(_inf != NULL) {
					_inf->read((char*)_buf, BUF_SZ);
					_buf_sz = _inf->gcount();
				} else if(_zIn != NULL) {
					_buf_sz = gzread(_zIn, (void *)_buf, BUF_SZ);
				} else if(_ins != NULL) {
					_ins->read((char*)_buf, BUF_SZ);
					_buf_sz = _ins->gcount();
#ifdef WITH_ZSTD
                                } else if (_zstdIn != NULL) {
					_buf_sz = zstdRead(_zstdIn, (void *)_buf, BUF_SZ);
#endif
                                } else {
					assert(_in != NULL);
					// TODO: consider an _unlocked function
					_buf_sz = fread(_buf, 1, BUF_SZ, _in);
                                }
                                _cur = 0;
				if(_buf_sz == 0) {
					// Exhausted, and we have nothing to return to the
					// caller
					_done = true;
					return -1;
				} else if(_buf_sz < BUF_SZ) {
					// Exhausted
					_done = true;
				}
			}
		}
		return (int)_buf[_cur];
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until we see a newline, EOF, or until 'len' characters have been
	 * read.
	 */
	size_t gets(char *buf, size_t len) {
		size_t stored = 0;
		while(true) {
			int c = get();
			if(c == -1) {
				// End-of-file
				buf[stored] = '\0';
				return stored;
			}
			if(stored == len-1 || isnewline(c)) {
				// End of string
				buf[stored] = '\0';
				// Skip over all end-of-line characters
				int pc = peek();
				while(isnewline(pc)) {
					get(); // discard
					pc = peek();
				}
				// Next get() will be after all newline characters
				return stored;
			}
			buf[stored++] = (char)c;
		}
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until we see a newline, EOF, or until 'len' characters have been
	 * read.
	 */
	size_t get(char *buf, size_t len) {
		size_t stored = 0;
		for(size_t i = 0; i < len; i++) {
			int c = get();
			if(c == -1) return i;
			buf[stored++] = (char)c;
		}
		return len;
	}

	static const size_t LASTN_BUF_SZ = 8 * 1024;

	/**
	 * Keep get()ing characters until a non-whitespace character (or
	 * -1) is reached, and return it.
	 */
	int getPastWhitespace() {
		int c;
		while(isspace(c = get()) && c != -1);
		return c;
	}

	/**
	 * Keep get()ing characters until a we've passed over the next
	 * string of newline characters (\r's and \n's) or -1 is reached,
	 * and return it.
	 */
	int getPastNewline() {
		int c = get();
		while(!isnewline(c) && c != -1) c = get();
		while(isnewline(c)) c = get();
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}

	/**
	 * Keep get()ing characters until a we've passed over the next
	 * string of newline characters (\r's and \n's) or -1 is reached,
	 * and return it.
	 */
	int peekPastNewline() {
		int c = peek();
		while(!isnewline(c) && c != -1) c = get();
		while(isnewline(c)) c = get();
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}

	/**
	 * Keep peek()ing then get()ing characters until the next return
	 * from peek() is just after the last newline of the line.
	 */
	int peekUptoNewline() {
		int c = peek();
		while(!isnewline(c) && c != -1) {
			get(); c = peek();
		}
		while(isnewline(c)) {
			get();
			c = peek();
		}
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}
	
	/**
	 * Parse a FASTA record.  Append name characters to 'name' and and append
	 * all sequence characters to 'seq'.  If gotCaret is true, assuming the
	 * file cursor has already moved just past the starting '>' character.
	 */
	template <typename TNameStr, typename TSeqStr>
	void parseFastaRecord(
		TNameStr& name,
		TSeqStr&  seq,
		bool      gotCaret = false)
	{
		int c;
		if(!gotCaret) {
			// Skip over caret and non-newline whitespace
			c = peek();
			while(isspace_notnl(c) || c == '>') { get(); c = peek(); }
		} else {
			// Skip over non-newline whitespace
			c = peek();
			while(isspace_notnl(c)) { get(); c = peek(); }
		}
		size_t namecur = 0, seqcur = 0;
		// c is the first character of the fasta name record, or is the first
		// newline character if the name record is empty
		while(!isnewline(c) && c != -1) {
			name[namecur++] = c; get(); c = peek();
		}
		// sequence consists of all the non-whitespace characters between here
		// and the next caret
		while(true) {
			// skip over whitespace
			while(isspace(c)) { get(); c = peek(); }
			// if we see caret or EOF, break
			if(c == '>' || c == -1) break;
			// append and continue
			seq[seqcur++] = c;
			get(); c = peek();
		}
	}

	/**
	 * Parse a FASTA record and return its length.  If gotCaret is true,
	 * assuming the file cursor has already moved just past the starting '>'
	 * character.
	 */
	void parseFastaRecordLength(
		size_t&   nameLen,
		size_t&   seqLen,
		bool      gotCaret = false)
	{
		int c;
		nameLen = seqLen = 0;
		if(!gotCaret) {
			// Skip over caret and non-newline whitespace
			c = peek();
			while(isspace_notnl(c) || c == '>') { get(); c = peek(); }
		} else {
			// Skip over non-newline whitespace
			c = peek();
			while(isspace_notnl(c)) { get(); c = peek(); }
		}
		// c is the first character of the fasta name record, or is the first
		// newline character if the name record is empty
		while(!isnewline(c) && c != -1) {
			nameLen++; get(); c = peek();
		}
		// sequence consists of all the non-whitespace characters between here
		// and the next caret
		while(true) {
			// skip over whitespace
			while(isspace(c)) { get(); c = peek(); }
			// if we see caret or EOF, break
			if(c == '>' || c == -1) break;
			// append and continue
			seqLen++;
			get(); c = peek();
		}
	}

	/**
	 * Reset to the beginning of the last-N-chars buffer.
	 */
	void resetLastN() {
		_lastn_cur = 0;
	}

	/**
	 * Copy the last several characters in the last-N-chars buffer
	 * (since the last reset) into the provided buffer.
	 */
	size_t copyLastN(char *buf) {
		memcpy(buf, _lastn_buf, _lastn_cur);
		return _lastn_cur;
	}

	/**
	 * Get const pointer to the last-N-chars buffer.
	 */
	const char *lastN() const {
		return _lastn_buf;
	}

	/**
	 * Get current size of the last-N-chars buffer.
	 */
	size_t lastNLen() const {
		return _lastn_cur;
	}

private:

	void init() {
		_in = NULL;
		_zIn = NULL;
		_inf = NULL;
		_ins = NULL;
#ifdef WITH_ZSTD
		_zstdIn = NULL;
#endif
		_cur = _buf_sz = BUF_SZ;
		_done = false;
		_lastn_cur = 0;
		// no need to clear _buf[]
	}

	static const size_t BUF_SZ = 256 * 1024;
	FILE     *_in;
	gzFile   _zIn;
#ifdef WITH_ZSTD
	zstdStrm *_zstdIn;
#endif
	std::ifstream *_inf;
	std::istream  *_ins;
	size_t    _cur;
	size_t    _buf_sz;
	bool      _done;
	uint8_t   _buf[BUF_SZ]; // (large) input buffer
	size_t    _lastn_cur;
	char      _lastn_buf[LASTN_BUF_SZ]; // buffer of the last N chars dispensed
};

/**
 * Wrapper for a buffered output stream that writes bitpairs.
 */
class BitpairOutFileBuf {
public:
	/**
	 * Open a new output stream to a file with given name.
	 */
	BitpairOutFileBuf(const char *in) : bpPtr_(0), cur_(0) {
		assert(in != NULL);
		out_ = fopen(in, "wb");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open bitpair-output file " << in << std::endl;
			throw 1;
		}
		memset(buf_, 0, BUF_SZ);
	}

	/**
	 * Write a single bitpair into the buf.  Flush the buffer if it's
	 * full.
	 */
	void write(int bp) {
		assert_lt(bp, 4);
		assert_geq(bp, 0);
		buf_[cur_] |= (bp << bpPtr_);
		if(bpPtr_ == 6) {
			bpPtr_ = 0;
			cur_++;
			if(cur_ == BUF_SZ) {
				// Flush the buffer
				if(!fwrite((const void *)buf_, BUF_SZ, 1, out_)) {
					std::cerr << "Error writing to the reference index file (.4.ebwt)" << std::endl;
					throw 1;
				}
				// Reset to beginning of the buffer
				cur_ = 0;
			}
			// Initialize next octet to 0
			buf_[cur_] = 0;
		} else {
			bpPtr_ += 2;
		}
	}

	/**
	 * Write any remaining bitpairs and then close the input
	 */
	void close() {
		if(cur_ > 0 || bpPtr_ > 0) {
			if(bpPtr_ == 0) cur_--;
			if(!fwrite((const void *)buf_, cur_ + 1, 1, out_)) {
				std::cerr << "Error writing to the reference index file (.4.ebwt)" << std::endl;
				throw 1;
			}
		}
		fclose(out_);
	}
private:
	static const size_t BUF_SZ = 128 * 1024;
	FILE    *out_;
	int      bpPtr_;
	size_t   cur_;
	char     buf_[BUF_SZ]; // (large) input buffer
};

/**
 * Wrapper for a buffered output stream that writes characters and
 * other data types.  This class is *not* synchronized; the caller is
 * responsible for synchronization.
 */
class OutFileBuf {

public:

	/**
	 * Open a new output stream to a file with given name.
	 */
	OutFileBuf(const std::string& out, bool binary = false) :
		name_(out.c_str()), out_(NULL), cur_(0), closed_(false),
		asyncData_(out_), asynct_(writeAsync, &asyncData_),
		buf1_(new char[BUF_SZ]), buf2_(new char[BUF_SZ]), cap1_(BUF_SZ), cap2_(BUF_SZ),
		buf_(buf1_), cap_(cap1_)
	{
		out_ = fopen(out.c_str(), binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out.c_str() << std::endl;
			throw 1;
		}
		if(setvbuf(out_, NULL, _IOFBF, 10* 1024* 1024)) 
			std::cerr << "Warning: Could not allocate the proper buffer size for output file stream. " << std::endl;
	}

	/**
	 * Open a new output stream to a file with given name.
	 */
	OutFileBuf(const char *out, bool binary = false) :
		name_(out), out_(NULL), cur_(0), closed_(false),
		asyncData_(out_), asynct_(writeAsync, &asyncData_),
		buf1_(new char[BUF_SZ]), buf2_(new char[BUF_SZ]), cap1_(BUF_SZ), cap2_(BUF_SZ),
		buf_(buf1_), cap_(cap1_)
	{
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
	}

	/**
	 * Open a new output stream to standard out.
	 */
	OutFileBuf() :
		name_("cout"), out_(stdout), cur_(0), closed_(false),
		asyncData_(out_), asynct_(writeAsync, &asyncData_),
		buf1_(new char[BUF_SZ]), buf2_(new char[BUF_SZ]), cap1_(BUF_SZ), cap2_(BUF_SZ),
		buf_(buf1_), cap_(cap1_) {}
	
	/**
	 * Close buffer when object is destroyed.
	 */
	~OutFileBuf() {
		close();
		asyncData_.notifyAbort();
		asynct_.join();
		delete[] buf2_;
		delete[] buf1_;
	}

	/**
	 * Open a new output stream to a file with given name.
	 */
	void setFile(const char *out, bool binary = false) {
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
		reset();
	}

	/**
	 * Write a single character into the write buffer and, if
	 * necessary, flush.
	 */
	void write(char c) {
		assert(!closed_);
		if(cur_ == cap_) {
			if(flushBlocking() && (cap_<MAX_BUF_SZ)) {
				// minimize blocking, increase buffer instead
				increaseBuffer(BUF_SZ);
			} else {
				flush();
			}
		}
		buf_[cur_++] = c;
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	void writeString(const std::string& s) {
		assert(!closed_);
		size_t slen = s.length();
		const char *cstr = s.c_str();
		writeChars(cstr, slen);
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	template<typename T>
	void writeString(const T& s) {
		assert(!closed_);
		size_t slen = s.length();
		const char *zbuf = s.toZBuf();
		writeChars(zbuf, slen);
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	void writeChars(const char * s, size_t len) {
		assert(!closed_);
		if(cur_ + len > cap_) {
			if(flushBlocking() && (cap_<MAX_BUF_SZ)) {
				// minimize blocking, increase buffer instead
				resizeBuffer(cur_ + len + BUF_SZ); // add a little spare
			} else {
				if(cur_ > 0) flush();
				// expand the buffer, if needed
				// so we can keep using async writes
				if(len >= cap_) resizeBufferNoCopy(len+BUF_SZ); // add a little spare
			}
		}
		memcpy(&buf_[cur_], s, len);
		cur_ += len;
		assert_leq(cur_, cap_);
	}

	/**
	 * Write a 0-terminated C string to the output stream.
	 */
	void writeChars(const char * s) {
		writeChars(s, strlen(s));
	}

	/**
	 * Write any remaining bitpairs and then close the input
	 */
	void close() {
		if(closed_) return;
		if(cur_ > 0) flush();
		asyncData_.waitIdle();
		closed_ = true;
		if(out_ != stdout) {
			fclose(out_);
		}
	}

	/**
	 * Reset so that the next write is as though it's the first.
	 */
	void reset() {
		cur_ = 0;
		closed_ = false;
	}

	void flush() {
		// there still could have been an outstanding async write
		asyncData_.waitIdle();
		// start the async write
		asyncData_.setBuf(buf_, cur_);
		// switch to the other buffer
		if(buf_==buf1_) {
			buf_ = buf2_;
			cap_ = cap2_;
		} else {
			buf_ = buf1_;
			cap_ = cap1_;
		}
		cur_ = 0;
	}

	/**
	 * Return true iff this stream is closed.
	 */
	bool closed() const {
		return closed_;
	}

	/**
	 * Return the filename.
	 */
	const char *name() {
		return name_;
	}

	bool flushBlocking() const {
		return asyncData_.buf!=NULL;
	}

private:
	size_t roundBufferSize(size_t len) const {
		// do exponential increases to reduce number of resizes during the process lifetime
		// but only up to MAX_BUF_SZ
		// Note: Cannot use std::min or max, as it will require MAX_BUF_SZ to have storage (before c++17)
		size_t dsize = cap_*2;
		if (dsize>MAX_BUF_SZ) dsize = MAX_BUF_SZ;
		if (len<dsize) len=dsize;
		// round up to multiple of BUF_SZ
		return ((len+BUF_SZ-1)/BUF_SZ)*BUF_SZ;
	}

	static char* _resizeSpecificBufferNoCopy(
						const size_t newCap, 
						char* &buf, size_t &cap) { // buffer to resize
		delete[] buf;
		buf = new char[newCap];
		cap = newCap;
		return buf;
	}



	// increase the current buffer to at least len
	void resizeBufferNoCopy(size_t len) {
		assert_eq(cur_, 0);
		// round up to multiple of BUF_SZ
		const size_t newCap = roundBufferSize(len);
		bool is1 = (buf_==buf1_);
		buf_ = _resizeSpecificBufferNoCopy(newCap,
						is1 ? buf1_ : buf2_,
						is1 ? cap1_ : cap2_);
		cap_ = newCap;
	}


	static char* _resizeSpecificBufferCopy(
						const size_t newCap,
						const size_t oldDataSize,  // number of bytes to preserve
						char* &buf, size_t &cap) { // buffer to resize
		const char* oldBuf = buf;
		buf = new char[newCap];
		memcpy(buf, oldBuf, oldDataSize);
		delete[] oldBuf;
		cap = newCap;
		return buf;
	}

	// increase the current buffer to at least len
	// copy over the content
	void resizeBuffer(size_t len) {
		if(cur_==0) {
			// nothing to copy, use the more efficient version
			resizeBufferNoCopy(len);
			return;
		}

		// round up to multiple of BUF_SZ
		size_t newCap = roundBufferSize(len);
		bool is1 = (buf_==buf1_);
		buf_ = _resizeSpecificBufferCopy(newCap, cur_,
						is1 ? buf1_ : buf2_,
						is1 ? cap1_ : cap2_);
		cap_ = newCap;
	}

	// increase the current buffer by delta
	// copy over the content
	void increaseBuffer(size_t delta) {
		size_t newCap = cap_+delta;
		resizeBuffer(newCap);
	}

	class AsyncData {
	public:
		bool abort; // the async thread will abort when this is set to true
		FILE*      &out;

		const char* buf;
		size_t      cur;

		std::mutex m;
		std::condition_variable cv;

		// m and cv default constructors are OK as-is
		AsyncData(FILE* &_out) : abort(false), out(_out), buf(NULL) {}

		void notifyAbort() {
			{
				std::lock_guard<std::mutex> lk(m);
				abort = true;
			}
			cv.notify_all();
		}

		void waitIdle() {
			std::unique_lock<std::mutex> lk(m);
			while(buf!=NULL) cv.wait(lk);
		}

		void setBuf(const char* _buf, size_t _cur) {
			{
				std::lock_guard<std::mutex> lk(m);
				buf = _buf;
				cur = _cur;
			}
			cv.notify_all();
		}

		// returns abort
		bool waitForBuf() {
			std::unique_lock<std::mutex> lk(m);
			while((buf==NULL)&&(!abort)) cv.wait(lk);
			return abort;
		}

		// returns abort
		bool writeComplete() {
			bool ret;
			{
				std::lock_guard<std::mutex> lk(m);
				buf = NULL;
				ret = abort;
			}
			cv.notify_all();
			return ret;
		}



	};

	static void writeAsync(AsyncData *asyncDataPtr) {
		AsyncData &asyncData = *asyncDataPtr;
		bool abort = false;
		while(!abort) {
			abort = asyncData.waitForBuf();
			if(abort) break;
			if(asyncData.cur != fwrite((const void *)asyncData.buf, 1, asyncData.cur, asyncData.out)) {
				if (errno == EPIPE) {
					exit(EXIT_SUCCESS);
				}
				std::cerr << "Error while flushing and closing output" << std::endl;
				throw 1;
			}
			abort = asyncData.writeComplete();
		}

	}

	static constexpr size_t BUF_SZ = 16ul * 1024ul;
	static constexpr size_t MAX_BUF_SZ = 16ul * 1024ul * 1024ul * 1024ul;

	const char *name_;
	FILE       *out_;
	size_t      cur_;  // how much of the buffer is currently filled
	bool        closed_;
	AsyncData   asyncData_;
	std::thread asynct_;
	char       *buf1_; // (large) output buffer
	char       *buf2_; // (large) output buffer
	size_t      cap1_; //capacity of buf1_
	size_t      cap2_; //capacity of buf2_
	char*       buf_; // points to one of the two buffers below
	size_t      cap_; // capacity of the pointed buffer
};

#endif /*ndef FILEBUF_H_*/
