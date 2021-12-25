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

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "bt2_idx.h"
#include <iomanip>

using namespace std;

///////////////////////////////////////////////////////////////////////
//
// Functions for reading and writing Ebwts
//
///////////////////////////////////////////////////////////////////////

/**
 * Read an Ebwt from file with given filename.
 */
void Ebwt::readIntoMemory(
	int color,
	int needEntireRev,
	bool loadSASamp,
	bool loadFtab,
	bool loadRstarts,
	bool justHeader,
	EbwtParams *params,
	bool mmSweep,
	bool loadNames,
	bool startVerbose)
{
#ifdef BOWTIE_MM
	char *mmFile[] = { NULL, NULL };
#endif
	if(_in1Str.length() > 0) {
		if(_verbose || startVerbose) {
			cerr << "  About to open input files: ";
			logTime(cerr);
		}
		// Initialize our primary and secondary input-stream fields
		if(_in1 != NULL) fclose(_in1);
		if(_verbose || startVerbose) cerr << "Opening \"" << _in1Str.c_str() << "\"" << endl;
		if((_in1 = fopen(_in1Str.c_str(), "rb")) == NULL) {
			cerr << "Could not open index file " << _in1Str.c_str() << endl;
		}
		if(loadSASamp) {
			if(_in2 != NULL) fclose(_in2);
			if(_verbose || startVerbose) cerr << "Opening \"" << _in2Str.c_str() << "\"" << endl;
			if((_in2 = fopen(_in2Str.c_str(), "rb")) == NULL) {
				cerr << "Could not open index file " << _in2Str.c_str() << endl;
			}
		}
		if(_verbose || startVerbose) {
			cerr << "  Finished opening input files: ";
			logTime(cerr);
		}

#ifdef BOWTIE_MM
		if(_useMm /*&& !justHeader*/) {
			const char *names[] = {_in1Str.c_str(), _in2Str.c_str()};
			int fds[] = { fileno(_in1), fileno(_in2) };
			for(int i = 0; i < (loadSASamp ? 2 : 1); i++) {
				if(_verbose || startVerbose) {
					cerr << "  Memory-mapping input file " << (i+1) << ": ";
					logTime(cerr);
				}
				struct stat sbuf;
				if (stat(names[i], &sbuf) == -1) {
					perror("stat");
					cerr << "Error: Could not stat index file " << names[i] << " prior to memory-mapping" << endl;
					throw 1;
				}
				mmFile[i] = (char*)mmap((void *)0, (size_t)sbuf.st_size,
										PROT_READ, MAP_SHARED, fds[(size_t)i], 0);
				if(mmFile[i] == (void *)(-1)) {
					perror("mmap");
					cerr << "Error: Could not memory-map the index file " << names[i] << endl;
					throw 1;
				}
				if(mmSweep) {
					int sum = 0;
					for(off_t j = 0; j < sbuf.st_size; j += 1024) {
						sum += (int) mmFile[i][j];
					}
					if(startVerbose) {
						cerr << "  Swept the memory-mapped ebwt index file 1; checksum: " << sum << ": ";
						logTime(cerr);
					}
				}
			}
			mmFile1_ = mmFile[0];
			mmFile2_ = loadSASamp ? mmFile[1] : NULL;
		}
#endif
	}
#ifdef BOWTIE_MM
	else if(_useMm && !justHeader) {
		mmFile[0] = mmFile1_;
		mmFile[1] = mmFile2_;
	}
	if(_useMm && !justHeader) {
		assert(mmFile[0] == mmFile1_);
		assert(mmFile[1] == mmFile2_);
	}
#endif

	if(_verbose || startVerbose) {
		cerr << "  Reading header: ";
		logTime(cerr);
	}

	// Read endianness hints from both streams
	uint64_t bytesRead = 0;
	_switchEndian = false;
	uint32_t one = readU<uint32_t>(_in1, _switchEndian); // 1st word of primary stream
	bytesRead += 4;
	if(loadSASamp) {
#ifndef NDEBUG
		assert_eq(one, readU<uint32_t>(_in2, _switchEndian)); // should match!
#else
		readU<uint32_t>(_in2, _switchEndian);
#endif
	}
	if(one != 1) {
		assert_eq((1u<<24), one);
		assert_eq(1, endianSwapU32(one));
		_switchEndian = true;
	}

	// Can't switch endianness and use memory-mapped files; in order to
	// support this, someone has to modify the file to switch
	// endiannesses appropriately, and we can't do this inside Bowtie
	// or we might be setting up a race condition with other processes.
	if(_switchEndian && _useMm) {
		cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
		throw 1;
	}

	// Reads header entries one by one from primary stream
	TIndexOffU len          = readU<TIndexOffU>(_in1, _switchEndian);
	bytesRead += OFF_SIZE;
	int32_t  lineRate     = readI<int32_t>(_in1, _switchEndian);
	bytesRead += 4;
	/*int32_t  linesPerSide =*/ readI<int32_t>(_in1, _switchEndian);
	bytesRead += 4;
	int32_t  offRate      = readI<int32_t>(_in1, _switchEndian);
	bytesRead += 4;
	// TODO: add isaRate to the actual file format (right now, the
	// user has to tell us whether there's an ISA sample and what the
	// sampling rate is.
	int32_t  ftabChars    = readI<int32_t>(_in1, _switchEndian);
	bytesRead += 4;
	// chunkRate was deprecated in an earlier version of Bowtie; now
	// we use it to hold flags.
	int32_t flags = readI<int32_t>(_in1, _switchEndian);
	bool entireRev = false;
	if(flags < 0 && (((-flags) & EBWT_COLOR) != 0)) {
		if(color != -1 && !color) {
			cerr << "Error: -C was not specified when running bowtie, but index is in colorspace.  If" << endl
			     << "your reads are in colorspace, please use the -C option.  If your reads are not" << endl
			     << "in colorspace, please use a normal index (one built without specifying -C to" << endl
			     << "bowtie-build)." << endl;
			throw 1;
		}
		color = 1;
	} else if(flags < 0) {
		if(color != -1 && color) {
			cerr << "Error: -C was specified when running bowtie, but index is not in colorspace.  If" << endl
			     << "your reads are in colorspace, please use a colorspace index (one built using" << endl
			     << "bowtie-build -C).  If your reads are not in colorspace, don't specify -C when" << endl
			     << "running bowtie." << endl;
			throw 1;
		}
		color = 0;
	}
	if(flags < 0 && (((-flags) & EBWT_ENTIRE_REV) == 0)) {
		if(needEntireRev != -1 && needEntireRev != 0) {
			cerr << "Error: This index is compatible with 0.* versions of Bowtie, but not with 2.*" << endl
			     << "versions.  Please build or download a version of the index that is compitble" << endl
				 << "with Bowtie 2.* (i.e. built with bowtie-build 2.* or later)" << endl;
			throw 1;
		}
	} else entireRev = true;
	bytesRead += 4;

	// Create a new EbwtParams from the entries read from primary stream
	EbwtParams *eh;
	bool deleteEh = false;
	if(params != NULL) {
		params->init(len, lineRate, offRate, ftabChars, color, entireRev);
		if(_verbose || startVerbose) params->print(cerr);
		eh = params;
	} else {
		eh = new EbwtParams(len, lineRate, offRate, ftabChars, color, entireRev);
		deleteEh = true;
	}

	// Set up overridden suffix-array-sample parameters
	TIndexOffU offsLen = eh->_offsLen;
	uint64_t offsSz = eh->_offsSz;
	TIndexOffU offRateDiff = 0;
	TIndexOffU offsLenSampled = offsLen;
	if(_overrideOffRate > offRate) {
		offRateDiff = _overrideOffRate - offRate;
	}
	if(offRateDiff > 0) {
		offsLenSampled >>= offRateDiff;
		if((offsLen & ~(OFF_MASK << offRateDiff)) != 0) {
			offsLenSampled++;
		}
	}

	// Can't override the offrate or isarate and use memory-mapped
	// files; ultimately, all processes need to copy the sparser sample
	// into their own memory spaces.
	if(_useMm && (offRateDiff)) {
		cerr << "Error: Can't use memory-mapped files when the offrate is overridden" << endl;
		throw 1;
	}

	// Read nPat from primary stream
	this->_nPat = readI<TIndexOffU>(_in1, _switchEndian);
	bytesRead += OFF_SIZE;
	_plen.reset();
	// Read plen from primary stream
	if(_useMm) {
#ifdef BOWTIE_MM
		_plen.init((TIndexOffU*)(mmFile[0] + bytesRead), _nPat, false);
		bytesRead += _nPat*OFF_SIZE;
		fseeko(_in1, _nPat*OFF_SIZE, SEEK_CUR);
#endif
	} else {
		try {
			if(_verbose || startVerbose) {
				cerr << "Reading plen (" << this->_nPat << "): ";
				logTime(cerr);
			}
			_plen.init(new TIndexOffU[_nPat], _nPat, true);
			if(_switchEndian) {
				for(TIndexOffU i = 0; i < this->_nPat; i++) {
					plen()[i] = readU<TIndexOffU>(_in1, _switchEndian);
				}
			} else {
				size_t r = MM_READ(_in1, (void*)(plen()), _nPat*OFF_SIZE);
				if(r != (size_t)(_nPat*OFF_SIZE)) {
					cerr << "Error reading _plen[] array: " << r << ", " << _nPat*OFF_SIZE << endl;
					throw 1;
				}
			}
		} catch(bad_alloc& e) {
			cerr << "Out of memory allocating plen[] in Ebwt::read()"
			<< " at " << __FILE__ << ":" << __LINE__ << endl;
			throw e;
		}
	}

	bool shmemLeader;

	// TODO: I'm not consistent on what "header" means.  Here I'm using
	// "header" to mean everything that would exist in memory if we
	// started to build the Ebwt but stopped short of the build*() step
	// (i.e. everything up to and including join()).
	if(justHeader) goto done;

	this->_nFrag = readU<TIndexOffU>(_in1, _switchEndian);
	bytesRead += OFF_SIZE;
	if(_verbose || startVerbose) {
		cerr << "Reading rstarts (" << this->_nFrag*3 << "): ";
		logTime(cerr);
	}
	// assert_geq(this->_nFrag, this->_nPat);
	_rstarts.reset();
	if(loadRstarts) {
		if(_useMm) {
#ifdef BOWTIE_MM
			_rstarts.init((TIndexOffU*)(mmFile[0] + bytesRead), _nFrag*3, false);
			bytesRead += this->_nFrag*OFF_SIZE*3;
			fseeko(_in1, this->_nFrag*OFF_SIZE*3, SEEK_CUR);
#endif
		} else {
			_rstarts.init(new TIndexOffU[_nFrag*3], _nFrag*3, true);
			if(_switchEndian) {
				for(TIndexOffU i = 0; i < this->_nFrag*3; i += 3) {
					// fragment starting position in joined reference
					// string, text id, and fragment offset within text
					this->rstarts()[i]   = readU<TIndexOffU>(_in1, _switchEndian);
					this->rstarts()[i+1] = readU<TIndexOffU>(_in1, _switchEndian);
					this->rstarts()[i+2] = readU<TIndexOffU>(_in1, _switchEndian);
				}
			} else {
				size_t r = MM_READ(_in1, (void *)rstarts(), this->_nFrag*OFF_SIZE*3);
				if(r != (size_t)(this->_nFrag*OFF_SIZE*3)) {
					cerr << "Error reading _rstarts[] array: " << r << ", " << (this->_nFrag*OFF_SIZE*3) << endl;
					throw 1;
				}
			}
		}
	} else {
		// Skip em
		assert(rstarts() == NULL);
		bytesRead += this->_nFrag*OFF_SIZE*3;
		fseeko(_in1, this->_nFrag*OFF_SIZE*3, SEEK_CUR);
	}

	_ebwt.reset();
	if(_useMm) {
#ifdef BOWTIE_MM
		_ebwt.init((uint8_t*)(mmFile[0] + bytesRead), eh->_ebwtTotLen, false);
		bytesRead += eh->_ebwtTotLen;
		fseek(_in1, eh->_ebwtTotLen, SEEK_CUR);
#endif
	} else {
		// Allocate ebwt (big allocation)
		if(_verbose || startVerbose) {
			cerr << "Reading ebwt (" << eh->_ebwtTotLen << "): ";
			logTime(cerr);
		}
		bool shmemLeader = true;
		if(useShmem_) {
			uint8_t *tmp = NULL;
			shmemLeader = ALLOC_SHARED_U8(
				(_in1Str + "[ebwt]"), eh->_ebwtTotLen, &tmp,
				"ebwt[]", (_verbose || startVerbose));
			assert(tmp != NULL);
			_ebwt.init(tmp, eh->_ebwtTotLen, false);
			if(_verbose || startVerbose) {
				cerr << "  shared-mem " << (shmemLeader ? "leader" : "follower") << endl;
			}
		} else {
			try {
				_ebwt.init(new uint8_t[eh->_ebwtTotLen], eh->_ebwtTotLen, true);
			} catch(bad_alloc& e) {
				cerr << "Out of memory allocating the ebwt[] array for the Bowtie index.  Please try" << endl
				<< "again on a computer with more memory." << endl;
				throw 1;
			}
		}
		if(shmemLeader) {
			// Read ebwt from primary stream
			uint64_t bytesLeft = eh->_ebwtTotLen;
			char *pebwt = (char*)this->ebwt();

			while (bytesLeft>0){
				size_t r = MM_READ(_in1, (void *)pebwt, bytesLeft);
				if(MM_IS_IO_ERR(_in1,r,bytesLeft)) {
					cerr << "Error reading _ebwt[] array: " << r << ", "
						 << bytesLeft << gLastIOErrMsg << endl;
					throw 1;
				} else if (r == 0) {
					cerr << "Error reading _ebwt[] array: no more data" << endl;
					throw 1;
				}
				pebwt += r;
				bytesLeft -= r;
			}
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) NOTIFY_SHARED(ebwt(), eh->_ebwtTotLen);
#endif
		} else {
			// Seek past the data and wait until master is finished
			fseeko(_in1, eh->_ebwtTotLen, SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) WAIT_SHARED(ebwt(), eh->_ebwtTotLen);
#endif
		}
	}

	// Read zOff from primary stream
	_zOff = readU<TIndexOffU>(_in1, _switchEndian);
	bytesRead += OFF_SIZE;
	assert_lt(_zOff, len);

	try {
		// Read fchr from primary stream
		if(_verbose || startVerbose) cerr << "Reading fchr (5)" << endl;
		_fchr.reset();
		if(_useMm) {
#ifdef BOWTIE_MM
			_fchr.init((TIndexOffU*)(mmFile[0] + bytesRead), 5, false);
			bytesRead += 5*OFF_SIZE;
			fseek(_in1, 5*OFF_SIZE, SEEK_CUR);
#endif
		} else {
			_fchr.init(new TIndexOffU[5], 5, true);
			for(int i = 0; i < 5; i++) {
				this->fchr()[i] = readU<TIndexOffU>(_in1, _switchEndian);
				assert_leq(this->fchr()[i], len);
				assert(i <= 0 || this->fchr()[i] >= this->fchr()[i-1]);
			}
		}
		assert_gt(this->fchr()[4], this->fchr()[0]);
		// Read ftab from primary stream
		if(_verbose || startVerbose) {
			if(loadFtab) {
				cerr << "Reading ftab (" << eh->_ftabLen << "): ";
				logTime(cerr);
			} else {
				cerr << "Skipping ftab (" << eh->_ftabLen << "): ";
			}
		}
		_ftab.reset();
		if(loadFtab) {
			if(_useMm) {
#ifdef BOWTIE_MM
				_ftab.init((TIndexOffU*)(mmFile[0] + bytesRead), eh->_ftabLen, false);
				bytesRead += eh->_ftabLen*OFF_SIZE;
				fseeko(_in1, eh->_ftabLen*OFF_SIZE, SEEK_CUR);
#endif
			} else {
				_ftab.init(new TIndexOffU[eh->_ftabLen], eh->_ftabLen, true);
				if(_switchEndian) {
					for(TIndexOffU i = 0; i < eh->_ftabLen; i++)
						this->ftab()[i] = readU<TIndexOffU>(_in1, _switchEndian);
				} else {
					size_t r = MM_READ(_in1, (void *)ftab(), eh->_ftabLen*OFF_SIZE);
					if(r != (size_t)(eh->_ftabLen*OFF_SIZE)) {
						cerr << "Error reading _ftab[] array: " << r << ", " << (eh->_ftabLen*OFF_SIZE) << endl;
						throw 1;
					}
				}
			}
			// Read etab from primary stream
			if(_verbose || startVerbose) {
				if(loadFtab) {
					cerr << "Reading eftab (" << eh->_eftabLen << "): ";
					logTime(cerr);
				} else {
					cerr << "Skipping eftab (" << eh->_eftabLen << "): ";
				}

			}
			_eftab.reset();
			if(_useMm) {
#ifdef BOWTIE_MM
				_eftab.init((TIndexOffU*)(mmFile[0] + bytesRead), eh->_eftabLen, false);
				bytesRead += eh->_eftabLen*OFF_SIZE;
				fseeko(_in1, eh->_eftabLen*OFF_SIZE, SEEK_CUR);
#endif
			} else {
				_eftab.init(new TIndexOffU[eh->_eftabLen], eh->_eftabLen, true);
				if(_switchEndian) {
					for(TIndexOffU i = 0; i < eh->_eftabLen; i++)
						this->eftab()[i] = readU<TIndexOffU>(_in1, _switchEndian);
				} else {
					size_t r = MM_READ(_in1, (void *)this->eftab(), eh->_eftabLen*OFF_SIZE);
					if(r != (size_t)(eh->_eftabLen*OFF_SIZE)) {
						cerr << "Error reading _eftab[] array: " << r << ", " << (eh->_eftabLen*OFF_SIZE) << endl;
						throw 1;
					}
				}
			}
			for(TIndexOffU i = 0; i < eh->_eftabLen; i++) {
				if(i > 0 && this->eftab()[i] > 0) {
					assert_geq(this->eftab()[i], this->eftab()[i-1]);
				} else if(i > 0 && this->eftab()[i-1] == 0) {
					assert_eq(0, this->eftab()[i]);
				}
			}
		} else {
			assert(ftab() == NULL);
			assert(eftab() == NULL);
			// Skip ftab
			bytesRead += eh->_ftabLen*OFF_SIZE;
			fseeko(_in1, eh->_ftabLen*OFF_SIZE, SEEK_CUR);
			// Skip eftab
			bytesRead += eh->_eftabLen*OFF_SIZE;
			fseeko(_in1, eh->_eftabLen*OFF_SIZE, SEEK_CUR);
		}
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating fchr[], ftab[] or eftab[] arrays for the Bowtie index." << endl
		<< "Please try again on a computer with more memory." << endl;
		throw 1;
	}

	// Read reference sequence names from primary index file (or not,
	// if --refidx is specified)
	if(loadNames) {
		while(true) {
			char c = '\0';
			if(MM_READ(_in1, (void *)(&c), (size_t)1) != (size_t)1) break;
			bytesRead++;
			if(c == '\0') break;
			else if(c == '\n') {
				this->_refnames.push_back("");
			} else {
				if(this->_refnames.size() == 0) {
					this->_refnames.push_back("");
				}
				this->_refnames.back().push_back(c);
			}
		}
	}

	_offs.reset();
	if(loadSASamp) {
		bytesRead = 4; // reset for secondary index file (already read 1-sentinel)

		shmemLeader = true;
		if(_verbose || startVerbose) {
			cerr << "Reading offs (" << offsLenSampled << std::setw(2) << OFF_SIZE*8 <<"-bit words): ";
			logTime(cerr);
		}

		if(!_useMm) {
			if(!useShmem_) {
				// Allocate offs_
				try {
					_offs.init(new TIndexOffU[offsLenSampled], offsLenSampled, true);
				} catch(bad_alloc& e) {
					cerr << "Out of memory allocating the offs[] array  for the Bowtie index." << endl
					<< "Please try again on a computer with more memory." << endl;
					throw 1;
				}
			} else {
				TIndexOffU *tmp = NULL;
				shmemLeader = ALLOC_SHARED_U(
					(_in2Str + "[offs]"), offsLenSampled*OFF_SIZE, &tmp,
					"offs", (_verbose || startVerbose));
				_offs.init((TIndexOffU*)tmp, offsLenSampled, false);
			}
		}

		if(_overrideOffRate < 32) {
			if(shmemLeader) {
				// Allocate offs (big allocation)
				if(_switchEndian || offRateDiff > 0) {
					assert(!_useMm);
					const TIndexOffU blockMaxSz = (2 * 1024 * 1024); // 2 MB block size
					const TIndexOffU blockMaxSzU = (blockMaxSz >> (OFF_SIZE/4 + 1)); // # U32s per block
					char *buf;
					try {
						buf = new char[blockMaxSz];
					} catch(std::bad_alloc& e) {
						cerr << "Error: Out of memory allocating part of _offs array: '" << e.what() << "'" << endl;
						throw e;
					}
					for(TIndexOffU i = 0; i < offsLen; i += blockMaxSzU) {
						TIndexOffU block = min<TIndexOffU>(blockMaxSzU, offsLen - i);
						size_t r = MM_READ(_in2, (void *)buf, block << (OFF_SIZE/4 + 1));
						if(r != (size_t)(block << (OFF_SIZE/4 + 1))) {
							cerr << "Error reading block of _offs[] array: " << r << ", " << (block << (OFF_SIZE/4 + 1)) << endl;
							throw 1;
						}
						TIndexOffU idx = i >> offRateDiff;
						for(TIndexOffU j = 0; j < block; j += (1 << offRateDiff)) {
							assert_lt(idx, offsLenSampled);
							this->offs()[idx] = ((TIndexOffU*)buf)[j];
							if(_switchEndian) {
								this->offs()[idx] = endianSwapU(this->offs()[idx]);
							}
							idx++;
						}
					}
					delete[] buf;
				} else {
					if(_useMm) {
#ifdef BOWTIE_MM
						_offs.init((TIndexOffU*)(mmFile[1] + bytesRead), offsLen, false);
						bytesRead += offsSz;
						fseeko(_in2, offsSz, SEEK_CUR);
#endif
					} else {
						// Workaround for small-index mode where MM_READ may
						// not be able to handle read amounts greater than 2^32
						// bytes.
						uint64_t bytesLeft = offsSz;
						char *offs = (char *)this->offs();

						while(bytesLeft > 0) {
							size_t r = MM_READ(_in2, (void*)offs, bytesLeft);
							if(MM_IS_IO_ERR(_in2,r,bytesLeft)) {
								cerr << "Error reading block of _offs[] array: "
								     << r << ", " << bytesLeft << gLastIOErrMsg << endl;
								throw 1;
							} else if (r == 0) {
								cerr << "Error reading block of _offs[] array: no more data" << endl;
								throw 1;
							}
							offs += r;
							bytesLeft -= r;
						}
					}
				}
#ifdef BOWTIE_SHARED_MEM
				if(useShmem_) NOTIFY_SHARED(offs(), offsLenSampled*OFF_SIZE);
#endif
			} else {
				// Not the shmem leader
				fseeko(_in2, offsLenSampled*OFF_SIZE, SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM
				if(useShmem_) WAIT_SHARED(offs(), offsLenSampled*OFF_SIZE);
#endif
			}
		}
	}

	this->postReadInit(*eh); // Initialize fields of Ebwt not read from file
	if(_verbose || startVerbose) print(cerr, *eh);

	// The fact that _ebwt and friends actually point to something
	// (other than NULL) now signals to other member functions that the
	// Ebwt is loaded into memory.

done: // Exit hatch for both justHeader and !justHeader

	// Be kind
	if(deleteEh) delete eh;
	if(_in1 != NULL) {
		rewind(_in1);
	}
	if(_in2 != NULL) {
		rewind(_in2);
	}
}

/**
 * Read reference names from an input stream 'in' for an Ebwt primary
 * file and store them in 'refnames'.
 */
void
readEbwtRefnames(FILE* fin, EList<string>& refnames) {
	// _in1 must already be open with the get cursor at the
	// beginning and no error flags set.
	assert(fin != NULL);
	assert_eq(ftello(fin), 0);

	// Read endianness hints from both streams
	bool _switchEndian = false;
	uint32_t one = readU<uint32_t>(fin, _switchEndian); // 1st word of primary stream
	if(one != 1) {
		assert_eq((1u<<24), one);
		_switchEndian = true;
	}

	// Reads header entries one by one from primary stream
	TIndexOffU len          = readU<TIndexOffU>(fin, _switchEndian);
	int32_t  lineRate     = readI<int32_t>(fin, _switchEndian);
	/*int32_t  linesPerSide =*/ readI<int32_t>(fin, _switchEndian);
	int32_t  offRate      = readI<int32_t>(fin, _switchEndian);
	int32_t  ftabChars    = readI<int32_t>(fin, _switchEndian);
	// BTL: chunkRate is now deprecated
	int32_t flags = readI<int32_t>(fin, _switchEndian);
	bool color = false;
	bool entireReverse = false;
	if(flags < 0) {
		color = (((-flags) & EBWT_COLOR) != 0);
		entireReverse = (((-flags) & EBWT_ENTIRE_REV) != 0);
	}

	// Create a new EbwtParams from the entries read from primary stream
	EbwtParams eh(len, lineRate, offRate, ftabChars, color, entireReverse);

	TIndexOffU nPat = readI<TIndexOffU>(fin, _switchEndian); // nPat
	fseeko(fin, nPat*OFF_SIZE, SEEK_CUR);

	// Skip rstarts
	TIndexOffU nFrag = readU<TIndexOffU>(fin, _switchEndian);
	fseeko(fin, nFrag*OFF_SIZE*3, SEEK_CUR);

	// Skip ebwt
	fseeko(fin, eh._ebwtTotLen, SEEK_CUR);

	// Skip zOff from primary stream
	readU<TIndexOffU>(fin, _switchEndian);

	// Skip fchr
	fseeko(fin, 5 * OFF_SIZE, SEEK_CUR);

	// Skip ftab
	fseeko(fin, eh._ftabLen*OFF_SIZE, SEEK_CUR);

	// Skip eftab
	fseeko(fin, eh._eftabLen*OFF_SIZE, SEEK_CUR);

	// Read reference sequence names from primary index file
	while(true) {
		char c = '\0';
		int read_value = 0;
		read_value = fgetc(fin);
		if(read_value == EOF) break;
		c = read_value;
		if(c == '\0') break;
		else if(c == '\n') {
			refnames.push_back("");
		} else {
			if(refnames.size() == 0) {
				refnames.push_back("");
			}
			refnames.back().push_back(c);
		}
	}
	if(refnames.back().empty()) {
		refnames.pop_back();
	}

	// Be kind
	fseeko(fin, 0, SEEK_SET);
	assert(ferror(fin) == 0);
}

/**
 * Read reference names from the index with basename 'in' and store
 * them in 'refnames'.
 */
void
readEbwtRefnames(const string& instr, EList<string>& refnames) {
    FILE* fin;
	// Initialize our primary and secondary input-stream fields
    fin = fopen((instr + ".1." + gEbwt_ext).c_str(),"rb");
	if(fin == NULL) {
		throw EbwtFileOpenException("Cannot open file " + instr);
	}
	assert_eq(ftello(fin), 0);
	readEbwtRefnames(fin, refnames);
    fclose(fin);
}

/**
 * Read just enough of the Ebwt's header to get its flags
 */
int32_t Ebwt::readFlags(const string& instr) {
	ifstream in;
	// Initialize our primary and secondary input-stream fields
	in.open((instr + ".1." + gEbwt_ext).c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		throw EbwtFileOpenException("Cannot open file " + instr);
	}
	assert(in.is_open());
	assert(in.good());
	bool _switchEndian = false;
	uint32_t one = readU<uint32_t>(in, _switchEndian); // 1st word of primary stream
	if(one != 1) {
		assert_eq((1u<<24), one);
		assert_eq(1, endianSwapU32(one));
		_switchEndian = true;
	}
	readU<TIndexOffU>(in, _switchEndian);
	readI<int32_t>(in, _switchEndian);
	readI<int32_t>(in, _switchEndian);
	readI<int32_t>(in, _switchEndian);
	readI<int32_t>(in, _switchEndian);
	int32_t flags = readI<int32_t>(in, _switchEndian);
	return flags;
}

/**
 * Read just enough of the Ebwt's header to determine whether it's
 * colorspace.
 */
bool
readEbwtColor(const string& instr) {
	int32_t flags = Ebwt::readFlags(instr);
	if(flags < 0 && (((-flags) & EBWT_COLOR) != 0)) {
		return true;
	} else {
		return false;
	}
}

/**
 * Read just enough of the Ebwt's header to determine whether it's
 * entirely reversed.
 */
bool
readEntireReverse(const string& instr) {
	int32_t flags = Ebwt::readFlags(instr);
	if(flags < 0 && (((-flags) & EBWT_ENTIRE_REV) != 0)) {
		return true;
	} else {
		return false;
	}
}

/**
 * Write an extended Burrows-Wheeler transform to a pair of output
 * streams.
 *
 * @param out1 output stream to primary file
 * @param out2 output stream to secondary file
 * @param be   write in big endian?
 */
void Ebwt::writeFromMemory(bool justHeader,
                           ostream& out1,
                           ostream& out2) const
{
	const EbwtParams& eh = this->_eh;
	assert(eh.repOk());
	assert(out1.good());
	assert(out2.good());

	// When building an Ebwt, these header parameters are known
	// "up-front", i.e., they can be written to disk immediately,
	// before we join() or buildToDisk()
	writeI<int32_t>(out1, 1, _switchEndian); // endian hint for priamry stream
	writeI<int32_t>(out2, 1, _switchEndian); // endian hint for secondary stream
	writeU<TIndexOffU>(out1, eh._len,          _switchEndian); // length of string (and bwt and suffix array)
	writeI<int32_t>(out1, eh._lineRate,     _switchEndian); // 2^lineRate = size in bytes of 1 line
	writeI<int32_t>(out1, 2,                _switchEndian); // not used
	writeI<int32_t>(out1, eh._offRate,      _switchEndian); // every 2^offRate chars is "marked"
	writeI<int32_t>(out1, eh._ftabChars,    _switchEndian); // number of 2-bit chars used to address ftab
	int32_t flags = 1;
	if(eh._color) flags |= EBWT_COLOR;
	if(eh._entireReverse) flags |= EBWT_ENTIRE_REV;
	writeI<int32_t>(out1, -flags, _switchEndian); // BTL: chunkRate is now deprecated

	if(!justHeader) {
		assert(rstarts() != NULL);
		assert(offs() != NULL);
		assert(ftab() != NULL);
		assert(eftab() != NULL);
		assert(isInMemory());
		// These Ebwt parameters are known after the inputs strings have
		// been joined() but before they have been built().  These can
		// written to the disk next and then discarded from memory.
		writeU<TIndexOffU>(out1, this->_nPat,      _switchEndian);
		for(TIndexOffU i = 0; i < this->_nPat; i++)
			writeU<TIndexOffU>(out1, this->plen()[i], _switchEndian);
		assert_geq(this->_nFrag, this->_nPat);
		writeU<TIndexOffU>(out1, this->_nFrag, _switchEndian);
		for(TIndexOffU i = 0; i < this->_nFrag*3; i++)
			writeU<TIndexOffU>(out1, this->rstarts()[i], _switchEndian);

		// These Ebwt parameters are discovered only as the Ebwt is being
		// built (in buildToDisk()).  Of these, only 'offs' and 'ebwt' are
		// terribly large.  'ebwt' is written to the primary file and then
		// discarded from memory as it is built; 'offs' is similarly
		// written to the secondary file and discarded.
		out1.write((const char *)this->ebwt(), eh._ebwtTotLen);
		writeU<TIndexOffU>(out1, this->zOff(), _switchEndian);
		TIndexOffU offsLen = eh._offsLen;
		for(TIndexOffU i = 0; i < offsLen; i++)
			writeU<TIndexOffU>(out2, this->offs()[i], _switchEndian);

		// 'fchr', 'ftab' and 'eftab' are not fully determined until the
		// loop is finished, so they are written to the primary file after
		// all of 'ebwt' has already been written and only then discarded
		// from memory.
		for(int i = 0; i < 5; i++)
			writeU<TIndexOffU>(out1, this->fchr()[i], _switchEndian);
		for(TIndexOffU i = 0; i < eh._ftabLen; i++)
			writeU<TIndexOffU>(out1, this->ftab()[i], _switchEndian);
		for(TIndexOffU i = 0; i < eh._eftabLen; i++)
			writeU<TIndexOffU>(out1, this->eftab()[i], _switchEndian);
	}
}

/**
 * Given a pair of strings representing output filenames, and assuming
 * this Ebwt object is currently in memory, write out this Ebwt to the
 * specified files.
 *
 * If sanity-checking is enabled, then once the streams have been
 * fully written and closed, we reopen them and read them into a
 * (hopefully) exact copy of this Ebwt.  We then assert that the
 * current Ebwt and the copy match in all of their fields.
 */
void Ebwt::writeFromMemory(bool justHeader,
                           const string& out1,
                           const string& out2) const
{
	ASSERT_ONLY(const EbwtParams& eh = this->_eh);
	assert(isInMemory());
	assert(eh.repOk());

	ofstream fout1(out1.c_str(), ios::binary);
	ofstream fout2(out2.c_str(), ios::binary);
	writeFromMemory(justHeader, fout1, fout2);
	fout1.close();
	fout2.close();

	// Read the file back in and assert that all components match
	if(_sanity) {
#if 0
		if(_verbose)
			cout << "Re-reading \"" << out1 << "\"/\"" << out2 << "\" for sanity check" << endl;
		Ebwt copy(out1, out2, _verbose, _sanity);
		assert(!isInMemory());
		copy.loadIntoMemory(eh._color ? 1 : 0, true, false, false);
		assert(isInMemory());
	    assert_eq(eh._lineRate,     copy.eh()._lineRate);
	    assert_eq(eh._offRate,      copy.eh()._offRate);
	    assert_eq(eh._ftabChars,    copy.eh()._ftabChars);
	    assert_eq(eh._len,          copy.eh()._len);
	    assert_eq(_zOff,             copy.zOff());
	    assert_eq(_zEbwtBpOff,       copy.zEbwtBpOff());
	    assert_eq(_zEbwtByteOff,     copy.zEbwtByteOff());
		assert_eq(_nPat,             copy.nPat());
		for(TIndexOffU i = 0; i < _nPat; i++)
			assert_eq(this->_plen[i], copy.plen()[i]);
		assert_eq(this->_nFrag, copy.nFrag());
		for(TIndexOffU i = 0; i < this->nFrag*3; i++) {
			assert_eq(this->_rstarts[i], copy.rstarts()[i]);
		}
		for(int i = 0; i < 5; i++)
			assert_eq(this->_fchr[i], copy.fchr()[i]);
		for(TIndexOffU i = 0; i < eh._ftabLen; i++)
			assert_eq(this->ftab()[i], copy.ftab()[i]);
		for(TIndexOffU i = 0; i < eh._eftabLen; i++)
			assert_eq(this->eftab()[i], copy.eftab()[i]);
		for(TIndexOffU i = 0; i < eh._offsLen; i++)
			assert_eq(this->_offs[i], copy.offs()[i]);
		for(TIndexOffU i = 0; i < eh._ebwtTotLen; i++)
			assert_eq(this->ebwt()[i], copy.ebwt()[i]);
		copy.sanityCheckAll();
		if(_verbose)
			cout << "Read-in check passed for \"" << out1 << "\"/\"" << out2 << "\"" << endl;
#endif
	}
}

/**
 * Write the rstarts array given the szs array for the reference.
 */
void Ebwt::szsToDisk(const EList<RefRecord>& szs, ostream& os, int reverse) {
	TIndexOffU seq = 0;
	TIndexOffU off = 0;
	TIndexOffU totlen = 0;
	for(unsigned int i = 0; i < szs.size(); i++) {
		if(szs[i].first) off = 0;
		off += szs[i].off;
		if(szs[i].first) seq++;
		if(szs[i].len == 0) continue;

		TIndexOffU seqm1 = seq-1;
		assert_lt(seqm1, _nPat);
		TIndexOffU fwoff = off;
		if(reverse == REF_READ_REVERSE) {
			// Invert pattern idxs
			seqm1 = _nPat - seqm1 - 1;
			// Invert pattern idxs
			assert_leq(off + szs[i].len, plen()[seqm1]);
			fwoff = plen()[seqm1] - (off + szs[i].len);
		}
		writeU<TIndexOffU>(os, totlen, _switchEndian); // offset from beginning of joined string
		writeU<TIndexOffU>(os, seqm1,  _switchEndian); // sequence id
		writeU<TIndexOffU>(os, fwoff,  _switchEndian); // offset into sequence
		totlen += szs[i].len;
		off += szs[i].len;
	}
}
