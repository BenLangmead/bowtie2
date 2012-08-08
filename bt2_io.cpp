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

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "bt2_idx.h"

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
	bool switchEndian; // dummy; caller doesn't care
#ifdef BOWTIE_MM
	char *mmFile[] = { NULL, NULL };
#endif
	if(_in1Str.length() > 0) {
		if(_verbose || startVerbose) {
			cerr << "  About to open input files: ";
			logTime(cerr);
		}
#ifdef BOWTIE_MM
		// Initialize our primary and secondary input-stream fields
		if(_in1 != -1) close(_in1);
		if(_verbose || startVerbose) {
			cerr << "Opening \"" << _in1Str << "\"" << endl;
		}
		if((_in1 = open(_in1Str.c_str(), O_RDONLY)) < 0) {
			cerr << "Could not open index file " << _in1Str << endl;
		}
		if(loadSASamp) {
			if(_in2 != -1) close(_in2);
			if(_verbose || startVerbose) {
				cerr << "Opening \"" << _in2Str << "\"" << endl;
			}
			if((_in2 = open(_in2Str.c_str(), O_RDONLY)) < 0) {
				cerr << "Could not open index file " << _in2Str << endl;
			}
		}
#else
		// Initialize our primary and secondary input-stream fields
		if(_in1 != NULL) fclose(_in1);
		if(_verbose || startVerbose) cerr << "Opening \"" << _in1Str << "\"" << endl;
		if((_in1 = fopen(_in1Str.c_str(), "rb")) == NULL) {
			cerr << "Could not open index file " << _in1Str << endl;
		}
		if(loadSASamp) {
			if(_in2 != NULL) fclose(_in2);
			if(_verbose || startVerbose) cerr << "Opening \"" << _in2Str << "\"" << endl;
			if((_in2 = fopen(_in2Str.c_str(), "rb")) == NULL) {
				cerr << "Could not open index file " << _in2Str << endl;
			}
		}
#endif
		if(_verbose || startVerbose) {
			cerr << "  Finished opening input files: ";
			logTime(cerr);
		}
		
#ifdef BOWTIE_MM
		if(_useMm /*&& !justHeader*/) {
			const char *names[] = {_in1Str.c_str(), _in2Str.c_str()};
			int fds[] = { _in1, _in2 };
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
	size_t bytesRead = 0;
	switchEndian = false;
	uint32_t one = readU32(_in1, switchEndian); // 1st word of primary stream
	bytesRead += 4;
	if(loadSASamp) {
#ifndef NDEBUG
		assert_eq(one, readU32(_in2, switchEndian)); // should match!
#else
		readU32(_in2, switchEndian);
#endif
	}
	if(one != 1) {
		assert_eq((1u<<24), one);
		assert_eq(1, endianSwapU32(one));
		switchEndian = true;
	}
	
	// Can't switch endianness and use memory-mapped files; in order to
	// support this, someone has to modify the file to switch
	// endiannesses appropriately, and we can't do this inside Bowtie
	// or we might be setting up a race condition with other processes.
	if(switchEndian && _useMm) {
		cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
		throw 1;
	}
	
	// Reads header entries one by one from primary stream
	uint32_t len          = readU32(_in1, switchEndian);
	bytesRead += 4;
	int32_t  lineRate     = readI32(_in1, switchEndian);
	bytesRead += 4;
	/*int32_t  linesPerSide =*/ readI32(_in1, switchEndian);
	bytesRead += 4;
	int32_t  offRate      = readI32(_in1, switchEndian);
	bytesRead += 4;
	// TODO: add isaRate to the actual file format (right now, the
	// user has to tell us whether there's an ISA sample and what the
	// sampling rate is.
	int32_t  ftabChars    = readI32(_in1, switchEndian);
	bytesRead += 4;
	// chunkRate was deprecated in an earlier version of Bowtie; now
	// we use it to hold flags.
	int32_t flags = readI32(_in1, switchEndian);
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
	uint32_t offsLen = eh->_offsLen;
	uint32_t offRateDiff = 0;
	uint32_t offsLenSampled = offsLen;
	if(_overrideOffRate > offRate) {
		offRateDiff = _overrideOffRate - offRate;
	}
	if(offRateDiff > 0) {
		offsLenSampled >>= offRateDiff;
		if((offsLen & ~(0xffffffff << offRateDiff)) != 0) {
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
	this->_nPat = readI32(_in1, switchEndian);
	bytesRead += 4;
	_plen.reset();
	// Read plen from primary stream
	if(_useMm) {
#ifdef BOWTIE_MM
		_plen.init((uint32_t*)(mmFile[0] + bytesRead), _nPat, false);
		bytesRead += _nPat*4;
		MM_SEEK(_in1, _nPat*4, SEEK_CUR);
#endif
	} else {
		try {
			if(_verbose || startVerbose) {
				cerr << "Reading plen (" << this->_nPat << "): ";
				logTime(cerr);
			}
			_plen.init(new uint32_t[_nPat], _nPat, true);
			if(switchEndian) {
				for(uint32_t i = 0; i < this->_nPat; i++) {
					plen()[i] = readU32(_in1, switchEndian);
				}
			} else {
				MM_READ_RET r = MM_READ(_in1, (void*)(plen()), _nPat*4);
				if(r != (MM_READ_RET)(_nPat*4)) {
					cerr << "Error reading _plen[] array: " << r << ", " << _nPat*4 << endl;
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
	
	this->_nFrag = readU32(_in1, switchEndian);
	bytesRead += 4;
	if(_verbose || startVerbose) {
		cerr << "Reading rstarts (" << this->_nFrag*3 << "): ";
		logTime(cerr);
	}
	assert_geq(this->_nFrag, this->_nPat);
	_rstarts.reset();
	if(loadRstarts) {
		if(_useMm) {
#ifdef BOWTIE_MM
			_rstarts.init((uint32_t*)(mmFile[0] + bytesRead), _nFrag*3, false);
			bytesRead += this->_nFrag*4*3;
			MM_SEEK(_in1, this->_nFrag*4*3, SEEK_CUR);
#endif
		} else {
			_rstarts.init(new uint32_t[_nFrag*3], _nFrag*3, true);
			if(switchEndian) {
				for(uint32_t i = 0; i < this->_nFrag*3; i += 3) {
					// fragment starting position in joined reference
					// string, text id, and fragment offset within text
					this->rstarts()[i]   = readU32(_in1, switchEndian);
					this->rstarts()[i+1] = readU32(_in1, switchEndian);
					this->rstarts()[i+2] = readU32(_in1, switchEndian);
				}
			} else {
				MM_READ_RET r = MM_READ(_in1, (void *)rstarts(), this->_nFrag*4*3);
				if(r != (MM_READ_RET)(this->_nFrag*4*3)) {
					cerr << "Error reading _rstarts[] array: " << r << ", " << (this->_nFrag*4*3) << endl;
					throw 1;
				}
			}
		}
	} else {
		// Skip em
		assert(rstarts() == NULL);
		bytesRead += this->_nFrag*4*3;
		MM_SEEK(_in1, this->_nFrag*4*3, SEEK_CUR);
	}
	
	_ebwt.reset();
	if(_useMm) {
#ifdef BOWTIE_MM
		_ebwt.init((uint8_t*)(mmFile[0] + bytesRead), eh->_ebwtTotLen, false);
		bytesRead += eh->_ebwtTotLen;
		MM_SEEK(_in1, eh->_ebwtTotLen, SEEK_CUR);
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
			MM_READ_RET r = MM_READ(_in1, (void *)this->ebwt(), eh->_ebwtTotLen);
			if(r != (MM_READ_RET)eh->_ebwtTotLen) {
				cerr << "Error reading _ebwt[] array: " << r << ", " << (eh->_ebwtTotLen) << endl;
				throw 1;
			}
			if(switchEndian) {
				uint8_t *side = this->ebwt();
				for(size_t i = 0; i < eh->_numSides; i++) {
					uint32_t *cums = reinterpret_cast<uint32_t*>(side + eh->_sideSz - 8);
					cums[0] = endianSwapU32(cums[0]);
					cums[1] = endianSwapU32(cums[1]);
					side += this->_eh._sideSz;
				}
			}
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) NOTIFY_SHARED(ebwt(), eh->_ebwtTotLen);
#endif
		} else {
			// Seek past the data and wait until master is finished
			MM_SEEK(_in1, eh->_ebwtTotLen, SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) WAIT_SHARED(ebwt(), eh->_ebwtTotLen);
#endif
		}
	}
	
	// Read zOff from primary stream
	_zOff = readU32(_in1, switchEndian);
	bytesRead += 4;
	assert_lt(_zOff, len);
	
	try {
		// Read fchr from primary stream
		if(_verbose || startVerbose) cerr << "Reading fchr (5)" << endl;
		_fchr.reset();
		if(_useMm) {
#ifdef BOWTIE_MM
			_fchr.init((uint32_t*)(mmFile[0] + bytesRead), 5, false);
			bytesRead += 5*4;
			MM_SEEK(_in1, 5*4, SEEK_CUR);
#endif
		} else {
			_fchr.init(new uint32_t[5], 5, true);
			for(int i = 0; i < 5; i++) {
				this->fchr()[i] = readU32(_in1, switchEndian);
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
				_ftab.init((uint32_t*)(mmFile[0] + bytesRead), eh->_ftabLen, false);
				bytesRead += eh->_ftabLen*4;
				MM_SEEK(_in1, eh->_ftabLen*4, SEEK_CUR);
#endif
			} else {
				_ftab.init(new uint32_t[eh->_ftabLen], eh->_ftabLen, true);
				if(switchEndian) {
					for(uint32_t i = 0; i < eh->_ftabLen; i++)
						this->ftab()[i] = readU32(_in1, switchEndian);
				} else {
					MM_READ_RET r = MM_READ(_in1, (void *)ftab(), eh->_ftabLen*4);
					if(r != (MM_READ_RET)(eh->_ftabLen*4)) {
						cerr << "Error reading _ftab[] array: " << r << ", " << (eh->_ftabLen*4) << endl;
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
				_eftab.init((uint32_t*)(mmFile[0] + bytesRead), eh->_eftabLen, false);
				bytesRead += eh->_eftabLen*4;
				MM_SEEK(_in1, eh->_eftabLen*4, SEEK_CUR);
#endif
			} else {
				_eftab.init(new uint32_t[eh->_eftabLen], eh->_eftabLen, true);
				if(switchEndian) {
					for(uint32_t i = 0; i < eh->_eftabLen; i++)
						this->eftab()[i] = readU32(_in1, switchEndian);
				} else {
					MM_READ_RET r = MM_READ(_in1, (void *)this->eftab(), eh->_eftabLen*4);
					if(r != (MM_READ_RET)(eh->_eftabLen*4)) {
						cerr << "Error reading _eftab[] array: " << r << ", " << (eh->_eftabLen*4) << endl;
						throw 1;
					}
				}
			}
			for(uint32_t i = 0; i < eh->_eftabLen; i++) {
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
			bytesRead += eh->_ftabLen*4;
			MM_SEEK(_in1, eh->_ftabLen*4, SEEK_CUR);
			// Skip eftab
			bytesRead += eh->_eftabLen*4;
			MM_SEEK(_in1, eh->_eftabLen*4, SEEK_CUR);
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
			if(MM_READ(_in1, (void *)(&c), (size_t)1) != (MM_READ_RET)1) break;
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
			cerr << "Reading offs (" << offsLenSampled << " 32-bit words): ";
			logTime(cerr);
		}
		
		if(!_useMm) {
			if(!useShmem_) {
				// Allocate offs_
				try {
					_offs.init(new uint32_t[offsLenSampled], offsLenSampled, true);
				} catch(bad_alloc& e) {
					cerr << "Out of memory allocating the offs[] array  for the Bowtie index." << endl
					<< "Please try again on a computer with more memory." << endl;
					throw 1;
				}
			} else {
				uint32_t *tmp = NULL;
				shmemLeader = ALLOC_SHARED_U32(
					(_in2Str + "[offs]"), offsLenSampled*4, &tmp,
					"offs", (_verbose || startVerbose));
				_offs.init((uint32_t*)tmp, offsLenSampled, false);
			}
		}
		
		if(_overrideOffRate < 32) {
			if(shmemLeader) {
				// Allocate offs (big allocation)
				if(switchEndian || offRateDiff > 0) {
					assert(!_useMm);
					const uint32_t blockMaxSz = (2 * 1024 * 1024); // 2 MB block size
					const uint32_t blockMaxSzU32 = (blockMaxSz >> 2); // # U32s per block
					char *buf;
					try {
						buf = new char[blockMaxSz];
					} catch(std::bad_alloc& e) {
						cerr << "Error: Out of memory allocating part of _offs array: '" << e.what() << "'" << endl;
						throw e;
					}
					for(uint32_t i = 0; i < offsLen; i += blockMaxSzU32) {
						uint32_t block = min<uint32_t>(blockMaxSzU32, offsLen - i);
						MM_READ_RET r = MM_READ(_in2, (void *)buf, block << 2);
						if(r != (MM_READ_RET)(block << 2)) {
							cerr << "Error reading block of _offs[] array: " << r << ", " << (block << 2) << endl;
							throw 1;
						}
						uint32_t idx = i >> offRateDiff;
						for(uint32_t j = 0; j < block; j += (1 << offRateDiff)) {
							assert_lt(idx, offsLenSampled);
							this->offs()[idx] = ((uint32_t*)buf)[j];
							if(switchEndian) {
								this->offs()[idx] = endianSwapU32(this->offs()[idx]);
							}
							idx++;
						}
					}
					delete[] buf;
				} else {
					if(_useMm) {
#ifdef BOWTIE_MM
						_offs.init((uint32_t*)(mmFile[1] + bytesRead), offsLen, false);
						bytesRead += (offsLen << 2);
						MM_SEEK(_in2, (offsLen << 2), SEEK_CUR);
#endif
					} else {
						// If any of the high two bits are set
						if((offsLen & 0xc0000000) != 0) {
							if(sizeof(char *) <= 4) {
								cerr << "Sanity error: sizeof(char *) <= 4 but offsLen is " << hex << offsLen << endl;
								throw 1;
							}
							// offsLen << 2 overflows, so do it in four reads
							char *offs = (char *)this->offs();
							for(int i = 0; i < 4; i++) {
								MM_READ_RET r = MM_READ(_in2, (void*)offs, offsLen);
								if(r != (MM_READ_RET)(offsLen)) {
									cerr << "Error reading block of _offs[] array: " << r << ", " << offsLen << endl;
									throw 1;
								}
								offs += offsLen;
							}
						} else {
							// Do it all in one read
							MM_READ_RET r = MM_READ(_in2, (void*)this->offs(), offsLen << 2);
							if(r != (MM_READ_RET)(offsLen << 2)) {
								cerr << "Error reading _offs[] array: " << r << ", " << (offsLen << 2) << endl;
								throw 1;
							}
						}
					}
				}
#if 0
				{
					ASSERT_ONLY(Bitset offsSeen(len+1));
					for(uint32_t i = 0; i < offsLenSampled; i++) {
						assert(!offsSeen.test(this->offs()[i]));
						ASSERT_ONLY(offsSeen.set(this->offs()[i]));
						assert_leq(this->offs()[i], len);
					}
				}
#endif
#ifdef BOWTIE_SHARED_MEM				
				if(useShmem_) NOTIFY_SHARED(offs(), offsLenSampled*4);
#endif
			} else {
				// Not the shmem leader
				MM_SEEK(_in2, offsLenSampled*4, SEEK_CUR);
#ifdef BOWTIE_SHARED_MEM				
				if(useShmem_) WAIT_SHARED(offs(), offsLenSampled*4);
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
#ifdef BOWTIE_MM
	MM_SEEK(_in1, 0, SEEK_SET);
	MM_SEEK(_in2, 0, SEEK_SET);
#else
	rewind(_in1); rewind(_in2);
#endif
}

/**
 * Read reference names from an input stream 'in' for an Ebwt primary
 * file and store them in 'refnames'.
 */
void
readEbwtRefnames(istream& in, EList<string>& refnames) {
	// _in1 must already be open with the get cursor at the
	// beginning and no error flags set.
	assert(in.good());
	assert_eq((streamoff)in.tellg(), ios::beg);
	
	// Read endianness hints from both streams
	bool switchEndian = false;
	uint32_t one = readU32(in, switchEndian); // 1st word of primary stream
	if(one != 1) {
		assert_eq((1u<<24), one);
		switchEndian = true;
	}
	
	// Reads header entries one by one from primary stream
	uint32_t len          = readU32(in, switchEndian);
	int32_t  lineRate     = readI32(in, switchEndian);
	/*int32_t  linesPerSide =*/ readI32(in, switchEndian);
	int32_t  offRate      = readI32(in, switchEndian);
	int32_t  ftabChars    = readI32(in, switchEndian);
	// BTL: chunkRate is now deprecated
	int32_t flags = readI32(in, switchEndian);
	bool color = false;
	bool entireReverse = false;
	if(flags < 0) {
		color = (((-flags) & EBWT_COLOR) != 0);
		entireReverse = (((-flags) & EBWT_ENTIRE_REV) != 0);
	}
	
	// Create a new EbwtParams from the entries read from primary stream
	EbwtParams eh(len, lineRate, offRate, ftabChars, color, entireReverse);
	
	uint32_t nPat = readI32(in, switchEndian); // nPat
	in.seekg(nPat*4, ios_base::cur); // skip plen
	
	// Skip rstarts
	uint32_t nFrag = readU32(in, switchEndian);
	in.seekg(nFrag*4*3, ios_base::cur);
	
	// Skip ebwt
	in.seekg(eh._ebwtTotLen, ios_base::cur);
	
	// Skip zOff from primary stream
	readU32(in, switchEndian);
	
	// Skip fchr
	in.seekg(5 * 4, ios_base::cur);
	
	// Skip ftab
	in.seekg(eh._ftabLen*4, ios_base::cur);
	
	// Skip eftab
	in.seekg(eh._eftabLen*4, ios_base::cur);
	
	// Read reference sequence names from primary index file
	while(true) {
		char c = '\0';
		in.read(&c, 1);
		if(in.eof()) break;
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
	in.clear(); in.seekg(0, ios::beg);
	assert(in.good());
}

/**
 * Read reference names from the index with basename 'in' and store
 * them in 'refnames'.
 */
void
readEbwtRefnames(const string& instr, EList<string>& refnames) {
	ifstream in;
	// Initialize our primary and secondary input-stream fields
	in.open((instr + ".1.bt2").c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		throw EbwtFileOpenException("Cannot open file " + instr);
	}
	assert(in.is_open());
	assert(in.good());
	assert_eq((streamoff)in.tellg(), ios::beg);
	readEbwtRefnames(in, refnames);
}

/**
 * Read just enough of the Ebwt's header to get its flags
 */
int32_t Ebwt::readFlags(const string& instr) {
	ifstream in;
	// Initialize our primary and secondary input-stream fields
	in.open((instr + ".1.bt2").c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		throw EbwtFileOpenException("Cannot open file " + instr);
	}
	assert(in.is_open());
	assert(in.good());
	bool switchEndian = false;
	uint32_t one = readU32(in, switchEndian); // 1st word of primary stream
	if(one != 1) {
		assert_eq((1u<<24), one);
		assert_eq(1, endianSwapU32(one));
		switchEndian = true;
	}
	readU32(in, switchEndian);
	readI32(in, switchEndian);
	readI32(in, switchEndian);
	readI32(in, switchEndian);
	readI32(in, switchEndian);
	int32_t flags = readI32(in, switchEndian);
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
	uint32_t be = this->toBe();
	assert(out1.good());
	assert(out2.good());
	
	// When building an Ebwt, these header parameters are known
	// "up-front", i.e., they can be written to disk immediately,
	// before we join() or buildToDisk()
	writeI32(out1, 1, be); // endian hint for priamry stream
	writeI32(out2, 1, be); // endian hint for secondary stream
	writeU32(out1, eh._len,          be); // length of string (and bwt and suffix array)
	writeI32(out1, eh._lineRate,     be); // 2^lineRate = size in bytes of 1 line
	writeI32(out1, 2,                be); // not used
	writeI32(out1, eh._offRate,      be); // every 2^offRate chars is "marked"
	writeI32(out1, eh._ftabChars,    be); // number of 2-bit chars used to address ftab
	int32_t flags = 1;
	if(eh._color) flags |= EBWT_COLOR;
	if(eh._entireReverse) flags |= EBWT_ENTIRE_REV;
	writeI32(out1, -flags, be); // BTL: chunkRate is now deprecated
	
	if(!justHeader) {
		assert(rstarts() != NULL);
		assert(offs() != NULL);
		assert(ftab() != NULL);
		assert(eftab() != NULL);
		assert(isInMemory());
		// These Ebwt parameters are known after the inputs strings have
		// been joined() but before they have been built().  These can
		// written to the disk next and then discarded from memory.
		writeU32(out1, this->_nPat,      be);
		for(uint32_t i = 0; i < this->_nPat; i++)
			writeU32(out1, this->plen()[i], be);
		assert_geq(this->_nFrag, this->_nPat);
		writeU32(out1, this->_nFrag, be);
		for(uint32_t i = 0; i < this->_nFrag*3; i++)
			writeU32(out1, this->rstarts()[i], be);
		
		// These Ebwt parameters are discovered only as the Ebwt is being
		// built (in buildToDisk()).  Of these, only 'offs' and 'ebwt' are
		// terribly large.  'ebwt' is written to the primary file and then
		// discarded from memory as it is built; 'offs' is similarly
		// written to the secondary file and discarded.
		out1.write((const char *)this->ebwt(), eh._ebwtTotLen);
		writeU32(out1, this->zOff(), be);
		uint32_t offsLen = eh._offsLen;
		for(uint32_t i = 0; i < offsLen; i++)
			writeU32(out2, this->offs()[i], be);
		
		// 'fchr', 'ftab' and 'eftab' are not fully determined until the
		// loop is finished, so they are written to the primary file after
		// all of 'ebwt' has already been written and only then discarded
		// from memory.
		for(int i = 0; i < 5; i++)
			writeU32(out1, this->fchr()[i], be);
		for(uint32_t i = 0; i < eh._ftabLen; i++)
			writeU32(out1, this->ftab()[i], be);
		for(uint32_t i = 0; i < eh._eftabLen; i++)
			writeU32(out1, this->eftab()[i], be);
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
		for(uint32_t i = 0; i < _nPat; i++)
			assert_eq(this->_plen[i], copy.plen()[i]);
		assert_eq(this->_nFrag, copy.nFrag());
		for(uint32_t i = 0; i < this->nFrag*3; i++) {
			assert_eq(this->_rstarts[i], copy.rstarts()[i]);
		}
		for(uint32_t i = 0; i < 5; i++)
			assert_eq(this->_fchr[i], copy.fchr()[i]);
		for(uint32_t i = 0; i < eh._ftabLen; i++)
			assert_eq(this->ftab()[i], copy.ftab()[i]);
		for(uint32_t i = 0; i < eh._eftabLen; i++)
			assert_eq(this->eftab()[i], copy.eftab()[i]);
		for(uint32_t i = 0; i < eh._offsLen; i++)
			assert_eq(this->_offs[i], copy.offs()[i]);
		for(uint32_t i = 0; i < eh._ebwtTotLen; i++)
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
	size_t seq = 0;
	uint32_t off = 0;
	uint32_t totlen = 0;
	for(unsigned int i = 0; i < szs.size(); i++) {
		if(szs[i].len == 0) continue;
		if(szs[i].first) off = 0;
		off += szs[i].off;
		if(szs[i].first && szs[i].len > 0) seq++;
		size_t seqm1 = seq-1;
		assert_lt(seqm1, _nPat);
		size_t fwoff = off;
		if(reverse == REF_READ_REVERSE) {
			// Invert pattern idxs
			seqm1 = _nPat - seqm1 - 1;
			// Invert pattern idxs
			assert_leq(off + szs[i].len, plen()[seqm1]);
			fwoff = plen()[seqm1] - (off + szs[i].len);
		}
		writeU32(os, totlen, this->toBe()); // offset from beginning of joined string
		writeU32(os, (uint32_t)seqm1,  this->toBe()); // sequence id
		writeU32(os, (uint32_t)fwoff,  this->toBe()); // offset into sequence
		totlen += szs[i].len;
		off += szs[i].len;
	}
}
