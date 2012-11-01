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

#include <string.h>
#include "reference.h"
#include "mem_ids.h"

using namespace std;

/**
 * Load from .3.bt2/.4.bt2 Bowtie index files.
 */
BitPairReference::BitPairReference(
	const string& in,
	bool color,
	bool sanity,
	EList<string>* infiles,
	EList<SString<char> >* origs,
	bool infilesSeq,
	bool useMm,
	bool useShmem,
	bool mmSweep,
	bool verbose,
	bool startVerbose) :
	buf_(NULL),
	sanityBuf_(NULL),
	loaded_(true),
	sanity_(sanity),
	useMm_(useMm),
	useShmem_(useShmem),
	verbose_(verbose)
{
	string s3 = in + ".3.bt2";
	string s4 = in + ".4.bt2";
	
#ifdef BOWTIE_MM
	int f3, f4;
	if((f3 = open(s3.c_str(), O_RDONLY)) < 0) {
		cerr << "Could not open reference-string index file " << s3 << " for reading." << endl;
		cerr << "This is most likely because your index was built with an older version" << endl
		<< "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
		<< "index (or download one from the Bowtie website) and try again." << endl;
		loaded_ = false;
		return;
	}
	if((f4 = open(s4.c_str(), O_RDONLY)) < 0) {
		cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
		loaded_ = false;
		return;
	}
	char *mmFile = NULL;
	if(useMm_) {
		if(verbose_ || startVerbose) {
			cerr << "  Memory-mapping reference index file " << s4 << ": ";
			logTime(cerr);
		}
		struct stat sbuf;
		if (stat(s4.c_str(), &sbuf) == -1) {
			perror("stat");
			cerr << "Error: Could not stat index file " << s4.c_str() << " prior to memory-mapping" << endl;
			throw 1;
		}
		mmFile = (char*)mmap((void *)0, (size_t)sbuf.st_size,
							 PROT_READ, MAP_SHARED, f4, 0);
		if(mmFile == (void *)(-1) || mmFile == NULL) {
			perror("mmap");
			cerr << "Error: Could not memory-map the index file " << s4.c_str() << endl;
			throw 1;
		}
		if(mmSweep) {
			int sum = 0;
			for(off_t i = 0; i < sbuf.st_size; i += 1024) {
				sum += (int) mmFile[i];
			}
			if(startVerbose) {
				cerr << "  Swept the memory-mapped ref index file; checksum: " << sum << ": ";
				logTime(cerr);
			}
		}
	}
#else
	FILE *f3, *f4;
	if((f3 = fopen(s3.c_str(), "rb")) == NULL) {
		cerr << "Could not open reference-string index file " << s3 << " for reading." << endl;
		cerr << "This is most likely because your index was built with an older version" << endl
		<< "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
		<< "index (or download one from the Bowtie website) and try again." << endl;
		loaded_ = false;
		return;
	}
	if((f4 = fopen(s4.c_str(), "rb"))  == NULL) {
		cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
		loaded_ = false;
		return;
	}
#endif
	
	// Read endianness sentinel, set 'swap'
	uint32_t one;
	bool swap = false;
	one = readU32(f3, swap);
	if(one != 1) {
		if(useMm_) {
			cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
			throw 1;
		}
		assert_eq(0x1000000, one);
		swap = true; // have to endian swap U32s
	}
	
	// Read # records
	uint32_t sz;
	sz = readU32(f3, swap);
	if(sz == 0) {
		cerr << "Error: number of reference records is 0 in " << s3 << endl;
		throw 1;
	}
	
	// Read records
	nrefs_ = 0;
	
	// Cumulative count of all unambiguous characters on a per-
	// stretch 8-bit alignment (i.e. count of bytes we need to
	// allocate in buf_)
	uint32_t cumsz = 0;
	uint32_t cumlen = 0;
	// For each unambiguous stretch...
	for(uint32_t i = 0; i < sz; i++) {
		recs_.push_back(RefRecord(f3, swap));
		if(recs_.back().first) {
			// This is the first record for this reference sequence (and the
			// last record for the one before)
			refRecOffs_.push_back((uint32_t)recs_.size()-1);
			// refOffs_ links each reference sequence with the total number of
			// unambiguous characters preceding it in the pasted reference
			refOffs_.push_back(cumsz);
			if(nrefs_ > 0) {
				// refLens_ links each reference sequence with the total number
				// of ambiguous and unambiguous characters in it.
				refLens_.push_back(cumlen);
			}
			cumlen = 0;
			nrefs_++;
		} else if(i == 0) {
			cerr << "First record in reference index file was not marked as "
			     << "'first'" << endl;
			throw 1;
		}
		cumsz += recs_.back().len;
		cumlen += recs_.back().off;
		cumlen += recs_.back().len;
	}
	if(verbose_ || startVerbose) {
		cerr << "Read " << nrefs_ << " reference strings from "
		     << sz << " records: ";
		logTime(cerr);
	}
	// Store a cap entry for the end of the last reference seq
	refRecOffs_.push_back((uint32_t)recs_.size());
	refOffs_.push_back(cumsz);
	refLens_.push_back(cumlen);
	bufSz_ = cumsz;
	assert_eq(nrefs_, refLens_.size());
	assert_eq(sz, recs_.size());
	MM_FILE_CLOSE(f3); // done with .3.bt2 file
	// Round cumsz up to nearest byte boundary
	if((cumsz & 3) != 0) {
		cumsz += (4 - (cumsz & 3));
	}
	bufAllocSz_ = cumsz >> 2;
	assert_eq(0, cumsz & 3); // should be rounded up to nearest 4
	if(useMm_) {
#ifdef BOWTIE_MM
		buf_ = (uint8_t*)mmFile;
		if(sanity_) {
			FILE *ftmp = fopen(s4.c_str(), "rb");
			sanityBuf_ = new uint8_t[cumsz >> 2];
			size_t ret = fread(sanityBuf_, 1, cumsz >> 2, ftmp);
			if(ret != (cumsz >> 2)) {
				cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4 << endl;
				throw 1;
			}
			fclose(ftmp);
			for(size_t i = 0; i < (cumsz >> 2); i++) {
				assert_eq(sanityBuf_[i], buf_[i]);
			}
		}
#else
		cerr << "Shouldn't be at " << __FILE__ << ":" << __LINE__ << " without BOWTIE_MM defined" << endl;
		throw 1;
#endif
	} else {
		bool shmemLeader = true;
		if(!useShmem_) {
			// Allocate a buffer to hold the reference string
			try {
				buf_ = new uint8_t[cumsz >> 2];
				if(buf_ == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				cerr << "Error: Ran out of memory allocating space for the bitpacked reference.  Please" << endl
				<< "re-run on a computer with more memory." << endl;
				throw 1;
			}
		} else {
			shmemLeader = ALLOC_SHARED_U8(
										  (s4 + "[ref]"), (cumsz >> 2), &buf_,
										  "ref", (verbose_ || startVerbose));
		}
		if(shmemLeader) {
			// Open the bitpair-encoded reference file
			FILE *f4 = fopen(s4.c_str(), "rb");
			if(f4 == NULL) {
				cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
				cerr << "This is most likely because your index was built with an older version" << endl
				<< "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
				<< "index (or download one from the Bowtie website) and try again." << endl;
				loaded_ = false;
				return;
			}
			// Read the whole thing in
			size_t ret = fread(buf_, 1, cumsz >> 2, f4);
			// Didn't read all of it?
			if(ret != (cumsz >> 2)) {
				cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4 << endl;
				throw 1;
			}
			// Make sure there's no more
			char c;
			ret = fread(&c, 1, 1, f4);
			assert_eq(0, ret); // should have failed
			fclose(f4);
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) NOTIFY_SHARED(buf_, (cumsz >> 2));
#endif
		} else {
#ifdef BOWTIE_SHARED_MEM
			if(useShmem_) WAIT_SHARED(buf_, (cumsz >> 2));
#endif
		}
	}
	
	// Populate byteToU32_
	bool big = currentlyBigEndian();
	for(int i = 0; i < 256; i++) {
		uint32_t word = 0;
		if(big) {
			word |= ((i >> 0) & 3) << 24;
			word |= ((i >> 2) & 3) << 16;
			word |= ((i >> 4) & 3) << 8;
			word |= ((i >> 6) & 3) << 0;
		} else {
			word |= ((i >> 0) & 3) << 0;
			word |= ((i >> 2) & 3) << 8;
			word |= ((i >> 4) & 3) << 16;
			word |= ((i >> 6) & 3) << 24;
		}
		byteToU32_[i] = word;
	}
	
#ifndef NDEBUG
	if(sanity_) {
		// Compare the sequence we just read from the compact index
		// file to the true reference sequence.
		EList<SString<char> > *os; // for holding references
		EList<SString<char> > osv(DEBUG_CAT); // for holding ref seqs
		EList<SString<char> > osn(DEBUG_CAT); // for holding ref names
		EList<size_t> osvLen(DEBUG_CAT); // for holding ref seq lens
		EList<size_t> osnLen(DEBUG_CAT); // for holding ref name lens
		SStringExpandable<uint32_t> tmp_destU32_;
		if(infiles != NULL) {
			if(infilesSeq) {
				for(size_t i = 0; i < infiles->size(); i++) {
					// Remove initial backslash; that's almost
					// certainly being used to protect the first
					// character of the sequence from getopts (e.g.,
					// when the first char is -)
					if((*infiles)[i].at(0) == '\\') {
						(*infiles)[i].erase(0, 1);
					}
					osv.push_back(SString<char>((*infiles)[i]));
				}
			} else {
				parseFastas(*infiles, osn, osnLen, osv, osvLen);
			}
			os = &osv;
		} else {
			assert(origs != NULL);
			os = origs;
		}
		
		// Go through the loaded reference files base-by-base and
		// sanity check against what we get by calling getBase and
		// getStretch
		for(size_t i = 0; i < os->size(); i++) {
			size_t olen = ((*os)[i]).length();
			size_t olenU32 = (olen + 12) / 4;
			uint32_t *buf = new uint32_t[olenU32];
			uint8_t *bufadj = (uint8_t*)buf;
			bufadj += getStretch(buf, i, 0, olen, tmp_destU32_);
			for(size_t j = 0; j < olen; j++) {
				assert_eq((int)(*os)[i][j], (int)bufadj[j]);
				assert_eq((int)(*os)[i][j], (int)getBase(i, j));
			}
			delete[] buf;
		}
	}
#endif
}

BitPairReference::~BitPairReference() {
	if(buf_ != NULL && !useMm_ && !useShmem_) delete[] buf_;
	if(sanityBuf_ != NULL) delete[] sanityBuf_;
}

/**
 * Return a single base of the reference.  Calling this repeatedly
 * is not an efficient way to retrieve bases from the reference;
 * use loadStretch() instead.
 *
 * This implementation scans linearly through the records for the
 * unambiguous stretches of the target reference sequence.  When
 * there are many records, binary search would be more appropriate.
 */
int BitPairReference::getBase(size_t tidx, size_t toff) const {
	uint32_t reci = refRecOffs_[tidx];   // first record for target reference sequence
	uint32_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
	assert_gt(recf, reci);
	uint32_t bufOff = refOffs_[tidx];
	uint32_t off = 0;
	// For all records pertaining to the target reference sequence...
	for(uint32_t i = reci; i < recf; i++) {
		assert_geq(toff, off);
		off += recs_[i].off;
		if(toff < off) {
			return 4;
		}
		assert_geq(toff, off);
		uint32_t recOff = off + recs_[i].len;
		if(toff < recOff) {
			toff -= off;
			bufOff += (uint32_t)toff;
			assert_lt(bufOff, bufSz_);
			const uint32_t bufElt = (bufOff) >> 2;
			const uint32_t shift = (bufOff & 3) << 1;
			return ((buf_[bufElt] >> shift) & 3);
		}
		bufOff += recs_[i].len;
		off = recOff;
		assert_geq(toff, off);
	} // end for loop over records
	return 4;
}

/**
 * Load a stretch of the reference string into memory at 'dest'.
 *
 * This implementation scans linearly through the records for the
 * unambiguous stretches of the target reference sequence.  When
 * there are many records, binary search would be more appropriate.
 */
int BitPairReference::getStretchNaive(
	uint32_t *destU32,
	size_t tidx,
	size_t toff,
	size_t count) const
{
	uint8_t *dest = (uint8_t*)destU32;
	uint32_t reci = refRecOffs_[tidx];   // first record for target reference sequence
	uint32_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
	assert_gt(recf, reci);
	uint32_t cur = 0;
	uint32_t bufOff = refOffs_[tidx];
	uint32_t off = 0;
	// For all records pertaining to the target reference sequence...
	for(uint32_t i = reci; i < recf; i++) {
		assert_geq(toff, off);
		off += recs_[i].off;
		for(; toff < off && count > 0; toff++) {
			dest[cur++] = 4;
			count--;
		}
		if(count == 0) break;
		assert_geq(toff, off);
		if(toff < off + recs_[i].len) {
			bufOff += (uint32_t)(toff - off); // move bufOff pointer forward
		} else {
			bufOff += recs_[i].len;
		}
		off += recs_[i].len;
		for(; toff < off && count > 0; toff++) {
			assert_lt(bufOff, bufSz_);
			const uint32_t bufElt = (bufOff) >> 2;
			const uint32_t shift = (bufOff & 3) << 1;
			dest[cur++] = (buf_[bufElt] >> shift) & 3;
			bufOff++;
			count--;
		}
		if(count == 0) break;
		assert_geq(toff, off);
	} // end for loop over records
	// In any chars are left after scanning all the records,
	// they must be ambiguous
	while(count > 0) {
		count--;
		dest[cur++] = 4;
	}
	assert_eq(0, count);
	return 0;
}

/**
 * Load a stretch of the reference string into memory at 'dest'.
 *
 * This implementation scans linearly through the records for the
 * unambiguous stretches of the target reference sequence.  When
 * there are many records, binary search would be more appropriate.
 */
int BitPairReference::getStretch(
	uint32_t *destU32,
	size_t tidx,
	size_t toff,
	size_t count
	ASSERT_ONLY(, SStringExpandable<uint32_t>& destU32_2)) const
{
	ASSERT_ONLY(size_t origCount = count);
	ASSERT_ONLY(size_t origToff = toff);
	if(count == 0) return 0;
	uint8_t *dest = (uint8_t*)destU32;
#ifndef NDEBUG
	destU32_2.clear();
	uint8_t *dest_2 = NULL;
	int off2;
	if((rand() % 10) == 0) {
		destU32_2.resize((origCount >> 2) + 2);
		off2 = getStretchNaive(destU32_2.wbuf(), tidx, origToff, origCount);
		dest_2 = ((uint8_t*)destU32_2.wbuf()) + off2;
	}
#endif
	destU32[0] = 0x04040404; // Add Ns, which we might end up using later
	uint32_t reci = refRecOffs_[tidx];   // first record for target reference sequence
	uint32_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
	assert_gt(recf, reci);
	uint32_t cur = 4; // keep a cushion of 4 bases at the beginning
	uint32_t bufOff = refOffs_[tidx];
	uint32_t off = 0;
	int offset = 4;
	bool firstStretch = true;
	// For all records pertaining to the target reference sequence...
	for(uint32_t i = reci; i < recf; i++) {
		ASSERT_ONLY(uint32_t origBufOff = bufOff);
		assert_geq(toff, off);
		off += recs_[i].off;
		assert_gt(count, 0);
		if(toff < off) {
			size_t cpycnt = min(off - toff, count);
			memset(&dest[cur], 4, cpycnt);
			count -= cpycnt;
			toff += cpycnt;
			cur += (uint32_t)cpycnt;
			if(count == 0) break;
		}
		assert_geq(toff, off);
		if(toff < off + recs_[i].len) {
			bufOff += (uint32_t)(toff - off); // move bufOff pointer forward
		} else {
			bufOff += recs_[i].len;
		}
		off += recs_[i].len;
		if(toff < off) {
			if(firstStretch) {
				if(toff + 8 < off && count > 8) {
					// We already added some Ns, so we have to do
					// a fixup at the beginning of the buffer so
					// that we can start clobbering at cur >> 2
					if(cur & 3) {
						offset -= (cur & 3);
					}
					uint32_t curU32 = cur >> 2;
					// Do the initial few bases
					if(bufOff & 3) {
						const uint32_t bufElt = (bufOff) >> 2;
						const int low2 = bufOff & 3;
						// Lots of cache misses on the following line
						destU32[curU32] = byteToU32_[buf_[bufElt]];
						for(int j = 0; j < low2; j++) {
							((char *)(&destU32[curU32]))[j] = 4;
						}
						curU32++;
						offset += low2;
						const int chars = 4 - low2;
						count -= chars;
						bufOff += chars;
						toff += chars;
					}
					assert_eq(0, bufOff & 3);
					uint32_t bufOffU32 = bufOff >> 2;
					uint32_t countLim = (uint32_t)count >> 2;
					uint32_t offLim = (uint32_t)((off - (toff + 4)) >> 2);
					uint32_t lim = min(countLim, offLim);
					// Do the fast thing for as far as possible
					for(uint32_t j = 0; j < lim; j++) {
						// Lots of cache misses on the following line
						destU32[curU32] = byteToU32_[buf_[bufOffU32++]];
#ifndef NDEBUG
						if(dest_2 != NULL) {
							assert_eq(dest[(curU32 << 2) + 0], dest_2[(curU32 << 2) - offset + 0]);
							assert_eq(dest[(curU32 << 2) + 1], dest_2[(curU32 << 2) - offset + 1]);
							assert_eq(dest[(curU32 << 2) + 2], dest_2[(curU32 << 2) - offset + 2]);
							assert_eq(dest[(curU32 << 2) + 3], dest_2[(curU32 << 2) - offset + 3]);
						}
#endif
						curU32++;
					}
					toff += (lim << 2);
					assert_leq(toff, off);
					assert_leq((lim << 2), count);
					count -= (lim << 2);
					bufOff = bufOffU32 << 2;
					cur = curU32 << 2;
				}
				// Do the slow thing for the rest
				for(; toff < off && count > 0; toff++) {
					assert_lt(bufOff, bufSz_);
					const uint32_t bufElt = (bufOff) >> 2;
					const uint32_t shift = (bufOff & 3) << 1;
					dest[cur++] = (buf_[bufElt] >> shift) & 3;
					bufOff++;
					count--;
				}
				firstStretch = false;
			} else {
				// Do the slow thing
				for(; toff < off && count > 0; toff++) {
					assert_lt(bufOff, bufSz_);
					const uint32_t bufElt = (bufOff) >> 2;
					const uint32_t shift = (bufOff & 3) << 1;
					dest[cur++] = (buf_[bufElt] >> shift) & 3;
					bufOff++;
					count--;
				}
			}
		}
		if(count == 0) break;
		assert_eq(recs_[i].len, bufOff - origBufOff);
		assert_geq(toff, off);
	} // end for loop over records
	// In any chars are left after scanning all the records,
	// they must be ambiguous
	while(count > 0) {
		count--;
		dest[cur++] = 4;
	}
	assert_eq(0, count);
	return offset;
}


/**
 * Parse the input fasta files, populating the szs list and writing the
 * .3.bt2 and .4.bt2 portions of the index as we go.
 */
pair<size_t, size_t>
BitPairReference::szsFromFasta(
	EList<FileBuf*>& is,
	const string& outfile,
	bool bigEndian,
	const RefReadInParams& refparams,
	EList<RefRecord>& szs,
	bool sanity)
{
	RefReadInParams parms = refparams;
	std::pair<size_t, size_t> sztot;
	if(!outfile.empty()) {
		string file3 = outfile + ".3.bt2";
		string file4 = outfile + ".4.bt2";
		// Open output stream for the '.3.bt2' file which will
		// hold the size records.
		ofstream fout3(file3.c_str(), ios::binary);
		if(!fout3.good()) {
			cerr << "Could not open index file for writing: \"" << file3 << "\"" << endl
				 << "Please make sure the directory exists and that permissions allow writing by" << endl
				 << "Bowtie." << endl;
			throw 1;
		}
		BitpairOutFileBuf bpout(file4.c_str());
		// Read in the sizes of all the unambiguous stretches of the genome
		// into a vector of RefRecords.  The input streams are reset once
		// it's done.
		writeU32(fout3, 1, bigEndian); // endianness sentinel
		bool color = parms.color;
		if(color) {
			parms.color = false;
			// Make sure the .3.bt2 and .4.bt2 files contain
			// nucleotides; not colors
			int numSeqs = 0;
			ASSERT_ONLY(std::pair<size_t, size_t> sztot2 =)
			fastaRefReadSizes(is, szs, parms, &bpout, numSeqs);
			parms.color = true;
			writeU32(fout3, (uint32_t)szs.size(), bigEndian); // write # records
			for(size_t i = 0; i < szs.size(); i++) {
				szs[i].write(fout3, bigEndian);
			}
			szs.clear();
			// Now read in the colorspace size records; these are
			// the ones that were indexed
			int numSeqs2 = 0;
			sztot = fastaRefReadSizes(is, szs, parms, NULL, numSeqs2);
			assert_eq(numSeqs, numSeqs2);
			assert_eq(sztot2.second, sztot.second + numSeqs);
		} else {
			int numSeqs = 0;
			sztot = fastaRefReadSizes(is, szs, parms, &bpout, numSeqs);
			writeU32(fout3, (uint32_t)szs.size(), bigEndian); // write # records
			for(size_t i = 0; i < szs.size(); i++) szs[i].write(fout3, bigEndian);
		}
		if(sztot.first == 0) {
			cerr << "Error: No unambiguous stretches of characters in the input.  Aborting..." << endl;
			throw 1;
		}
		assert_gt(sztot.first, 0);
		assert_gt(sztot.second, 0);
		bpout.close();
		fout3.close();
	} else {
		// Read in the sizes of all the unambiguous stretches of the
		// genome into a vector of RefRecords
		int numSeqs = 0;
		sztot = fastaRefReadSizes(is, szs, parms, NULL, numSeqs);
#ifndef NDEBUG
		if(parms.color) {
			parms.color = false;
			EList<RefRecord> szs2(EBWTB_CAT);
			int numSeqs2 = 0;
			ASSERT_ONLY(std::pair<size_t, size_t> sztot2 =)
			fastaRefReadSizes(is, szs2, parms, NULL, numSeqs2);
			assert_eq(numSeqs, numSeqs2);
			// One less color than base
			assert_geq(sztot2.second, sztot.second + numSeqs);
			parms.color = true;
		}
#endif
	}
	return sztot;
}
