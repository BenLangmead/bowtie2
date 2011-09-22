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

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <stdexcept>
#include <fcntl.h>
#include <sys/stat.h>
#include <utility>
#ifdef BOWTIE_MM
#include <sys/mman.h>
#include <sys/shm.h>
#endif
#include "endian_swap.h"
#include "ref_read.h"
#include "sequence_io.h"
#include "mm.h"
#include "shmem.h"
#include "timer.h"
#include "sstring.h"

/**
 * Concrete reference representation that bulk-loads the reference from
 * the bit-pair-compacted binary file and stores it in memory also in
 * bit-pair-compacted format.  The user may request reference
 * characters either on a per-character bases or by "stretch" using
 * getBase(...) and getStretch(...) respectively.
 *
 * Most of the complexity in this class is due to the fact that we want
 * to represent references with ambiguous (non-A/C/G/T) characters but
 * we don't want to use more than two bits per base.  This means we
 * need a way to encode the ambiguous stretches of the reference in a
 * way that is external to the bitpair sequence.  To accomplish this,
 * we use the RefRecords vector, which is stored in the .3.ebwt index
 * file.  The bitpairs themselves are stored in the .4.ebwt index file.
 *
 * Once it has been loaded, a BitPairReference is read-only, and is
 * safe for many threads to access at once.
 */
class BitPairReference {

public:
	/**
	 * Load from .3.ebwt/.4.ebwt Bowtie index files.
	 */
	BitPairReference(
		const string& in,
		bool color,
		bool sanity = false,
		EList<string>* infiles = NULL,
		EList<SString<char> >* origs = NULL,
		bool infilesSeq = false,
		bool useMm = false,
		bool useShmem = false,
		bool mmSweep = false,
		bool verbose = false,
		bool startVerbose = false);

	~BitPairReference();

	/**
	 * Return a single base of the reference.  Calling this repeatedly
	 * is not an efficient way to retrieve bases from the reference;
	 * use loadStretch() instead.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getBase(size_t tidx, size_t toff) const;

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getStretchNaive(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count) const;

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getStretch(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count
		ASSERT_ONLY(, SStringExpandable<uint32_t>& destU32_2)) const;

	/**
	 * Return the number of reference sequences.
	 */
	uint32_t numRefs() const {
		return nrefs_;
	}

	/**
	 * Return the approximate length of a reference sequence (it might leave
	 * off some Ns on the end).
	 *
	 * TODO: Is it still true that it might leave off Ns?
	 */
	uint32_t approxLen(uint32_t elt) const {
		assert_lt(elt, nrefs_);
		return refLens_[elt];
	}

	/**
	 * Return true iff buf_ and all the vectors are populated.
	 */
	bool loaded() const {
		return loaded_;
	}
	
	/**
	 * Given a reference sequence id, return its offset into the pasted
	 * reference string; i.e., return the number of unambiguous nucleotides
	 * preceding it.
	 */
	uint32_t pastedOffset(uint32_t idx) const {
		return refOffs_[idx];
	}

	/**
	 * Parse the input fasta files, populating the szs list and writing the
	 * .3.ebwt and .4.ebwt portions of the index as we go.
	 */
	static std::pair<size_t, size_t>
	szsFromFasta(
		EList<FileBuf*>& is,
		const string& outfile,
		bool bigEndian,
		const RefReadInParams& refparams,
		EList<RefRecord>& szs,
		bool sanity);
	
protected:

	uint32_t byteToU32_[256];

	EList<RefRecord> recs_;       /// records describing unambiguous stretches
	EList<uint32_t>  refLens_;    /// approx lens of ref seqs (excludes trailing ambig chars)
	EList<uint32_t>  refOffs_;    /// buf_ begin offsets per ref seq
	EList<uint32_t>  refRecOffs_; /// record begin/end offsets per ref seq
	uint8_t *buf_;      /// the whole reference as a big bitpacked byte array
	uint8_t *sanityBuf_;/// for sanity-checking buf_
	uint32_t bufSz_;    /// size of buf_
	uint32_t bufAllocSz_;
	uint32_t nrefs_;    /// the number of reference sequences
	bool     loaded_;   /// whether it's loaded
	bool     sanity_;   /// do sanity checking
	bool     useMm_;    /// load the reference as a memory-mapped file
	bool     useShmem_; /// load the reference into shared memory
	bool     verbose_;
	ASSERT_ONLY(SStringExpandable<uint32_t> tmp_destU32_);
};

#endif
