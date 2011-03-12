/*
 *  sim.h
 */

#include <iostream>
#include <string>
#include "ds.h"
#include "sstring.h"
#include "random_source.h"

/**
 * Encapsulates a procedure for randomly generating a reference sequence.
 */
class SimRef {

	friend class SimRefBuilder;

public:
	
	SimRef() : refs_(), names_(), totlen_(0) { }
	
	/**
	 * Copy given SimRef into this one, using operator= to make deep
	 * copies of constituent lists.
	 */
	SimRef& operator=(const SimRef& o) {
		refs_ = o.refs_;
		names_ = o.names_;
		totlen_ = o.totlen_;
		return (*this);
	}
	
	/**
	 * Reset SimRef so that it doesn't contain any sequences.
	 */
	void reset() {
		refs_.clear();
		names_.clear();
		totlen_ = 0;
		assert(empty());
	}

	/**
	 * Return true iff there are no reference sequences in this SimRef.
	 */
	bool empty() const {
		return totlen_ == 0;
	}

	/**
	 * Check if the SimRef is internally consistent.
	 */
	bool repOk() const {
		assert_eq(refs_.size(), names_.size());
#ifndef NDEBUG
		size_t tot = 0;
		for(size_t i = 0; i < refs_.size(); i++) {
			assert_gt(refs_[i].length(), 0);
			tot += refs_[i].length();
		}
		assert_eq(totlen_, tot);
#endif
		return true;
	}
	
	/**
	 * Extract a random substring from one of the reference sequences
	 * and return it.
	 *
	 * len:            length of substring to extract
	 * watson:         true -> possibly generate from Watson strand
	 * crick:          true -> possibly generate from Crick strand
	 */
	void randSubstr(
		RandomSource& rnd,  // pseudo-random generator
		BTDnaMask& dst,     // put sampled substring here
		size_t len,         // length of substring to extract
		bool watson = true, // true -> possibly extract from Watson strand
		bool crick = true)  // true -> possibly extract from Crick strand
	{
		assert(!empty());
		assert_gt(totlen_, 0);
		//
		// First, select which reference sequence to extract from
		//
		uint32_t rndoff = rnd.nextU32() % totlen_;
		size_t tot = totlen_;
		ASSERT_ONLY(bool sampled = false);
		for(size_t i = 0; i < refs_.size(); i++) {
			if(rndoff >= tot) {
				// Get substring from reference 'i'
				ASSERT_ONLY(sampled = true);
				if(len <= refs_[i].length()) {
					refs_[i].randSubstr(rnd, dst, len, watson, crick);
				} else {
					// Fill 'dst' with the entire reference string
					dst = refs_[i];
				}
				break;
			}
			tot -= refs_[i].length();
		}
		assert(sampled);
	}
	
	/**
	 * Given ostream, dump the reference strings in FASTA format to the
	 * output stream.
	 */
	void dumpToFasta(std::ostream& os) {
		for(size_t i = 0; i < refs_.size(); i++) {
			os << '>' << names_[i] << '\n' << refs_[i].toZBuf() << '\n';
		}
	}

protected:

	// Reference sequences - expandable strings allowing IUPAC chars
	EList<BTDnaMask> refs_;
	// Names of reference sequences
	EList<std::string> names_;
	// Total length of all reference sequences
	size_t totlen_;
};

typedef std::pair<size_t, size_t> TSzPair;

/**
 * Markov model random reference generator.
 */
class SimRefBuilder {

public:

	/**
	 * Generates a random set of reference with number of references and
	 * length of references chosen uniformly randomly from the given
	 * ranges.
	 */
	void build(
		RandomSource& rnd,
		SimRef& sr,
		TSzPair numref,
		TSzPair reflen);

protected:

	/**
	 * Generate a random DNA sequence with length chosen uniformly randomly
	 * from the given range.  Return the length.
	 */
	uint32_t buildSeq(
		RandomSource& rnd,
		BTDnaMask& seq,
		TSzPair lenrange);

};

/**
 * Encapsulates a simulated read.
 */
class SimRead {

public:

protected:

	BTDnaMask seq_;
	

};

/**
 * Encapsulates a sequencing error model.  Provides a routine that,
 * given a read sequence, will generate a corresponding quality value
 * string and mutates the sequence accordingly.
 */
class SimErrorModel {

public:

	SimErrorModel() { }
	
	/**
	 * Apply error model to read.
	 */
	void apply(SimRead& r) {
	}
};

/**
 * Random read generator given a reference.
 */
class SimReadBuilder {

public:

	/**
	 * Generates a random set of unpaied reads from the given reference
	 * with number of reads and length of reads chosen uniformly
	 * randomly from the given ranges.
	 */
	void build(
		RandomSource& rnd,        // pseudo-random generator
		const SimRef& sr,         // reference
		EList<BTDnaMask>& reads,  // append reads here
		TSzPair numreads,         // # reads
		bool    samelen,          // reads all the same length?
		TSzPair readlen);         // range of read lengths

	/**
	 * Generates a random set of paired-end reads from the given
	 * reference with number of pairs, length of fragments, and length
	 * of mates chosen uniformly randomly from the given ranges.
	 */
	void buildPairs(
		RandomSource&     rnd,      // pseudo-random generator
		const SimRef&     sr,       // reference
		EList<BTDnaMask>& mate1s,   // append reads here
		EList<BTDnaMask>& mate2s,   // append reads here
		TSzPair           numreads, // # reads
		bool              samelen,  // reads all the same length?
		TSzPair           readlen); // range of read lengths

protected:

	/**
	 * Generate a random DNA sequence with length chosen uniformly randomly
	 * from the given range.  Return the length.
	 */
	uint32_t buildFragment(
		RandomSource& rnd,        // pseudo-random generator
		const SimRef& sr,         // reference
		BTDnaMask&    read,       // append reads here
		size_t        readlen);   // length of read to generate

};

/**
 * Encapsultes methods and data for driving the random testing of
 * Bowtie2.
 */
class SimDriver {
public:
	
	void run(
		size_t refs,
		TSzPair readsPerRef,
		uint32_t seed)
	{
		RandomSource rnd(seed);
		for(size_t i = 0; i < refs; i++) {
			SimRef sr;
			SimRefBuilder sb;
			sb.build(rnd, sr, make_pair<size_t>(20, 30), make_pair<size_t>(1, 5000));
		}
	}
};

