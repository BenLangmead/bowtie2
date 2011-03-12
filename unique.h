/*
 * mapq.h
 *
 * Encapsulates objects and routines for determining whether and to
 * what extend the best alignment for a read is "unique."  In the
 * simplest scenario, uniqueness is determined by whether we found only
 * one alignment.  More complex scenarios might assign a uniqueness
 * score based that is a function of (of a summarized version of): all
 * the alignments found and their scores.
 *
 * Since mapping quality is related to uniqueness, objects and
 * routings for calculating mapping quality are also included here.
 */

#ifndef UNIQUE_H_
#define UNIQUE_H_

#include "aligner_result.h"

typedef int64_t TMapq;

/**
 * Class that returns yes-or-no answers to the question of whether a 
 */
class Uniqueness {
public:

	/**
	 * Given an AlnSetSumm, determine if the best alignment is "unique"
	 * according to some definition.
	 */
	static bool bestIsUnique(const AlnSetSumm& s) {
		assert(!s.empty());
		return !VALID_AL_SCORE(s.secbest());
	}
};

/**
 * Collection of routines for calculating mapping quality.
 */
class Mapq {
public:
	virtual ~Mapq() { }
	virtual TMapq mapq(const AlnSetSumm& s) const = 0;
};

class BlatMapq : public Mapq {
public:
	virtual ~BlatMapq() { }
	/**
	 * Given an AlnSetSumm, return a mapping quality calculated in a
	 * manner similar to Blat.
	 */
	virtual TMapq mapq(const AlnSetSumm& s) const {
		return 10; // TODO
	}
};

class BwaMapq : public Mapq {
public:
	virtual ~BwaMapq() { }
	/**
	 * Given an AlnSetSumm, return a mapping quality calculated in a
	 * manner similar to BWA.
	 */
	virtual TMapq mapq(const AlnSetSumm& s) const {
		return 20; // TODO
	}
};

class BowtieMapq : public Mapq {
public:
	virtual ~BowtieMapq() { }
	/**
	 * Given an AlnSetSumm, return a mapping quality calculated.
	 */
	virtual TMapq mapq(const AlnSetSumm& s) const {
		if(VALID_AL_SCORE(s.secbest())) {
			return 0;
		} else {
			return 40;
		}
	}
};

#endif /*ndef UNIQUE_H_*/
