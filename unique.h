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
#include "simple_func.h"

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
	static bool bestIsUnique(
		const AlnSetSumm& s,
		const AlnFlags& flags,
		bool mate1,
		size_t rdlen)
	{
		assert(!s.empty());
		return !VALID_AL_SCORE(s.secbest(mate1));
	}
};

/**
 * Collection of routines for calculating mapping quality.
 */
class Mapq {

public:

	virtual ~Mapq() { }
	
	virtual TMapq mapq(
		const AlnSetSumm& s,
		const AlnFlags& flags,
		bool mate1,
		size_t rdlen) const = 0;
};

/**
 * TODO: Do BowtieMapq on a per-thread basis prior to the mutex'ed output
 * function.
 */
class BowtieMapq : public Mapq {

public:

	BowtieMapq(
		const SimpleFunc& scoreMin,
		float mapqDiffcoeff,
		float mapqHorizon,
		float mapqMax) :
		scoreMin_(scoreMin),
		mapqDiffcoeff_(mapqDiffcoeff),
		mapqHorizon_(mapqHorizon),
		mapqMax_(mapqMax)
	{ }

	virtual ~BowtieMapq() { }

	/**
	 * Given an AlnSetSumm, return a mapping quality calculated.
	 */
	virtual TMapq mapq(
		const AlnSetSumm& s,
		const AlnFlags& flags,
		bool mate1,
		size_t rdlen) const
	{
		EList<TAlScore>& scores = const_cast<EList<TAlScore>&>(scores_);
		if(VALID_AL_SCORE(s.secbest(mate1))) {
			TAlScore minsc = scoreMin_.f<TAlScore>((double)rdlen);
			scores.clear();
			scores.push_back(s.best(mate1).score());
			scores.push_back(s.secbest(mate1).score());
			scores.push_back((TAlScore)(mapqHorizon_ * minsc + 0.5f));
			return calc();
		} else {
			if(!flags.canMax() && !s.exhausted(mate1)) {
				return 255;
			} else {
				TAlScore minsc = scoreMin_.f<TAlScore>((double)rdlen);
				scores.clear();
				scores.push_back(s.best(mate1).score());
				scores.push_back((TAlScore)(mapqHorizon_ * minsc + 0.5f));
				return calc();
			}
		}
	}

protected:

	/**
	 * Given a collection of alignment scores, calculate a mapping quality.
	 */
	TMapq calc() const {
		assert_geq(scores_.size(), 2);
		TMapq sc1 = std::abs(scores_[0]);
		TMapq sc2 = std::abs(scores_[1]);
		float scaledDiff = fabs((float)(sc2 - sc1) * mapqDiffcoeff_);
		return min<TMapq>((TMapq)scaledDiff, (TMapq)mapqMax_);
	}

	SimpleFunc scoreMin_;
	float mapqDiffcoeff_;
	float mapqHorizon_;
	float mapqMax_;
	EList<TAlScore> scores_;
};

#endif /*ndef UNIQUE_H_*/
