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

#include <string>
#include "aligner_result.h"
#include "simple_func.h"
#include "util.h"
#include "scoring.h"

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
		size_t rdlen,
		char *inps)
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
		size_t rdlen,
		char *inps) const = 0;
};

/**
 * TODO: Do BowtieMapq on a per-thread basis prior to the mutex'ed output
 * function.
 *
 * topCoeff :: top_coeff
 * botCoeff :: bot_coeff
 * mx :: mapqMax
 * horiz :: mapqHorizon (sort of)
 *
 * sc1 <- tab$sc1
 * sc2 <- tab$sc2
 * mapq <- rep(mx, length(sc1))
 * diff_top <- ifelse(sc1 != best & sc2 != best, abs(best - abs(pmax(sc1, sc2))), 0)
 * mapq <- mapq - diff_top * top_coeff
 * diff_bot <- ifelse(sc2 != horiz, abs(abs(sc2) - abs(horiz)), 0)
 * mapq <- mapq - diff_bot * bot_coeff
 * mapq <- round(pmax(0, pmin(mx, mapq)))
 * tab$mapq <- mapq
 */
class BowtieMapq : public Mapq {

public:

	BowtieMapq(
		const SimpleFunc& scoreMin,
		float topCoeff,
		float botCoeff,
		float mapqMax,
		const Scoring& sc) :
		scoreMin_(scoreMin),
		topCoeff_(topCoeff),
		botCoeff_(botCoeff),
		mapqMax_(mapqMax),
		sc_(sc)
	{
		assert_gt(topCoeff_, 0);
		assert_gt(botCoeff_, 0);
	}

	virtual ~BowtieMapq() { }

	/**
	 * Given an AlnSetSumm, return a mapping quality calculated.
	 */
	virtual TMapq mapq2(
		const AlnSetSumm& s,
		const AlnFlags& flags,
		bool mate1,
		size_t rdlen,
		char *inps)     // put string representation of inputs here
		const
	{
		bool hasSecbest = VALID_AL_SCORE(s.secbest(mate1));
		if(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest) {
			return 255;
		}
		TAlScore scPer = (TAlScore)sc_.perfectScore(rdlen);
		TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
		float mapq = mapqMax_;
		TAlScore best = s.best(mate1).score();
		assert_leq(best, scPer);
		if(best < scPer) {
			mapq -= (float)(scPer - best) * topCoeff_;
		}
		TAlScore secbest = scMin-1;
		if(hasSecbest) {
			secbest = s.secbest(mate1).score();
			assert_geq(secbest, scMin);
			mapq -= (float)(secbest - scMin) * botCoeff_;
			TAlScore diff = abs(abs(best) - abs(secbest));
			if(diff < mapqMax_) {
				mapq -= (mapqMax_ - diff);
			}
		}
		TMapq ret = (TMapq)(std::max(0.0f, std::min(mapqMax_, mapq + 0.5f)));
		if(inps != NULL) {
			inps = itoa10<TAlScore>(best, inps);
			*inps++ = ',';
			inps = itoa10<TAlScore>(secbest, inps);
			*inps++ = ',';
			inps = itoa10<TMapq>(ret, inps);
		}
		return ret;
	}

	/**
	 * Given an AlnSetSumm, return a mapping quality calculated.
	 */
	virtual TMapq mapq(
		const AlnSetSumm& s,
		const AlnFlags& flags,
		bool mate1,
		size_t rdlen,
		char *inps)     // put string representation of inputs here
		const
	{
		bool hasSecbest = VALID_AL_SCORE(s.secbest(mate1));
		if(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest) {
			return 255;
		}
		TAlScore scPer = (TAlScore)sc_.perfectScore(rdlen);
		TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
		TAlScore secbest = scMin-1;
		TAlScore diff = (scPer - scMin);
		float sixth_2 = (float)scPer - diff * 0.1666f * 2; 
		float sixth_3 = (float)scPer - diff * 0.1666f * 3;
		TMapq ret = 0;
		TAlScore best = s.best(mate1).score();
		if(!hasSecbest) {
			// Top third?
			if(best >= sixth_2) {
				ret = 37;
			}
			// Otherwise in top half?
			else if(best >= sixth_3) {
				ret = 25;
			}
			// Otherwise has no second-best?
			else {
				ret = 8;
			}
		} else {
			secbest = s.secbest(mate1).score();
			TAlScore bestdiff = abs(abs(best)-abs(secbest));
			if(bestdiff >= diff * 0.1666f * 5) {
				ret = 5;
			} else if(bestdiff >= diff * 0.1666f * 4) {
				ret = 4;
			} else if(bestdiff >= diff * 0.1666f * 3) {
				ret = 3;
			} else if(bestdiff >= diff * 0.1666f * 2) {
				ret = 2;
			} else if(bestdiff >= diff * 0.1666f * 1) {
				ret = 1;
			}
		}
		if(inps != NULL) {
			inps = itoa10<TAlScore>(best, inps);
			*inps++ = ',';
			inps = itoa10<TAlScore>(secbest, inps);
			*inps++ = ',';
			inps = itoa10<TMapq>(ret, inps);
		}
		return ret;
	}

protected:

	SimpleFunc      scoreMin_;
	float           topCoeff_;
	float           botCoeff_;
	float           mapqMax_;
	const Scoring&  sc_;
};

#endif /*ndef UNIQUE_H_*/
