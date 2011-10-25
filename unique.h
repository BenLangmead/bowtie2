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

/*
 * unique.h
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
 * V2 of the MAPQ calculator
 */
class BowtieMapq2 : public Mapq {

public:

	BowtieMapq2(
		const SimpleFunc& scoreMin,
		const Scoring& sc) :
		scoreMin_(scoreMin),
		sc_(sc)
	{ }

	virtual ~BowtieMapq2() { }

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
		TMapq ret = 0;
		TAlScore best = s.best(mate1).score();
		if(!hasSecbest) {
			// Top third?
			if(best >= diff * 0.9f) {
				ret = 47;
			} else if(best >= diff * 0.8f) {
				ret = 45;
			} else if(best >= diff * 0.7f) {
				ret = 43;
			} else if(best >= diff * 0.6f) {
				ret = 41;
			} else if(best >= diff * 0.5f) {
				ret = 39;
			} else if(best >= diff * 0.4f) {
				ret = 37;
			} else if(best >= diff * 0.3f) {
				ret = 35;
			} else if(best >= diff * 0.2f) {
				ret = 33;
			} else if(best >= diff * 0.1f) {
				ret = 31;
			} else {
				ret = 29;
			}
		} else {
			secbest = s.secbest(mate1).score();
			TAlScore bestdiff = abs(abs(best)-abs(secbest));
			if(bestdiff >= diff * 0.9f) {
				ret = 19;
			} else if(bestdiff >= diff * 0.8f) {
				ret = 17;
			} else if(bestdiff >= diff * 0.7f) {
				ret = 15;
			} else if(bestdiff >= diff * 0.6f) {
				ret = 13;
			} else if(bestdiff >= diff * 0.5f) {
				ret = 11;
			} else if(bestdiff >= diff * 0.4f) {
				ret = 9;
			} else if(bestdiff >= diff * 0.3f) {
				ret = 7;
			} else if(bestdiff >= diff * 0.2f) {
				ret = 5;
			} else if(bestdiff >= diff * 0.1f) {
				ret = 3;
			} else {
				ret = 1;
			}
		}
		if(flags.alignedConcordant()) {
			ret = (TMapq)(ret * 1.3f);
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
	const Scoring&  sc_;
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
		const Scoring& sc) :
		scoreMin_(scoreMin),
		sc_(sc)
	{ }

	virtual ~BowtieMapq() { }

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
				ret = 10;
			}
		} else {
			secbest = s.secbest(mate1).score();
			TAlScore bestdiff = abs(abs(best)-abs(secbest));
			if(bestdiff >= diff * 0.1666f * 5) {
				ret = 6;
			} else if(bestdiff >= diff * 0.1666f * 4) {
				ret = 5;
			} else if(bestdiff >= diff * 0.1666f * 3) {
				ret = 4;
			} else if(bestdiff >= diff * 0.1666f * 2) {
				ret = 3;
			} else if(bestdiff >= diff * 0.1666f * 1) {
				ret = 2;
			} else {
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
	const Scoring&  sc_;
};

/**
 * Create and return new MAPQ calculating object.
 */
static inline Mapq *new_mapq(
	int version,
	const SimpleFunc& scoreMin,
	const Scoring& sc)
{
	if(version == 2) {
		return new BowtieMapq2(scoreMin, sc);
	} else {
		return new BowtieMapq(scoreMin, sc);
	}
}

#endif /*ndef UNIQUE_H_*/
