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
		size_t ordlen,
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
		size_t ordlen,
		char *inps) const = 0;
};

extern const TMapq unp_nosec_perf;
extern const TMapq unp_nosec[11];
extern const TMapq unp_sec_perf[11];
extern const TMapq unp_sec[11][11];
extern const TMapq pair_nosec_perf;

/**
 * V3 of the MAPQ calculator
 */
class BowtieMapq3 : public Mapq {

public:

	BowtieMapq3(
		const SimpleFunc& scoreMin,
		const Scoring& sc) :
		scoreMin_(scoreMin),
		sc_(sc)
	{ }

	virtual ~BowtieMapq3() { }

	/**
	 * Given an AlnSetSumm, return a mapping quality calculated.
	 */
	virtual TMapq mapq(
		const AlnSetSumm& s,
		const AlnFlags& flags,
		bool mate1,
		size_t rdlen,
		size_t ordlen,
		char *inps)     // put string representation of inputs here
		const
	{
		if(s.paired()) {
			return pair_nosec_perf;
		} else {
			bool hasSecbest = VALID_AL_SCORE(s.secbest(mate1));
			if(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest) {
				return 255;
			}
			TAlScore scMax = (TAlScore)sc_.perfectScore(rdlen);
			TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
			assert_geq(scMax, scMin);
			TAlScore best  = scMax - s.best(mate1).score(); // best score (lower=better)
			size_t best_bin = (size_t)((double)best * (10.0 / (double)(scMax - scMin)) + 0.5);
			assert_geq(best_bin, 0);
			assert_lt(best_bin, 11);
			if(hasSecbest) {
				assert_geq(s.best(mate1).score(), s.secbest(mate1).score());
				size_t diff = s.best(mate1).score() - s.secbest(mate1).score();
				size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(scMax - scMin)) + 0.5);
				assert_geq(diff_bin, 0);
				assert_lt(diff_bin, 11);
				// A valid second-best alignment was found
				if(best == scMax) {
					// Best alignment has perfect score
					return unp_sec_perf[best_bin];
				} else {
					// Best alignment has less than perfect score
					return unp_sec[diff_bin][best_bin];
				}
			} else {
				// No valid second-best alignment was found
				if(best == scMax) {
					// Best alignment has perfect score
					return unp_nosec_perf;
				} else {
					// Best alignment has less than perfect score
					return unp_nosec[best_bin];
				}
			}
		}
	}

protected:

	SimpleFunc      scoreMin_;
	const Scoring&  sc_;
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
		const AlnFlags&   flags,
		bool              mate1,
		size_t            rdlen,
		size_t            ordlen,
		char             *inps)   // put string representation of inputs here
		const
	{
		bool hasSecbest = s.paired() ?
			VALID_AL_SCORE(s.secbestPaired()) :
			VALID_AL_SCORE(s.secbest(mate1));
		if(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest) {
			return 255;
		}
		TAlScore scPer = (TAlScore)sc_.perfectScore(rdlen);
		if(s.paired()) {
			scPer += (TAlScore)sc_.perfectScore(ordlen);
		}
		TAlScore scMin = scoreMin_.f<TAlScore>((float)rdlen);
		if(s.paired()) {
			scMin += scoreMin_.f<TAlScore>((float)ordlen);
		}
		TAlScore secbest = scMin-1;
		TAlScore diff = (scPer - scMin);  // scores can vary by up to this much
		TMapq ret = 0;
		TAlScore best = s.paired() ?
			s.bestPaired().score() : s.best(mate1).score();
		// best score but normalized so that 0 = worst valid score
		TAlScore bestOver = best - scMin;
		if(sc_.monotone) {
			// End-to-end alignment
			if(!hasSecbest) {
				if     (bestOver >= diff * (double)0.8f) ret = 42;
				else if(bestOver >= diff * (double)0.7f) ret = 40;
				else if(bestOver >= diff * (double)0.6f) ret = 24;
				else if(bestOver >= diff * (double)0.5f) ret = 23;
				else if(bestOver >= diff * (double)0.4f) ret = 8;
				else if(bestOver >= diff * (double)0.3f) ret = 3;
				else                                     ret = 0;
			} else {
				secbest = s.paired() ?
					s.secbestPaired().score() : s.secbest(mate1).score();
				TAlScore bestdiff = abs(abs(best)-abs(secbest));
				if(bestdiff >= diff * (double)0.9f) {
					if(bestOver == diff) {
						ret = 39;
					} else {
						ret = 33;
					}
				} else if(bestdiff >= diff * (double)0.8f) {
					if(bestOver == diff) {
						ret = 38;
					} else {
						ret = 27;
					}
				} else if(bestdiff >= diff * (double)0.7f) {
					if(bestOver == diff) {
						ret = 37;
					} else {
						ret = 26;
					}
				} else if(bestdiff >= diff * (double)0.6f) {
					if(bestOver == diff) {
						ret = 36;
					} else {
						ret = 22;
					}
				} else if(bestdiff >= diff * (double)0.5f) {
					// Top third is still pretty good
					if       (bestOver == diff) {
						ret = 35;
					} else if(bestOver >= diff * (double)0.84f) {
						ret = 25;
					} else if(bestOver >= diff * (double)0.68f) {
						ret = 16;
					} else {
						ret = 5;
					}
				} else if(bestdiff >= diff * (double)0.4f) {
					// Top third is still pretty good
					if       (bestOver == diff) {
						ret = 34;
					} else if(bestOver >= diff * (double)0.84f) {
						ret = 21;
					} else if(bestOver >= diff * (double)0.68f) {
						ret = 14;
					} else {
						ret = 4;
					}
				} else if(bestdiff >= diff * (double)0.3f) {
					// Top third is still pretty good
					if       (bestOver == diff) {
						ret = 32;
					} else if(bestOver >= diff * (double)0.88f) {
						ret = 18;
					} else if(bestOver >= diff * (double)0.67f) {
						ret = 15;
					} else {
						ret = 3;
					}
				} else if(bestdiff >= diff * (double)0.2f) {
					// Top third is still pretty good
					if       (bestOver == diff) {
						ret = 31;
					} else if(bestOver >= diff * (double)0.88f) {
						ret = 17;
					} else if(bestOver >= diff * (double)0.67f) {
						ret = 11;
					} else {
						ret = 0;
					}
				} else if(bestdiff >= diff * (double)0.1f) {
					// Top third is still pretty good
					if       (bestOver == diff) {
						ret = 30;
					} else if(bestOver >= diff * (double)0.88f) {
						ret = 12;
					} else if(bestOver >= diff * (double)0.67f) {
						ret = 7;
					} else {
						ret = 0;
					}
				} else if(bestdiff > 0) {
					// Top third is still pretty good
					if(bestOver >= diff * (double)0.67f) {
						ret = 6;
					} else {
						ret = 2;
					}
				} else {
					assert_eq(bestdiff, 0);
					// Top third is still pretty good
					if(bestOver >= diff * (double)0.67f) {
						ret = 1;
					} else {
						ret = 0;
					}
				}
			}
		} else {
			// Local alignment
			if(!hasSecbest) {
				if     (bestOver >= diff * (double)0.8f) ret = 44;
				else if(bestOver >= diff * (double)0.7f) ret = 42;
				else if(bestOver >= diff * (double)0.6f) ret = 41;
				else if(bestOver >= diff * (double)0.5f) ret = 36;
				else if(bestOver >= diff * (double)0.4f) ret = 28;
				else if(bestOver >= diff * (double)0.3f) ret = 24;
				else                                     ret = 22;
			} else {
				secbest = s.paired() ?
					s.secbestPaired().score() : s.secbest(mate1).score();
				TAlScore bestdiff = abs(abs(best)-abs(secbest));
				if     (bestdiff >= diff * (double)0.9f) ret = 40;
				else if(bestdiff >= diff * (double)0.8f) ret = 39;
				else if(bestdiff >= diff * (double)0.7f) ret = 38;
				else if(bestdiff >= diff * (double)0.6f) ret = 37;
				else if(bestdiff >= diff * (double)0.5f) {
					if     (bestOver == diff)                 ret = 35;
					else if(bestOver >= diff * (double)0.50f) ret = 25;
					else                                      ret = 20;
				} else if(bestdiff >= diff * (double)0.4f) {
					if     (bestOver == diff)                 ret = 34;
					else if(bestOver >= diff * (double)0.50f) ret = 21;
					else                                      ret = 19;
				} else if(bestdiff >= diff * (double)0.3f) {
					if     (bestOver == diff)                 ret = 33;
					else if(bestOver >= diff * (double)0.5f)  ret = 18;
					else                                      ret = 16;
				} else if(bestdiff >= diff * (double)0.2f) {
					if     (bestOver == diff)                 ret = 32;
					else if(bestOver >= diff * (double)0.5f)  ret = 17;
					else                                      ret = 12;
				} else if(bestdiff >= diff * (double)0.1f) {
					if     (bestOver == diff)                 ret = 31;
					else if(bestOver >= diff * (double)0.5f)  ret = 14;
					else                                      ret = 9;
				} else if(bestdiff > 0) {
					if(bestOver >= diff * (double)0.5f)       ret = 11;
					else                                      ret = 2;
				} else {
					assert_eq(bestdiff, 0);
					if(bestOver >= diff * (double)0.5f)       ret = 1;
					else                                      ret = 0;
				}
			}
		}
		//if(flags.alignedConcordant()) {
		//	ret = (TMapq)((double)ret * 1.15);
		//}
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
		size_t ordlen,
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
		float sixth_2 = (float)(scPer - diff * (double)0.1666f * 2); 
		float sixth_3 = (float)(scPer - diff * (double)0.1666f * 3);
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
			if(bestdiff >= diff * 0.1666 * 5) {
				ret = 6;
			} else if(bestdiff >= diff * 0.1666 * 4) {
				ret = 5;
			} else if(bestdiff >= diff * 0.1666 * 3) {
				ret = 4;
			} else if(bestdiff >= diff * 0.1666 * 2) {
				ret = 3;
			} else if(bestdiff >= diff * 0.1666 * 1) {
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
	if(version == 3) {
		return new BowtieMapq3(scoreMin, sc);
	} else if(version == 2) {
		return new BowtieMapq2(scoreMin, sc);
	} else {
		return new BowtieMapq(scoreMin, sc);
	}
}

#endif /*ndef UNIQUE_H_*/
