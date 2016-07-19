/*
 * Copyright 2016, Ben Langmead <langmea@cs.jhu.edu>
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

#ifndef HINTS_H_
#define HINTS_H_

#include <string>
#include "ds.h"
#include "read.h"

struct SeedHit {
	Interval refival;
	size_t rd5primeOff;
};

/** Return true iff read name indicates hints are present */
extern bool has_hint(const Read& r);

/**
 * Parse hints out of the read name and into the given list of SeedHits.
 */
extern void parse_hints(const Read& r,
                        EList<SeedHit>& hints,
                        const EMap<std::string, TRefId>& refidMap);

#endif /* HINTS_H_ */
