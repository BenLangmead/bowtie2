/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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

#ifndef BANDED_H_
#define BANDED_H_

#include "sse_util.h"

/**
 * Use SSE instructions to quickly find stretches with lots of matches, then
 * resolve alignments.
 */
class BandedSseAligner {

public:

	void init(
		int    *q,      // query, maskized
		size_t  qi,     // query start
		size_t  qf,     // query end
		int    *r,      // reference, maskized
		size_t  ri,     // reference start
		size_t  rf)     // reference end
	{
		
	}
	
	void nextAlignment() {
	}

protected:

	EList_m128i mat_;
};

#endif
