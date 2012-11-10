/*
 * Copyright 2012, Ben Langmead <blangmea@jhsph.edu>
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

#include "aligner_driver.h"

void AlignerDriverRootSelector::select(
	const DescentQuery& q,
	const DescentQuery* qo,
	EList<DescentConfig>& confs,
	EList<DescentRoot> roots)
{
	// Calculate interval length for both mates
	int interval = rootIval_.f<int>((double)q.length());
	if(qo != NULL) {
		// Boost interval length by 20% for paired-end reads
		interval = (int)(interval * 1.2 + 0.5);
	}
	float pri = 0.0f;
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		// Put down left-to-right roots w/r/t forward and reverse-complement reads
		{
			bool first = true;
			size_t i = 0;
			while(first || (i + landing_ <= q.length())) {
				confs.expand();
				confs.back().cons.init(descCons_);
				roots.expand();
				roots.back().init(
					i,          // offset from 5' end
					true,       // left-to-right?
					fw,         // fw?
					q.length(), // query length
					pri);       // root priority
				i += interval;
				first = false;
			}
		}
		// Put down right-to-left roots w/r/t forward and reverse-complement reads
		{
			bool first = true;
			size_t i = 0;
			while(first || (i + landing_ <= q.length())) {
				confs.expand();
				confs.back().cons.init(descCons_);
				roots.expand();
				roots.back().init(
					q.length() - i - 1, // offset from 5' end
					false,              // left-to-right?
					fw,                 // fw?
					q.length(),         // query length
					pri);               // root priority
				i += interval;
				first = false;
			}
		}
	}
}
