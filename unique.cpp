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

#include "unique.h"

using namespace std;

// There is no valid second-best alignment and the best alignment has a
// perfect score.
const TMapq unp_nosec_perf = 44;

// There is no valid second-best alignment.  We stratify the alignment
// score of the best alignment into 10 bins.
const TMapq unp_nosec[11] = {
	43, 42, 41, 36, 32, 27, 20, 11, 4, 1, 0
};

// The best alignment has a perfect score, and we stratify the distance
// between best and second-best alignment scores into 10 bins.
const TMapq unp_sec_perf[11] = {
	2, 16, 23, 30, 31, 32, 34, 36, 38, 40, 42
};

// The best alignment has a non-perfect score, and we stratify both by best
// alignment score (specifically, the maximum score minus the best "best")
// and by the distance between the best and second-best alignment scores
// ("difference").  Each is stratified into 10 bins.  Each row is a
// difference (smaller elts = smaller differences) and each column is a
// best score (smaller elts = higher best alignment scores).
const TMapq unp_sec[11][11] = {
	{  2,  2,  2,  1,  1, 0, 0, 0, 0, 0, 0},
	{ 20, 14,  7,  3,  2, 1, 0, 0, 0, 0, 0},
	{ 20, 16, 10,  6,  3, 1, 0, 0, 0, 0, 0},
	{ 20, 17, 13,  9,  3, 1, 1, 0, 0, 0, 0},
	{ 21, 19, 15,  9,  5, 2, 2, 0, 0, 0, 0},
	{ 22, 21, 16, 11, 10, 5, 0, 0, 0, 0, 0},
	{ 23, 22, 19, 16, 11, 0, 0, 0, 0, 0, 0},
	{ 24, 25, 21, 30,  0, 0, 0, 0, 0, 0, 0},
	{ 30, 26, 29,  0,  0, 0, 0, 0, 0, 0, 0},
	{ 30, 27,  0,  0,  0, 0, 0, 0, 0, 0, 0},
	{ 30,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
};

//
// Paired mapping quality:
//

// There is no valid second-best alignment and the best alignment has a
// perfect score.
const TMapq pair_nosec_perf = 44;
