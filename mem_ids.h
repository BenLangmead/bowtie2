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

// For holding index data
#define EBWT_CAT  ((int) 1)
// For holding index-building data
#define EBWTB_CAT ((int) 2)
// For holding cache data
#define CA_CAT    ((int) 3)
// For holding group-walk-left bookkeeping data
#define GW_CAT    ((int) 4)
// For holding alignment bookkeeping data
#define AL_CAT    ((int) 5)
// For holding dynamic programming bookkeeping data
#define DP_CAT    ((int) 6)
// For holding alignment results and other hit objects
#define RES_CAT   ((int) 7)
#define MISC_CAT  ((int) 9)
#define DEBUG_CAT ((int)10)
