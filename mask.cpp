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

#include "mask.h"

// 5-bit pop count
int alts5[32] = {
	 0, 1, 1, 2, 1, 2, 2, 3,
	 1, 2, 2, 3, 2, 3, 3, 4,
	 1, 2, 2, 3, 2, 3, 3, 4,
	 2, 3, 3, 4, 3, 4, 4, 5
};

// Index of lowest set bit
int firsts5[32] = {
	-1, 0, 1, 0, 2, 0, 1, 0,
	 3, 0, 1, 0, 2, 0, 1, 0,
	 4, 0, 1, 0, 2, 0, 1, 0,
	 3, 0, 1, 0, 2, 0, 1, 0
};
