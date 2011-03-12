/*
 *  mask.cpp
 *  bowtie-beta1
 *
 *  Created by Benjamin Langmead on 2/8/11.
 *  Copyright 2011 Johns Hopkins University. All rights reserved.
 *
 */

#include "mask.h"

// 5-bit pop count
int alts5[32] = {
	-1, 1, 1, 2, 1, 2, 2, 3,
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
