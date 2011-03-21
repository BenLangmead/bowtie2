//
//  multikey_qsort.cpp
//

#include "multikey_qsort.h"

// 5 64-element buckets for bucket-sorting A, C, G, T, $
uint32_t bkts[4][4 * 1024 * 1024];
