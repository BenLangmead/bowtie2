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

#include <limits>
#include "limit.h"

uint8_t  MIN_U8  = std::numeric_limits<uint8_t>::min();
uint8_t  MAX_U8  = std::numeric_limits<uint8_t>::max();
uint16_t MIN_U16 = std::numeric_limits<uint16_t>::min();
uint16_t MAX_U16 = std::numeric_limits<uint16_t>::max();
uint32_t MIN_U32 = std::numeric_limits<uint32_t>::min();
uint32_t MAX_U32 = std::numeric_limits<uint32_t>::max();
uint64_t MIN_U64 = std::numeric_limits<uint64_t>::min();
uint64_t MAX_U64 = std::numeric_limits<uint64_t>::max();
size_t   MIN_SIZE_T = std::numeric_limits<size_t>::min();
size_t   MAX_SIZE_T = std::numeric_limits<size_t>::max();

int      MIN_I   = std::numeric_limits<int>::min();
int      MAX_I   = std::numeric_limits<int>::max();
int8_t   MIN_I8  = std::numeric_limits<int8_t>::min();
int8_t   MAX_I8  = std::numeric_limits<int8_t>::max();
int16_t  MIN_I16 = std::numeric_limits<int16_t>::min();
int16_t  MAX_I16 = std::numeric_limits<int16_t>::max();
int32_t  MIN_I32 = std::numeric_limits<int32_t>::min();
int32_t  MAX_I32 = std::numeric_limits<int32_t>::max();
int64_t  MIN_I64 = std::numeric_limits<int64_t>::min();
int64_t  MAX_I64 = std::numeric_limits<int64_t>::max();
