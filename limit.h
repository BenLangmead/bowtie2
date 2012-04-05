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

#ifndef LIMIT_H_
#define LIMIT_H_

#include <stdint.h>
#include <cstring>

extern uint8_t  MIN_U8;
extern uint8_t  MAX_U8;
extern uint16_t MIN_U16;
extern uint16_t MAX_U16;
extern uint32_t MIN_U32;
extern uint32_t MAX_U32;
extern uint64_t MIN_U64;
extern uint64_t MAX_U64;
extern size_t   MIN_SIZE_T;
extern size_t   MAX_SIZE_T;

extern int     MIN_I;
extern int     MAX_I;
extern int8_t  MIN_I8;
extern int8_t  MAX_I8;
extern int16_t MIN_I16;
extern int16_t MAX_I16;
extern int32_t MIN_I32;
extern int32_t MAX_I32;
extern int64_t MIN_I64;
extern int64_t MAX_I64;

#endif
