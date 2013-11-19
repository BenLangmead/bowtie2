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


#ifndef BOWTIE_INDEX_TYPES_H
#define	BOWTIE_INDEX_TYPES_H

#ifdef BOWTIE_64BIT_INDEX
#define OFF_MASK 0xffffffffffffffff
#define OFF_LEN_MASK 0xc000000000000000
#define LS_SIZE 0x100000000000000
#define OFF_SIZE 8

typedef uint64_t TIndexOffU;
typedef int64_t TIndexOff;
    
#else
#define OFF_MASK 0xffffffff
#define OFF_LEN_MASK 0xc0000000
#define LS_SIZE 0x10000000
#define OFF_SIZE 4

typedef uint32_t TIndexOffU;
typedef int TIndexOff;

#endif /* BOWTIE_64BIT_INDEX */

extern const std::string gEbwt_ext;

#endif	/* BOWTIE_INDEX_TYPES_H */

