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

#ifndef MM_H_
#define MM_H_

/**
 * mm.h:
 *
 * Defines that make it easier to handle files in the two different MM
 * contexts: i.e. on Linux and Mac where MM is supported and POSIX I/O
 * functions work as expected, and on Windows where MM is not supported
 * and where there isn't POSIX I/O,
 */

#ifdef BOWTIE_MM
#define MM_FILE_CLOSE(x) if(x > 3) { close(x); }
#define MM_READ_RET ssize_t
#define MM_READ read
#define MM_SEEK lseek
#define MM_FILE int
#define MM_FILE_INIT -1
#else
#define MM_FILE_CLOSE(x) if(x != NULL) { fclose(x); }
#define MM_READ_RET size_t
#define MM_READ(file, dest, sz) fread(dest, 1, sz, file)
#define MM_SEEK fseek
#define MM_FILE FILE*
#define MM_FILE_INIT NULL
#endif

#endif /* MM_H_ */
