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

#ifdef BOWTIE_SHARED_MEM

#include <iostream>
#include <string>
#include <unistd.h>
#include <sys/shm.h>
#include <errno.h>
#include "shmem.h"

using namespace std;

/**
 * Notify other users of a shared-memory chunk that the leader has
 * finished initializing it.
 */
void notifySharedMem(void *mem, size_t len) {
	((volatile uint32_t*)((char*)mem + len))[0] = SHMEM_INIT;
}

/**
 * Wait until the leader of a shared-memory chunk has finished
 * initializing it.
 */
void waitSharedMem(void *mem, size_t len) {
	while(((volatile uint32_t*)((char*)mem + len))[0] != SHMEM_INIT) {
		sleep(1);
	}
}

#endif
