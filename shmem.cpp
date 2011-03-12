/*
 * shmem.cpp
 *
 *  Created on: August 13, 2009
 *      Author: Ben Langmead
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
