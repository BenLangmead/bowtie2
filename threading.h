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

#ifndef THREADING_H_
#define THREADING_H_

#include <iostream>
#include "tinythread.h"
#include "fast_mutex.h"

#ifdef NO_SPINLOCK
#   define MUTEX_T tthread::mutex
#else
#  	define MUTEX_T tthread::fast_mutex
#endif /* NO_SPINLOCK */


/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:
    ThreadSafe(MUTEX_T* ptr_mutex, bool locked = true) {
		if(locked) {
		    this->ptr_mutex = ptr_mutex;
		    ptr_mutex->lock();
		}
		else
		    this->ptr_mutex = NULL;
	}

	~ThreadSafe() {
	    if (ptr_mutex != NULL)
	        ptr_mutex->unlock();
	}
    
private:
	MUTEX_T *ptr_mutex;
};

#endif
