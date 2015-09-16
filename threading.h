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

#ifdef WITH_TBB
# include <tbb/mutex.h>
# include <tbb/spin_mutex.h>
# ifdef WITH_AFFINITY
#  include <cstdlib>
#  include <sched.h>
#  include <tbb/task_group.h>
#  include <tbb/task_scheduler_observer.h>
#  include <tbb/atomic.h>
#  include <tbb/task_scheduler_init.h>
# endif
#else
# include "tinythread.h"
# include "fast_mutex.h"
#endif

#ifdef NO_SPINLOCK
# ifdef WITH_TBB
#   define MUTEX_T tbb::mutex
# else
#   define MUTEX_T tthread::mutex
# endif
#else
# ifdef WITH_TBB
#  	define MUTEX_T tbb::spin_mutex
# else
#  	define MUTEX_T tthread::fast_mutex
# endif
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

#ifdef WITH_TBB
#ifdef WITH_AFFINITY
//ripped entirely from;
//https://software.intel.com/en-us/blogs/2013/10/31/applying-intel-threading-building-blocks-observers-for-thread-affinity-on-intel
class concurrency_tracker: public tbb::task_scheduler_observer {
    tbb::atomic<int> num_threads;
public:
    concurrency_tracker() : num_threads() { observe(true); }
    /*override*/ void on_scheduler_entry( bool ) { ++num_threads; }
    /*override*/ void on_scheduler_exit( bool ) { --num_threads; }

    int get_concurrency() { return num_threads; }
};

class pinning_observer: public tbb::task_scheduler_observer {
    cpu_set_t *mask;
    int ncpus;

    const int pinning_step;
    tbb::atomic<int> thread_index;
public:
    pinning_observer( int pinning_step=1 ) : pinning_step(pinning_step), thread_index() {
        for ( ncpus = sizeof(cpu_set_t)/CHAR_BIT; ncpus < 16*1024 /* some reasonable limit */; ncpus <<= 1 ) {
            mask = CPU_ALLOC( ncpus );
            if ( !mask ) break;
            const size_t size = CPU_ALLOC_SIZE( ncpus );
            CPU_ZERO_S( size, mask );
            const int err = sched_getaffinity( 0, size, mask );
            if ( !err ) break;

            CPU_FREE( mask );
            mask = NULL;
            if ( errno != EINVAL )  break;
        }
        if ( !mask )
            std::cout << "Warning: Failed to obtain process affinity mask. Thread affinitization is disabled." << std::endl;
    }

/*override*/ void on_scheduler_entry( bool ) {
    if ( !mask ) return;

    const size_t size = CPU_ALLOC_SIZE( ncpus );
    const int num_cpus = CPU_COUNT_S( size, mask );
    int thr_idx =
//cwilks: we're one interface version lower than what
//is required for task arena (7000 vs. 7001)
#if USE_TASK_ARENA_CURRENT_SLOT
        tbb::task_arena::current_slot();
#else
        thread_index++;
#endif
#if __MIC__
    thr_idx += 1; // To avoid logical thread zero for the master thread on Intel(R) Xeon Phi(tm)
#endif
    thr_idx %= num_cpus; // To limit unique number in [0; num_cpus-1] range

        // Place threads with specified step
        int cpu_idx = 0;
        for ( int i = 0, offset = 0; i<thr_idx; ++i ) {
            cpu_idx += pinning_step;
            if ( cpu_idx >= num_cpus )
                cpu_idx = ++offset;
        }

        // Find index of 'cpu_idx'-th bit equal to 1
        int mapped_idx = -1;
        while ( cpu_idx >= 0 ) {
            if ( CPU_ISSET_S( ++mapped_idx, size, mask ) )
                --cpu_idx;
        }

        cpu_set_t *target_mask = CPU_ALLOC( ncpus );
        CPU_ZERO_S( size, target_mask );
        CPU_SET_S( mapped_idx, size, target_mask );
        const int err = sched_setaffinity( 0, size, target_mask );

        //std::cout << "Just set affinity for thread " << thr_idx << "\n";
        if ( err ) {
            std::cout << "Failed to set thread affinity!n";
            exit( EXIT_FAILURE );
        }
#if LOG_PINNING
        else {
            std::stringstream ss;
            ss << "Set thread affinity: Thread " << thr_idx << ": CPU " << mapped_idx << std::endl;
            std::cerr << ss.str();
        }
#endif
        CPU_FREE( target_mask );
    }

    ~pinning_observer() {
        if ( mask )
            CPU_FREE( mask );
    }
};

#endif
#endif

#endif
