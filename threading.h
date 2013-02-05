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
// TODO: TTR
//#ifdef BOWTIE_PTHREADS
//#   include <pthread.h>
//#endif

#ifdef USE_SPINLOCK
#  	include "fast_mutex.h"
#  	define MUTEX_T tthread::fast_mutex
//    typedef MUTEX_T tthread::fast_mutex;
//#  	define MUTEX_INIT(l)  // no init required
//#  	define MUTEX_LOCK(l) (l).Enter() // use lock_guard<>
//#  	define MUTEX_UNLOCK(l) (l).Leave() // this will happen after lock_guard<> destruction.
#else
#  	ifdef BOWTIE_PTHREADS
#       define MUTEX_T tthread::mutex
//        typedef MUTEX_T tthread::mutex;
//#       define MUTEX_INIT(l) pthread_mutex_init(&l, NULL)
//#       define MUTEX_LOCK(l) pthread_mutex_lock(&l)
//#       define MUTEX_UNLOCK(l) pthread_mutex_unlock(&l)
#  	else
#       define MUTEX_T tthread::mutex
        //define MUTEX_T int
//#       define MUTEX_INIT(l) l = 0
//#       define MUTEX_LOCK(l) l = 1
//#       define MUTEX_UNLOCK(l) l = 0
#  endif /* BOWTIE_PTHREADS */
#endif /* USE_SPINLOCK */

// TODO: Are these really needed ?
//#ifdef BOWTIE_PTHREADS
//#    define THREAD_T               pthread_t
//#    define THREAD_CREATE(t, f, v) createThread(t, f, v)
//#    define THREAD_JOIN(t)         joinThread(t)
//#    define THREAD_YIELD()  0
//#else
//#    define THREAD_T        int
//#    define THREAD_CREATE(t, f, v) (*f)(v)
//#    define THREAD_JOIN(thr)
//#    define THREAD_YIELD()	0
//#endif /* USE_PTHREADS */
// end TODO

/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:
//	ThreadSafe(MUTEX_T* lock, bool really = true) {
    ThreadSafe(MUTEX_T* ptr_mutex, bool locked = true) {
		if(locked) {
		    // lock_ = lock;
		    // MUTEX_LOCK(*lock_);
		    this->ptr_mutex = ptr_mutex;
		    ptr_mutex->lock();
		}
		else
		    this->ptr_mutex = NULL;

//		else
//			lock_ = NULL;
	}
	~ThreadSafe() {
		//if(lock_ != NULL) delete lock_;
	    if (ptr_mutex != NULL)
	        ptr_mutex->unlock();
	}
private:
	//MUTEX_T mutex;
	MUTEX_T *ptr_mutex;
	//MUTEX_T *lock_;
};

// Not sure if we still need these.
// We can join and create threads using the tinythreads thread object.
//
//#ifdef BOWTIE_PTHREADS
//static inline void joinThread(pthread_t th) {
//	int ret, *tmp;
//	if((ret = pthread_join(th, (void**)(int**)&tmp)) != 0) {
//		std::cerr << "Error: pthread_join returned non-zero status: "
//		          << ret << std::endl;
//		throw 1;
//	}
//}
//static inline void createThread(pthread_t* th,
//                                void *(*start_routine) (void *),
//                                void *arg)
//{
//	int ret;
//	pthread_attr_t pt_attr;
//	pthread_attr_init(&pt_attr);
//	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
//	pthread_attr_setstacksize(&pt_attr, 2 << 20);
//	if((ret = pthread_create(th, &pt_attr, start_routine, arg)) != 0) {
//		std::cerr << "Error: pthread_create returned non-zero status: "
//		          << ret << std::endl;
//		throw 1;
//	}
//}
//#endif

#endif
