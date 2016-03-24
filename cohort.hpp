#if WITH_TBB && WITH_AFFINITY && WITH_COHORTLOCK

#ifndef COHORT_H_
#define COHORT_H_

#include <iostream>
#include <cstdio>
#include <assert.h>
#include <tbb/mutex.h>
#include <tbb/atomic.h>
#include "cpu_numa_info.h"

#if WITH_QUEUELOCK
# include <tbb/queuing_mutex.h>
# define MUTEX_G tbb::mutex
# define MUTEX_L tbb::queuing_mutex
#else
# include "tkt.hpp"
# include "ptl.hpp"
# define MUTEX_G PTLLock
# define MUTEX_L TKTLock
#endif

const int STARVATION_LIMIT=100;

class LocalLock
{
public:
	LocalLock();
	~LocalLock();
	void set_id(int id);
	void lock();
	void unlock();
	uint64_t fetch_counter();
private:
	int id;
	MUTEX_L* local_lock;
#if WITH_QUEUELOCK
	MUTEX_L::scoped_lock* last_scoped_lock;
#endif
	tbb::atomic<uint64_t> local_counter;
};

class CohortLock
{
public:
	CohortLock();
	~CohortLock();
	void reset_lock(int nthreads);
	int determine_numa_idx();
	void lock(int numa_idx);
	void lock();
	void unlock(int numa_idx);
	void unlock();
private:
	uint64_t num_numa_nodes;
	int starvation_limit;
	volatile int* starvation_counters;
	int* lock_called;
	int* lock_released;
	volatile bool* own_global;
	volatile int lockers_numa_idx;
	MUTEX_G* global_lock;
	LocalLock* local_locks;
};

#endif

#endif
