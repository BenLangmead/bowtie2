#if WITH_TBB && WITH_AFFINITY && WITH_COHORTLOCK

#ifndef COHORT_H_
#define COHORT_H_

#include <numa.h>
#include <iostream>
#include <cstdio>
#include <assert.h>
#include <tbb/mutex.h>
#include <tbb/queuing_mutex.h>
#include <tbb/atomic.h>
#include "cpu_numa_info.h"

#define MUTEX_G tbb::mutex
#define MUTEX_L tbb::queuing_mutex

const int STARVATION_LIMIT=100;

class LocalLock
{
public:
	LocalLock();
	LocalLock(int id);
	~LocalLock();
	void lock();
	void unlock();
	uint64_t fetch_counter();
private:
	int id;
	MUTEX_L* local_lock;
	MUTEX_L::scoped_lock* last_scoped_lock;
	tbb::atomic<uint64_t> local_counter;
};

class CohortLock
{
public:
	CohortLock();
	~CohortLock();
	int determine_numa_idx();
	void lock(int numa_idx);
	void lock();
	void unlock(int numa_idx);
	void unlock();
private:
	uint64_t num_numa_nodes;
	int starvation_limit;
	int* starvation_counters;
	bool* own_global;
	MUTEX_G* global_lock;
	LocalLock* local_locks;
};

#endif

#endif
