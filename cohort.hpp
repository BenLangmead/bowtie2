#ifdef WITH_TBB
#ifdef WITH_COHORT

#ifndef COHORT_H_
#define COHORT_H_

#include <iostream>
#include <tbb/mutex.h>
#include <tbb/queuing_mutex.h>
#include <tbb/atomic.h>

#define MUTEX_G tbb::mutex
#define MUTEX_L tbb::queuing_mutex

class LocalLock
{
public:
	LocalLock();
	~LocalLock();
	//MUTEX_L::scoped_lock* lock();
	void lock();
	void unlock();
	uint64_t fetch_counter();
private:
	MUTEX_L* local_lock;
	MUTEX_L::scoped_lock* last_scoped_lock;
	tbb::atomic<uint64_t> local_counter;
};

class CohortLock
{
public:
	CohortLock(uint64_t num_numa_nodes);
	~CohortLock();
private:
	uint64_t num_numa_nodes;
	MUTEX_G* global_lock;
	LocalLock* local_locks;
};

#endif

#endif
#endif
