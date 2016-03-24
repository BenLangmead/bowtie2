#if WITH_TBB && WITH_AFFINITY && WITH_COHORTLOCK

#include "cohort.hpp"

	LocalLock::LocalLock()
	{
		this->id=0;
		local_lock = new MUTEX_L();
		local_counter = 0;
#if WITH_QUEUELOCK
		last_scoped_lock = NULL;
#endif
	}

	LocalLock::~LocalLock()
	{
		//printf("LocalLock destructor b4 %d\n",id);
		/*if(last_scoped_lock)
			delete last_scoped_lock;*/
		if(local_lock)
			delete local_lock;
		//printf("LocalLock destructor af %d\n",id);
	}

	void LocalLock::set_id(int id)
	{
		this->id=id;
	}

	void LocalLock::lock()
	{
		//update local counter to implement the "alone?" 
		//functionality of the local part of the cohort property
		local_counter.fetch_and_increment();
		//actually attempt to acquire lock
		assert(local_lock!=NULL);
#if WITH_QUEUELOCK
		MUTEX_L::scoped_lock* temp_ = new MUTEX_L::scoped_lock(*(local_lock));
		//now that lock is acquired update the pointer to the current held lock
		//for unlocking and destructing
		last_scoped_lock = temp_;
#else
		local_lock->lock();
#endif
		if(this->fetch_counter() > 0)
			local_counter.fetch_and_add(-1);
	}

	void LocalLock::unlock()
	{
#if WITH_QUEUELOCK
		if(last_scoped_lock)
			delete last_scoped_lock;
#else
		local_lock->unlock();
#endif
	}

	uint64_t LocalLock::fetch_counter()
	{
		return local_counter;
	}	

	CohortLock::CohortLock()
	{
		assert(numa_available()!=-1);
		num_numa_nodes = 1024;
		//printf("num numa nodes %d\n",num_numa_nodes);
		//this->starvation_limit = starvation_limit;
		starvation_counters = new int [num_numa_nodes]();
		own_global = new bool [num_numa_nodes]();
		local_locks = new LocalLock [num_numa_nodes];
		lock_called = new int [num_numa_nodes];
		lock_released = new int [num_numa_nodes];
		uint64_t i = 0;
		for(i=0;i<num_numa_nodes;i++)
		{
			local_locks[i].set_id(i+1);
		}
		global_lock = new MUTEX_G();
	}
	
	CohortLock::~CohortLock()
	{
		/*uint64_t i=0;
		for(i=0;i<num_numa_nodes;i++)
		{
			printf("thread: %p lock calls/releases for %d numa node: %d/%d\n",this,i,lock_called[i],lock_released[i]);
		}*/
		delete[] starvation_counters;
		delete[] own_global;
		delete[] local_locks;
		delete[] lock_called;
		delete[] lock_released;
		delete global_lock;
	}

	void CohortLock::reset_lock(int nthreads)
	{
#if not WITH_QUEUELOCK
		global_lock->reset_lock(nthreads);
#endif
	}

	int CohortLock::determine_numa_idx()
	{
		//TODO: figure out numa idx
		//uint64_t numa_idx = 0;
		//uint64_t numa_idx = numa_node_of_cpu(cpu_idx);
		//TODO: check this, use perferred for now, not sure if this is OK
		//int numa_idx = numa_preferred();
		int numa_idx = -1;
		int cpu_idx = -1;
		//assume we're running WITH_AFFINITY
		get_cpu_and_node_(cpu_idx, numa_idx);
		return numa_idx;
	}

	void CohortLock::lock()
	{
		assert(lockers_numa_idx == -1);
		int idx = this->determine_numa_idx();
		this->lock(idx);
		lockers_numa_idx = idx;
	}

	void CohortLock::lock(int numa_idx)
	{
		// get the local lock
		local_locks[numa_idx].lock();
		if(!own_global[numa_idx])
		{
			// now try for global
			global_lock->lock();
		}
		starvation_counters[numa_idx]++;
		own_global[numa_idx]=true;
	}
	
	void CohortLock::unlock()
	{
		assert(lockers_numa_idx != -1);
		int idx = lockers_numa_idx;
		lockers_numa_idx = -1;
		this->unlock(idx);
	}

	void CohortLock::unlock(int numa_idx)
	{
		//possible race condition, but shouldn't hurt us as then 
		//lock contention is presumably low per NUMA node
		if(local_locks[numa_idx].fetch_counter() == 0 
			|| starvation_counters[numa_idx] > STARVATION_LIMIT)
		{
			//relinquish global lock
			//printf("thread: giving up global lock %d\n",numa_idx);
			lock_released[numa_idx]++;
			global_lock->unlock();
			//reset NUMA node specific vars
			starvation_counters[numa_idx]=0;
			own_global[numa_idx]=false;
		}
		local_locks[numa_idx].unlock();
	}	

#endif
