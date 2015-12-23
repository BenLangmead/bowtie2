#ifdef WITH_TBB
#ifdef WITH_COHORT

#include "cohort.hpp"

	LocalLock::LocalLock()
	{
		local_lock = new MUTEX_L();
		local_counter = 0;
		last_scoped_lock = NULL;
	}

	LocalLock::~LocalLock()
	{
		if(last_scoped_lock)
			delete last_scoped_lock;
		if(local_lock)
			delete local_lock;
	}

	void LocalLock::lock()
	{
		//update local counter to implement the "alone?" 
		//functionality of the local part of the cohort property
		local_counter.fetch_and_increment();
		//actually attempt to acquire lock
		MUTEX_L::scoped_lock* temp_ = new MUTEX_L::scoped_lock(*(local_lock));
		//now that lock is acquired update the pointer to the current held lock
		//for unlocking and destructing
		last_scoped_lock = temp_;
		if(this->fetch_counter() > 0)
			local_counter.fetch_and_add(-1);
	}

	void LocalLock::unlock()
	{
		if(last_scoped_lock)
			delete last_scoped_lock;
	}

	uint64_t LocalLock::fetch_counter()
	{
		return local_counter;
	}	

	CohortLock::CohortLock(uint64_t num_numa_nodes,int starvation_limit)
	{
		this->num_numa_nodes = num_numa_nodes;
		this->starvation_limit = starvation_limit;
		starvation_counters = new int [num_numa_nodes]();
		own_global = new bool [num_numa_nodes]();
		local_locks = new LocalLock [num_numa_nodes];
		global_lock = new MUTEX_G();
	}
	
	CohortLock::~CohortLock()
	{
		delete[] starvation_counters;
		delete[] own_global;
		delete[] local_locks;
		delete global_lock;
	}

	void CohortLock::lock()
	{
		//TODO: figure out numa idx
		uint64_t numa_idx = 0;
		//get the local lock
		local_locks[numa_idx].lock();
		if(!own_global[numa_idx])
		{
			//now try for global
			global_lock->lock();
		}
		starvation_counters[numa_idx]++;
		own_global[numa_idx]=true;
	}
	
	void CohortLock::unlock()
	{
		//TODO: figure out numa idx
		uint64_t numa_idx = 0;
		//possible race condition, but shouldn't hurt us as then 
		//lock contention is presumably low per NUMA node
		if(local_locks[numa_idx].fetch_counter() == 0 
			|| starvation_counters[numa_idx] > starvation_limit)
		{
			//relinquish global lock
			global_lock->unlock();
			//reset NUMA node specific vars
			starvation_counters[numa_idx]=0;
			own_global[numa_idx]=false;
		}
		local_locks[numa_idx].unlock();
	}	

#endif
#endif
