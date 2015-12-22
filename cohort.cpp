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

	//MUTEX_L::scoped_lock* LocalLock::lock()
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

	CohortLock::CohortLock(uint64_t num_numa_nodes)
	{
		this->num_numa_nodes = num_numa_nodes;
		local_locks = new LocalLock [num_numa_nodes];
		global_lock = new MUTEX_G();
	}

	CohortLock::~CohortLock()
	{
		delete[] local_locks;
		delete global_lock;
	}

#endif
#endif
