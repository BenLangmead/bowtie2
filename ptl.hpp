#ifndef PTLLOCK_H_
#define PTLLOCK_H_

#include <tbb/atomic.h>
#include <assert.h>


//same as number of threads unless that's not-appropriate
const int DEFAULT_NUM_PARTITIONS=120;

//from https://www.quora.com/How-does-an-MCS-lock-work
//and Lock Cohorting: A General Technique for Designing NUMA Locks
//(DAVID DICE and VIRENDRA J. MARATHE, Oracle Labs NIR SHAVIT, MIT)
//Partition Lock
class PTLLock
{
public:
	PTLLock()
	{
		grant=0;
		request=0;
		this->reset_lock(DEFAULT_NUM_PARTITIONS);
	}

	PTLLock(int num_partitions_)
	{
		grant=0;
		request=0;
		this->reset_lock(num_partitions_);
	}
		
	~PTLLock()
	{
		delete[] partitions;
	}

	void reset_lock(int num_partitions_)
	{
		assert(num_partitions_>=1);
		this->num_partitions = num_partitions_;
		if(grant >= 1)
		{
			delete[] partitions;
		}
		//offset by 1 since we're starting grant/request at 1 due to 
		//defauls to 0
		partitions = new uint64_t [num_partitions_+1]();
		request=1;
		grant=1;
		partitions[request % num_partitions_]=grant;
	}

	void lock()
	{
		uint64_t my_request = request.fetch_and_increment();
		uint64_t my_request_modded = my_request % num_partitions;
		while(partitions[my_request_modded] != my_request){}
		grant = my_request;
	}

	void unlock()
	{
		uint64_t my_grant = grant + 1;
		uint64_t my_grant_modded = my_grant % num_partitions;
		partitions[my_grant_modded] = my_grant;
	}
private:
	tbb::atomic<uint64_t> request;
	uint64_t grant;
	//use a large integer to avoid false sharing
	volatile uint64_t* partitions;
	int num_partitions;
};

#endif
