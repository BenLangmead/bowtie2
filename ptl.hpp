#ifndef TKTLOCK_H_
#define TKTLOCK_H_

#include <tbb/atomic.h>
#include <assert.h>

//from https://www.quora.com/How-does-an-MCS-lock-work
//Partition Lock
class PTLLock
{
public:
	PTLLock()
	{
		PTLLock(1);
	}

	PTLLock(int num_partitions)
	{
		this->reset_lock(num_partitions);
	}
		
	~PTLLock()
	{
		delete[] partitions;
	}

	void reset_lock(int num_partitions)
	{
		assert(num_partitions>=1);
		this->num_partitions = num_partitions;
		//offset by 1 since we're starting grant/request at 1 due to 
		//defauls to 0
		partitions = new uint64_t [num_partitions+1]();
		request=1;
		grant=1;
		partitions[request % num_partitions]=grant;
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
	uint64_t* partitions;
	int num_partitions;
};

#endif
