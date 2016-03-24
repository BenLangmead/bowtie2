#ifndef TKTLOCK_H_
#define TKTLOCK_H_

#include <tbb/atomic.h>

//from https://www.quora.com/How-does-an-MCS-lock-work
//and Lock Cohorting: A General Technique for Designing NUMA Locks
//(DAVID DICE and VIRENDRA J. MARATHE, Oracle Labs NIR SHAVIT, MIT)
class TKTLock
{
public:
	TKTLock()
	{
		//printf("TKTLock initialized\n");
		request=0;
		grant=0;
	}

	void lock()
	{
		//printf("lock %p b4 request:%u grant:%u\n",this,request,grant);
		//volatile uint64_t my_request = request.fetch_and_increment();
		volatile uint64_t my_request = request.fetch_and_increment();
		//request.fetch_and_increment();
		//printf("lock %p mid request:%u my_request:%d grant:%u\n",this,request,my_request,grant);
		while(my_request != grant)
		{
			//printf("lock %p mid request:%u my_request:%d grant:%u\r",this,request,my_request,grant);
		}
		//printf("lock %p af request:%u grant:%u\n",this,request,grant);
	}

	void unlock()
	{
		//printf("unlock %p b4 request:%u grant:%u\n",this,request,grant);
		grant++;
		//printf("unlock %p af request:%u grant:%u\n",this,request,grant);
	}	

	//here for compatibility with PTLLock
	void reset_lock(int void_){}

private:
	tbb::atomic<uint64_t> request;
	volatile uint64_t grant;
};

#endif
