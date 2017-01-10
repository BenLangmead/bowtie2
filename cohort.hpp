#ifndef COHORT_H_
#define COHORT_H_

#include <iostream>
#include <assert.h>
#include <tbb/atomic.h>
#include "tkt.hpp"
#include "ptl.hpp"


class CohortLock {

public:
  CohortLock() : lockers_numa_idx(-1) {
    starvation_counters = new int[MAX_NODES]();
    own_global = new bool[MAX_NODES]();
    local_locks = new TKTLock[MAX_NODES];
  }

  ~CohortLock() {
    delete[] starvation_counters;
    delete[] own_global;
    delete[] local_locks;
  }

  void lock();
  void unlock();
private:
  static const int STARVATION_LIMIT = 100;
  static const int MAX_NODES = 128;
  volatile int*  starvation_counters; // 1 per node
  volatile bool* own_global;          // 1 per node
  volatile int   lockers_numa_idx;    // 1 per node
  TKTLock* local_locks;               // 1 per node
  PTLLock global_lock;
};

#endif
