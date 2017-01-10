#include <cstdio>
#include <numa.h>
#include "cohort.hpp"
#include "cpu_numa_info.h"

static int determine_numa_idx() {
  int numa_idx = -1;
  int cpu_idx = -1;
  get_cpu_and_node_(cpu_idx, numa_idx);
  return numa_idx;
}

void CohortLock::lock() {
  const int numa_idx = determine_numa_idx();
  local_locks[numa_idx].lock();
  if(!own_global[numa_idx]) {
      global_lock.lock();
  }
  starvation_counters[numa_idx]++;
  own_global[numa_idx] = true;
  lockers_numa_idx = numa_idx;
}

void CohortLock::unlock() {
  assert(lockers_numa_idx != -1);
  int numa_idx = lockers_numa_idx;
  lockers_numa_idx = -1; 
  if(local_locks[numa_idx].q_length() == 1 ||
     starvation_counters[numa_idx] > STARVATION_LIMIT)
  {
    // Either no one is waiting or we've held it too long
    global_lock.unlock();
    starvation_counters[numa_idx] = 0;
    own_global[numa_idx] = false;
  }
  local_locks[numa_idx].unlock();
}
