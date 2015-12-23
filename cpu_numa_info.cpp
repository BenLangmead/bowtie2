#include "cpu_numa_info.h"

//adapted from http://stackoverflow.com/questions/7259363/measuring-numa-non-uniform-memory-access-no-observable-asymmetry-why
int bind_to_node(int thread_idx)
{
	//int numcpus = numa_num_task_cpus()
	//bitmask* bm = numa_bitmask_alloc(numcpus);
	int num_nodes = numa_max_possible_node();
	int numa_node = thread_idx % num_nodes;
	bitmask* nm = numa_allocate_nodemask();
	numa_bitmask_setbit(nm,numa_node);
	//void numa_bind(struct bitmask *nodemask); sets node AND memory
	//this one only sets node (not memeory)
	//numa_run_on_node_mask(struct bitmask *nodemask);	
	numa_run_on_node_mask(nm);
	numa_free_nodemask(nm);
	return numa_node;
	//numa_bitmask_free(bm);	
}

/// Based on http://stackoverflow.com/questions/16862620/numa-get-current-node-core
void get_cpu_and_node_(int& cpu, int& node) {
	unsigned long a,d,c;
	__asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
	node = (c & 0xFFF000)>>12;
	cpu = c & 0xFFF;
}
