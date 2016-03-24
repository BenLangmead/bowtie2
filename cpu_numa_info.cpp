#include "cpu_numa_info.h"

/// Based on http://stackoverflow.com/questions/16862620/numa-get-current-node-core
void get_cpu_and_node_(int& cpu, int& node) {
	unsigned long a,d,c;
	__asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
	node = (c & 0xFFF000)>>12;
	cpu = c & 0xFFF;
}
