#ifndef CPU_AND_NODE_H_
#define CPU_AND_NODE_H_

#include <numa.h>

extern void get_cpu_and_node_(int& cpu, int& node);
extern int bind_to_node(int thread_idx);

#endif
