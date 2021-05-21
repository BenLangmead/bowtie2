#include "bt2_locks.h"

void mcs_lock::lock(mcs_node &node) {
	node.next = nullptr;

	mcs_node *pred = q.exchange(&node, std::memory_order_acq_rel);
	if (pred) {
		node.unlocked = false;
		pred->next.store(&node, std::memory_order_release);
		spin_while_eq(node.unlocked, false);
	}
	node.unlocked.load(std::memory_order_acquire);
}

void mcs_lock::unlock(mcs_node &node) {
	if (!node.next.load(std::memory_order_acquire)) {
		mcs_node *node_ptr = &node;
		if (q.compare_exchange_strong(node_ptr,
					      (mcs_node *)nullptr,
					      std::memory_order_release))
			return;
		spin_while_eq(node.next, (mcs_node *)nullptr);
	}
	node.next.load(std::memory_order_acquire)->unlocked.store(true, std::memory_order_release);
}


void spin_lock::lock() {
	cpu_backoff backoff;
	while (flag.test_and_set(std::memory_order_acquire))
		backoff.pause();
}

void spin_lock::unlock() {
	flag.clear(std::memory_order_release);
}
