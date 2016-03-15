/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BT2_SEARCH_H_
#define BT2_SEARCH_H_
#ifdef WITH_TBB

#include <tbb/tbb.h>
#include <tbb/task_group.h>

class multiseedSearchWorker {
	int tid;

public:
	multiseedSearchWorker(const multiseedSearchWorker& W): tid(W.tid) {};
	multiseedSearchWorker(int id):tid(id) {};
	void operator()() const;
};

class multiseedSearchWorker_2p5 {
	int tid;

public:
	multiseedSearchWorker_2p5(const multiseedSearchWorker_2p5& W): tid(W.tid) {};
	multiseedSearchWorker_2p5(int id):tid(id) {};
	void operator()() const;

};

#endif /* WITH_TBB */
#endif /* BT2_SEARCH_H_ */
