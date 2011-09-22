/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
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

#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include "threading.h"

class SyncLogger {
public:
	SyncLogger() {
		MUTEX_INIT(lock_);
	}

	void msg(const char *s) {
		MUTEX_LOCK(lock_);
		std::cout << s << std::endl;
		MUTEX_UNLOCK(lock_);
	}

	void msg(const std::string& s) {
		MUTEX_LOCK(lock_);
		std::cout << s << std::endl;
		MUTEX_UNLOCK(lock_);
	}

private:
	MUTEX_T lock_;
};

extern SyncLogger glog;

#endif /*LOG_H_*/
