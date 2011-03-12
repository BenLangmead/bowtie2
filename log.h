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
