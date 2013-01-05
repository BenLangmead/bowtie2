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

#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

/**
 * Use time() call to keep track of elapsed time between creation and
 * destruction.  If verbose is true, Timer will print a message showing
 * elapsed time to the given output stream upon destruction.
 */
class Timer {
public:
	Timer(ostream& out = cout, const char *msg = "", bool verbose = true) :
		_t(time(0)), _out(out), _msg(msg), _verbose(verbose) { }

	/// Optionally print message
	~Timer() {
		if(_verbose) write(_out);
	}
	
	/// Return elapsed time since Timer object was created
	time_t elapsed() const {
		return time(0) - _t;
	}
	
	void write(ostream& out) {
		time_t passed = elapsed();
		// Print the message supplied at construction time followed
		// by time elapsed formatted HH:MM:SS 
		time_t hours   = (passed / 60) / 60;
		time_t minutes = (passed / 60) % 60;
		time_t seconds = (passed % 60);
		std::ostringstream oss;
		oss << _msg << setfill ('0') << setw (2) << hours << ":"
		           << setfill ('0') << setw (2) << minutes << ":"
		           << setfill ('0') << setw (2) << seconds << endl;
		out << oss.str().c_str();
	}
	
private:
	time_t      _t;
	ostream&    _out;
	const char *_msg;
	bool        _verbose;
};

static inline void logTime(std::ostream& os, bool nl = true) {
	struct tm *current;
	time_t now;
	time(&now);
	current = localtime(&now);
	std::ostringstream oss;
	oss << setfill('0') << setw(2)
	    << current->tm_hour << ":"
	    << setfill('0') << setw(2)
	    << current->tm_min << ":"
	    << setfill('0') << setw(2)
	    << current->tm_sec;
	if(nl) oss << std::endl;
	os << oss.str().c_str();
}

#endif /*TIMER_H_*/
