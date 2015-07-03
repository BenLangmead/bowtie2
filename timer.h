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
#ifdef USE_FINE_TIMER
#include <sys/time.h>
#endif

using namespace std;

#ifdef USE_FINE_TIMER

/* Subtract the ‘struct timeval’ values X and Y, storing the result in RESULT.
 Return 1 if the difference is negative, otherwise 0.
 Borrowed from: https://www.gnu.org/software/libc/manual/html_node/Elapsed-Time.html
 */
static inline bool timeval_subtract(timeval& result, const timeval& xin, const timeval& yin) {
	/* Perform the carry for the later subtraction by updating y. */
	timeval x = xin;
	timeval y = yin;
	if (x.tv_usec < y.tv_usec) {
		int nsec = (y.tv_usec - x.tv_usec) / 1000000 + 1;
		y.tv_usec -= 1000000 * nsec;
		y.tv_sec += nsec;
	}
	if (x.tv_usec - y.tv_usec > 1000000) {
		int nsec = (x.tv_usec - y.tv_usec) / 1000000;
		y.tv_usec += 1000000 * nsec;
		y.tv_sec -= nsec;
	}

	/* Compute the time remaining to wait. tv_usec is certainly positive. */
	result.tv_sec = x.tv_sec - y.tv_sec;
	result.tv_usec = x.tv_usec - y.tv_usec;

	/* Return 1 if result is negative. */
	return x.tv_sec < y.tv_sec;
}

/**
 * Use gettimeofday() call to keep track of elapsed time between creation and
 * destruction.  If verbose is true, Timer will print a message showing
 * elapsed time to the given output stream upon destruction.
 */
class Timer {
public:
	Timer(ostream& out = cout, const char *msg = "", bool verbose = true) :
		_t(), _out(out), _msg(msg), _verbose(verbose)
	{
		gettimeofday(&_t, NULL);
	}

	/// Optionally print message
	~Timer() {
		if(_verbose) write(_out);
	}

	/// Return elapsed time since Timer object was created
	time_t elapsed() const {
		timeval f;
		gettimeofday(&f, NULL);
		return f.tv_sec - _t.tv_sec;
	}

	void write(ostream& out) {
		timeval f;
		gettimeofday(&f, NULL);
		timeval diff;
		timeval_subtract(diff, f, _t);
		time_t hours   = (diff.tv_sec / 60) / 60;
		time_t minutes = (diff.tv_sec / 60) % 60;
		time_t seconds = (diff.tv_sec % 60);
		time_t milliseconds = (diff.tv_usec / 1000);
		std::ostringstream oss;
		oss << _msg << setfill ('0') << setw (2) << hours << ":"
			<< setfill ('0') << setw (2) << minutes << ":"
			<< setfill ('0') << setw (2) << seconds << "."
			<< setfill ('0') << setw (3) << milliseconds << endl;
		out << oss.str().c_str();
	}

private:
	timeval     _t;
	ostream&    _out;
	const char *_msg;
	bool        _verbose;
};
#else
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
#endif // USE_FINE_TIMER

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
