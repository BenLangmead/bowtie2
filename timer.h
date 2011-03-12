#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
#include <iostream>
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
		unsigned int hours   = (passed / 60) / 60;
		unsigned int minutes = (passed / 60) % 60;
		unsigned int seconds = (passed % 60);
		out << _msg << setfill ('0') << setw (2) << hours << ":"
		            << setfill ('0') << setw (2) << minutes << ":"
		            << setfill ('0') << setw (2) << seconds << endl;
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
	os << setfill('0') << setw(2)
	    << current->tm_hour << ":"
	    << setfill('0') << setw(2)
	    << current->tm_min << ":"
	    << setfill('0') << setw(2)
	    << current->tm_sec;
	if(nl) os << std::endl;
}

#endif /*TIMER_H_*/
