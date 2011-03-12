/*
 * ref_coord.cpp
 */

#include "ref_coord.h"
#include <iostream>

using namespace std;

ostream& operator<<(ostream& out, const Interval& c) {
	out << c.upstream() << "+" << c.len();
	return out;
}

ostream& operator<<(ostream& out, const Coord& c) {
	out << c.ref() << ":" << c.off();
	return out;
}
