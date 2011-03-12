/*
 * refmap.cpp
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#include <stdexcept>
#include "refmap.h"
#include "assert_helpers.h"

using namespace std;

/**
 * Given a refid,offset pair in the index space, transform it into the
 * reference coordinate space according to the reference mappings
 * provided by the user.
 */
void ReferenceMap::map(U32Pair& h) const {
	if(h.first >= map_.size()) {
		cerr << "Could not find a reference-map entry for reference "
				  << h.first << " in map file \"" << fname_ << "\""
				  << endl;
		throw 1;
	}
	h.second += map_[h.first].second;
	h.first = map_[h.first].first;
}

/**
 * Parse a reference-map file.
 */
void ReferenceMap::parse() {
	ifstream in(fname_);
	if(!in.good() || !in.is_open()) {
		cerr << "Could not open reference map file " << fname_ << endl;
		throw 1;
	}
	int c;
	while((c = in.peek()) != EOF) {
		if(c == '>') {
			// This appears to be a name line
			in.get(); // chop off the initial '>'
			uint32_t off;
			in >> off;
			in.get(); // chop off tab
			char buf[1024];
			in.getline(buf, 1023);
			if(parseNames_) {
				if(names_.size() <= off) names_.resize(off+1);
				names_[off] = string(buf);
			}
			continue;
		}
		uint32_t id, off;
		in >> id >> off;
		map_.resize(map_.size()+1);
		map_.back().first = id;
		map_.back().second = off;
		while(isspace(in.peek())) in.get();
	}
	assert_eq(EOF, c);
	in.close();
}
