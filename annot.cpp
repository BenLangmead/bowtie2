/*
 * annot.cpp
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#include <stdexcept>
#include "annot.h"

using namespace std;

/**
 * Parse an annotation-map file.
 */
void AnnotationMap::parse() {
	ifstream in(fname_);
	if(!in.good() && in.is_open()) {
		cerr << "Could not open annotation file " << fname_ << endl;
		throw 1;
	}
	while(in.peek() != EOF) {
		U32Pair pos;
		CharPair an;
		in >> pos.first >> pos.second >> an.first >> an.second;
		map_[pos] = an;
		while(isspace(in.peek())) in.get();
	}
	in.close();
}
