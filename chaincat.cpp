/*
 * chaincat.cpp
 *
 * Print out a chained-hit file.
 *
 *  Created on: Jul 31, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <stdexcept>
#include "hit_set.h"
#include "filebuf.h"

using namespace std;

int main(int argc, char **argv) {
	try {
		if(argc <= 1) {
			cerr << "Error: must specify chain file as first argument" << endl;
			return 1;
		}
		FILE *in = fopen(argv[1], "rb");
		if(in == NULL) {
			cerr << "Could not open " << argv[1] << endl;
			return 1;
		}
		FileBuf fb(in);
		while(!fb.eof()) {
			HitSet s(fb);
			s.reportUpTo(cout);
		}
		fb.close();
	} catch(std::exception& e) {
		return 1;
	}
}
