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

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "tokenize.h"
#include "ds.h"
#include "mem_ids.h"

using namespace std;

extern "C" {
	int bowtie_build(int argc, const char **argv);
}

/**
 * bowtie-build main function.  It is placed in a separate source file
 * to make it slightly easier to compile as a library.
 *
 * If the user specifies -A <file> as the first two arguments, main
 * will interpret that file as having one set of command-line arguments
 * per line, and will dispatch each batch of arguments one at a time to
 * bowtie-build.
 */
int main(int argc, const char **argv) {
	if(argc > 2 && strcmp(argv[1], "-A") == 0) {
		const char *file = argv[2];
		ifstream in;
		in.open(file);
		char buf[4096];
		int lastret = -1;
		while(in.getline(buf, 4095)) {
			EList<string> args(MISC_CAT);
			args.push_back(string(argv[0]));
			tokenize(buf, " \t", args);
			const char **myargs = (const char**)malloc(sizeof(char*)*args.size());
			for(size_t i = 0; i < args.size(); i++) {
				myargs[i] = args[i].c_str();
			}
			if(args.size() == 1) continue;
			lastret = bowtie_build((int)args.size(), myargs);
			free(myargs);
		}
		if(lastret == -1) {
			cerr << "Warning: No arg strings parsed from " << file << endl;
			return 0;
		}
		return lastret;
	} else {
		return bowtie_build(argc, argv);
	}
}
