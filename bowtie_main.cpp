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

#ifdef ENABLE_x86_64_v3
#include <unistd.h>
#endif

using namespace std;

extern "C" {
	int bowtie(int argc, const char **argv);
}

#ifdef ENABLE_x86_64_v3
void check_x86_64_v3(int argc, const char **argv) {
	if (__builtin_cpu_supports ("x86-64-v3") && (argc<126)) {
		const char* new_argv[128]; // should always be enough, but above check enforces it, too
		const char * org_path = argv[0];
		// Append -v256 to the original path 
		const int fn_len = strlen(org_path);
		char *new_path = (char*) malloc(fn_len+16);
		memcpy(new_path,org_path,fn_len);
		strncpy(new_path+fn_len,"-v256",15);

		for (int i=1; i<=argc; i++) new_argv[i] = argv[i]; // all but first the same, also copy final NULL
		new_argv[0] = new_path;
		// now replace the executable with new variant
		// assuming the executable exists... execvp will gracefully fail else, which is OK
		execvp(new_argv[0], (char *const *)new_argv);
		// we should never get out of the above call
		fprintf(stderr,"[WARNING] Failed to launch x86-64-v3 version, staying with default\n");
	}
}
#endif

/**
 * Bowtie main function.  It is placed in a separate source file to
 * make it slightly easier to compile Bowtie as a library.
 *
 * If the user specifies -A <file> as the first two arguments, main
 * will interpret that file as having one set of command-line arguments
 * per line, and will dispatch each batch of arguments one at a time to
 * bowtie.
 */
int main(int argc, const char **argv) {
#ifdef ENABLE_x86_64_v3
	check_x86_64_v3(argc, argv);
#endif
	if(argc > 2 && strcmp(argv[1], "-A") == 0) {
		const char *file = argv[2];
		ifstream in;
			in.open(file);
		char buf[4096];
		int lastret = -1;
		while(in.getline(buf, 4095)) {
			EList<string> args;
			args.push_back(string(argv[0]));
			tokenize(buf, " \t", args);
			const char **myargs = (const char**)malloc(sizeof(char*)*args.size());
			for(size_t i = 0; i < args.size(); i++) {
				myargs[i] = args[i].c_str();
			}
			if(args.size() == 1) continue;
			lastret = bowtie((int)args.size(), myargs);
			free(myargs);
		}
		if(lastret == -1) {
			cerr << "Warning: No arg strings parsed from " << file << endl;
			return 0;
		}
		return lastret;
	} else {
		return bowtie(argc, argv);
	}
}
