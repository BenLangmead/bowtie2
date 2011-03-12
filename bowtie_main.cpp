/*
 * bowtie_main.cpp
 *
 *  Created on: Sep 10, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "tokenize.h"
#include "ds.h"

using namespace std;

extern "C" {
	int bowtie(int argc, const char **argv);
}

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
			lastret = bowtie(args.size(), myargs);
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
