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
#include <sstream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include "tokenize.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <readline/readline.h>
#include <readline/history.h>

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

static volatile sig_atomic_t done = false;

static const char *options[] = {
"--al",                 "--al-conc",               "--dpad",            "--end-to-end",
"--fast",               "--fast-local",            "--fr",              "--gbar",
"--ignore-quals",       "--int-quals",             "--local",           "--ma",
"--met",                "--met-file",              "--met-stderr",      "--mm",
"--mp",                 "--n-ceil",                "--no-1mm-upfront",  "--no-contain",
"--no-discordant",      "--no-dovetail",           "--no-head",         "--no-mixed",
"--no-overlap",         "--no-sq",                 "--no-unal",         "--nofw",
"--non-deterministic",  "--norc",                  "--np",              "--omit-sec-seq",
"--phred33",            "--phred64",               "--qc-filter",       "--qseq",
"--quiet",              "--rdg",                   "--reorder",         "--rfg",
"--rg",                 "--rg-id",                 "--score-min",       "--seed",
"--sensitive",          "--sensitive-local",       "--un",              "--un-conc",
"--un-gz",              "--version",               "--very-fast",       "--very-fast-local",
"--very-sensitive",     "--very-sensitive-local",  "-3",                "-5",
"-D",                   "-I",                      "-L",                "-N",
"-R",                   "-X",                      "-a",                "-c",
"-f",                   "-h",                      "-i",                "-k",
"-p",                   "-q",                      "-r",                "-s",
"-t",                   "-u",                      "-1",                "-2",
"-S",                   "-U",                      "--all",             "--ff",
"--help",               "--maxins",                "--minins",          "--rf",
"--skip",               "--threads",               "--time",            "--trim3",
"--trim5",              "--upto",                  NULL
};


static bool isdirectory(const char *path) {
	struct stat statbuf;
	if(stat(path, &statbuf) != 0) {
		perror("stat");
		return true;
	}
	return S_ISDIR(statbuf.st_mode);
}

static char *optgen(const char *text, int state) {
	static int list_index, len;
	const char *name = NULL;

	if (!state) {
		list_index = 0;
		len = (int)strlen(text);
	}

	name = rl_filename_completion_function(text, state);
	if (name != NULL) {
		rl_completion_append_character = isdirectory(name) ? '/': ' ';
		return strdup(name);
	}

	if (text[0] == '-') {
		while ((name = options[list_index++])) {
			if (strncmp(name, text, len) == 0) {
				return strdup(name);
			}
		}
	}
	return NULL;
}

static char **optcomplete(const char *text, int start, int end) {
	rl_attempted_completion_over = 1;
	return rl_completion_matches(text, optgen);
}

static void rlinit() {
	rl_attempted_completion_function = optcomplete;
}

static void handler(int sig) {
	done = true;
}

static int _getline(istream *in, char *buf, size_t len) {
	if (in == &cin) {
		char *input = readline("bowtie2> ");
		if (!input) {
			buf[0] = '\0';
		}
		else {
			strncpy(buf, input, len-1);
			add_history(buf);
			free(input);
		}
	}
	else {
		in->getline(buf, len-1);
		if (!in->good())
			buf[0] = '\0';
	}
	buf[len-1] = '\0';
	return (int)strlen(buf);
}

static int createfifo(vector<char *>& v) {
	const char *pattern = "/tmp/btfifo.XX";
	char *fn = (char *)malloc(strlen(pattern)+1);
	memset(fn, 0, strlen(pattern)+1);
	strcpy(fn, pattern);
	mktemp(fn);

	if (mkfifo(fn, S_IRUSR | S_IWUSR) == -1) {
		perror("mkfifo");
		return 1;
	}
	v.push_back(fn);
	return 0; 
}

void printargs(const char **args, size_t size) {
	copy(args, args+size, ostream_iterator<string>(cout, " "));
	cout << endl;
}

bool called_from_wrapper(int argc, const char **argv) {
	if (argc > 2 && strcmp(argv[1], "--wrapper") == 0
			&& strcmp(argv[2], "basic-0") == 0)
		return true;
	return false;
}

int main(int argc, const char **argv) {
	int offset = called_from_wrapper(argc, argv) ? 3 : 1;

	if(argc > offset + 1 && strcmp(argv[offset], "-A") == 0) {
		const char *file = argv[offset+1];
		ifstream in;
		istream *inptr = &in;
		if (strcmp(file, "-") == 0) {
			inptr = &cin;
		}
		else {
			in.open(file);
		}
		char buf[4096];
		int lastret = -1;

		rlinit();
		while(_getline(inptr, buf, 4096)) {
			done = false;
			vector<string> args;
			args.push_back(string(argv[0]));

			if (offset > 1) {
				args.push_back(string(argv[1]));
				args.push_back(string(argv[2]));
			}

			tokenize(buf, " \t", args);
			const char **myargs = (const char**)malloc(sizeof(char*)*args.size());
			vector<char *> fifonames;
			int sam_outfile_pos = -1;

			for(size_t i = 0; i < args.size(); i++) {
				if (args[i] == "_") {
					if (i > 0 && args[i-1] == "-S") {
						sam_outfile_pos = (int)i;
					}
					else {
						createfifo(fifonames);
						args[i] = fifonames.back();
					}
				}
				myargs[i] = args[i].c_str();
			}
			if(args.size() == 1) continue;

			if (fifonames.size() > 0) {
				struct sigaction sa;
				sigemptyset(&sa.sa_mask);
				sa.sa_flags = 0;
				sa.sa_handler = handler;

				if (sigaction(SIGINT, &sa, NULL) == -1
						|| signal(SIGPIPE, SIG_IGN) == SIG_ERR) {
					perror("sigaction");
					exit(EXIT_FAILURE);
				}

				struct pollfd *pollfds = (struct pollfd *)calloc(fifonames.size(), sizeof(struct pollfd));
				if (pollfds == NULL) {
					perror("calloc");
					exit(EXIT_FAILURE);
				}

				for (size_t i = 0; i < fifonames.size(); i++) {
					pollfds[i].fd = open(fifonames[i], O_NONBLOCK | O_EXCL);
					pollfds[i].events = POLLIN;
				}

				printargs(myargs, args.size());

				for (int count = 0;; count++) {
					size_t r = poll(pollfds, (nfds_t)fifonames.size(), -1);
					// wait until all fifos are ready
					if (r != fifonames.size()) {
						if (done)
							break;
						continue;
					}
					if (sam_outfile_pos >= 0) {
						ostringstream os;
						os << "out" << getpid() << "_" << setfill('0') << setw(3) << count << ".sam"; 
						myargs[sam_outfile_pos] = os.str().c_str();
						printargs(myargs, args.size());
					}

					lastret = bowtie((int)args.size(), myargs);

					// replace the args shuffled by getopt
					for (size_t i = 0; i < args.size(); i++) {
						myargs[i] = args[i].c_str();
					}
				}

				for (size_t i = 0; i < fifonames.size(); i++) {
					if (close(pollfds[i].fd))
						perror("close");
					if (remove(fifonames[i]))
						perror("remove");
					free(fifonames[i]);
				}
				free(pollfds);
			}
			else {
				lastret = bowtie((int)args.size(), myargs);
			}
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
