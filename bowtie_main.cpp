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
#ifdef _MSC_VER
#include <intrin.h>
#include <immintrin.h>
#endif
#endif

using namespace std;

extern "C" {
	int bowtie(int argc, const char **argv);
}

#ifdef ENABLE_x86_64_v3
// The -v256 binaries are compiled with -march=x86-64-v3: x86-64-v2 (SSE3,
// SSSE3, SSE4.1, SSE4.2, POPCNT, CMPXCHG16B, LAHF/SAHF) plus AVX, AVX2, BMI1,
// BMI2, F16C, FMA, LZCNT, MOVBE and XSAVE.  The OS must also save the XMM and
// YMM registers (OSXSAVE, then XCR0 bits 1 and 2), as __builtin_cpu_supports
// ("x86-64-v3") requires; that builtin is not available on every compiler.
static const unsigned X86_64_V3_LEAF1_ECX = // CPUID leaf 1, ECX
	(1u << 0) | (1u << 9) | (1u << 12) | (1u << 13) | // SSE3, SSSE3, FMA, CMPXCHG16B
	(1u << 19) | (1u << 20) | (1u << 22) | (1u << 23) | // SSE4.1, SSE4.2, MOVBE, POPCNT
	(1u << 26) | (1u << 27) | (1u << 28) | (1u << 29); // XSAVE, OSXSAVE, AVX, F16C
static const unsigned X86_64_V3_LEAF7_EBX = // CPUID leaf 7 sub-leaf 0, EBX
	(1u << 3) | (1u << 5) | (1u << 8); // BMI1, AVX2, BMI2
static const unsigned X86_64_V3_EXT1_ECX = // CPUID leaf 0x80000001, ECX
	(1u << 0) | (1u << 5); // LAHF/SAHF, LZCNT

// A leaf the CPU does not report is passed as 0.
static bool x86_64_v3_from_cpuid(unsigned leaf1_ecx, unsigned leaf7_ebx,
                                 unsigned ext1_ecx, unsigned long long xcr0) {
	return (leaf1_ecx & X86_64_V3_LEAF1_ECX) == X86_64_V3_LEAF1_ECX &&
	       (leaf7_ebx & X86_64_V3_LEAF7_EBX) == X86_64_V3_LEAF7_EBX &&
	       (ext1_ecx & X86_64_V3_EXT1_ECX) == X86_64_V3_EXT1_ECX &&
	       (xcr0 & 6) == 6;
}

static void x86_cpuid(unsigned leaf, unsigned subleaf, unsigned regs[4]) {
#ifdef _MSC_VER
	int r[4];
	__cpuidex(r, (int)leaf, (int)subleaf);
	for (int i = 0; i < 4; i++) regs[i] = (unsigned)r[i];
#else
	asm volatile("cpuid" : "=a" (regs[0]), "=b" (regs[1]), "=c" (regs[2]), "=d" (regs[3])
	                     : "a" (leaf), "c" (subleaf));
#endif
}

static bool x86_64_v3_supported(void) {
	unsigned r[4];
	x86_cpuid(0, 0, r);
	const unsigned max_leaf = r[0];
	x86_cpuid(0x80000000u, 0, r);
	const unsigned max_ext_leaf = r[0];
	unsigned leaf1_ecx = 0, leaf7_ebx = 0, ext1_ecx = 0;
	unsigned long long xcr0 = 0;
	if (max_leaf >= 1) { x86_cpuid(1, 0, r); leaf1_ecx = r[2]; }
	if (max_leaf >= 7) { x86_cpuid(7, 0, r); leaf7_ebx = r[1]; }
	if (max_ext_leaf >= 0x80000001u) { x86_cpuid(0x80000001u, 0, r); ext1_ecx = r[2]; }
	if ((leaf1_ecx >> 27) & 1) { // XGETBV faults unless the OS set OSXSAVE
#ifdef _MSC_VER
		xcr0 = _xgetbv(0);
#else
		unsigned lo, hi;
		asm volatile("xgetbv" : "=a" (lo), "=d" (hi) : "c" (0));
		xcr0 = ((unsigned long long)hi << 32) | lo;
#endif
	}
	return x86_64_v3_from_cpuid(leaf1_ecx, leaf7_ebx, ext1_ecx, xcr0);
}

void check_x86_64_v3(int argc, const char **argv) {
	if (x86_64_v3_supported() && (argc<126)) {
		const char* new_argv[128]; // should always be enough, but above check enforces it, too
		const char * org_path = argv[0];
		// Append -v256 to the original path 
		const int fn_len = strlen(org_path);
		char *new_path = (char*) malloc(fn_len+16);
		memcpy(new_path,org_path,fn_len);
		strncpy(new_path+fn_len,"-v256",15);

		for (int i=1; i<=argc; i++) new_argv[i] = argv[i]; // all but first the same, also copy final NULL
		// keep argv[0], so --version and the @PG line name the binary that was run,
		// whichever variant does the work
		new_argv[0] = org_path;
		// now replace the executable with new variant
		// assuming the executable exists... execvp will gracefully fail else, which is OK
		execvp(new_path, (char *const *)new_argv);
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
