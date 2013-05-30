/*
 * Copyright 2013, Ben Langmead <langmea@cs.jhu.edu>
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

#include <getopt.h>
#include "assert_helpers.h"
#include "ds.h"
#include "simple_func.h"
#include "aligner_seed_policy.h"
#include "scoring.h"
#include "opts.h"
#include "aligner_sw.h"

using namespace std;

int gVerbose;               // be talkative
int gQuiet;                 // print nothing but the alignments
static int sanityCheck;     // enable expensive sanity checks
static int seed;            // srandom() seed
static bool showVersion;    // just print version and quit?
static uint32_t qUpto;      // max # of queries to read
static int nthreads;        // number of pthreads operating concurrently
static bool useSpinlock;    // false -> don't use of spinlocks even if they're #defines
static uint32_t skipReads;  // # reads/read pairs to skip
int gGapBarrier;            // # diags on top/bot only to be entered diagonally
static int bonusMatchType;  // how to reward matches
static int bonusMatch;      // constant reward if bonusMatchType=constant
static int penMmcType;      // how to penalize mismatches
static int penMmcMax;       // max mm penalty
static int penMmcMin;       // min mm penalty
static int penNType;        // how to penalize Ns in the read
static int penN;            // constant if N pelanty is a constant
static bool penNCatPair;    // concatenate mates before N filtering?
static bool localAlign;     // do local alignment in DP steps
static int   penRdGapConst;   // constant cost of extending a gap in the read
static int   penRfGapConst;   // constant cost of extending a gap in the reference
static int   penRdGapLinear;  // coeff of linear term for cost of gap extension in read
static int   penRfGapLinear;  // coeff of linear term for cost of gap extension in ref
static SimpleFunc scoreMin;   // minimum valid score as function of read len
static SimpleFunc nCeil;      // max # Ns allowed as function of read len
static SimpleFunc msIval;     // interval between seeds as function of read len
static bool enable8;          // use 8-bit SSE where possible?
static size_t cminlen;        // longer reads use checkpointing
static size_t cpow2;          // checkpoint interval log2
static bool doTri;            // do triangular mini-fills?
static bool ignoreQuals;      // all mms incur same penalty, regardless of qual
static EList<string> queries; // list of query files
static string outfile;        // write output to this file

static void resetOptions() {
	gVerbose                = 0;
	gQuiet					= false;
	sanityCheck				= 0;  // enable expensive sanity checks
	seed					= 0; // srandom() seed
	showVersion				= false; // just print version and quit?
	qUpto					= 0xffffffff; // max # of queries to read
	nthreads				= 1;     // number of pthreads operating concurrently
	useSpinlock				= true;  // false -> don't use of spinlocks even if they're #defines
	skipReads				= 0;     // # reads/read pairs to skip
	gGapBarrier				= 4;     // disallow gaps within this many chars of either end of alignment
	bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch      = DEFAULT_MATCH_BONUS;
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmcMax       = DEFAULT_MM_PENALTY_MAX;
	penMmcMin       = DEFAULT_MM_PENALTY_MIN;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	penNCatPair     = DEFAULT_N_CAT_PAIR; // concatenate mates before N filtering?
	localAlign      = false;     // do local alignment in DP steps
	penRdGapConst   = DEFAULT_READ_GAP_CONST;
	penRfGapConst   = DEFAULT_REF_GAP_CONST;
	penRdGapLinear  = DEFAULT_READ_GAP_LINEAR;
	penRfGapLinear  = DEFAULT_REF_GAP_LINEAR;
	scoreMin.init  (SIMPLE_FUNC_LINEAR, DEFAULT_MIN_CONST,   DEFAULT_MIN_LINEAR);
	nCeil.init     (SIMPLE_FUNC_LINEAR, 0.0f, std::numeric_limits<double>::max(), 2.0f, 0.1f);
	msIval.init    (SIMPLE_FUNC_LINEAR, 1.0f, std::numeric_limits<double>::max(), DEFAULT_IVAL_B, DEFAULT_IVAL_A);
	enable8            = true;  // use 8-bit SSE where possible?
	cminlen            = 2000;  // longer reads use checkpointing
	cpow2              = 4;     // checkpoint interval log2
	doTri              = false; // do triangular mini-fills?
	ignoreQuals = false;     // all mms incur same penalty, regardless of qual
	queries.clear();         // list of query files
	outfile.clear();         // write output to this file
}

static const char *short_options = "u:hp:P:S:";

static struct option long_options[] = {
	{(char*)"verbose",          no_argument,       0, ARG_VERBOSE},
	{(char*)"quiet",            no_argument,       0, ARG_QUIET},
	{(char*)"sanity",           no_argument,       0, ARG_SANITY},
	{(char*)"qupto",            required_argument, 0, 'u'},
	{(char*)"upto",             required_argument, 0, 'u'},
	{(char*)"version",          no_argument,       0, ARG_VERSION},
	{(char*)"help",             no_argument,       0, 'h'},
	{(char*)"threads",          required_argument, 0, 'p'},
	{(char*)"usage",            no_argument,       0, ARG_USAGE},
	{(char*)"gbar",             required_argument, 0, ARG_GAP_BAR},
	{(char*)"policy",           required_argument, 0, ARG_ALIGN_POLICY},
	{(char*)"454",              no_argument,       0, ARG_NOISY_HPOLY},
	{(char*)"ion-torrent",      no_argument,       0, ARG_NOISY_HPOLY},
	{(char*)"local",            no_argument,       0, ARG_LOCAL},
	{(char*)"end-to-end",       no_argument,       0, ARG_END_TO_END},
	{(char*)"sse8",             no_argument,       0, ARG_SSE8},
	{(char*)"no-sse8",          no_argument,       0, ARG_SSE8_NO},
	{(char*)"ma",               required_argument, 0, ARG_SCORE_MA},
	{(char*)"mp",               required_argument, 0, ARG_SCORE_MMP},
	{(char*)"np",               required_argument, 0, ARG_SCORE_NP},
	{(char*)"rdg",              required_argument, 0, ARG_SCORE_RDG},
	{(char*)"rfg",              required_argument, 0, ARG_SCORE_RFG},
	{(char*)"score-min",        required_argument, 0, ARG_SCORE_MIN},
	{(char*)"min-score",        required_argument, 0, ARG_SCORE_MIN},
	{(char*)"n-ceil",           required_argument, 0, ARG_N_CEIL},
	{(char*)"ignore-quals",     no_argument,       0, ARG_IGNORE_QUALS},
	{(char*)"output",           required_argument, 0, 'S'},
	{(char*)"cp-min",           required_argument, 0, ARG_CP_MIN},
	{(char*)"cp-ival",          required_argument, 0, ARG_CP_IVAL},
	{(char*)"tri",              no_argument,       0, ARG_TRI},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Bowtie 2 dynamic programming engine, by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)" << endl;
	string tool_name = "bowtie2-dp";
	out << "Usage: " << endl
	    << "  " << tool_name.c_str() << " [options]* <in> <out>" << endl
	    << endl
	    <<     "  <in>           File with DP input problems (default: stdin)" << endl
	    <<     "  <out>          File with DP output solutions (default: stdout)" << endl
		<< endl
	    << "Options (defaults in parentheses):" << endl
		<< endl
	    << " Input:" << endl
	    << "  -s/--skip <int>    skip the first <int> problems in the input (none)" << endl
	    << "  -u/--upto <int>    stop after first <int> problems (no limit)" << endl
		<< endl
	    << " Alignment:" << endl
		<< "  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)" << endl
		<< "  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)" << endl
		<< "  --ignore-quals     treat all quality values as 30 on Phred scale (off)" << endl
		<< endl
		<< "  --end-to-end       entire read must align; no clipping (on)" << endl
		<< "   OR" << endl
		<< "  --local            local alignment; ends might be soft clipped (off)" << endl
		<< endl
	    << " Scoring:" << endl
		<< "  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) " << endl
		<< "  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)" << endl
		<< "  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)" << endl
		<< "  --rdg <int>,<int>  read gap open, extend penalties (5,3)" << endl
		<< "  --rfg <int>,<int>  reference gap open, extend penalties (5,3)" << endl
		<< "  --score-min <func> min acceptable alignment score w/r/t read length" << endl
		<< "                     (G,20,8 for local, L,-0.6,-0.6 for end-to-end)" << endl
	    << "  --quiet            print nothing to stderr except serious errors" << endl
		<< endl
	    << " Performance:" << endl
	    << "  -p/--threads <int> number of alignment threads to launch (1)" << endl
		<< endl
	    << " Other:" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print this usage message" << endl
	    ;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, int upper, const char *errmsg, const char *arg) {
	long l;
	char *endPtr= NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower || l > upper) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Upper is maximum int by default.
 */
static int parseInt(int lower, const char *errmsg, const char *arg) {
	return parseInt(lower, std::numeric_limits<int>::max(), errmsg, arg);
}

/**
 * Parse a T string 'str'.
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

static int parseFuncType(const std::string& otype) {
	string type = otype;
	if(type == "C" || type == "Constant") {
		return SIMPLE_FUNC_CONST;
	} else if(type == "L" || type == "Linear") {
		return SIMPLE_FUNC_LINEAR;
	} else if(type == "S" || type == "Sqrt") {
		return SIMPLE_FUNC_SQRT;
	} else if(type == "G" || type == "Log") {
		return SIMPLE_FUNC_LOG;
	}
	std::cerr << "Error: Bad function type '" << otype.c_str()
	          << "'.  Should be C (constant), L (linear), "
	          << "S (square root) or G (natural log)." << std::endl;
	throw 1;
}

#define PARSE_FUNC(fv) { \
	if(args.size() >= 1) { \
		fv.setType(parseFuncType(args[0])); \
	} \
	if(args.size() >= 2) { \
		double co; \
		istringstream tmpss(args[1]); \
		tmpss >> co; \
		fv.setConst(co); \
	} \
	if(args.size() >= 3) { \
		double ce; \
		istringstream tmpss(args[2]); \
		tmpss >> ce; \
		fv.setCoeff(ce); \
	} \
	if(args.size() >= 4) { \
		double mn; \
		istringstream tmpss(args[3]); \
		tmpss >> mn; \
		fv.setMin(mn); \
	} \
	if(args.size() >= 5) { \
		double mx; \
		istringstream tmpss(args[4]); \
		tmpss >> mx; \
		fv.setMin(mx); \
	} \
}

/**
 * TODO: Argument parsing is very, very flawed.  The biggest problem is that
 * there are two separate worlds of arguments, the ones set via polstr, and
 * the ones set directly in variables.  This makes for nasty interactions,
 * e.g., with the -M option being resolved at an awkward time relative to
 * the -k and -a options.
 */
static void parseOption(int next_option, const char *arg) {
	switch (next_option) {
		case 's':
			skipReads = (uint32_t)parseInt(0, "-s arg must be positive", arg);
			break;
		case ARG_GAP_BAR:
			gGapBarrier = parseInt(1, "--gbar must be no less than 1", arg);
			break;
		case 'u':
			qUpto = (uint32_t)parseInt(1, "-u/--qupto arg must be at least 1", arg);
			break;
		case 'p':
			nthreads = parseInt(1, "-p/--threads arg must be at least 1", arg);
			break;
		case 'h': printUsage(cout); throw 0; break;
		case ARG_USAGE: printUsage(cout); throw 0; break;
		case ARG_VERBOSE: gVerbose = 1; break;
		case ARG_QUIET: gQuiet = true; break;
		case ARG_SANITY: sanityCheck = true; break;
		case ARG_CP_MIN:
			cminlen = parse<size_t>(arg);
			break;
		case ARG_CP_IVAL:
			cpow2 = parse<size_t>(arg);
			break;
		case ARG_TRI:
			doTri = true;
			break;
		case ARG_LOCAL: localAlign = true; break;
		case ARG_END_TO_END: localAlign = false; break;
		case ARG_SSE8: enable8 = true; break;
		case ARG_SSE8_NO: enable8 = false; break;
		case ARG_IGNORE_QUALS: ignoreQuals = true; break;
		case ARG_N_CEIL: {
			// Split argument by comma
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3) {
				cerr << "Error: expected 3 or fewer comma-separated "
					 << "arguments to --n-ceil option, got "
					 << args.size() << endl;
				throw 1;
			}
			if(args.size() == 0) {
				cerr << "Error: expected at least one argument to --n-ceil option" << endl;
				throw 1;
			}
			PARSE_FUNC(nCeil);
			break;
		}
		case ARG_SCORE_MA: {
			// Split argument by comma
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() != 1) {
				cerr << "Error parsing --ma; RHS must have 1 token" << endl;
				assert(false); throw 1;
			}
			string tmp = args[0];
			istringstream tmpss(tmp);
			tmpss >> bonusMatch;
			break;
		}
		case ARG_SCORE_MMP: {
			// Split argument by comma
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3) {
				cerr << "Error parsing --mmp "
				     << "; RHS must have at most 3 tokens" << endl;
				assert(false); throw 1;
			}
			if(args[0][0] == 'C') {
				string tmp = args[0].substr(1);
				// Parse constant penalty
				istringstream tmpss(tmp);
				tmpss >> penMmcMax;
				penMmcMin = penMmcMax;
				// Parse constant penalty
				penMmcType = COST_MODEL_CONSTANT;
			} else if(args[0][0] == 'Q') {
				if(args.size() >= 2) {
					string tmp = args[1];
					istringstream tmpss(tmp);
					tmpss >> penMmcMax;
				} else {
					penMmcMax = DEFAULT_MM_PENALTY_MAX;
				}
				if(args.size() >= 3) {
					string tmp = args[2];
					istringstream tmpss(tmp);
					tmpss >> penMmcMin;
				} else {
					penMmcMin = DEFAULT_MM_PENALTY_MIN;
				}
				if(penMmcMin > penMmcMax) {
					cerr << "Error: Maximum mismatch penalty (" << penMmcMax
					     << ") is less than minimum penalty (" << penMmcMin
						 << endl;
					throw 1;
				}
				// Set type to =quality
				penMmcType = COST_MODEL_QUAL;
			} else if(args[0][0] == 'R') {
				// Set type to=Maq-quality
				penMmcType = COST_MODEL_ROUNDED_QUAL;
			} else {
				cerr << "Error parsing --mmp "
				     << "; RHS must start with C, Q or R" << endl;
				assert(false); throw 1;
			}
			break;
		}
		case ARG_SCORE_NP: {
			// Split argument by comma
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() != 1) {
				cerr << "Error parsing --np "
				     << "; RHS must have 1 token" << endl;
				assert(false); throw 1;
			}
			if(args[0][0] == 'C') {
				string tmp = args[0].substr(1);
				// Parse constant penalty
				istringstream tmpss(tmp);
				tmpss >> penN;
				// Parse constant penalty
				penNType = COST_MODEL_CONSTANT;
			} else if(args[0][0] == 'Q') {
				// Set type to =quality
				penNType = COST_MODEL_QUAL;
			} else if(args[0][0] == 'R') {
				// Set type to=Maq-quality
				penNType = COST_MODEL_ROUNDED_QUAL;
			} else {
				cerr << "Error parsing --np "
				     << "; RHS must start with C, Q or R" << endl;
				assert(false); throw 1;
			}
		}
		case ARG_SCORE_RDG: {
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() >= 1) {
				istringstream tmpss(args[0]);
				tmpss >> penRdGapConst;
			} else {
				penRdGapConst = DEFAULT_READ_GAP_CONST;
			}
			if(args.size() >= 2) {
				istringstream tmpss(args[1]);
				tmpss >> penRdGapLinear;
			} else {
				penRdGapLinear = DEFAULT_READ_GAP_LINEAR;
			}
		}
		case ARG_SCORE_RFG: {
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() >= 1) {
				istringstream tmpss(args[0]);
				tmpss >> penRfGapConst;
			} else {
				penRfGapConst = DEFAULT_REF_GAP_CONST;
			}
			if(args.size() >= 2) {
				istringstream tmpss(args[1]);
				tmpss >> penRfGapLinear;
			} else {
				penRfGapLinear = DEFAULT_REF_GAP_LINEAR;
			}
		}
		case ARG_SCORE_MIN: {
			EList<string> args;
			tokenize(arg, ",", args);
			if(args.size() > 3 && args.size() == 0) {
				cerr << "Error: expected 3 or fewer comma-separated "
					 << "arguments to --n-ceil option, got "
					 << args.size() << endl;
				throw 1;
			}
			PARSE_FUNC(scoreMin);
			break;
		}
		case 'S': outfile = arg; break;
		case 'U': {
			EList<string> args;
			tokenize(arg, ",", args);
			for(size_t i = 0; i < args.size(); i++) {
				queries.push_back(args[i]);
			}
			break;
		}
		case ARG_VERSION: showVersion = 1; break;
		default:
			printUsage(cerr);
			throw 1;
	}
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	while(true) {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		const char * arg = optarg;
		if(next_option == EOF) {
			break;
		}
		parseOption(next_option, arg);
	}
	// If both -s and -u are used, we need to adjust qUpto accordingly
	// since it uses rdid to know if we've reached the -u limit (and
	// rdids are all shifted up by skipReads characters)
	if(qUpto + skipReads > qUpto) {
		qUpto += skipReads;
	}
	if(gGapBarrier < 1) {
		cerr << "Warning: --gbar was set less than 1 (=" << gGapBarrier
		     << "); setting to 1 instead" << endl;
		gGapBarrier = 1;
	}
#ifndef NDEBUG
	if(!gQuiet) {
		cerr << "Warning: Running in debug mode.  Please use debug mode only "
			 << "for diagnosing errors, and not for typical use of Bowtie 2."
			 << endl;
	}
#endif
}

struct DpProblem {
	void reset() {
		ref.clear();
	}

	TRefId   refidx;
	TRefOff  reflen;
	TAlScore minsc;
	BTString ref;
	bool     fw;
	DPRect   rect;
	bool     aligned;
	TAlScore score;
};

class DpLogReader {

public:

	DpLogReader() { }
	
	~DpLogReader() { reset(); }
	
	void init(const string& fn) {
		reset();
		fn_ = fn;
		ih_.open(fn_.c_str());
		ih_.sync_with_stdio(false);
	}
	
	void reset() {
		if(ih_.is_open()) {
			ih_.close();
		}
	}
	
	bool nextRead(
		BTDnaString& seq,
		BTString& qual,
		EList<DpProblem>& refs)
	{
		if(done()) {
			return false;
		}
		ln_.clear();
		getline(ih_, ln_);
		while(ln_.empty() && ih_.good()) {
			getline(ih_, ln_);
		}
		if(ln_.empty() && !ih_.good()) {
			return false;
		}
		EList<string> buf;
		tokenize(ln_, '\t', buf);
		assert_gt(buf.size(), 2);
		seq.install(buf[0].c_str(), true);
		qual = buf[1];
		for(size_t i = 2; i < buf.size(); i++) {
			refs.expand();
			istringstream is(buf[i]);
			char comma, tmp;
			// ref idx
			is >> refs.back().refidx;
			is >> comma; assert_eq(',', comma);
			// ref length
			is >> refs.back().reflen;
			is >> comma; assert_eq(',', comma);
			// minimum score
			is >> refs.back().minsc;
			is >> comma; assert_eq(',', comma);
			// read orientation
			is >> tmp;
			assert(tmp == '-' || tmp == '+');
			refs.back().fw = (tmp == '+');
			is >> comma; assert_eq(',', comma);
			// DP rectangle
			is >> refs.back().rect.refl;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.refr;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.refl_pretrim;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.refr_pretrim;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.triml;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.trimr;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.corel;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.corer;
			is >> comma; assert_eq(',', comma);
			is >> refs.back().rect.maxgap;
			is >> comma; assert_eq(',', comma);
			// reference string
			string ref;
			while(true) {
				char c;
				is >> c;
				if(c == ',') break;
				ref.push_back(c);
			}
			refs.back().ref.install(ref.c_str());
			for(size_t i = 0; i < ref.length(); i++) {
				int m = asc2dnamask[(int)refs.back().ref[i]];
				if(m == 15) {
					m = 16; // N
				}
				refs.back().ref.set(m, i);
			}
			// whether the DP alignment was successful
			int aligned;
			is >> aligned;
			refs.back().aligned = (aligned == 1);
			// alignment score
			is >> comma; assert_eq(',', comma);
			is >> refs.back().score;
		}
		return true;
	}
	
	bool done() const {
		return !ih_.good();
	}

protected:

	string   fn_; // file name
	ifstream ih_; // file handle
	string   ln_; // line buffer
};

int main(int argc, const char **argv) {
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();
		parseOptions(argc, argv);
		if(showVersion) {
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
				 << sizeof(int)
				 << ", " << sizeof(long) << ", " << sizeof(long long)
				 << ", " << sizeof(void *) << ", " << sizeof(size_t)
				 << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}
		while(optind < argc) {
			queries.push_back(argv[optind++]);
		}
		{
			// Optionally summarize
			if(gVerbose) {
				cout << "DP inputs:" << endl;
				for(size_t i = 0; i < queries.size(); i++) {
					cout << "  " << queries[i].c_str() << endl;
				}
				cout << "Output file: \"" << outfile.c_str() << "\"" << endl;
				cout << "Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
			#ifdef NDEBUG
				cout << "Assertions: disabled" << endl;
			#else
				cout << "Assertions: enabled" << endl;
			#endif
			}
		}
		// Do stuff
		SwAligner sw(NULL);
		DpLogReader logrd;
		Scoring sc(
			bonusMatch,     // constant reward for match
			penMmcType,     // how to penalize mismatches
			penMmcMax,      // max mm pelanty
			penMmcMin,      // min mm pelanty
			scoreMin,       // min score as function of read len
			nCeil,          // max # Ns as function of read len
			penNType,       // how to penalize Ns in the read
			penN,           // constant if N pelanty is a constant
			penNCatPair,    // whether to concat mates before N filtering
			penRdGapConst,  // constant coeff for read gap cost
			penRfGapConst,  // constant coeff for ref gap cost
			penRdGapLinear, // linear coeff for read gap cost
			penRfGapLinear, // linear coeff for ref gap cost
			gGapBarrier);   // # rows at top/bot only entered diagonally
		RandomSource rnd(seed);
		{
			Timer tim(std::cerr, "Alignment ", true);
			BTDnaString seq, seqrc;
			BTString qual, qualrc;
			EList<DpProblem> probs;
			size_t qid = 0;
			size_t totnuc = 0, totcup = 0;
			for(size_t i = 0; i < queries.size(); i++) {
				logrd.init(queries[i]);
				while(logrd.nextRead(seq, qual, probs)) {
					totnuc += seq.length();
					seqrc = seq;
					seqrc.reverseComp();
					qualrc = qual;
					qualrc.reverse();
					//cerr << "Initing read with " << probs.size() << " problems" << endl;
					sw.initRead(seq, seqrc, qual, qualrc, 0, seq.length(), sc);
					// Calculate minimum score
					bool extend = true;
					for(size_t j = 0; j < probs.size(); j++) {
						sw.initRef(
							probs[j].fw,
							probs[j].refidx,
							probs[j].rect,
							const_cast<char *>(probs[j].ref.toZBuf()),
							0,
							probs[j].ref.length(),
							probs[j].reflen,
							sc,
							probs[j].minsc,
							enable8,
							cminlen,
							cpow2,
							doTri,
							extend);
						// Now fill the dynamic programming matrix and return true iff
						// there is at least one valid alignment
						TAlScore bestCell = std::numeric_limits<TAlScore>::min();
						ASSERT_ONLY(bool aligned =) sw.align(bestCell);
						assert(aligned == probs[j].aligned);
						assert(!aligned || bestCell == probs[j].score);
						totcup += (seq.length() * probs[j].ref.length());
					}
					seq.clear();  seqrc.clear();
					qual.clear(); qualrc.clear();
					probs.clear();
					qid++;
				}
			}
			size_t el = (size_t)tim.elapsed();
			double cups = 0.0;
			double totnucps = 0.0;
			double readps = 0.0;
			if(el > 0) {
				cups = totcup / (double)el;
				totnucps = totnuc / (double)el;
				readps = qid / (double)el;
			}
			cerr << qid << " reads" << endl;
			cerr << std::setprecision(4) << "  " << readps << " reads per second" << endl;
			cerr << totnuc << " nucleotides" << endl;
			cerr << std::setprecision(4) << "  " << totnucps << " nucleotides per second" << endl;
			cerr << totcup << " cell updates" << endl;
			cerr << std::setprecision(4) << "  " << cups << " cell updates per second (CUPS)" << endl;
		}
		return 0;
	} catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal Bowtie 2 exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
}

