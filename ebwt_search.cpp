#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <getopt.h>
#include <pthread.h>
#include <math.h>
#include <utility>
#include "alphabet.h"
#include "assert_helpers.h"
#include "endian_swap.h"
#include "ebwt.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "hit.h"
#include "aln_sink.h"
#include "pat.h"
#include "bitset.h"
#include "threading.h"
#include "range_cache.h"
#include "refmap.h"
#include "annot.h"
#include "ds.h"
#include "aligner.h"
#include "aligner_0mm.h"
#include "aligner_1mm.h"
#include "aligner_23mm.h"
#include "aligner_seed_mm.h"
#include "aligner_metrics.h"
#include "sam.h"
#include "aligner_seed.h"
#include "aligner_seed_policy.h"
#include "aligner_sw.h"
#include "aligner_sw_driver.h"
#include "aligner_counters.h"
#include "aligner_cache.h"
#include "util.h"
#include "pe.h"
#ifdef CHUD_PROFILING
#include <CHUD/CHUD.h>
#endif

using namespace std;

static EList<string> mates1;  // mated reads (first mate)
static EList<string> mates2;  // mated reads (second mate)
static EList<string> mates12; // mated reads (1st/2nd interleaved in 1 file)
static string adjustedEbwtFileBase;
int gVerbose;      // be talkative
static bool startVerbose; // be talkative at startup
int gQuiet;        // print nothing but the alignments
static int sanityCheck;   // enable expensive sanity checks
static int format;        // default read format is FASTQ
static string origString; // reference text, or filename(s)
static int seed;          // srandom() seed
static int timing;        // whether to report basic timing data
static int metricsIval;   // interval between alignment metrics messages (0 = no messages)
static string metricsFile;// output file to put alignment metrics in
static bool metricsStderr;// output file to put alignment metrics in
static bool allHits;      // for multihits, report just one
bool gRangeMode;    // report BWT ranges instead of ref locs
static int showVersion;   // just print version and quit?
static int ipause;        // pause before maching?
static uint32_t qUpto;    // max # of queries to read
int gTrim5;         // amount to trim from 5' end
int gTrim3;         // amount to trim from 3' end
static int reportOpps;    // whether to report # of other mappings
static int offRate;       // keep default offRate
static int mismatches;    // allow 0 mismatches by default
static char *patDumpfile; // filename to dump patterns to
static bool solexaQuals;  // quality strings are solexa quals, not phred, and subtract 64 (not 33)
static bool phred64Quals; // quality chars are phred, but must subtract 64 (not 33)
static bool integerQuals; // quality strings are space-separated strings of integers, not ASCII
static int maqLike;       // do maq-like searching
static int seedLen;       // seed length (changed in Maq 0.6.4 from 24)
static int seedMms;       // # mismatches allowed in seed (maq's -n)
static int qualThresh;    // max qual-weighted hamming dist (maq's -e)
static int maxBtsBetter;  // max # backtracks allowed in half-and-half mode
static int maxBts;        // max # backtracks allowed in half-and-half mode
static int nthreads;      // number of pthreads operating concurrently
static output_types outType;  // style of output
static int numRandomReads;    // # random reads (see Random*PatternSource in pat.h)
static int lenRandomReads;    // len of random reads (see Random*PatternSource in pat.h)
static bool noRefNames;       // true -> print reference indexes; not names
static string dumpAlBase;     // basename of same-format files to dump aligned reads to
static string dumpUnalBase;   // basename of same-format files to dump unaligned reads to
static string dumpMaxBase;    // basename of same-format files to dump reads with more than -m valid alignments to
static uint32_t khits;  // number of hits per read; >1 is much slower
static uint32_t mhits;  // don't report any hits if there are > mhits
static bool strata;     // true -> don't stop at stratum boundaries
static int partitionSz; // output a partitioning key in first field
bool gNoMaqRound; // true -> don't round quals to nearest 10 like maq
static bool useSpinlock;  // false -> don't use of spinlocks even if they're #defines
static bool fileParallel; // separate threads read separate input files in parallel
static bool useShmem;     // use shared memory to hold the index
static bool useMm;        // use memory-mapped files to hold the index
static bool mmSweep;      // sweep through memory-mapped files immediately after mapping
int gMinInsert;     // minimum insert size
int gMaxInsert;     // maximum insert size
bool gMate1fw;           // -1 mate aligns in fw orientation on fw strand
bool gMate2fw;           // -2 mate aligns in rc orientation on fw strand
bool mateFwSet;
bool gLocalAlign;      // use local alignment when searching for mate
bool gFlippedMatesOK;  // allow mates to be in wrong order
bool gDovetailMatesOK; // allow one mate to extend off the end of the other
bool gContainMatesOK;  // allow one mate to contain the other in PE alignment
bool gOlapMatesOK;     // allow mates to overlap in PE alignment
bool gExpandToFrag;    // incr max frag length to =larger mate len if necessary
bool gReportDiscordant; // find and report discordant paired-end alignments
bool gReportMixed;      // find and report unpaired alignments for paired reads
static uint32_t mixedThresh;   // threshold for when to switch to paired-end mixed mode (see aligner.h)
static uint32_t mixedAttemptLim; // number of attempts to make in "mixed mode" before giving up on orientation
static bool dontReconcileMates;  // suppress pairwise all-versus-all way of resolving mates
static uint32_t cacheLimit;      // ranges w/ size > limit will be cached
static uint32_t cacheSize;       // # words per range cache
static int offBase;              // offsets are 0-based by default, but configurable
static bool tryHard;             // set very high maxBts, mixedAttemptLim
static uint32_t skipReads;       // # reads/read pairs to skip
bool gNofw; // don't align fw orientation of read
bool gNorc; // don't align rc orientation of read
bool gStrandFix;  // attempt to fix strand bias
static bool stats; // print performance stats
static int chunkPoolMegabytes;    // max MB to dedicate to best-first search frames per thread
static int chunkSz;    // size of single chunk disbursed by ChunkPool
static bool chunkVerbose; // have chunk allocator output status messages?
bool gReportSe;
static const char * refMapFile;  // file containing a map from index coordinates to another coordinate system
static const char * annotMapFile;  // file containing a map from reference coordinates to annotations
static uint32_t fastaContLen;
static uint32_t fastaContFreq;
static int randomNum;    // number of randomly-generated reads
static int randomLength; // length of randomly-generated reads
static bool hadoopOut; // print Hadoop status and summary messages
static bool fuzzy;
static bool fullRef;
static bool samNoHead; // don't print any header lines in SAM output
static bool samNoSQ;   // don't print @SQ header lines
bool gColor;     // true -> inputs are colorspace
bool gColorExEnds; // true -> nucleotides on either end of decoded cspace alignment should be excluded
bool gReportOverhangs; // false -> filter out alignments that fall off the end of a reference sequence
static string rgs; // SAM outputs for @RG header line
int gSnpPhred; // probability of SNP, for scoring colorspace alignments
static EList<bool> suppressOuts; // output fields to suppress
static bool sampleMax; // whether to report a random alignment when maxed-out via -m/-M
static int defaultMapq; // default mapping quality to print in SAM mode
bool gColorSeq; // true -> show colorspace alignments as colors, not decoded bases
bool gColorEdit; // true -> show edits as colors, not decoded bases
bool gColorQual; // true -> show colorspace qualities as original quals, not decoded quals
static bool printPlaceholders; // true -> print records for maxed-out, unaligned reads
static bool printFlags; // true -> print alignment flags
static bool printCost; // true -> print stratum and cost
static bool printParams; // true -> print parameters like seed spacing, cost ceiling
bool gShowSeed;
bool gGaps; // true -> allow gapped alignments
uint32_t gInsOpen;   // insertion open penalty
uint32_t gDelOpen;   // deletion open penalty
uint32_t gInsExtend; // insertion gap extension penalty
uint32_t gDelExtend; // deletion gap extension penalty
int gGapBarrier;
int gAllowRedundant;
static EList<string> qualities;
static EList<string> qualities1;
static EList<string> qualities2;
static bool doMultiseed; // true -> use multiseed aligner

static string seedAlignmentPolicyString; // temporary holder for policy string
static bool  msNoCache;      // true -> disable local cache
static int   penMmcType;     // how to penalize mismatches
static int   penMmc;         // constant if mm pelanty is a constant
static int   penSnp;         // penalty for nucleotide mismatches in decoded colorspace als
static int   penNType;       // how to penalize Ns in the read
static int   penN;           // constant if N pelanty is a constant
static bool  penNCatPair;    // concatenate mates before N filtering?
static bool  noisyHpolymer;  // set to true if gap penalties should be reduced to be consistent with a sequencer that under- and overcalls homopolymers
static int   penRdGapConst;   // constant cost of extending a gap in the read
static int   penRfGapConst;   // constant cost of extending a gap in the reference
static int   penRdGapLinear;  // coeff of linear term for cost of gap extension in read
static int   penRfGapLinear;  // coeff of linear term for cost of gap extension in ref
static float costCeilConst;  // constant factor in cost ceiling w/r/t read length
static float costCeilLinear; // coeff of linear term in cost ceiling w/r/t read length
static float nCeilConst;     // constant factor in N ceiling w/r/t read length
static float nCeilLinear;    // coeff of linear term in N ceiling w/r/t read length
static int   multiseedMms;   // mismatches permitted in a multiseed seed
static int   multiseedLen;   // length of multiseed seeds
static int   multiseedPeriod;// space between multiseed seeds
static int   multiseedIvalType; // formula for seed spacing
static float multiseedIvalA; // linear coefficient in formula
static float multiseedIvalB; // constant coefficient in formula
static string saCountersFn;  // filename to dump per-read SeedAligner counters to
static string saActionsFn;   // filename to dump all alignment actions to
static string saHitsFn;      // filename to dump all seed hits to
static uint32_t seedCacheLocalMB;   // # MB to use for non-shared seed alignment cacheing
static uint32_t seedCacheCurrentMB; // # MB to use for current-read seed hit cacheing
static int maxposs;          // maximum number of positions to consider per read
static int maxrows;          // maximum number of rows to consider per position
static bool seedSummaryOnly; // print summary information about seed hits, not alignments

static void resetOptions() {
	mates1.clear();
	mates2.clear();
	mates12.clear();
	adjustedEbwtFileBase	= "";
	gVerbose                = 0;
	startVerbose			= 0;
	gQuiet					= false;
	sanityCheck				= 0;  // enable expensive sanity checks
	format					= FASTQ; // default read format is FASTQ
	origString				= ""; // reference text, or filename(s)
	seed					= 0; // srandom() seed
	timing					= 0; // whether to report basic timing data
	metricsIval				= 0; // interval between alignment metrics messages (0 = no messages)
	metricsFile             = ""; // output file to put alignment metrics in
	metricsStderr           = false; // print metrics to stderr (in addition to --metrics-file if it's specified
	allHits					= false; // for multihits, report just one
	gRangeMode				= false; // report BWT ranges instead of ref locs
	showVersion				= 0; // just print version and quit?
	ipause					= 0; // pause before maching?
	qUpto					= 0xffffffff; // max # of queries to read
	gTrim5					= 0; // amount to trim from 5' end
	gTrim3					= 0; // amount to trim from 3' end
	reportOpps				= 0; // whether to report # of other mappings
	offRate					= -1; // keep default offRate
	mismatches				= 0; // allow 0 mismatches by default
	patDumpfile				= NULL; // filename to dump patterns to
	solexaQuals				= false; // quality strings are solexa quals, not phred, and subtract 64 (not 33)
	phred64Quals			= false; // quality chars are phred, but must subtract 64 (not 33)
	integerQuals			= false; // quality strings are space-separated strings of integers, not ASCII
	maqLike					= 1;   // do maq-like searching
	seedLen					= 28;  // seed length (changed in Maq 0.6.4 from 24)
	seedMms					= 2;   // # mismatches allowed in seed (maq's -n)
	qualThresh				= 70;  // max qual-weighted hamming dist (maq's -e)
	maxBtsBetter			= 125; // max # backtracks allowed in half-and-half mode
	maxBts					= 800; // max # backtracks allowed in half-and-half mode
	nthreads				= 1;     // number of pthreads operating concurrently
	outType					= OUTPUT_FULL;  // style of output
	numRandomReads			= 50000000; // # random reads (see Random*PatternSource in pat.h)
	lenRandomReads			= 35;    // len of random reads (see Random*PatternSource in pat.h)
	noRefNames				= false; // true -> print reference indexes; not names
	dumpAlBase				= "";    // basename of same-format files to dump aligned reads to
	dumpUnalBase			= "";    // basename of same-format files to dump unaligned reads to
	dumpMaxBase				= "";    // basename of same-format files to dump reads with more than -m valid alignments to
	khits					= 1;     // number of hits per read; >1 is much slower
	mhits					= 0xffffffff; // don't report any hits if there are > mhits
	strata					= false; // true -> don't stop at stratum boundaries
	partitionSz				= 0;     // output a partitioning key in first field
	gNoMaqRound				= false; // true -> don't round quals to nearest 10 like maq
	useSpinlock				= true;  // false -> don't use of spinlocks even if they're #defines
	fileParallel			= false; // separate threads read separate input files in parallel
	useShmem				= false; // use shared memory to hold the index
	useMm					= false; // use memory-mapped files to hold the index
	mmSweep					= false; // sweep through memory-mapped files immediately after mapping
	gMinInsert				= 0;     // minimum insert size
	gMaxInsert				= 250;   // maximum insert size
	gMate1fw				= true;  // -1 mate aligns in fw orientation on fw strand
	gMate2fw				= false; // -2 mate aligns in rc orientation on fw strand
	mateFwSet				= false; // true -> user set mate1fw/mate2fw with --ff/--fr/--rf
	gLocalAlign             = false; // use local alignment when searching for mate
	gFlippedMatesOK         = false; // allow mates to be in wrong order
	gDovetailMatesOK        = true;  // allow one mate to extend off the end of the other
	gContainMatesOK         = true;  // allow one mate to contain the other in PE alignment
	gOlapMatesOK            = true;  // allow mates to overlap in PE alignment
	gExpandToFrag           = true;  // incr max frag length to =larger mate len if necessary
	gReportDiscordant       = true;  // find and report discordant paired-end alignments
	gReportMixed            = true;  // find and report unpaired alignments for paired reads
	mixedThresh				= 4;     // threshold for when to switch to paired-end mixed mode (see aligner.h)
	mixedAttemptLim			= 100;   // number of attempts to make in "mixed mode" before giving up on orientation
	dontReconcileMates		= true;  // suppress pairwise all-versus-all way of resolving mates
	cacheLimit				= 5;     // ranges w/ size > limit will be cached
	cacheSize				= 0;     // # words per range cache
	offBase					= 0;     // offsets are 0-based by default, but configurable
	tryHard					= false; // set very high maxBts, mixedAttemptLim
	skipReads				= 0;     // # reads/read pairs to skip
	gNofw					= false; // don't align fw orientation of read
	gNorc					= false; // don't align rc orientation of read
	gStrandFix				= true;  // attempt to fix strand bias
	stats					= false; // print performance stats
	chunkPoolMegabytes		= 32;    // max MB to dedicate to best-first search frames per thread
	chunkSz					= 256;   // size of single chunk disbursed by ChunkPool
	chunkVerbose			= false; // have chunk allocator output status messages?
	gReportSe				= false;
	refMapFile				= NULL;  // file containing a map from index coordinates to another coordinate system
	annotMapFile			= NULL;  // file containing a map from reference coordinates to annotations
	fastaContLen			= 0;
	fastaContFreq			= 0;
    randomNum               = 100;   // number of randomly-generated reads
    randomLength            = 35;    // length of randomly-generated reads
	hadoopOut				= false; // print Hadoop status and summary messages
	fuzzy					= false; // reads will have alternate basecalls w/ qualities
	fullRef					= false; // print entire reference name instead of just up to 1st space
	samNoHead				= false; // don't print any header lines in SAM output
	samNoSQ					= false; // don't print @SQ header lines
	gColor					= false; // don't align in colorspace by default
	gColorExEnds			= true;  // true -> nucleotides on either end of decoded cspace alignment should be excluded
	gReportOverhangs        = false; // false -> filter out alignments that fall off the end of a reference sequence
	rgs						= "";    // SAM outputs for @RG header line
	gSnpPhred				= 30;    // probability of SNP, for scoring colorspace alignments
	suppressOuts.clear();            // output fields to suppress
	suppressOuts.resize(64);
	suppressOuts.fill(false);
	sampleMax				= false;
	defaultMapq				= 255;
	gColorSeq				= false; // true -> show colorspace alignments as colors, not decoded bases
	gColorEdit				= false; // true -> show edits as colors, not decoded bases
	gColorQual				= false; // true -> show colorspace qualities as original quals, not decoded quals
	printPlaceholders       = false; // true -> print records for maxed-out, unaligned reads
	printFlags              = false; // true -> print alignment flags
	printCost				= false; // true -> print cost and stratum
	printParams				= false; // true -> print parameters regarding seeding, ceilings
	gShowSeed				= false; // true -> print per-read pseudo-random seed
	gGaps					= false; // true -> allow gaps
	gInsOpen				= 100;   // insertion open penalty
	gDelOpen				= 100;   // deletion open penalty
	gInsExtend				= 40;    // insertion gap extension penalty
	gDelExtend				= 40;    // deletion gap extension penalty
	gGapBarrier				= 4;     // disallow gaps within this many chars of either end of alignment
	gAllowRedundant			= 0;     // 1 -> allow alignments with 1 anchor in common, 2 -> allow alignments with both anchors in common
	qualities.clear();
	qualities1.clear();
	qualities2.clear();
#ifdef BOWTIE2
	doMultiseed     = true; // true -> use multiseed aligner
#else
	doMultiseed     = false; // false -> use normal bowtie1 aligner
#endif
	seedAlignmentPolicyString.clear();
	msNoCache       = false; // true -> disable local cache
	penMmcType      = DEFAULT_MM_PENALTY_TYPE;
	penMmc          = DEFAULT_MM_PENALTY;
	penSnp          = DEFAULT_SNP_PENALTY;
	penNType        = DEFAULT_N_PENALTY_TYPE;
	penN            = DEFAULT_N_PENALTY;
	penNCatPair     = DEFAULT_N_CAT_PAIR; // concatenate mates before N filtering?
	noisyHpolymer   = false;
	penRdGapConst    = DEFAULT_READ_EXTEND_CONST;
	penRfGapConst    = DEFAULT_REF_EXTEND_CONST;
	penRdGapLinear   = DEFAULT_READ_EXTEND_LINEAR;
	penRfGapLinear   = DEFAULT_REF_EXTEND_LINEAR;
	costCeilConst   = DEFAULT_CEIL_CONST;
	costCeilLinear  = DEFAULT_CEIL_LINEAR;
	nCeilConst      = 2.0f; // constant factor in N ceiling w/r/t read length
	nCeilLinear     = 0.1f; // coeff of linear term in N ceiling w/r/t read length
	multiseedMms    = DEFAULT_SEEDMMS;
	multiseedLen    = DEFAULT_SEEDLEN;
	// Note multiseedPeriod is re-instantiated for each read, since it
	// may ultimately depend on the read length.
	multiseedPeriod = DEFAULT_SEEDPERIOD;
	multiseedIvalType = DEFAULT_IVAL; // formula for seed spacing
	multiseedIvalA  = DEFAULT_IVAL_A; // linear coefficient in formula
	multiseedIvalB  = DEFAULT_IVAL_B; // constant coefficient in formula
	saCountersFn.clear();    // filename to dump per-read SeedAligner counters to
	saActionsFn.clear();     // filename to dump all alignment actions to
	saHitsFn.clear();        // filename to dump all seed hits to
	seedCacheLocalMB   = 32; // # MB to use for non-shared seed alignment cacheing
	seedCacheCurrentMB = 16; // # MB to use for current-read seed hit cacheing
	maxposs            = 25; // maximum number of positions to consider per read
	maxrows            = 10; // maximum number of rows to consider per position
	seedSummaryOnly    = false; // print summary information about seed hits, not alignments
}

static const char *short_options = "fF:qbzhcu:rv:s:aP:t3:5:o:e:n:l:w:p:k:m:M:1:2:I:X:x:B:ySCgO:E:Q:";

enum {
	ARG_ORIG = 256,
	ARG_SEED,
	ARG_DUMP_PATS,
	ARG_RANGE,
	ARG_CONCISE,
	ARG_SOLEXA_QUALS,
	ARG_MAXBTS,
	ARG_VERBOSE,
	ARG_STARTVERBOSE,
	ARG_QUIET,
	ARG_RANDOM_READS,
	ARG_NOOUT,
	ARG_FAST,
	ARG_METRIC_IVAL,
	ARG_METRIC_FILE,
	ARG_METRIC_STDERR,
	ARG_AL,
	ARG_UN,
	ARG_MAXDUMP,
	ARG_REFIDX,
	ARG_SANITY,
	ARG_OLDBEST,
	ARG_BETTER,
	ARG_BEST,
	ARG_PARTITION,
	ARG_INTEGER_QUALS,
	ARG_NOMAQROUND,
	ARG_USE_SPINLOCK,
	ARG_FILEPAR,
	ARG_SHMEM,
	ARG_MM,
	ARG_MMSWEEP,
	ARG_FF,
	ARG_FR,
	ARG_RF,
	ARG_MIXED_ATTEMPTS,
	ARG_NO_MIXED,
	ARG_NO_DISCORDANT,
	ARG_NO_RECONCILE,
	ARG_CACHE_LIM,
	ARG_CACHE_SZ,
	ARG_NO_FW,
	ARG_NO_RC,
	ARG_SKIP,
	ARG_STRAND_FIX,
	ARG_STATS,
	ARG_ONETWO,
	ARG_PHRED64,
	ARG_PHRED33,
	ARG_CHUNKMBS,
	ARG_CHUNKSZ,
	ARG_CHUNKVERBOSE,
	ARG_STRATA,
	ARG_PEV2,
	ARG_CHAINOUT,
	ARG_CHAININ,
	ARG_REFMAP,
	ARG_ANNOTMAP,
	ARG_REPORTSE,
	ARG_HADOOPOUT,
	ARG_FUZZY,
	ARG_FULLREF,
	ARG_USAGE,
	ARG_SNPPHRED,
	ARG_SNPFRAC,
	ARG_SAM_NOHEAD,
	ARG_SAM_NOSQ,
	ARG_SAM_RG,
	ARG_SUPPRESS_FIELDS,
	ARG_DEFAULT_MAPQ,
	ARG_COLOR_SEQ,
	ARG_COLOR_EDIT,
	ARG_COLOR_QUAL,
	ARG_PRINT_PLACEHOLDERS,
	ARG_PRINT_FLAGS,
	ARG_COST,
	ARG_PRINT_PARAMS,
	ARG_COLOR_KEEP_ENDS,
	ARG_SHOWSEED,
	ARG_GAP_BAR,
	ARG_REDUNDANTS,
	ARG_QUALS1,
	ARG_QUALS2,
	ARG_QSEQ,
	ARG_SA_DUMP,
	ARG_SEED_SUMM,
	ARG_OVERHANG,
	ARG_NO_CACHE,
	ARG_NOISY_HPOLY
};

static struct option long_options[] = {
	{(char*)"verbose",      no_argument,       0,            ARG_VERBOSE},
	{(char*)"startverbose", no_argument,       0,            ARG_STARTVERBOSE},
	{(char*)"quiet",        no_argument,       0,            ARG_QUIET},
	{(char*)"sanity",       no_argument,       0,            ARG_SANITY},
	{(char*)"pause",        no_argument,       &ipause,      1},
	{(char*)"orig",         required_argument, 0,            ARG_ORIG},
	{(char*)"all",          no_argument,       0,            'a'},
	{(char*)"concise",      no_argument,       0,            ARG_CONCISE},
	{(char*)"noout",        no_argument,       0,            ARG_NOOUT},
	{(char*)"solexa-quals", no_argument,       0,            ARG_SOLEXA_QUALS},
	{(char*)"integer-quals",no_argument,       0,            ARG_INTEGER_QUALS},
	{(char*)"metrics",      required_argument, 0,            ARG_METRIC_IVAL},
	{(char*)"metrics-file", required_argument, 0,            ARG_METRIC_FILE},
	{(char*)"metrics-stderr",no_argument,      0,            ARG_METRIC_STDERR},
	{(char*)"time",         no_argument,       0,            't'},
	{(char*)"trim3",        required_argument, 0,            '3'},
	{(char*)"trim5",        required_argument, 0,            '5'},
	{(char*)"seed",         required_argument, 0,            ARG_SEED},
	{(char*)"qupto",        required_argument, 0,            'u'},
	{(char*)"al",           required_argument, 0,            ARG_AL},
	{(char*)"un",           required_argument, 0,            ARG_UN},
	{(char*)"max",          required_argument, 0,            ARG_MAXDUMP},
	{(char*)"offrate",      required_argument, 0,            'o'},
	{(char*)"reportopps",   no_argument,       &reportOpps,  1},
	{(char*)"version",      no_argument,       &showVersion, 1},
	{(char*)"dumppats",     required_argument, 0,            ARG_DUMP_PATS},
	{(char*)"maqerr",       required_argument, 0,            'e'},
	{(char*)"seedlen",      required_argument, 0,            'l'},
	{(char*)"seedmms",      required_argument, 0,            'n'},
	{(char*)"filepar",      no_argument,       0,            ARG_FILEPAR},
	{(char*)"help",         no_argument,       0,            'h'},
	{(char*)"threads",      required_argument, 0,            'p'},
	{(char*)"khits",        required_argument, 0,            'k'},
	{(char*)"mhits",        required_argument, 0,            'm'},
	{(char*)"minins",       required_argument, 0,            'I'},
	{(char*)"maxins",       required_argument, 0,            'X'},
	{(char*)"quals",        required_argument, 0,            'Q'},
	{(char*)"Q1",           required_argument, 0,            ARG_QUALS1},
	{(char*)"Q2",           required_argument, 0,            ARG_QUALS2},
	{(char*)"strata",       no_argument,       0,            ARG_STRATA},
	{(char*)"nomaqround",   no_argument,       0,            ARG_NOMAQROUND},
	{(char*)"refidx",       no_argument,       0,            ARG_REFIDX},
	{(char*)"range",        no_argument,       0,            ARG_RANGE},
	{(char*)"maxbts",       required_argument, 0,            ARG_MAXBTS},
	{(char*)"randread",     no_argument,       0,            ARG_RANDOM_READS},
	{(char*)"partition",    required_argument, 0,            ARG_PARTITION},
	{(char*)"nospin",       no_argument,       0,            ARG_USE_SPINLOCK},
	{(char*)"ff",           no_argument,       0,            ARG_FF},
	{(char*)"fr",           no_argument,       0,            ARG_FR},
	{(char*)"rf",           no_argument,       0,            ARG_RF},
	{(char*)"mixthresh",    required_argument, 0,            'x'},
	{(char*)"pairtries",    required_argument, 0,            ARG_MIXED_ATTEMPTS},
	{(char*)"noreconcile",  no_argument,       0,            ARG_NO_RECONCILE},
	{(char*)"cachelim",     required_argument, 0,            ARG_CACHE_LIM},
	{(char*)"cachesz",      required_argument, 0,            ARG_CACHE_SZ},
	{(char*)"nofw",         no_argument,       0,            ARG_NO_FW},
	{(char*)"norc",         no_argument,       0,            ARG_NO_RC},
	{(char*)"offbase",      required_argument, 0,            'B'},
	{(char*)"tryhard",      no_argument,       0,            'y'},
	{(char*)"skip",         required_argument, 0,            's'},
	{(char*)"strandfix",    no_argument,       0,            ARG_STRAND_FIX},
	{(char*)"stats",        no_argument,       0,            ARG_STATS},
	{(char*)"12",           required_argument, 0,            ARG_ONETWO},
	{(char*)"phred33-quals", no_argument,      0,            ARG_PHRED33},
	{(char*)"phred64-quals", no_argument,      0,            ARG_PHRED64},
	{(char*)"solexa1.3-quals", no_argument,    0,            ARG_PHRED64},
	{(char*)"chunkmbs",     required_argument, 0,            ARG_CHUNKMBS},
	{(char*)"chunksz",      required_argument, 0,            ARG_CHUNKSZ},
	{(char*)"chunkverbose", no_argument,       0,            ARG_CHUNKVERBOSE},
	{(char*)"mm",           no_argument,       0,            ARG_MM},
	{(char*)"shmem",        no_argument,       0,            ARG_SHMEM},
	{(char*)"mmsweep",      no_argument,       0,            ARG_MMSWEEP},
	{(char*)"pev2",         no_argument,       0,            ARG_PEV2},
	{(char*)"chainout",     no_argument,       0,            ARG_CHAINOUT},
	{(char*)"chainin",      no_argument,       0,            ARG_CHAININ},
	{(char*)"refmap",       required_argument, 0,            ARG_REFMAP},
	{(char*)"annotmap",     required_argument, 0,            ARG_ANNOTMAP},
	{(char*)"reportse",     no_argument,       0,            ARG_REPORTSE},
	{(char*)"hadoopout",    no_argument,       0,            ARG_HADOOPOUT},
	{(char*)"fuzzy",        no_argument,       0,            ARG_FUZZY},
	{(char*)"fullref",      no_argument,       0,            ARG_FULLREF},
	{(char*)"usage",        no_argument,       0,            ARG_USAGE},
	{(char*)"sam",          no_argument,       0,            'S'},
	{(char*)"gaps",         no_argument,       0,            'g'},
	{(char*)"sam-nohead",   no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-nosq",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"sam-noSQ",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"color",        no_argument,       0,            'C'},
	{(char*)"sam-RG",       required_argument, 0,            ARG_SAM_RG},
	{(char*)"snpphred",     required_argument, 0,            ARG_SNPPHRED},
	{(char*)"snpfrac",      required_argument, 0,            ARG_SNPFRAC},
	{(char*)"suppress",     required_argument, 0,            ARG_SUPPRESS_FIELDS},
	{(char*)"mapq",         required_argument, 0,            ARG_DEFAULT_MAPQ},
	{(char*)"col-cseq",     no_argument,       0,            ARG_COLOR_SEQ},
	{(char*)"col-cqual",    no_argument,       0,            ARG_COLOR_QUAL},
	{(char*)"col-cedit",    no_argument,       0,            ARG_COLOR_EDIT},
	{(char*)"col-keepends", no_argument,       0,            ARG_COLOR_KEEP_ENDS},
	{(char*)"print-placeholders", no_argument, 0,            ARG_PRINT_PLACEHOLDERS},
	{(char*)"print-flags",  no_argument,       0,            ARG_PRINT_FLAGS},
	{(char*)"cost",         no_argument,       0,            ARG_COST},
	{(char*)"print-params", no_argument,       0,            ARG_PRINT_PARAMS},
	{(char*)"showseed",     no_argument,       0,            ARG_SHOWSEED},
	{(char*)"gbar",         required_argument, 0,            ARG_GAP_BAR},
	{(char*)"gopen",        required_argument, 0,            'O'},
	{(char*)"gextend",      required_argument, 0,            'E'},
	{(char*)"redundants",   required_argument, 0,            ARG_REDUNDANTS},
	{(char*)"redundant",    required_argument, 0,            ARG_REDUNDANTS},
	{(char*)"qseq",         no_argument,       0,            ARG_QSEQ},
	{(char*)"align-policy", required_argument, 0,            'P'},
	{(char*)"seed-summ",    no_argument,       0,            ARG_SEED_SUMM},
	{(char*)"seed-summary", no_argument,       0,            ARG_SEED_SUMM},
	{(char*)"overhang",     no_argument,       0,            ARG_OVERHANG},
	{(char*)"no-cache",     no_argument,       0,            ARG_NO_CACHE},
	{(char*)"454",          no_argument,       0,            ARG_NOISY_HPOLY},
	{(char*)"ion-torrent",  no_argument,       0,            ARG_NOISY_HPOLY},
	{(char*)"no-mixed",     no_argument,       0,            ARG_NO_MIXED},
	{(char*)"no-discordant",no_argument,       0,            ARG_NO_DISCORDANT},
	{(char*)0, 0, 0, 0} // terminator
};

#ifdef BOWTIE2
/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: " << endl
	    << "  bowtie2 [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]" << endl
	    << endl
	    << "  <m1>    Comma-separated list of files containing upstream mates (or the" << endl
	    << "          sequences themselves, if -c is set) paired with mates in <m2>" << endl
	    << "  <m2>    Comma-separated list of files containing downstream mates (or the" << endl
	    << "          sequences themselves if -c is set) paired with mates in <m1>" << endl
	    << "  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be" << endl
	    << "          a mixture of paired and unpaired.  Specify \"-\" for stdin." << endl
	    << "  <s>     Comma-separated list of files containing unpaired reads, or the" << endl
	    << "          sequences themselves, if -c is set.  Specify \"-\" for stdin." << endl
	    << "  <hit>   File to write hits to (default: stdout)" << endl
	    << "Key missing features:" << endl
	    << "  - No paired-end support" << endl
	    << "  - No colorspace support" << endl
	    << "  - No SAM output support" << endl
	    << "Input:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    << "  --qseq             query input files are in Illumina's qseq format" << endl
	    << "  -c                 query sequences given on cmd line (as <mates>, <singles>)" << endl
	    << "  -C                 reads and index are in colorspace" << endl
	    << "  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C" << endl
	    << "  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input" << endl
	    << "  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
	    << "  --phred33-quals    input quals are Phred+33 (default)" << endl
	    << "  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)" << endl
	    << "  --solexa-quals     input quals are from GA Pipeline ver. < 1.3" << endl
	    << "  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3" << endl
	    << "  --integer-quals    qualities are given as space-separated integers (not ASCII)" << endl
	    << "Alignment:" << endl
	    << "  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)" << endl
	    << "  --gNofw/--gNorc    do not align to forward/reverse-complement reference strand" << endl
		<< "Paired-end:" << endl
	    << "  -I/--minins <int>  minimum fragment length (default: 0)" << endl
	    << "  -X/--maxins <int>  maximum fragment length (default: 250)" << endl
	    << "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)" << endl
		<< "  --no-mixed         report only paired alns, ignore unpaired" << endl
		<< "  --no-discordant    report only concordant paired-end alns, ignore discordant" << endl
	    << "Reporting:" << endl
	    << "  -k <int>           report up to <int> good alignments per read (default: 1)" << endl
	    << "  -a/--all           report all alignments per read (much slower than low -k)" << endl
	    << "  -m <int>           suppress all alignments if > <int> exist (def: no limit)" << endl
	    << "  -M <int>           like -m, but reports 1 random hit (MAPQ=0)" << endl
	    << "Output:" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
	    << "  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)" << endl
	    << "  --quiet            print nothing but the alignments" << endl
	    << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --al <fname>       write aligned reads/pairs to file(s) <fname>" << endl
	    << "  --un <fname>       write unaligned reads/pairs to file(s) <fname>" << endl
	    << "  --max <fname>      write reads/pairs over -m limit to file(s) <fname>" << endl
	    << "  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output" << endl
	    << "  --fullref          write entire ref name (default: only up to 1st space)" << endl
	    << "Performance:" << endl
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
#ifdef BOWTIE_PTHREADS
	    << "  -p/--threads <int> number of alignment threads to launch (default: 1)" << endl
#endif
#ifdef BOWTIE_MM
	    << "  --mm               use memory-mapped I/O for index; many 'bowtie's can share" << endl
#endif
#ifdef BOWTIE_SHARED_MEM
	    << "  --shmem            use shared mem for index; many 'bowtie's can share" << endl
#endif
	    << "Other:" << endl
	    << "  --seed <int>       seed for random number generator" << endl
	    << "  --verbose          verbose output (for debugging)" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print this usage message" << endl
	    ;
}
#else
/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: " << endl
	    << "  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]" << endl
	    << endl
	    << "  <m1>    Comma-separated list of files containing upstream mates (or the" << endl
	    << "          sequences themselves, if -c is set) paired with mates in <m2>" << endl
	    << "  <m2>    Comma-separated list of files containing downstream mates (or the" << endl
	    << "          sequences themselves if -c is set) paired with mates in <m1>" << endl
	    << "  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be" << endl
	    << "          a mixture of paired and unpaired.  Specify \"-\" for stdin." << endl
	    << "  <s>     Comma-separated list of files containing unpaired reads, or the" << endl
	    << "          sequences themselves, if -c is set.  Specify \"-\" for stdin." << endl
	    << "  <hit>   File to write hits to (default: stdout)" << endl
	    << "Input:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    << "  --qseq             query input files are in Illumina's qseq format" << endl
	    << "  -c                 query sequences given on cmd line (as <mates>, <singles>)" << endl
	    << "  -C                 reads and index are in colorspace" << endl
	    << "  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C" << endl
	    << "  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input" << endl
	    << "  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
	    << "  --phred33-quals    input quals are Phred+33 (default)" << endl
	    << "  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)" << endl
	    << "  --solexa-quals     input quals are from GA Pipeline ver. < 1.3" << endl
	    << "  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3" << endl
	    << "  --integer-quals    qualities are given as space-separated integers (not ASCII)" << endl
	    << "Alignment:" << endl
	    << "  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities" << endl
	    << "    or" << endl
	    << "  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)" << endl
	    << "  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)" << endl
	    << "  -l/--seedlen <int> seed length for -n (default: 28)" << endl
	    << "  -g/--gaps          allow gapped alignment" << endl
	    << "  --gbar <int>       disallow gaps within <int> chars of either end of read" << endl
	    << "  -O/--gopen <int> gap open penalty (default: 100)" << endl
	    << "  -E/--gextend <int> gap extension penalty (default: 40)" << endl
	    << "  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)" << endl
	    << "  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)" << endl
	    << "  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)" << endl
	    << "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)" << endl
	    << "  --gNofw/--gNorc      do not align to forward/reverse-complement reference strand" << endl
	    << "  --maxbts <int>     max # backtracks for -n 2/3 (default: 800)" << endl
	    << "  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)" << endl
	    << "  -y/--tryhard       try hard to find valid alignments, at the expense of speed" << endl
	    << "  --chunkmbs <int>   max megabytes of RAM for search frames (def: 32)" << endl
	    << "Reporting:" << endl
	    << "  -k <int>           report up to <int> good alignments per read (default: 1)" << endl
	    << "  -a/--all           report all alignments per read (much slower than low -k)" << endl
	    << "  -m <int>           suppress all alignments if > <int> exist (def: no limit)" << endl
	    << "  -M <int>           like -m, but reports 1 random hit (MAPQ=0)" << endl
	    << "  --strata           hits in sub-optimal strata aren't reported" << endl
	    << "Output:" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
	    << "  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)" << endl
	    << "  --quiet            print nothing but the alignments" << endl
	    << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --al <fname>       write aligned reads/pairs to file(s) <fname>" << endl
	    << "  --un <fname>       write unaligned reads/pairs to file(s) <fname>" << endl
	    << "  --max <fname>      write reads/pairs over -m limit to file(s) <fname>" << endl
	    << "  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output" << endl
	    << "  --fullref          write entire ref name (default: only up to 1st space)" << endl
	    << "Colorspace:" << endl
	    << "  --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)" << endl
	    << "     or" << endl
	    << "  --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred" << endl
	    << "  --col-cseq         print aligned colorspace seqs as colors, not decoded bases" << endl
	    << "  --col-cqual        print original colorspace quals, not decoded quals" << endl
	    << "  --col-keepends     keep nucleotides at extreme ends of decoded alignment" << endl
	    << "SAM:" << endl
	    << "  -S/--sam           write hits in SAM format" << endl
	    << "  --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments" << endl
	    << "  --sam-nohead       supppress header lines (starting with @) for SAM output" << endl
	    << "  --sam-nosq         supppress @SQ header lines for SAM output" << endl
	    << "  --sam-RG <text>    add <text> (usually \"lab=value\") to @RG line of SAM header" << endl
	    << "Performance:" << endl
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
#ifdef BOWTIE_PTHREADS
	    << "  -p/--threads <int> number of alignment threads to launch (default: 1)" << endl
#endif
#ifdef BOWTIE_MM
	    << "  --mm               use memory-mapped I/O for index; many 'bowtie's can share" << endl
#endif
#ifdef BOWTIE_SHARED_MEM
	    << "  --shmem            use shared mem for index; many 'bowtie's can share" << endl
#endif
	    << "Other:" << endl
	    << "  --seed <int>       seed for random number generator" << endl
	    << "  --verbose          verbose output (for debugging)" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print this usage message" << endl
	    ;
}
#endif

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
 * Parse from optarg by default.
 */
static int parseInt(int lower, const char *errmsg) {
	return parseInt(lower, INT_MAX, errmsg, optarg);
}

/**
 * Upper is INT_MAX by default.
 */
static int parseInt(int lower, const char *errmsg, const char *arg) {
	return parseInt(lower, INT_MAX, errmsg, arg);
}

/**
 * Upper is INT_MAX, parse from optarg by default.
 */
static int parseInt(int lower, int upper, const char *errmsg) {
	return parseInt(lower, upper, errmsg, optarg);
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

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
pair<T, T> parsePair(const char *str, char delim) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	pair<T, T> ret;
	ret.first = parse<T>(ss[0].c_str());
	ret.second = parse<T>(ss[1].c_str());
	return ret;
}

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
void parseTuple(const char *str, char delim, EList<T>& ret) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	for(size_t i = 0; i < ss.size(); i++) {
		ret.push_back(parse<T>(ss[i].c_str()));
	}
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	if(startVerbose) { cerr << "Parsing options: "; logTime(cerr, true); }
	do {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		switch (next_option) {
			case '1': tokenize(optarg, ",", mates1); break;
			case '2': tokenize(optarg, ",", mates2); break;
			case ARG_ONETWO: tokenize(optarg, ",", mates12); format = TAB_MATE; break;
			case 'f': format = FASTA; break;
			case 'g': gGaps = true; break;
			case 'F': {
				format = FASTA_CONT;
				pair<uint32_t, uint32_t> p = parsePair<uint32_t>(optarg, ',');
				fastaContLen = p.first;
				fastaContFreq = p.second;
				break;
			}
			case 'q': format = FASTQ; break;
			case 'r': format = RAW; break;
			case 'c': format = CMDLINE; break;
			case ARG_QSEQ: format = QSEQ; break;
			case 'C': gColor = true; break;
			case ARG_CHAININ: format = INPUT_CHAIN; break;
			case 'I':
				gMinInsert = parseInt(0, "-I arg must be positive");
				break;
			case 'X':
				gMaxInsert = parseInt(1, "-X arg must be at least 1");
				break;
			case ARG_NO_DISCORDANT:
				gReportDiscordant = false;
				break;
			case ARG_NO_MIXED:
				gReportMixed = false;
				break;
			case 's':
				skipReads = (uint32_t)parseInt(0, "-s arg must be positive");
				break;
			case ARG_FF: gMate1fw = true;  gMate2fw = true;  mateFwSet = true; break;
			case ARG_RF: gMate1fw = false; gMate2fw = true;  mateFwSet = true; break;
			case ARG_FR: gMate1fw = true;  gMate2fw = false; mateFwSet = true; break;
			case ARG_RANDOM_READS: format = RANDOM; break;
			case ARG_RANGE: gRangeMode = true; break;
			case ARG_CONCISE: outType = OUTPUT_CONCISE; break;
			case ARG_CHAINOUT: outType = OUTPUT_CHAIN; break;
			case 'S': outType = OUTPUT_SAM; break;
			case ARG_NOOUT: outType = OUTPUT_NONE; break;
			case ARG_REFMAP: refMapFile = optarg; break;
			case ARG_ANNOTMAP: annotMapFile = optarg; break;
			case ARG_USE_SPINLOCK: useSpinlock = false; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_COLOR_SEQ: gColorSeq = true; break;
			case ARG_COLOR_EDIT: gColorEdit = true; break;
			case ARG_COLOR_QUAL: gColorQual = true; break;
			case ARG_SHOWSEED: gShowSeed = true; break;
			case ARG_SEED_SUMM: seedSummaryOnly = true; break;
			case ARG_SUPPRESS_FIELDS: {
				EList<string> supp;
				tokenize(optarg, ",", supp);
				for(size_t i = 0; i < supp.size(); i++) {
					int ii = parseInt(1, "--suppress arg must be at least 1", supp[i].c_str());
					suppressOuts[ii-1] = true;
				}
				break;
			}
			case ARG_MM: {
#ifdef BOWTIE_MM
				useMm = true;
				break;
#else
				cerr << "Memory-mapped I/O mode is disabled because bowtie was not compiled with" << endl
				     << "BOWTIE_MM defined.  Memory-mapped I/O is not supported under Windows.  If you" << endl
				     << "would like to use memory-mapped I/O on a platform that supports it, please" << endl
				     << "refrain from specifying BOWTIE_MM=0 when compiling Bowtie." << endl;
				throw 1;
#endif
			}
			case ARG_MMSWEEP: mmSweep = true; break;
			case ARG_HADOOPOUT: hadoopOut = true; break;
			case ARG_AL: dumpAlBase = optarg; break;
			case ARG_UN: dumpUnalBase = optarg; break;
			case ARG_MAXDUMP: dumpMaxBase = optarg; break;
			case ARG_SOLEXA_QUALS: solexaQuals = true; break;
			case ARG_INTEGER_QUALS: integerQuals = true; break;
			case ARG_PHRED64: phred64Quals = true; break;
			case ARG_PHRED33: solexaQuals = false; phred64Quals = false; break;
			case ARG_NOMAQROUND: gNoMaqRound = true; break;
			case ARG_COLOR_KEEP_ENDS: gColorExEnds = false; break;
			case ARG_OVERHANG: gReportOverhangs = true; break;
			case ARG_NO_CACHE: msNoCache = true; break;
			case ARG_SNPPHRED: gSnpPhred = parseInt(0, "--snpphred must be at least 0"); break;
			case ARG_SNPFRAC: {
				double p = parse<double>(optarg);
				if(p <= 0.0) {
					cerr << "Error: --snpfrac parameter must be > 0.0" << endl;
					throw 1;
				}
				p = (log10(p) * -10);
				gSnpPhred = (int)(p + 0.5);
				if(gSnpPhred < 10)
				cout << "gSnpPhred: " << gSnpPhred << endl;
				break;
			}
			case 'z': {
				cerr << "Error: -z/--phased mode is no longer supported" << endl;
				throw 1;
			}
			case ARG_REFIDX: noRefNames = true; break;
			case ARG_FUZZY: fuzzy = true; break;
			case ARG_REPORTSE: gReportSe = true; break;
			case ARG_FULLREF: fullRef = true; break;
			case 'B':
				offBase = parseInt(-999999, "-B/--offbase cannot be a large negative number");
				break;
			case ARG_GAP_BAR:
				gGapBarrier = parseInt(1, "--gbar must be no less than 1");
				break;
			case 'O':
				gInsOpen = gDelOpen = parseInt(1, "-O/--gopen must be at least 1");
				break;
			case 'E':
				gInsExtend = gDelExtend = parseInt(1, "-E/--gextend must be at least 1");
				break;
			case ARG_REDUNDANTS:
				gAllowRedundant = parseInt(0, "--redundant must be at least 0");
				break;
			case ARG_SEED:
				seed = parseInt(0, "--seed arg must be at least 0");
				break;
			case 'u':
				qUpto = (uint32_t)parseInt(1, "-u/--qupto arg must be at least 1");
				break;
			case 'k':
				khits = (uint32_t)parseInt(1, "-k arg must be at least 1");
				break;
			case 'Q':
				tokenize(optarg, ",", qualities);
				integerQuals = true;
				break;
			case ARG_QUALS1:
				tokenize(optarg, ",", qualities1);
				integerQuals = true;
				break;
			case ARG_QUALS2:
				tokenize(optarg, ",", qualities2);
				integerQuals = true;
				break;
			case 'M':
				sampleMax = true;
			case 'm':
				mhits = (uint32_t)parseInt(1, "-m arg must be at least 1");
				break;
			case 'x':
				mixedThresh = (uint32_t)parseInt(0, "-x arg must be at least 0");
				break;
			case ARG_MIXED_ATTEMPTS:
				mixedAttemptLim = (uint32_t)parseInt(1, "--mixatt arg must be at least 1");
				break;
			case ARG_CACHE_LIM:
				cacheLimit = (uint32_t)parseInt(1, "--cachelim arg must be at least 1");
				break;
			case ARG_CACHE_SZ:
				cacheSize = (uint32_t)parseInt(1, "--cachesz arg must be at least 1");
				cacheSize *= (1024 * 1024); // convert from MB to B
				break;
			case ARG_NO_RECONCILE:
				dontReconcileMates = true;
				break;
			case 'p':
#ifndef BOWTIE_PTHREADS
				cerr << "-p/--threads is disabled because bowtie was not compiled with pthreads support" << endl;
				throw 1;
#endif
				nthreads = parseInt(1, "-p/--threads arg must be at least 1");
				break;
			case ARG_FILEPAR:
#ifndef BOWTIE_PTHREADS
				cerr << "--filepar is disabled because bowtie was not compiled with pthreads support" << endl;
				throw 1;
#endif
				fileParallel = true;
				break;
			case 'v':
				maqLike = 0;
				mismatches = parseInt(0, 3, "-v arg must be at least 0 and at most 3");
				break;
			case '3': gTrim3 = parseInt(0, "-3/--trim3 arg must be at least 0"); break;
			case '5': gTrim5 = parseInt(0, "-5/--trim5 arg must be at least 0"); break;
			case 'o': offRate = parseInt(1, "-o/--offrate arg must be at least 1"); break;
			case 'e': qualThresh = parseInt(1, "-e/--err arg must be at least 1"); break;
			case 'n': seedMms = parseInt(0, 3, "-n/--seedmms arg must be at least 0 and at most 3"); maqLike = 1; break;
			case 'l': seedLen = parseInt(5, "-l/--seedlen arg must be at least 5"); break;
			case 'h': printUsage(cout); throw 0; break;
			case ARG_USAGE: printUsage(cout); throw 0; break;
			case 'a': allHits = true; break;
			case 'y': tryHard = true; break;
			case ARG_CHUNKMBS: chunkPoolMegabytes = parseInt(1, "--chunkmbs arg must be at least 1"); break;
			case ARG_CHUNKSZ: chunkSz = parseInt(1, "--chunksz arg must be at least 1"); break;
			case ARG_CHUNKVERBOSE: chunkVerbose = true; break;
			case ARG_STRATA: strata = true; break;
			case ARG_VERBOSE: gVerbose = 1; break;
			case ARG_STARTVERBOSE: startVerbose = true; break;
			case ARG_QUIET: gQuiet = true; break;
			case ARG_SANITY: sanityCheck = true; break;
			case 't': timing = true; break;
			case ARG_METRIC_IVAL: {
#ifdef BOWTIE_PTHREADS
				metricsIval = parseInt(1, "--metrics arg must be at least 1");
#else
				cerr << "Must compile with BOWTIE_PTHREADS to use --metrics" << endl;
				throw 1;
#endif
				break;
			}
			case ARG_METRIC_FILE: metricsFile = optarg; break;
			case ARG_METRIC_STDERR: metricsStderr = true; break;
			case ARG_NO_FW: gNofw = true; break;
			case ARG_NO_RC: gNorc = true; break;
			case ARG_STATS: stats = true; break;
			case ARG_SAM_NOHEAD: samNoHead = true; break;
			case ARG_SAM_NOSQ: samNoSQ = true; break;
			case ARG_SAM_RG: {
				if(!rgs.empty()) rgs += '\t';
				rgs += optarg;
				break;
			}
			case ARG_PRINT_PLACEHOLDERS: printPlaceholders = true; break;
			case ARG_PRINT_FLAGS: printFlags = true; break;
			case ARG_COST: printCost = true; break;
			case ARG_PRINT_PARAMS: printParams = true; break;
			case ARG_DEFAULT_MAPQ:
				defaultMapq = parseInt(0, "--mapq must be positive");
				break;
			case ARG_MAXBTS: {
				maxBts  = parseInt(0, "--maxbts must be positive");
				maxBtsBetter = maxBts;
				break;
			}
			case ARG_DUMP_PATS: patDumpfile = optarg; break;
			case ARG_STRAND_FIX: gStrandFix = true; break;
			case ARG_PARTITION: partitionSz = parse<int>(optarg); break;
			case ARG_ORIG:
				if(optarg == NULL || strlen(optarg) == 0) {
					cerr << "--orig arg must be followed by a string" << endl;
					printUsage(cerr);
					throw 1;
				}
				origString = optarg;
				break;
			case ARG_NOISY_HPOLY: noisyHpolymer = true; break;
			case 'P': seedAlignmentPolicyString = optarg; break;
			case ARG_SA_DUMP: {
				EList<string> fns;
				parseTuple<string>(optarg, ',', fns);
				if(fns.size() >= 1) saHitsFn     = fns[0];
				if(fns.size() >= 2) saCountersFn = fns[1];
				if(fns.size() >= 3) saActionsFn  = fns[2];
				break;
			}
			case -1: break; /* Done with options. */
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
		}
	} while(next_option != -1);
	bool paired = mates1.size() > 0 || mates2.size() > 0 || mates12.size() > 0;
	if(gRangeMode) {
		// Tell the Ebwt loader to ignore the suffix-array portion of
		// the index.  We don't need it because the user isn't asking
		// for bowtie to report reference positions (just matrix
		// ranges).
		offRate = 32;
	}
	if(!seedAlignmentPolicyString.empty()) {
		SeedAlignmentPolicy::parseString(
			seedAlignmentPolicyString,
			noisyHpolymer,
			penMmcType,
			penMmc,
			penSnp,
			penNType,
			penN,
			penRdGapConst,
			penRfGapConst,
			penRdGapLinear,
			penRfGapLinear,
			costCeilConst,
			costCeilLinear,
			nCeilConst,
			nCeilLinear,
			penNCatPair,
			multiseedMms,
			multiseedLen,
			multiseedPeriod,
			multiseedIvalType,
			multiseedIvalA,
			multiseedIvalB);
	} else {
		penMmcType        = DEFAULT_MM_PENALTY_TYPE;
		penMmc            = DEFAULT_MM_PENALTY;
		penSnp            = DEFAULT_SNP_PENALTY;
		penNType          = DEFAULT_N_PENALTY_TYPE;
		penN              = DEFAULT_N_PENALTY;
		penRdGapConst     = noisyHpolymer ? DEFAULT_READ_EXTEND_CONST_BADHPOLY
		                                  : DEFAULT_READ_EXTEND_CONST;
		penRfGapConst     = noisyHpolymer ? DEFAULT_REF_EXTEND_CONST_BADHPOLY
		                                  : DEFAULT_REF_EXTEND_CONST;
		penRdGapLinear    = noisyHpolymer ? DEFAULT_READ_EXTEND_LINEAR_BADHPOLY
		                                  : DEFAULT_READ_EXTEND_LINEAR;
		penRfGapLinear     = noisyHpolymer ? DEFAULT_REF_EXTEND_LINEAR_BADHPOLY
		                                  : DEFAULT_REF_EXTEND_LINEAR;
		costCeilConst     = DEFAULT_CEIL_CONST;
		costCeilLinear    = DEFAULT_CEIL_LINEAR;
		nCeilConst        = DEFAULT_N_CEIL_CONST;
		nCeilLinear       = DEFAULT_N_CEIL_LINEAR;
		multiseedMms      = DEFAULT_SEEDMMS;
		multiseedLen      = DEFAULT_SEEDLEN;
		multiseedPeriod   = DEFAULT_SEEDPERIOD;
		multiseedIvalType = DEFAULT_IVAL;
		multiseedIvalA    = DEFAULT_IVAL_A;
		multiseedIvalB    = DEFAULT_IVAL_B;
	}
	if(mates1.size() != mates2.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
		     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
		     << "sequences must be specified with -1 and -2." << endl;
		throw 1;
	}
	if(qualities.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with -Q but -f was not" << endl
		     << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities.size() && !gColor) {
		cerr << "Error: one or more quality files were specified with -Q but -C was not" << endl
		     << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q1 but -f was not" << endl
		     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() && !gColor) {
		cerr << "Error: one or more quality files were specified with --Q1 but -C was not" << endl
		     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities2.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q2 but -f was not" << endl
		     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities2.size() && !gColor) {
		cerr << "Error: one or more quality files were specified with --Q2 but -C was not" << endl
		     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() > 0 && mates1.size() != qualities1.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << qualities1.size() << endl
		     << "quality files were specified with --Q1.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -1 and --Q1." << endl;
		throw 1;
	}
	if(qualities2.size() > 0 && mates2.size() != qualities2.size()) {
		cerr << "Error: " << mates2.size() << " mate files/sequences were specified with -2, but " << qualities2.size() << endl
		     << "quality files were specified with --Q2.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -2 and --Q2." << endl;
		throw 1;
	}
	// Check for duplicate mate input files
	if(format != CMDLINE) {
		for(size_t i = 0; i < mates1.size(); i++) {
			for(size_t j = 0; j < mates2.size(); j++) {
				if(mates1[i] == mates2[j] && !gQuiet) {
					cerr << "Warning: Same mate file \"" << mates1[i] << "\" appears as argument to both -1 and -2" << endl;
				}
			}
		}
	}
	if(tryHard) {
		// Increase backtracking limit to huge number
		maxBts = maxBtsBetter = INT_MAX;
		// Increase number of paired-end scan attempts to huge number
		mixedAttemptLim = UINT_MAX;
	}
	if(strata && !allHits && khits == 1 && mhits == 0xffffffff) {
		cerr << "--strata has no effect unless combined with -k, -m or -a" << endl;
		throw 1;
	}
	// If both -s and -u are used, we need to adjust qUpto accordingly
	// since it uses patid to know if we've reached the -u limit (and
	// patids are all shifted up by skipReads characters)
	if(qUpto + skipReads > qUpto) {
		qUpto += skipReads;
	}
	if(useShmem && useMm && !gQuiet) {
		cerr << "Warning: --shmem overrides --mm..." << endl;
		useMm = false;
	}
	if(gSnpPhred <= 10 && gColor && !gQuiet) {
		cerr << "Warning: the colorspace SNP penalty (--snpphred) is very low: " << gSnpPhred << endl;
	}
	if(format == INPUT_CHAIN) {
		bool error = false;
		if(paired) {
			cerr << "Error: --chainin cannot be combined with paired-end alignment; aborting..." << endl;
			error = true;
		}
		if(error) throw 1;
	}
	if(outType == OUTPUT_CHAIN) {
		bool error = false;
		if(paired) {
			cerr << "Error: --chainout cannot be combined with paired-end alignment; aborting..." << endl;
			error = true;
		}
		if(error) throw 1;
	}
	if(!mateFwSet) {
		if(gColor) {
			// Set colorspace default (--ff)
			gMate1fw = true;
			gMate2fw = true;
		} else {
			// Set nucleotide space default (--fr)
			gMate1fw = true;
			gMate2fw = false;
		}
	}
	size_t fieldsSuppressed = 0;
	for(size_t i = 0; i < suppressOuts.size(); i++) {
		if(suppressOuts[i]) fieldsSuppressed++;
	}
	if(outType != OUTPUT_FULL && fieldsSuppressed > 0 && !gQuiet) {
		cerr << "Warning: Ignoring --suppress because output type is not default." << endl;
		cerr << "         --suppress is only available for the default output type." << endl;
		suppressOuts.fill(false);
	}
	if(gInsExtend > gInsOpen) {
		cerr << "Warning: Read gap extension penalty (" << gInsExtend << ") is greater than gap open penalty (" << gInsOpen << ")." << endl;
		cerr << "Adjusting extension penalty to be equal to open penalty." << endl;
		gInsExtend = gInsOpen;
	}
	if(gDelExtend > gDelOpen) {
		cerr << "Warning: Reference gap extension penalty (" << gDelExtend << ") is greater than gap open penalty (" << gDelOpen << ")." << endl;
		cerr << "Adjusting extension penalty to be equal to open penalty." << endl;
		gDelExtend = gDelOpen;
	}
	if(gGapBarrier < 1) {
		cerr << "Error: --gbar must be >= 1; --gbar is " << gGapBarrier << endl;
	}
	if(gColor && gColorExEnds) {
		gGapBarrier++;
	}
}

static const char *argv0 = NULL;

#ifdef CHUD_PROFILING
#define CHUD_START() chudStartRemotePerfMonitor("Bowtie");
#define CHUD_STOP()  chudStopRemotePerfMonitor();
#else
#define CHUD_START()
#define CHUD_STOP()
#endif

/// Create a PatternSourcePerThread for the current thread according
/// to the global params and return a pointer to it
static PatternSourcePerThreadFactory*
createPatsrcFactory(PairedPatternSource& _patsrc, int tid) {
	PatternSourcePerThreadFactory *patsrcFact;
	patsrcFact = new WrappedPatternSourcePerThreadFactory(_patsrc);
	assert(patsrcFact != NULL);
	return patsrcFact;
}

/**
 * Allocate a HitSinkPerThreadFactory on the heap according to the
 * global params and return a pointer to it.
 */
static HitSinkPerThreadFactory*
createSinkFactory(HitSink& _sink) {
	HitSinkPerThreadFactory *sink = NULL;
	if(format == INPUT_CHAIN) {
		sink = new ChainingHitSinkPerThreadFactory(_sink, allHits ? 0xffffffff : khits, mhits, strata);
	} else if(!strata) {
		// Unstratified
		if(!allHits) {
			// First N good; "good" inherently ignores strata
			sink = new NGoodHitSinkPerThreadFactory(_sink, khits, mhits);
		} else {
			// All hits, spanning strata
			sink = new AllHitSinkPerThreadFactory(_sink, mhits);
		}
	} else {
		// Stratified
		if(!allHits) {
			// Buffer best hits, assuming they're arriving in best-
			// to-worst order
			sink = new NBestFirstStratHitSinkPerThreadFactory(_sink, khits, mhits);
		} else {
			// Buffer best hits, assuming they're arriving in best-
			// to-worst order
			sink = new NBestFirstStratHitSinkPerThreadFactory(_sink, 0xffffffff/2, mhits);
		}
	}
	assert(sink != NULL);
	return sink;
}


static PairedPatternSource*   exactSearch_patsrc;
static HitSink*               exactSearch_sink;
//static ACHitSink*             exactSearch_acsink;
static Ebwt*                  exactSearch_ebwt;
static EList<SString<char> >* exactSearch_os;
static BitPairReference*      exactSearch_refs;

/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
static void *exactSearchWorkerStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt& ebwt                   = *exactSearch_ebwt;
	EList<SString<char> >& os    = *exactSearch_os;
	BitPairReference* refs       =  exactSearch_refs;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);

	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);
	UnpairedExactAlignerV1Factory alSEfact(
			ebwt,
			NULL,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs,
			os,
			seed);
	PairedExactAlignerV1Factory alPEfact(
			ebwt,
			NULL,
			_sink,
			*sinkFact,
			dontReconcileMates,
			mhits,       // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs, os,
			seed);
	{
		MixedMultiAligner multi(
				1,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) pthread_exit(NULL);
#endif
	return NULL;
}

#ifdef BOWTIE_PTHREADS
#define PTHREAD_ATTRS (PTHREAD_CREATE_JOINABLE | PTHREAD_CREATE_DETACHED)
#endif

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static void exactSearch(PairedPatternSource& _patsrc,
                        HitSink& _sink,
                        Ebwt& ebwt,
                        EList<SString<char> >& os)
{
	exactSearch_patsrc = &_patsrc;
	exactSearch_sink   = &_sink;
	exactSearch_ebwt   = &ebwt;
	exactSearch_os     = &os;

	assert(!ebwt.isInMemory());
	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwt.loadIntoMemory(
			gColor ? 1 : 0, // we expect index to be colorspace? (for sanity)
			0,              // it's exact search, so we don't care how reverse is constructed
			true,           // load the SA sample
			true,           // load the ftab
			true,           // load the rstarts
			!noRefNames,    
			startVerbose);
	}

	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	// Do we need to load the 2-bit-encoded reference string?  Only if:
	//  1. We're doing Aho-Corasick to find seed hits
	//  2. We're doing colorspace (and we need nuc ref for decoding)
	//  3. We're doing paried-end (and we need ref to scan for opposite mate)
	if(gColor || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(
			adjustedEbwtFileBase, gColor, sanityCheck, NULL, &os, false, useMm,
			useShmem, mmSweep, gVerbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	exactSearch_refs   = refs;

#ifdef BOWTIE_PTHREADS
	AutoArray<pthread_t> threads(nthreads-1);
	AutoArray<int> tids(nthreads-1);
#endif
	CHUD_START();
	{
		Timer _t(cerr, "Time for 0-mismatch search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			createThread(&threads[i],
						 exactSearchWorkerStateful,
						 (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		exactSearchWorkerStateful((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) joinThread(threads[i]);
#endif
	}
	if(refs != NULL) delete refs;
}

/**
 * Search through a pair of Ebwt indexes, one for the forward direction
 * and one for the backward direction, for exact end-to-end hits and 1-
 * mismatch end-to-end hits.  In my experience, this is slightly faster
 * than Maq (default) mode with the -n 1 option.
 *
 * Forward Ebwt (ebwtFw) is already loaded into memory and backward
 * Ebwt (ebwtBw) is not loaded into memory.
 */
static PairedPatternSource*           mismatchSearch_patsrc;
static HitSink*                       mismatchSearch_sink;
static Ebwt*            mismatchSearch_ebwtFw;
static Ebwt*            mismatchSearch_ebwtBw;
static EList<SString<char> >*         mismatchSearch_os;
static BitPairReference*              mismatchSearch_refs;

/**
 * A statefulness-aware worker driver.  Uses Unpaired/Paired1mmAlignerV1.
 */
static void *mismatchSearchWorkerFullStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc = *mismatchSearch_patsrc;
	HitSink&               _sink   = *mismatchSearch_sink;
	Ebwt&    ebwtFw  = *mismatchSearch_ebwtFw;
	Ebwt&    ebwtBw  = *mismatchSearch_ebwtBw;
	EList<SString<char> >& os      = *mismatchSearch_os;
	BitPairReference*      refs    =  mismatchSearch_refs;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);

	Unpaired1mmAlignerV1Factory alSEfact(
			ebwtFw,
			&ebwtBw,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs,
			os,
			seed);
	Paired1mmAlignerV1Factory alPEfact(
			ebwtFw,
			&ebwtBw,
			_sink,
			*sinkFact,
			dontReconcileMates,
			mhits,     // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs, os,
			seed);
	{
		MixedMultiAligner multi(
				1,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) pthread_exit(NULL);
#endif
	return NULL;
}

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static void mismatchSearchFull(PairedPatternSource& _patsrc,
                               HitSink& _sink,
                               Ebwt& ebwtFw,
                               Ebwt& ebwtBw,
                               EList<SString<char> >& os)
{
	mismatchSearch_patsrc       = &_patsrc;
	mismatchSearch_sink         = &_sink;
	mismatchSearch_ebwtFw       = &ebwtFw;
	mismatchSearch_ebwtBw       = &ebwtBw;
	mismatchSearch_os           = &os;

	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // it's not bidirectional search, so we don't care how reverse is constructed
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
	}
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // it's not bidirectional search, so we don't care how reverse is constructed
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(gColor || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, gColor, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, gVerbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	mismatchSearch_refs = refs;

#ifdef BOWTIE_PTHREADS
	// Allocate structures for threads
	AutoArray<pthread_t> threads(nthreads-1);
	AutoArray<int> tids(nthreads-1);
#endif
    CHUD_START();
    {
		Timer _t(cerr, "Time for 1-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			createThread(&threads[i], mismatchSearchWorkerFullStateful, (void *)&tids[i]);
		}
#endif
		// Go to town
		int tmp = 0;
		mismatchSearchWorkerFullStateful((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) joinThread(threads[i]);
#endif
	}
	if(refs != NULL) delete refs;
}

static PairedPatternSource*           twoOrThreeMismatchSearch_patsrc;
static HitSink*                       twoOrThreeMismatchSearch_sink;
static Ebwt*            twoOrThreeMismatchSearch_ebwtFw;
static Ebwt*            twoOrThreeMismatchSearch_ebwtBw;
static EList<SString<char> >*         twoOrThreeMismatchSearch_os;
static bool                           twoOrThreeMismatchSearch_two;
static BitPairReference*              twoOrThreeMismatchSearch_refs;

/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
static void *twoOrThreeMismatchSearchWorkerStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc = *twoOrThreeMismatchSearch_patsrc;
	HitSink&               _sink   = *twoOrThreeMismatchSearch_sink;
	Ebwt&    ebwtFw  = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt&    ebwtBw  = *twoOrThreeMismatchSearch_ebwtBw;
	EList<SString<char> >& os      = *twoOrThreeMismatchSearch_os;
	BitPairReference*      refs    =  twoOrThreeMismatchSearch_refs;
	static bool            two     =  twoOrThreeMismatchSearch_two;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);

	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);
	Unpaired23mmAlignerV1Factory alSEfact(
			ebwtFw,
			&ebwtBw,
			two,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs,
			os,
			seed);
	Paired23mmAlignerV1Factory alPEfact(
			ebwtFw,
			&ebwtBw,
			two,
			_sink,
			*sinkFact,
			dontReconcileMates,
			mhits,       // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs, os,
			seed);
	{
		MixedMultiAligner multi(
				1,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) pthread_exit(NULL);
#endif
	return NULL;
}

static void twoOrThreeMismatchSearchFull(
		PairedPatternSource& _patsrc,   /// pattern source
		HitSink& _sink,                 /// hit sink
		Ebwt& ebwtFw,                   /// index of original text
		Ebwt& ebwtBw,                   /// index of mirror text
		EList<SString<char> >& os,      /// text strings, if available (empty otherwise)
		bool two = true)                /// true -> 2, false -> 3
{
	// Global initialization
	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // it's not bidirectional search, so we don't care how reverse is constructed
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
	}
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // it's not bidirectional search, so we don't care how reverse is constructed
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(gColor || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, gColor, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, gVerbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	twoOrThreeMismatchSearch_refs     = refs;
	twoOrThreeMismatchSearch_patsrc   = &_patsrc;
	twoOrThreeMismatchSearch_sink     = &_sink;
	twoOrThreeMismatchSearch_ebwtFw   = &ebwtFw;
	twoOrThreeMismatchSearch_ebwtBw   = &ebwtBw;
	twoOrThreeMismatchSearch_os       = &os;
	twoOrThreeMismatchSearch_two      = two;

#ifdef BOWTIE_PTHREADS
	AutoArray<pthread_t> threads(nthreads-1);
	AutoArray<int> tids(nthreads-1);
#endif
	CHUD_START();
	{
		Timer _t(cerr, "End-to-end 2/3-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			createThread(&threads[i], twoOrThreeMismatchSearchWorkerStateful, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		twoOrThreeMismatchSearchWorkerStateful((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) joinThread(threads[i]);
#endif
	}
	if(refs != NULL) delete refs;
	return;
}

static PairedPatternSource*     multiseed_patsrc;
static Ebwt*                    multiseed_ebwtFw;
static Ebwt*                    multiseed_ebwtBw;
static Penalties*               multiseed_pens;
static EList<Seed>*             multiseed_seeds;
static BitPairReference*        multiseed_refs;
static AlignmentCache*          multiseed_sc; // seed cache
static AlnSink*                 multiseed_msink;
static OutFileBuf*              multiseed_metricsOfb;

static EList<ReadCounterSink*>* multiseed_readCounterSink;
static EList<SeedHitSink*>*     multiseed_seedHitSink;
static EList<SeedCounterSink*>* multiseed_seedCounterSink;
static EList<SeedActionSink*>*  multiseed_seedActionSink;
static EList<SwCounterSink*>*   multiseed_swCounterSink;
static EList<SwActionSink*>*    multiseed_swActionSink;

/**
 * Metrics for measuring the work done by the outer read alignment
 * loop.
 */
struct OuterLoopMetrics {

	OuterLoopMetrics() { reset(); MUTEX_INIT(lock); }

	/**
	 * Set all counters to 0.
	 */
	void reset() {
		reads = bases = srreads = srbases =
		freads = fbases = ureads = ubases = 0;
	}

	/**
	 * Sum the counters in m in with the conters in this object.  This
	 * is the only safe way to update an OuterLoopMetrics that's shared
	 * by multiple threads.
	 */
	void merge(
		const OuterLoopMetrics& m,
		bool getLock = false)
	{
		ThreadSafe ts(&lock, getLock);
		reads += m.reads;
		bases += m.bases;
		srreads += m.srreads;
		srbases += m.srbases;
		freads += m.freads;
		fbases += m.fbases;
		ureads += m.ureads;
		ubases += m.ubases;
	}

	uint64_t reads;   // total reads
	uint64_t bases;   // total bases
	uint64_t srreads; // same-read reads
	uint64_t srbases; // same-read bases
	uint64_t freads;  // filtered reads
	uint64_t fbases;  // filtered bases
	uint64_t ureads;  // unfiltered reads
	uint64_t ubases;  // unfiltered bases
	MUTEX_T lock;
};

/**
 * Collection of all relevant performance metrics when aligning in
 * multiseed mode.
 */
struct PerfMetrics {

	PerfMetrics() : first(true) { reset(); }

	/**
	 * Set all counters to 0.
	 */
	void reset() {
		olm.reset();
		sdm.reset();
		wlm.reset();
		swm.reset();
		rpm.reset();
		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmu.reset();
		rpmu.reset();
	}

	/**
	 * Merge a set of specific metrics into this object.
	 */
	void merge(
		const OuterLoopMetrics *ol,
		const SeedSearchMetrics *sd,
		const WalkMetrics *wl,
		const SwMetrics *sw,
		const ReportingMetrics *rm,
		bool getLock)
	{
		ThreadSafe ts(&lock, getLock);
		if(ol != NULL) {
			olmu.merge(*ol, false);
		}
		if(sd != NULL) {
			sdmu.merge(*sd, false);
		}
		if(wl != NULL) {
			wlmu.merge(*wl, false);
		}
		if(sw != NULL) {
			swmu.merge(*sw, false);
		}
		if(rm != NULL) {
			rpmu.merge(*rm, false);
		}
	}

	/**
	 * Reports a matrix of results, incl. column labels, to an
	 * OutFileBuf.  Optionally also send results to stderr (unbuffered).
	 */
	void reportInterval(
		OutFileBuf* o,      // file to send output to
		bool metricsStderr, // additionally output to stderr?
		bool sync = true)
	{
		ThreadSafe ts(&lock, sync);
		time_t curtime = time(0);
		char buf[1024];
		if(first) {
			const char *str =
				/*  1 */ "Time"           "\t"
				/*  2 */ "Read"           "\t"
				/*  3 */ "Base"           "\t"
				/*  4 */ "SameRead"       "\t"
				/*  5 */ "SameReadBase"   "\t"
				/*  6 */ "UnfilteredRead" "\t"
				/*  7 */ "UnfilteredBase" "\t"
				/*  8 */ "Aligned"        "\t"
				/*  9 */ "Unaligned"      "\t"
				/* 10 */ "Maxed"          "\t"
				/* 11 */ "SeedSearch"     "\t"
				/* 12 */ "IntraSCacheHit" "\t"
				/* 13 */ "InterSCacheHit" "\t"
				/* 14 */ "OutOfMemory"    "\t"
				/* 15 */ "AlBWOp"         "\t"
				/* 16 */ "AlBWBranch"     "\t"
				/* 17 */ "ResBWOp"        "\t"
				/* 18 */ "ResBWBranch"    "\t"
				/* 19 */ "ResResolve"     "\t"
				/* 20 */ "ResReport"      "\t"
				/* 21 */ "RedundantSHit"  "\t"
				/* 22 */ "DynProg"        "\t"
				/* 23 */ "DynProgRow"     "\t"
				/* 24 */ "DynProgRowSkip" "\t"
				/* 25 */ "DynProgSucc"    "\t"
				/* 26 */ "DynProgFail"    "\t"
				/* 27 */ "DynProgCell"    "\t"
				/* 28 */ "DynProgBt"      "\t"    
				/* 29 */ "MemPeak"        "\t"
				/* 30 */ "UncatMemPeak"   "\t" // 0
				/* 31 */ "EbwtMemPeak"    "\t" // EBWT_CAT
				/* 32 */ "CacheMemPeak"   "\t" // CA_CAT
				/* 33 */ "ResolveMemPeak" "\t" // GW_CAT
				/* 34 */ "AlignMemPeak"   "\t" // AL_CAT
				/* 35 */ "DPMemPeak"      "\t" // SW_CAT
				/* 36 */ "MiscMemPeak"    "\t" // MISC_CAT
				/* 37 */ "DebugMemPeak"   "\t" // DEBUG_CAT
				"\n";
			if(o != NULL) o->writeChars(str);
			if(metricsStderr) cerr << str;
			first = false;
		}
		// 1. Current time in secs
		itoa10<time_t>(curtime, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 2. Reads
		itoa10<uint64_t>(olmu.reads, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 3. Bases
		itoa10<uint64_t>(olmu.bases, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 4. Same-read reads
		itoa10<uint64_t>(olmu.srreads, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 5. Same-read bases
		itoa10<uint64_t>(olmu.srbases, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 6. Unfiltered reads
		itoa10<uint64_t>(olmu.ureads, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 7. Unfiltered bases
		itoa10<uint64_t>(olmu.ubases, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 8. Aligned reads
		itoa10<uint64_t>(rpmu.al, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 9. Unaligned reads
		itoa10<uint64_t>(rpmu.unal, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 10. Non-uniquely aligning reads
		itoa10<uint64_t>(rpmu.max, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		
		// 11. Seed searches
		itoa10<uint64_t>(sdmu.seedsearch, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 12. Hits in 'current' cache
		itoa10<uint64_t>(sdmu.intrahit, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 13. Hits in 'local' cache
		itoa10<uint64_t>(sdmu.interhit, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 14. Out of memory
		itoa10<uint64_t>(sdmu.ooms, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 15. Burrows-Wheeler ops in aligner
		itoa10<uint64_t>(sdmu.bwops, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 16. Burrows-Wheeler branches (edits) in aligner
		itoa10<uint64_t>(sdmu.bweds, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 17. Burrows-Wheeler ops in resolver
		itoa10<uint64_t>(wlmu.bwops, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 18. Burrows-Wheeler branches in resolver
		itoa10<uint64_t>(wlmu.branches, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 19. Offset resolutions
		itoa10<uint64_t>(wlmu.resolutions, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 20. Offset reports
		itoa10<uint64_t>(wlmu.reports, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 21. Redundant seed hit
		itoa10<uint64_t>(swmu.rshit, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 22. Dynamic programming problems
		itoa10<uint64_t>(swmu.sws, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }

		// 23. Dynamic programming rows completed
		itoa10<uint64_t>(swmu.swrows, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 24. Dynamic programming rows skipped
		itoa10<uint64_t>(swmu.swskiprows, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 25. Dynamic programming successes
		itoa10<uint64_t>(swmu.swsucc, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 26. Dynamic programming failures
		itoa10<uint64_t>(swmu.swfail, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 27. Dynamic programming CUPs
		itoa10<uint64_t>(swmu.swcups, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 28. Dynamic programming backtraces
		itoa10<uint64_t>(swmu.swbts, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 29. Overall memory peak
		itoa10<size_t>(gMemTally.peak() >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 30. Uncategorized memory peak
		itoa10<size_t>(gMemTally.peak(0) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 31. Ebwt memory peak
		itoa10<size_t>(gMemTally.peak(EBWT_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 32. Cache memory peak
		itoa10<size_t>(gMemTally.peak(CA_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 33. Resolver memory peak
		itoa10<size_t>(gMemTally.peak(GW_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 34. Seed aligner memory peak
		itoa10<size_t>(gMemTally.peak(AL_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 35. Dynamic programming aligner memory peak
		itoa10<size_t>(gMemTally.peak(SW_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 36. Miscellaneous memory peak
		itoa10<size_t>(gMemTally.peak(MISC_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf << '\t';
		if(o != NULL) { o->writeChars(buf); o->write('\t'); }
		// 37. Debug memory peak
		itoa10<size_t>(gMemTally.peak(DEBUG_CAT) >> 20, buf);
		if(metricsStderr) cerr << buf;
		if(o != NULL) { o->writeChars(buf); }

		if(o != NULL) { o->write('\n'); }
		if(metricsStderr) cerr << endl;
		finishReport();
	}
	
	void finishReport() {
		olm.merge(olmu, false);
		sdm.merge(sdmu, false);
		wlm.merge(wlmu, false);
		swm.merge(swmu, false);
		rpm.merge(rpmu, false);
		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmu.reset();
		rpmu.reset();
	}

	// Total over the whole job
	OuterLoopMetrics  olm;   // overall metrics
	SeedSearchMetrics sdm;   // metrics related to seed alignment
	WalkMetrics       wlm;   // metrics related to walking left (i.e. resolving reference offsets)
	SwMetrics         swm;   // metrics related to dynamic prog alignment
	ReportingMetrics  rpm;   // metrics related to reporting

	// Just since the last update
	OuterLoopMetrics  olmu;  // overall metrics
	SeedSearchMetrics sdmu;  // metrics related to seed alignment
	WalkMetrics       wlmu;  // metrics related to walking left (i.e. resolving reference offsets)
	SwMetrics         swmu;  // metrics related to dynamic prog alignment
	ReportingMetrics  rpmu;  // metrics related to reporting

	MUTEX_T           lock;  // lock for when one ob
	bool              first; // yet to print first line?
	time_t            lastElapsed; // used in reportInterval to measure time since last call
};

static PerfMetrics metrics;

// Cyclic rotations
#define ROTL(n, x) (((x) << (n)) | ((x) >> (32-n)))
#define ROTR(n, x) (((x) >> (n)) | ((x) << (32-n)))

/**
 * Called once per thread.  Sets up per-thread pointers to the shared global
 * data structures, creates per-thread structures, then enters the alignment
 * loop.  The general flow of the alignment loop is:
 *
 * - If it's been a while and we're the master thread, report some alignment
 *   metrics
 * - Get the next read/pair
 * - Check if this read/pair is identical to the previous
 *   + If identical, check whether we can skip any or all alignment stages.  If
 *     we can skip all stages, report the result immediately and move to next
 *     read/pair
 *   + If not identical, continue
 * - 
 */
static void* multiseedSearchWorker(void *vp) {
	int tid = *((int*)vp);
	assert(multiseed_ebwtFw != NULL);
	assert(multiseedMms == 0 || multiseed_ebwtBw != NULL);
	PairedPatternSource&    patsrc   = *multiseed_patsrc;
	const Ebwt&             ebwtFw   = *multiseed_ebwtFw;
	const Ebwt&             ebwtBw   = *multiseed_ebwtBw;
	const Penalties&        pens     = *multiseed_pens;
	const EList<Seed>&      seeds    = *multiseed_seeds;
	const BitPairReference& ref      = *multiseed_refs;
	AlignmentCache&         scShared = *multiseed_sc;
	AlnSink&                msink    = *multiseed_msink;
	OutFileBuf*             metricsOfb = multiseed_metricsOfb;

	// Sinks: these are so that we can print tables encoding counts for
	// events of interest on a per-read, per-seed, per-join, or per-SW
	// level.  These in turn can be used to diagnose performance
	// problems, or generally characterize performance.
	
	// readCounterSink: for each read, keep a record of counts related
	// to all phases of alignment: seed alignment, joining, SW
	EList<ReadCounterSink*>& readCounterSink = *multiseed_readCounterSink;

	// seedHitSink: for each seed hit, keep a record of counts
	EList<SeedHitSink*>&     seedHitSink     = *multiseed_seedHitSink;
	// seedCounterSink: for each seed, keep a record of counts related
	// to the seed alignment
	EList<SeedCounterSink*>& seedCounterSink = *multiseed_seedCounterSink;
	// seedActionSink: keep a detailed record of seed alignment actions
	// taken per read
	EList<SeedActionSink*>&  seedActionSink  = *multiseed_seedActionSink;

	// swCounterSink: for each SW, keep a record of counts related to
	// the SW alignment
	EList<SwCounterSink*>&   swCounterSink   = *multiseed_swCounterSink;
	// swActionSink: keep a detailed record of SW alignment actions
	// taken per read
	EList<SwActionSink*>&    swActionSink    = *multiseed_swActionSink;
	
	//const BitPairReference& refs   = *multiseed_refs;
	auto_ptr<PatternSourcePerThreadFactory> patsrcFact(createPatsrcFactory(patsrc, tid));
	auto_ptr<PatternSourcePerThread> ps(patsrcFact->create());
	
	// Thread-local cache for seed alignments
	AlignmentCache scLocal(msNoCache ? 0 : (seedCacheLocalMB * 1024 * 1024), false);
	// Thread-local cache for current seed alignments
	AlignmentCache scCurrent(seedCacheCurrentMB * 1024 * 1024, false);
	
	// Interfaces for alignment and seed caches
	AlignmentCacheIface sc(
		&scCurrent,
		msNoCache ? NULL : &scLocal,
		msNoCache ? NULL : &scShared);
	
	// Instantiate an object for holding reporting-related parameters.
	ReportingParams rp(
		(allHits ? std::numeric_limits<THitInt>::max() : khits), // -k
		mhits,             // -m/-M
		0,                 // penalty gap (not used now)
		sampleMax,         // true -> -M was specified, otherwise assume -m
		gReportDiscordant, // report discordang paired-end alignments?
		gReportMixed);     // report unpaired alignments for paired reads?
	
	// Make a per-thread wrapper for the global MHitSink object.
	AlnSinkWrap msinkwrap(msink, rp);

	SeedAligner al;
	SwDriver sd;
	SwAligner sw, osw;
	SwParams pa;
	SeedResults shs[2];
	QVal *qv;
	GroupWalk gws[2];
	OuterLoopMetrics olm;
	SeedSearchMetrics sdm;
	WalkMetrics wlm;
	SwMetrics swm;
	ReportingMetrics rpm;
	RandomSource rnd;

	int pepolFlag;
	if(gMate1fw && gMate2fw) {
		pepolFlag = PE_POLICY_FF;
	} else if(gMate1fw && !gMate2fw) {
		pepolFlag = PE_POLICY_FR;
	} else if(!gMate1fw && gMate2fw) {
		pepolFlag = PE_POLICY_RF;
	} else {
		pepolFlag = PE_POLICY_RR;
	}
	assert_geq(gMaxInsert, gMinInsert);
	assert_geq(gMinInsert, 0);
	PairedEndPolicy pepol(
		pepolFlag,
		gMaxInsert,
		gMinInsert,
		gLocalAlign,
		gFlippedMatesOK,
		gDovetailMatesOK,
		gContainMatesOK,
		gOlapMatesOK,
		gExpandToFrag);
	
	// Used by thread with threadid == 1 to measure time elapsed
	time_t iTime = time(0);

	int mergei = 0;
	int mergeival = 16;
	while(true) {
		ps->nextReadPair();
		TReadId patid = ps->patid();
		if(patid < qUpto && !ps->empty()) {
			// Align this read/pair
			bool retry = true;
			if(metricsIval > 0 &&
			   (metricsOfb != NULL || metricsStderr) &&
			   ++mergei == mergeival)
			{
				// Update global metrics, in a synchronized manner if needed
				metrics.merge(&olm, &sdm, &wlm, &swm, &rpm, nthreads > 1);
				olm.reset();
				sdm.reset();
				wlm.reset();
				swm.reset();
				rpm.reset();
				mergei = 0;
				// Check if a progress message should be printed
				if(tid == 0) {
					// Only thread 1 prints progress messages
					time_t curTime = time(0);
					if(curTime - iTime >= metricsIval) {
						metrics.reportInterval(metricsOfb, metricsStderr, true);
						iTime = curTime;
					}
				}
			}
			while(retry) {
				qv = NULL;
				retry = false;
				assert_eq(ps->bufa().color, gColor);
				sc.nextRead();
				olm.reads++;
				assert(!sc.aligning());
				// NB: read may be either unpaired or paired-end at this point
				bool pair = ps->paired();
				const size_t rdlen1 = ps->bufa().length();
				const size_t rdlen2 = pair ? ps->bufb().length() : 0;
				olm.bases += (rdlen1 + rdlen2);
				// Check if read is identical to previous read
				rnd.init(ROTL(ps->bufa().seed, 5));
				int skipStages = msinkwrap.nextRead(
					&ps->bufa(),
					pair ? &ps->bufb() : NULL,
					patid,
					pens.qualitiesMatter());
				assert(msinkwrap.inited());
				if(skipStages == -1) {
					// Read or pair is identical to previous.  Re-report from
					// the msinkwrap immediately.
					olm.srreads++;
					olm.srbases += (rdlen1 + rdlen2);
					rnd.init(ROTL(ps->bufa().seed, 20));
					msinkwrap.finishRead(
						NULL,                 // seed results for mate 1
						NULL,                 // seed results for mate 2
						rnd,                  // pseudo-random generator
						rpm,                  // reporting metrics
						!seedSummaryOnly,     // suppress seed summaries?
						seedSummaryOnly);     // suppress alignments?
					break; // next read
				}
				bool nfilt[2];
				pens.nFilterPair(
					&ps->bufa().patFw,
					pair ? &ps->bufb().patFw : NULL,
					nfilt[0],
					nfilt[1]);
				// Re-set 'pair': true only if both mates passed the filter
				pair = nfilt[0] && nfilt[1];
				size_t rdlens[2]   = { rdlen1, rdlen2 };
				const Read* rds[2] = { &ps->bufa(), &ps->bufb() };
				// For each mate...
				// TODO: Examine mates in a random order instead of #1 then #2
				assert(msinkwrap.empty());
				sd.nextRead(pair);
				for(size_t mate = 0; mate < 2; mate++) {
					if(!nfilt[mate]) {
						// Mate was rejected by N filter
						olm.freads++;               // reads filtered out
						olm.fbases += rdlens[mate]; // bases filtered out
						continue; // on to next mate
					}
					if(msinkwrap.state().doneWithMate(mate == 0)) {
						// Done with this mate
						continue;
					}
					// Passed N filter
					assert_gt(rds[mate]->length(), 0);
					assert(!msinkwrap.maxed());
					assert(msinkwrap.repOk());
					olm.ureads++;               // reads passing filter
					olm.ubases += rdlens[mate]; // bases passing filter
					QKey qkr(rds[mate]->patFw);
					rnd.init(ROTL(rds[mate]->seed, 10));
					// Seed search
					shs[mate].clear(); // clear seed hits
					assert(shs[mate].repOk(&sc.current()));
					int interval = multiseedPeriod;
					if(interval == -1) {
						// Calculate the seed interval as a
						// function of the read length
						float x = (float)rdlens[mate];
						if(multiseedIvalType == SEED_IVAL_SQUARE_ROOT) {
							x = pow(x, 0.5f);
						} else if(multiseedIvalType == SEED_IVAL_CUBE_ROOT) {
							x = pow(x, 0.333333f);
						}
						interval = (int)((multiseedIvalA * x) + multiseedIvalB + 0.5f);
						if(interval < 1) interval = 1;
					}
					assert_geq(interval, 1);
					// Instantiate the seeds
					std::pair<int, int> inst = al.instantiateSeeds(
						seeds,                // seed templates
						interval,             // interval between seeds
						*rds[mate],           // read
						pens,                 // edit penalties
						nCeilConst,           // ceiling on # Ns w/r/t read length, constant coeff
						nCeilLinear,          // ceiling on # Ns w/r/t read length, linear coeff
						sc,                   // seed cache
						shs[mate],            // seed hits
						sdm);                 // metrics
					assert(shs[mate].repOk(&sc.current()));
					if(inst.first + inst.second == 0) {
						continue; // on to next mate
					}
					// Align seeds
					al.searchAllSeeds(
						seeds,            // search seeds
						&ebwtFw,          // BWT index
						&ebwtBw,          // BWT' index
						*rds[mate],       // read
						pens,             // edit penalties
						sc,               // alignment cache
						shs[mate],        // store seed hits here
						sdm,              // metrics
						&readCounterSink, // send counter summary for each read to this sink
						&seedHitSink,     // send seed hits to this sink
						&seedCounterSink, // send counter summary for each seed to this sink
						&seedActionSink); // send search action list for each read to this sink
					assert_eq(0, sdm.ooms);
					assert(shs[mate].repOk(&sc.current()));
					if(!seedSummaryOnly) {
						// If there aren't any seed hits...
						if(shs[mate].empty()) {
							continue; // on to the next mate
						}
						// Sort seed hits into ranks
						shs[mate].sort();
						// Calculate the penalty ceiling for the read
						int penceil = Constraint::instantiate(
							rdlens[mate],
							costCeilConst,
							costCeilLinear);
						int nceil = (int)pens.nCeil(rdlens[mate]);
						assert_geq(penceil, 0);
						bool done = false;
						if(pair) {
							int openceil = Constraint::instantiate(
								rdlens[mate ^ 1],
								costCeilConst,
								costCeilLinear);
							int onceil = (int)pens.nCeil(rdlens[mate ^ 1]);
							// Paired-end dynamic programming driver
							done = sd.extendSeedsPaired(
								*rds[mate],     // mate to align as anchor
								*rds[mate ^ 1], // mate to align as opp.
								mate == 0,      // anchor is mate 1?
								gColor,         // colorspace?
								shs[mate],      // seed hits for anchor
								ebwtFw,         // bowtie index
								ref,            // packed reference strings
								gws[mate],      // walk left for anchor
								sw,             // dyn prog aligner, anchor
								osw,            // dyn prog aligner, opposite
								pa,             // parameters for DP
								pens,           // penalties for edits
								pepol,          // paired-end policy
								multiseedMms,   // # mms allowed in a seed
								multiseedLen,   // length of a seed
								interval,       // interval between seeds
								penceil,        // penalty ceil for anchor
								openceil,       // penalty ceil for opp.
								nceil,          // N ceil for anchor
								onceil,         // N ceil for opp.
								maxposs,        // max off/orient combos
								maxrows,        // max offset resolutions
								sc,             // seed alignment cache
								rnd,            // pseudo-random source
								wlm,            // group walk left metrics
								swm,            // dynamic prog metrics
								rpm,            // reporting metrics
								&msinkwrap,     // for organizing hits
								true,           // seek mate immediately
								true,           // report hits once found
								gReportDiscordant,// look for discordant alns?
								gReportMixed,   // look for unpaired alns?
								&swCounterSink, // send counter info here
								&swActionSink); // send action info here
							
							// Might be done, but just with this mate
						} else {
							// Unpaired dynamic programming driver
							done = sd.extendSeeds(
								*rds[mate],     // read
								gColor,         // colorspace?
								shs[mate],      // seed hits
								ebwtFw,         // bowtie index
								ref,            // packed reference strings
								gws[mate],      // group walk left
								sw,             // dynamic prog aligner
								pa,             // parameters for DP
								pens,           // penalties for edits
								multiseedMms,   // # mms allowed in a seed
								multiseedLen,   // length of a seed
								interval,       // interval between seeds
								penceil,        // penalty ceiling
								nceil,          // N ceil for anchor
								maxposs,        // max off/orient combos
								maxrows,        // max offset resolutions
								sc,             // seed alignment cache
								rnd,            // pseudo-random source
								wlm,            // group walk left metrics
								swm,            // dynamic prog metrics
								rpm,            // reporting metrics
								&msinkwrap,     // for organizing hits
								true,           // report hits once found
								&swCounterSink, // send counter info here
								&swActionSink); // send action info here
						}
						// Are we done with this read/pair?
						if(done) {
							break; // ...break out of the loop over mates
						}
					} // if(!seedSummaryOnly)
				} // for(size_t mate = 0; mate < 2; mate++)
				
				// Commit and report paired-end/unpaired alignments
				uint32_t seed = rds[0]->seed ^ rds[1]->seed;
				rnd.init(ROTL(seed, 20));
				msinkwrap.finishRead(
					&shs[0],              // seed results for mate 1
					&shs[1],              // seed results for mate 2
					rnd,                  // pseudo-random generator
					rpm,                  // reporting metrics
					!seedSummaryOnly,     // suppress seed summaries?
					seedSummaryOnly);     // suppress alignments?
				assert(!retry || msinkwrap.empty());
				
				// Add local metrics to global metrics in a synchronized fashion
				gws[0].reset();  // reset the group walk-left object
				gws[1].reset();  // reset the group walk-left object
			} // while(retry)
		} // if(patid < qUpto && !ps->empty())
		else {
			break;
		}
	} // while(true)
	
	// One last metrics merge
	metrics.merge(&olm, &sdm, &wlm, &swm, &rpm, nthreads > 1);
#ifdef BOWTIE_PTHREADS
	if(tid > 0) { pthread_exit(NULL); }
#endif
	return NULL;
}

/**
 * Called once per alignment job.  Sets up global pointers to the
 * shared global data structures, creates per-thread structures, then
 * enters the search loop.
 */
static void multiseedSearch(
	Penalties& pens,
	EList<Seed>& seeds,
	PairedPatternSource& patsrc,  // pattern source
	AlnSink& msink,             // hit sink
	Ebwt& ebwtFw,                 // index of original text
	Ebwt& ebwtBw,                 // index of mirror text
	EList<SeedHitSink*>& seedHitSink,
	EList<SeedCounterSink*>& seedCounterSink,
	EList<SeedActionSink*>&  seedActionSink,
	OutFileBuf *metricsOfb)
{
	multiseed_patsrc = &patsrc;
	multiseed_msink  = &msink;
	multiseed_ebwtFw = &ebwtFw;
	multiseed_ebwtBw = &ebwtBw;
	multiseed_pens   = &pens;
	multiseed_seeds  = &seeds;
	multiseed_metricsOfb      = metricsOfb;
	multiseed_seedHitSink     = &seedHitSink;
	multiseed_seedCounterSink = &seedCounterSink;
	multiseed_seedActionSink  = &seedActionSink;
	Timer *_t = new Timer(cerr, "Time loading reference: ", timing);
	auto_ptr<BitPairReference> refs(
		new BitPairReference(
			adjustedEbwtFileBase,
			gColor,
			sanityCheck,
			NULL,
			NULL,
			false,
			useMm,
			useShmem,
			mmSweep,
			gVerbose,
			startVerbose)
	);
	delete _t;
	if(!refs->loaded()) throw 1;
	multiseed_refs = refs.get();
#ifdef BOWTIE_PTHREADS
	AutoArray<pthread_t> threads(nthreads-1);
	AutoArray<int> tids(nthreads-1);
#endif
	{
		// Load the other half of the index into memory
		assert(!ebwtFw.isInMemory());
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(
			gColor ? 1 : 0, // colorspace?
			-1, // not the reverse index
			true,         // load SA samp? (yes, need forward index's SA samp)
			true,         // load ftab (in forward index)
			true,         // load rstarts (in forward index)
			!noRefNames,  // load names?
			startVerbose);
	}
	if(multiseedMms > 0) {
		// Load the other half of the index into memory
		assert(!ebwtBw.isInMemory());
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(
			gColor ? 1 : 0, // colorspace?
			// It's bidirectional search, so we need the reverse to be
			// constructed as the reverse of the concatenated strings.
			1,
			false,        // don't load SA samp in reverse index
			true,         // yes, need ftab in reverse index
			false,        // don't load rstarts in reverse index
			!noRefNames,  // load names?
			startVerbose);
	}
	CHUD_START();
	// Start the metrics thread
	//MetricsThread mett(cerr, metrics, (time_t)metricsIval);
	//if(metricsIval > 0) mett.run();
	{
		Timer _t(cerr, "Multiseed full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			// Thread IDs start at 1
			tids[i] = i+1;
			createThread(&threads[i], multiseedSearchWorker, (void*)&tids[i]);
		}
#endif
		int tmp = 0;
		multiseedSearchWorker((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) joinThread(threads[i]);
#endif
	}
	if(metricsIval > 0 && (metricsOfb != NULL || metricsStderr)) {
		metrics.reportInterval(metricsOfb, metricsStderr);
	}
	//if(metricsIval > 0) { mett.kill(); mett.join(); }
}


static PairedPatternSource*   seededQualSearch_patsrc;
static HitSink*               seededQualSearch_sink;
static Ebwt*                  seededQualSearch_ebwtFw;
static Ebwt*                  seededQualSearch_ebwtBw;
static EList<SString<char> >* seededQualSearch_os;
static int                    seededQualSearch_qualCutoff;
static BitPairReference*      seededQualSearch_refs;

static void* seededQualSearchWorkerFullStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource&     _patsrc    = *seededQualSearch_patsrc;
	HitSink&                 _sink      = *seededQualSearch_sink;
	Ebwt&                    ebwtFw     = *seededQualSearch_ebwtFw;
	Ebwt&                    ebwtBw     = *seededQualSearch_ebwtBw;
	EList<SString<char> >&   os         = *seededQualSearch_os;
	int                      qualCutoff = seededQualSearch_qualCutoff;
	BitPairReference*        refs       = seededQualSearch_refs;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);

	AlignerMetrics *metrics = NULL;
	if(stats) {
		metrics = new AlignerMetrics();
	}
	UnpairedSeedAlignerFactory alSEfact(
			ebwtFw,
			&ebwtBw,
			seedMms,
			seedLen,
			qualCutoff,
			maxBts,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs,
			os,
			seed,
			metrics);
	PairedSeedAlignerFactory alPEfact(
			ebwtFw,
			&ebwtBw,
			seedMms,
			seedLen,
			qualCutoff,
			maxBts,
			_sink,
			*sinkFact,
			dontReconcileMates,
			mhits,       // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs,
			os,
			seed);
	{
		MixedMultiAligner multi(
				1,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}
	if(metrics != NULL) {
		metrics->printSummary();
		delete metrics;
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) { pthread_exit(NULL); }
#endif
	return NULL;
}

/**
 * Search for a good alignments for each read using criteria that
 * correspond somewhat faithfully to Maq's.  Search is aided by a pair
 * of Ebwt indexes, one for the original references, and one for the
 * transpose of the references.  Neither index should be loaded upon
 * entry to this function.
 *
 * Like Maq, we treat the first 24 base pairs of the read (those
 * closest to the 5' end) differently from the remainder of the read.
 * We call the first 24 base pairs the "seed."
 */
static void seededQualCutoffSearchFull(
        int seedLen,                    /// length of seed (not a maq option)
        int qualCutoff,                 /// maximum sum of mismatch qualities
                                        /// like maq map's -e option
                                        /// default: 70
        int seedMms,                    /// max # mismatches allowed in seed
                                        /// (like maq map's -n option)
                                        /// Can only be 1 or 2, default: 1
        PairedPatternSource& _patsrc,   /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt& ebwtFw,                   /// index of original text
        Ebwt& ebwtBw,                   /// index of mirror text
        EList<SString<char> >& os)      /// text strings, if available (empty otherwise)
{
	// Global intialization
	assert_leq(seedMms, 3);

	seededQualSearch_patsrc   = &_patsrc;
	seededQualSearch_sink     = &_sink;
	seededQualSearch_ebwtFw   = &ebwtFw;
	seededQualSearch_ebwtBw   = &ebwtBw;
	seededQualSearch_os       = &os;
	seededQualSearch_qualCutoff = qualCutoff;

	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(gColor || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, gColor, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, gVerbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	seededQualSearch_refs = refs;

#ifdef BOWTIE_PTHREADS
	AutoArray<pthread_t> threads(nthreads-1);
	AutoArray<int> tids(nthreads-1);
#endif

	{
		// Load the other half of the index into memory
		assert(!ebwtFw.isInMemory());
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // it's not bidirectional search, so we don't care how reverse is constructed
			true, // load SA sample?
			true, // load ftab?
			true, // load rstarts?
			!noRefNames,
			startVerbose);
	}
	{
		// Load the other half of the index into memory
		assert(!ebwtBw.isInMemory());
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // it's not bidirectional search, so we don't care how reverse is constructed
			true, // load SA sample?
			true, // load ftab?
			true, // load rstarts?
			!noRefNames,
			startVerbose);
	}
	CHUD_START();
	{
		// Phase 1: Consider cases 1R and 2R
		Timer _t(cerr, "Seeded quality full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			createThread(&threads[i],
			             seededQualSearchWorkerFullStateful,
			             (void*)&tids[i]);
		}
#endif
		int tmp = 0;
		seededQualSearchWorkerFullStateful((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) joinThread(threads[i]);
#endif
	}
	if(refs != NULL) {
		delete refs;
	}
	ebwtBw.evictFromMemory();
}


#define PASS_DUMP_FILES dumpAlBase, dumpUnalBase, dumpMaxBase

static string argstr;

template<typename TStr>
static void driver(
	const char * type,
	const string& ebwtFileBase,
	const string& query,
	const EList<string>& queries,
	const EList<string>& qualities,
	const string& outfile)
{
	if(gVerbose || startVerbose)  {
		cerr << "Entered driver(): "; logTime(cerr, true);
	}
	// Vector of the reference sequences; used for sanity-checking
	EList<SString<char> > names, os;
	EList<size_t> nameLens, seqLens;
	// Read reference sequences from the command-line or from a FASTA file
	if(!origString.empty()) {
		// Read fasta file(s)
		EList<string> origFiles;
		tokenize(origString, ",", origFiles);
		parseFastas(origFiles, names, nameLens, os, seqLens);
	}
	PatternParams pp(
		format,        // file format
		fileParallel,  // true -> wrap files with separate PairedPatternSources
		seed,          // pseudo-random seed
		useSpinlock,   // use spin locks instead of pthreads
		patDumpfile,   // file name of file to dump to
		solexaQuals,   // true -> qualities are on solexa64 scale
		phred64Quals,  // true -> qualities are on phred64 scale
		integerQuals,  // true -> qualities are space-separated numbers
		fuzzy,         // true -> try to parse fuzzy fastq
		randomNum,     // number of randomly-generated reads to make
		randomLength,  // length of randomly-generated reads
		fastaContLen,  // length of sampled reads for FastaContinuous...
		fastaContFreq, // frequency of sampled reads for FastaContinuous...
		skipReads      // skip the first 'skip' patterns
	);
	if(gVerbose || startVerbose) {
		cerr << "Creating PatternSource: "; logTime(cerr, true);
	}
	PairedPatternSource *patsrc = PairedPatternSource::setupPatternSources(
		queries,     // singles, from argv
		mates1,      // mate1's, from -1 arg
		mates2,      // mate2's, from -2 arg
		mates12,     // both mates on each line, from --12 arg
		qualities,   // qualities associated with singles
		qualities1,  // qualities associated with m1
		qualities2,  // qualities associated with m2
		pp,          // read read-in parameters
		gVerbose || startVerbose); // be talkative

	// Open hit output file
	if(gVerbose || startVerbose) {
		cerr << "Opening hit output file: "; logTime(cerr, true);
	}
	OutFileBuf *fout;
	if(!outfile.empty()) {
		fout = new OutFileBuf(outfile.c_str(), false);
	} else {
		if(outType == OUTPUT_CHAIN) {
			cerr << "Error: Must specify an output file when output mode is --chain" << endl;
			throw 1;
		}
		fout = new OutFileBuf();
	}
	ReferenceMap* rmap = NULL;
	if(refMapFile != NULL) {
		if(gVerbose || startVerbose) {
			cerr << "About to load in a reference map file with name "
			     << refMapFile << ": "; logTime(cerr, true);
		}
		rmap = new ReferenceMap(refMapFile, !noRefNames);
	}
	AnnotationMap* amap = NULL;
	if(annotMapFile != NULL) {
		if(gVerbose || startVerbose) {
			cerr << "About to load in an annotation map file with name "
			     << annotMapFile << ": "; logTime(cerr, true);
		}
		amap = new AnnotationMap(annotMapFile);
	}
	// Initialize Ebwt object and read in header
	if(gVerbose || startVerbose) {
		cerr << "About to initialize fw Ebwt: "; logTime(cerr, true);
	}
	adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, gVerbose);
	Ebwt ebwt(
		adjustedEbwtFileBase,
	    gColor,   // index is colorspace
		-1,       // fw index
	    true,     // index is for the forward direction
	    /* overriding: */ offRate,
	    useMm,    // whether to use memory-mapped files
	    useShmem, // whether to use shared memory
	    mmSweep,  // sweep memory-mapped files
	    !noRefNames, // load names?
		true,        // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
	    rmap,     // reference map, or NULL if none is needed
	    gVerbose, // whether to be talkative
	    startVerbose, // talkative during initialization
	    false /*passMemExc*/,
	    sanityCheck);
	Ebwt* ebwtBw = NULL;
	// We need the mirror index if mismatches are allowed
	if((!doMultiseed && (mismatches > 0 || maqLike)) ||
	   (doMultiseed && multiseedMms > 0))
	{
		if(gVerbose || startVerbose) {
			cerr << "About to initialize rev Ebwt: "; logTime(cerr, true);
		}
		ebwtBw = new Ebwt(
			adjustedEbwtFileBase + ".rev",
			gColor,  // index is colorspace
			1,       // TODO: maybe not
		    false, // index is for the reverse direction
		    /* overriding: */ offRate,
		    useMm,    // whether to use memory-mapped files
		    useShmem, // whether to use shared memory
		    mmSweep,  // sweep memory-mapped files
		    !noRefNames, // load names?
			true,        // load SA sample?
			true,        // load ftab?
			true,        // load rstarts?
		    rmap,        // reference map, or NULL if none is needed
		    gVerbose,    // whether to be talkative
		    startVerbose, // talkative during initialization
		    false /*passMemExc*/,
		    sanityCheck);
	}
	if(sanityCheck && !os.empty()) {
		// Sanity check number of patterns and pattern lengths in Ebwt
		// against original strings
		assert_eq(os.size(), ebwt.nPat());
		for(size_t i = 0; i < os.size(); i++) {
			assert_eq(os[i].length(), ebwt.plen()[i] + (gColor ? 1 : 0));
		}
	}
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		ebwt.loadIntoMemory(
			gColor ? 1 : 0,
			-1, // fw index
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
		ebwt.checkOrigs(os, gColor, false);
		ebwt.evictFromMemory();
	}
	{
		Timer _t(cerr, "Time searching: ", timing);
		if(gVerbose || startVerbose) {
			cerr << "Creating ReadSink, HitSink: "; logTime(cerr, true);
		}
		ReadSink readSink(
			dumpAlBase,   // basename of same-format files to dump aligned reads to
			dumpUnalBase, // basename of same-format files to dump unaligned reads to
			dumpMaxBase,  // basename of same-format files to dump reads with more than -m valid alignments to
			format == TAB_MATE); // both mates to same file?
		// Set up hit sink; if sanityCheck && !os.empty() is true,
		// then instruct the sink to "retain" hits in a vector in
		// memory so that we can easily sanity check them later on
		HitSink *sink;
		AlnSink *mssink = NULL;
		auto_ptr<Mapq> bmapq(new BowtieMapq());
		const EList<string>* refnames = &ebwt.refnames();
		if(noRefNames) refnames = NULL;
		switch(outType) {
			case OUTPUT_FULL: {
				sink = new VerboseHitSink(
					fout,         // pointer to output object
					&readSink,    // dump reads according to alignment status
					refnames,     // reference names
					offBase,      // add this to 0-based offsets before printing
					gColorSeq,    // whether to print color sequence
					gColorQual,   // whether to print color qualities
					printCost,    // whether to print penalty
					suppressOuts, // suppress alignment columns
					rmap,         // if != NULL, ReferenceMap to use
					amap,         // if != NULL, AnnotationMap to use
					fullRef,      // true -> print whole reference name, not just up to first whitespace
					partitionSz); // size of partition, so we can check for straddling alignments
				mssink = new AlnSinkVerbose(
					fout,         // initial output stream
					suppressOuts, // suppress alignment columns
					&readSink,    // read sink
					*bmapq.get(), // mapping quality calculator
					false,        // delete output stream objects upon destruction
					refnames,     // reference names
					gQuiet,       // don't print alignment summary at end   
					offBase,      // add this to 0-based offsets before printing
					gColorSeq,    // colorspace: print color seq instead of decoded nucs
					gColorQual,   // colorspace: print color quals instead of decoded quals
					gColorExEnds, // whether to exclude end positions from decoded colorspace alignments
					printPlaceholders, // print maxs/unals
					printFlags,   // print alignment flags
					printCost,    // print penalty in extra column
					printParams,  // print alignment parameters in extra column
					rmap,         // if != NULL, ReferenceMap to use
					fullRef,      // print entire reference name including whitespace
					partitionSz); // size of partition, so we can check for straddling alignments
				break;
			}
			case OUTPUT_SAM: {
				SAMHitSink *sam = new SAMHitSink(
					fout,
					&readSink,           // dump reads according to alignment status
					refnames,            // reference names
					1,
					rmap,
					amap,
					fullRef,
					sampleMax,
					defaultMapq);
				if(!samNoHead) {
					EList<string> refnames;
					if(!samNoSQ) {
						readEbwtRefnames(adjustedEbwtFileBase, refnames);
					}
					sam->appendHeaders(
							sam->out(0), ebwt.nPat(),
							refnames, samNoSQ, rmap,
							ebwt.plen(), fullRef,
							argstr.c_str(),
							rgs.empty() ? NULL : rgs.c_str());
				}
				sink = sam;
				break;
			}
			case OUTPUT_CHAIN: {
				throw 1;
				sink = new ChainingHitSink(
					fout,
					&readSink,           // dump reads according to alignment status
					refnames,            // reference names
					strata,
					amap);
				break;
			}
			case OUTPUT_NONE: {
				sink = new StubHitSink();
				break;
			}
			default:
				cerr << "Invalid output type: " << outType << endl;
				throw 1;
		}
		if(gVerbose || startVerbose) {
			cerr << "Dispatching to search driver: "; logTime(cerr, true);
		}
		if(doMultiseed) {
			// Bowtie 2
			// Set up penalities
			Penalties pens(
				0,             // reward for a match
				penMmcType,    // how to penalize mismatches
				penMmc,        // constant if mm pelanty is a constant
				penSnp,        // penalty for nucleotide mismatch in decoded colorspace als
				nCeilConst,    // constant coeff for penalty ceil w/r/t read length
				nCeilLinear,   // linear coeff for penalty ceil w/r/t read length
				penNType,      // how to penalize Ns in the read
				penN,          // constant if N pelanty is a constant
				penNCatPair,   // whether to concat mates before N filtering
				penRdGapConst,  // constant cost of extending a gap in the read
				penRfGapConst,  // constant cost of extending a gap in the reference
				penRdGapLinear, // coefficient of linear term for cost of gap extension in read
				penRfGapLinear  // coefficient of linear term for cost of gap extension in ref
			);
			// Set up global constraint
			Constraint gc = Constraint::penaltyFuncBased(
				costCeilConst,
				costCeilLinear
			);
			// Set up seeds
			EList<Seed> seeds;
			Seed::mmSeeds(
				multiseedMms,    // max # mms allowed in a multiseed seed
				multiseedLen,    // length of a multiseed seed (scales down if read is shorter)
				seeds,           // seeds
				gc               // global constraint
			);
			// Set up listeners for alignment progress
			EList<SeedHitSink*> seedHitSink;
			EList<SeedCounterSink*> seedCounterSink;
			EList<SeedActionSink*> seedActionSink;
			ofstream *hitsOf = NULL, *cntsOf = NULL, *actionsOf = NULL;
			if(!saHitsFn.empty()) {
				hitsOf = new ofstream(saHitsFn.c_str());
				if(!hitsOf->is_open()) {
					cerr << "Error: Unable to open seed hit dump file " << saHitsFn << endl;
					throw 1;
				}
				seedHitSink.push_back(new StreamTabSeedHitSink(*hitsOf));
			}
			if(!saCountersFn.empty()) {
				cntsOf = new ofstream(saCountersFn.c_str());
				if(!cntsOf->is_open()) {
					cerr << "Error: Unable to open seed counter dump file " << saCountersFn << endl;
					throw 1;
				}
				seedCounterSink.push_back(new StreamTabSeedCounterSink(*cntsOf));
			}
			if(!saActionsFn.empty()) {
				actionsOf = new ofstream(saActionsFn.c_str());
				if(!actionsOf->is_open()) {
					cerr << "Error: Unable to open seed action dump file " << saActionsFn << endl;
					throw 1;
				}
				seedActionSink.push_back(new StreamTabSeedActionSink(*actionsOf));
			}
			OutFileBuf *metricsOfb = NULL;
			if(!metricsFile.empty() && metricsIval > 0) {
				metricsOfb = new OutFileBuf(metricsFile);
			}
			// Do the search for all input reads
			assert(patsrc != NULL);
			assert(mssink != NULL);
			multiseedSearch(
				pens,    // edit penalties
				seeds,   // seeds
				*patsrc, // pattern source
				*mssink, // hit sink
				ebwt,    // BWT
				*ebwtBw, // BWT'
				seedHitSink,
				seedCounterSink,
				seedActionSink,
				metricsOfb
			);
			for(size_t i = 0; i < seedHitSink.size();     i++) delete seedHitSink[i];
			for(size_t i = 0; i < seedCounterSink.size(); i++) delete seedCounterSink[i];
			for(size_t i = 0; i < seedActionSink.size();  i++) delete seedActionSink[i];
			if(hitsOf != NULL) delete hitsOf;
			if(cntsOf != NULL) delete cntsOf;
			if(actionsOf != NULL) delete actionsOf;
			if(metricsOfb != NULL) delete metricsOfb;
		} else {
			// Bowtie 1
			if(maqLike) {
				seededQualCutoffSearchFull(seedLen,
										   qualThresh,
										   seedMms,
										   *patsrc,
										   *sink,
										   ebwt,    // forward index
										   *ebwtBw, // mirror index (not optional)
										   os);     // references, if available
			}
			else if(mismatches > 0) {
				if(mismatches == 1) {
					assert(ebwtBw != NULL);
					mismatchSearchFull(*patsrc, *sink, ebwt, *ebwtBw, os);
				} else if(mismatches == 2 || mismatches == 3) {
					twoOrThreeMismatchSearchFull(*patsrc, *sink, ebwt, *ebwtBw, os, mismatches == 2);
				} else {
					cerr << "Error: " << mismatches << " is not a supported number of mismatches" << endl;
					throw 1;
				}
			} else {
				// Search without mismatches
				// Note that --fast doesn't make a difference here because
				// we're only loading half of the index anyway
				exactSearch(*patsrc, *sink, ebwt, os);
			}
		}
		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}
		if(ebwtBw != NULL) {
			delete ebwtBw;
		}
		if(!gQuiet) {
			if(doMultiseed) {
				mssink->finish(hadoopOut, sampleMax); // end the hits section of the hit file
			} else {
				sink->finish(hadoopOut, sampleMax); // end the hits section of the hit file
			}
		}
		delete patsrc;
		delete sink;
		delete mssink;
		delete amap;
		delete rmap;
		if(fout != NULL) delete fout;
	}
}

// C++ name mangling is disabled for the bowtie() function to make it
// easier to use Bowtie as a library.
extern "C" {

/**
 * Main bowtie entry function.  Parses argc/argv style command-line
 * options, sets global configuration variables, and calls the driver()
 * function.
 */
int bowtie(int argc, const char **argv) {
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();
		for(int i = 0; i < argc; i++) {
			argstr += argv[i];
			if(i < argc-1) argstr += " ";
		}
		string ebwtFile;  // read serialized Ebwt from this file
		string query;   // read query string(s) from this file
		EList<string> queries(0);
		string outfile; // write query results to this file
		if(startVerbose) { cerr << "Entered main(): "; logTime(cerr, true); }
		parseOptions(argc, argv);
		argv0 = argv[0];
		if(showVersion) {
			cout << argv0 << " version " << BOWTIE_VERSION << endl;
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Built on " << BUILD_HOST << endl;
			cout << BUILD_TIME << endl;
			cout << "Compiler: " << COMPILER_VERSION << endl;
			cout << "Options: " << COMPILER_OPTIONS << endl;
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
				 << sizeof(int)
				 << ", " << sizeof(long) << ", " << sizeof(long long)
				 << ", " << sizeof(void *) << ", " << sizeof(size_t)
				 << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}
	#ifdef CHUD_PROFILING
		chudInitialize();
		chudAcquireRemoteAccess();
	#endif
		{
			Timer _t(cerr, "Overall time: ", timing);
			if(startVerbose) {
				cerr << "Parsing index and read arguments: "; logTime(cerr, true);
			}

			// Get index basename
			if(optind >= argc) {
				cerr << "No index, query, or output file specified!" << endl;
				printUsage(cerr);
				return 1;
			}
			ebwtFile = argv[optind++];

			// Get query filename
			if(optind >= argc) {
				if(mates1.size() > 0 || mates12.size() > 0) {
					query = "";
				} else {
					cerr << "No query or output file specified!" << endl;
					printUsage(cerr);
					return 1;
				}
			} else if (mates1.size() == 0 && mates12.size() == 0) {
				query = argv[optind++];
				// Tokenize the list of query files
				tokenize(query, ",", queries);
				if(queries.size() < 1) {
					cerr << "Tokenized query file list was empty!" << endl;
					printUsage(cerr);
					return 1;
				}
			}

			// Get output filename
			if(optind < argc) {
				outfile = argv[optind++];
			}

			// Extra parametesr?
			if(optind < argc) {
				cerr << "Extra parameter(s) specified: ";
				for(int i = optind; i < argc; i++) {
					cerr << "\"" << argv[i] << "\"";
					if(i < argc-1) cerr << ", ";
				}
				cerr << endl;
				if(mates1.size() > 0) {
					cerr << "Note that if <mates> files are specified using -1/-2, a <singles> file cannot" << endl
						 << "also be specified.  Please run bowtie separately for mates and singles." << endl;
				}
				throw 1;
			}

			// Optionally summarize
			if(gVerbose) {
				cout << "Input ebwt file: \"" << ebwtFile << "\"" << endl;
				cout << "Query inputs (DNA, " << file_format_names[format] << "):" << endl;
				for(size_t i = 0; i < queries.size(); i++) {
					cout << "  " << queries[i] << endl;
				}
				cout << "Quality inputs:" << endl;
				for(size_t i = 0; i < qualities.size(); i++) {
					cout << "  " << qualities[i] << endl;
				}
				cout << "Output file: \"" << outfile << "\"" << endl;
				cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
				cout << "Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
			#ifdef NDEBUG
				cout << "Assertions: disabled" << endl;
			#else
				cout << "Assertions: enabled" << endl;
			#endif
			}
			if(ipause) {
				cout << "Press key to continue..." << endl;
				getchar();
			}
			driver<SString<char> >("DNA", ebwtFile, query, queries, qualities, outfile);
			CHUD_STOP();
		}
#ifdef CHUD_PROFILING
		chudReleaseRemoteAccess();
#endif
		return 0;
	} catch(exception& e) {
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
} // bowtie()
} // extern "C"
