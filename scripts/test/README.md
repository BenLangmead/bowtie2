Bowtie 2 tests
==============

Many of the tests use the debug versions of the binaries, so first do `make allall` in the root.

#### Simple test cases

A collection of small tests with hand-crafted input and output checks are collected in `scripts/test/simple_tests.pl`, with `scripts/test/simple_tests.sh` as a driver script.  Uses Perl module `DNA.pm`.

From root:

    sh scripts/test/simple_tests.sh

#### Regression tests

They use the Python3 infrastructure in the `scripts/test` subdirectory, including `bt2face.py`, `dataface.py` and `btdata.py`.

From root:

    python3 scripts/test/regressions.py --verbose

Val Antonescu originally set these up.

#### Big index test

Builds an index consisting of both human and mouse genomes, pushing the genome size above the 2^32 limit, and necessitating a "big" 64-bit index.  This takes a lot of time and RAM.

From root:

    python3 scripts/test/large_idx.py --verbose

#### Randomized tests

`scripts/sim` contains Perl infrastructure for running randomized tests.  These scripts generate random reference genomes and reads, and Bowtie 2 programs are run in various combinations with cross-checks to ensure, for example, that release-binary output matches debug-binary output and that single-threaded output matches multithreaded output.  The non-core `Math::Random` Perl module is needed.  You can use CPAN to install it.

As of 10/1/2016, the unit tests in `AlignmentCheck.pm` don't all succeed.  Also, random tests fail for reasons that seem related to the test infrastructure, e.g., not specifying the command-line arguments properly.

To run the unit tests for the infrastructure, from `scripts/sim`:

    sh unit.sh

To run the actual suite of random tests in parallel on 8 cores, from the root:

    sh scripts/sim/run.sh 8
