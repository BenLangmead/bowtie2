Bowtie 2 tests
==============

Many of the tests use the debug versions of the binaries, so first do `make allall` in the root.

#### Simple test cases

A collection of small tests with hand-crafted input and output checks are collected in `scripts/test/simple_tests.pl`, with `scripts/test/simple_tests.sh` as a driver script.  Uses Perl module `DNA.pm`.

From root:

    sh scripts/test/simple_tests.sh

#### Regression tests

They use the Python infrastructure in the `scripts/test` subdirectory, including `bt2face.py`, `dataface.py` and `btdata.py`.

From root:

    python scripts/test/regressions.py --verbose

Val Antonescu originally set these up.

####
