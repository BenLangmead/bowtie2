#!/bin/sh

#  simple_tests.sh

make $* bowtie2 bowtie2-debug bowtie2-build bowtie2-build-debug && \
perl scripts/test/simple_tests.pl \
	--bowtie2=./bowtie2-debug \
	--bowtie2-build=./bowtie2-build-debug
