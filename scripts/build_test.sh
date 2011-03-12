#!/bin/sh

#
# Run this from the bowtie directory to sanity-check the bowtie index
# builder
#

make bowtie-build-debug
if ./bowtie-build-debug -s -v genomes/NC_008253.fna .build_test ; then
	echo Build test PASSED
else
	echo Build test FAILED
fi
rm -f .tmp*.ebwt
