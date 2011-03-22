#!/bin/sh

#  run.sh

make bowtie2 bowtie2-debug bowtie2-build bowtie2-build-debug
perl scripts/sim/run.pl \
	--bowtie2=./bowtie2 \
	--bowtie2-debug=./bowtie2-debug \
	--bowtie2-build=./bowtie2-build \
	--bowtie2-build-debug=./bowtie2-build-debug
