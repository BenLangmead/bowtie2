#!/bin/sh

# Spawn a bunch of parallel processes running random_bowtie_tests.sh
# and dumping results into text files.  Run from the bowtie dir.

NUM=$1
if [ -z "$NUM" ] ; then
	NUM=4
fi

make bowtie bowtie-debug bowtie-build-debug bowtie-inspect-debug
if [ "$?" != "0" ] ; then
	echo "Error during build"
	exit 1
fi

echo > .randpids
while [ $NUM -gt 0 ] ; do
	echo "Spawning: sh scripts/random_bowtie_tests.sh $NUM > .randstdout.$NUM 2> .randstderr.$NUM &"
	sh scripts/random_bowtie_tests.sh $NUM > .randstdout.$NUM 2> .randstderr.$NUM &
	echo $! >> .randpids
	NUM=`expr $NUM - 1`
done

for p in `cat .randpids` ; do
	echo "Waiting for pid $p"
	wait $p
done
