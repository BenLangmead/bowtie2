#!/bin/sh

for al in bt2 bt2loc bwam ; do
    for level in high low ; do
	awk '{print "@"$1; print $10; print "+"; print $11}' \
	    misscored_${al}_${level}.sam > \
	    misscored_${al}_${level}.fastq
    done
done

