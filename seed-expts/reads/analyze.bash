#!/bin/bash

for i in *.sam ; do
	sed 's/.*AS:i://' < ${i} \
		| awk '{t[$1] += 1; tot += $1; n += 1} END {print tot/n ; for(k in t) {print k,t[k]}}' \
		> $i.stats
done
