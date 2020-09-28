#!/bin/bash

URL="https://genome-idx.s3.amazonaws.com/bt"
INDEXES="GRCh38_noalt_as GRCm38 BDGP6 Ash1v1.7 Dmel_A4_1.0"

for idx in ${INDEXES} ; do
	if [[ ! -f ${idx}/${idx}.1.bt2 ]] ; then
		wget ${URL}/${idx}.zip && unzip -o ${idx}.zip && rm -f ${idx}.zip
	fi
done
