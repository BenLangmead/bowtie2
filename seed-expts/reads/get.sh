#!/bin/bash

URL="https://bt2-bench.s3.amazonaws.com/reads/vargas_annotated"
READS="fly_A4_local.sam
fly_A4_semiglobal.sam
fly_BDGP6_local.sam
fly_BDGP6_semiglobal.sam
human_Ash1_local.sam
human_Ash1_semiglobal.sam
human_GRCh38_local.sam
human_GRCh38_semiglobal.sam
mouse_GRCm38_local.sam
mouse_GRCm38_semiglobal.sam"

for rds in ${READS} ; do
	if [[ ! -f ${rds} ]] ; then
		wget "${URL}/${rds}"
	fi
done
