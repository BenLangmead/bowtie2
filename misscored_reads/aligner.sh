#!/usr/bin/env bash

set -ex

for al in bt2 bt2loc ; do
  for level in high low ; do
    READS_SUMM="${al}_${level}.reads.csv"
    CMDS_SUMM="${al}_${level}.cmds.csv"
    if [[ ! -f ${READS_SUMM} || ! -f ${CMDS_SUMM} ]] ; then
      python vassess_aligner.py \
        cmds_${al}.txt ${al}_${level}.fastq \
        ${READS_SUMM} ${CMDS_SUMM}
    else
      echo "*** Found ${READS_SUMM} ${CMDS_SUMM} so skipping... ***"
    fi
  done
done
