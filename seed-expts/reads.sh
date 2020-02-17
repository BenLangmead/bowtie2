#!/usr/bin/env bash

set -ex

for al in bt2 bt2loc ; do
  for level in high low ; do
    python vassess_reads.py vargas_misscored_${al}_${level}.sam > ${al}_${level}.fastq
  done
done
