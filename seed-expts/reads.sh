#!/usr/bin/env bash

set -ex

for mode in local semiglobal ; do
  for genome in fly_A4 fly_BDGP6 human_Ash1 human_GRCh38 mouse_GRCm38 ; do
    python vassess_reads.py reads/${genome}_${mode}.sam > reads/${genome}_${mode}.fastq
  done
done
