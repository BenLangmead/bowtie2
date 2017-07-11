#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import itertools

fq1 = 'example/reads/reads_1.fq'
fq2 = 'example/reads/reads_2.fq'

def parse(readlen):
    reads1, reads2 = [], []
    for fq, reads in [(fq1, reads1), (fq2, reads2)]:
        with open(fq) as fh:
            while True:
                nm = fh.readline()
                if len(nm) == 0:
                    break
                nm = bytes(nm.rstrip())
                seq = bytes(fh.readline().rstrip()[:readlen])
                nm2 = bytes(fh.readline().rstrip())
                qual = bytes(fh.readline().rstrip()[:readlen])
                reads.append([nm, seq, nm2, qual, len(nm) + len(seq) + len(nm2) + len(qual) + 4])
    return reads1, reads2

for block_sz, read_len in itertools.product([0, 400, 1600, 6000, 12000], [30, 50, 100]):
    reads1, reads2 = parse(read_len)
    reads_per_block = int(block_sz / (read_len * 4))
    lab = '%d_%d' % (block_sz, read_len)
    ofn1, ofn2, osamfn = 'reads_%s_1.fq' % lab, 'reads_%s_2.fq' % lab, 'reads_%s.sam' % lab
    with open(ofn1, 'wb') as ofh1:
        with open(ofn2, 'wb') as ofh2:
            if block_sz == 0:
                for rd in reads1:
                    ofh1.write(b'\n'.join(rd[:4]) + b'\n')
                for rd in reads2:
                    ofh2.write(b'\n'.join(rd[:4]) + b'\n')
            else:
                block_i = 0
                while block_i < len(reads1):
                    block_reads1 = reads1[block_i:block_i+reads_per_block]
                    block_reads2 = reads2[block_i:block_i+reads_per_block]
                    if block_i+reads_per_block <= len(reads1):  # not the last
                        block_len1 = sum(map(lambda x: x[-1], block_reads1))
                        block_len2 = sum(map(lambda x: x[-1], block_reads2))
                        assert block_len1 < block_sz, (block_len1, block_sz, read_len, reads_per_block)
                        assert block_len2 < block_sz, (block_len1, block_sz, read_len, reads_per_block)
                        block_reads1[-1][0] += b' ' * (block_sz - block_len1)
                        block_reads2[-1][0] += b' ' * (block_sz - block_len2)
                    for rd in block_reads1:
                        ofh1.write(b'\n'.join(rd[:4]) + b'\n')
                    for rd in block_reads2:
                        ofh2.write(b'\n'.join(rd[:4]) + b'\n')
                    block_i += reads_per_block
    cmd = './bowtie2 --no-hd --reads-per-block %d --block-bytes %d -x example/index/lambda_virus -1 %s -2 %s -S %s' % (reads_per_block, block_sz, ofn1, ofn2, osamfn)
    print(cmd, file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError('Non-zero return value (%d) from bowtie2' % ret)
