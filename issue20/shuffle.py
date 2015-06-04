#!/usr/bin/env python

# Shuffle the hard way, black box style
import os
import sys
from random import randrange

        
if __name__ == "__main__":
    pair_1_fname = sys.argv[1]
    pair_2_fname = sys.argv[2]

    fs1 = open(pair_1_fname, "r")
    fs2 = open(pair_2_fname, "r")
    fo1 = open(pair_1_fname+".out", "w")
    fo2 = open(pair_2_fname+".out", "w")
    # Gulp everything into memory
    all_lines_1 = fs1.readlines()
    all_lines_2 = fs2.readlines()

    times_shuffle = 100000
    for x in xrange(times_shuffle):
        records = len(all_lines_1)/4
        src_rec = randrange(0, records)*4

        for i in range(4):
            fo1.write(all_lines_1[src_rec + i])
            fo2.write(all_lines_2[src_rec + i])
        
        del(all_lines_1[src_rec:src_rec + 4])
        del(all_lines_2[src_rec:src_rec + 4])
            
        if len(all_lines_1) == 0:
            break
        
    fo1.writelines(all_lines_1)
    fo2.writelines(all_lines_2)
    fs1.close()
    fs2.close()
    fo1.close()
    fo2.close()

