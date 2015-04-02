#!/usr/bin/env python

import glob
import sys
    
    
def getEachRowFor(file_name, row=1):
    threads = file_name.split("_", 1)[0]
    with open(file_name, "r") as fin:
        while True:
            for i in range(1, row):
                buff = fin.readline()
                if not buff:
                    return
            yield (threads, fin.readline())


def get_seconds(val):
    time_in_seconds = 0.0
    time_raw_vals = val.split(":")
    coef = 1
    for v in reversed(time_raw_vals):
        time_in_seconds += coef * float(v)
        coef *= 60
    return time_in_seconds


def translate_time(t_raw_tuple):
    s = t_raw_tuple[1].split(",")
    rez = [t_raw_tuple[0], []]

    for i in range(0, 3):
        t_sec = get_seconds(s[i])
        rez[1].append(t_sec)
    return rez


def compute_stats(thread_list):
    n = 0
    sums = [0, 0, 0]
    for thread in thread_list:
        n += 1
        sums = [thread[1][x]+sums[x] for x in range(0, 3)]
    return [thread_list[0][0], [sums[x]/n for x in range(0, 3)]]
    
    
file_mask = "*" + sys.argv[1]  # _newHW_Imr_def_sample.dat
rezult = list()

for fname in glob.glob(file_mask):
    file_handler = getEachRowFor(fname, row=8)  # 17 for pairs, 8 for non-paired
    per_thread_rez = list()
    for pair in file_handler:
        per_thread_rez.append(translate_time(pair))
    rezult.append(compute_stats(per_thread_rez))

for sample in rezult:
    print "%s,%.2f,%.2f" % (sample[0], sample[1][0] + sample[1][1], sample[1][2])

