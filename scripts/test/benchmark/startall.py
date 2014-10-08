#!/usr/bin/env python

import os
import inspect
import logging
import benchmarks as bm
from optparse import OptionParser
    
            
def parse_args():
    usage = " %prog [options] \n\n"
    usage += "This is a long and resource consuming performance test.\n"
    usage += "It is meant to be run inside a git repo. It will run the bowtie tool\n"
    usage += "suite for various cases and record the performance for the current HEAD\n"
    usage += "(default) or any other desired commit.\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--commit", 
                    action="store",type="string",dest="commit", default=None,
                    help="For what commit to record performance metrics. Default HEAD.")
    parser.add_option("-i", "--benchmark-id", 
                    action="store",type="string",dest="benchmark_id", default=None,
                    help="What name/id to use for this recording. Default to branch name.")
    parser.add_option("-v", "--verbose", 
                    action="store_true",dest="verbose", default=False,
                    help="Print more info about each step.")

    (options, args) = parser.parse_args()
    return options
    


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.ERROR)    
    options = parse_args()

    if options.verbose:
        logging.getLogger().setLevel(level=logging.INFO)
   
    src_file_path  = os.path.realpath(inspect.getsourcefile(parse_args))
    curr_path      = os.path.dirname(src_file_path)
    bw2_subdir     = 'bowtie2'
    
    benchmarks_subdir = os.path.join(curr_path,"data","conf") 
    download_subdir   = os.path.join(curr_path,"data","downloads")
    out_dir           = os.path.join(curr_path,"data","out")
    
    i = curr_path.find(bw2_subdir)
    bt2_path = curr_path[:i+len(bw2_subdir)] 

    batch_benchmarks = bm.Benchmarks(benchmarks_dir = benchmarks_subdir, 
                                     data_dir = download_subdir,
                                     output_dir = out_dir,
                                     bin_dir = bt2_path,
                                     commit = options.commit,
                                     benchmark_id = options.benchmark_id)

    for set in batch_benchmarks:
        if not set.input_data_loaded:
            set.load()
        set.run();





