#!/usr/bin/env python3
"""
   Runs benchmark sets specified in JSON files.

Example:   
   # Runs a test set specified in simple.json and
   # record this run under id "test".
   
   run.py -t simple.json -i test 
   
   # runs same test using another bowtie binaries from 
   # directory bowtie-devel.
   
   run.py -t simple.json -i test -b ~/bowtie-devel
   
   # Runs all *.json benchmark test sets from ~/bt2_benchmarks
   # recording them under the output directory "records/all_tests". 
   
   run.py -s ~/bt2_benchmarks -i all_tests -o records

"""

import os
import logging
import benchmarks as bm
from optparse import OptionParser


def parse_args():
    usage = " %prog [options] \n\n"
    usage += "Runs specified bowtie2 benchmarking tests.\n"
    usage += "Benchmarks are defined using JSON test sets files. For each test\n"
    usage += "it is defined what input data prerequisites are required, commands\n"
    usage += "used to generate them in case they are missing, what command to test\n"
    usage += "and what metrics to record."
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--benchmark-test",
                      action="store", type="string", dest="benchmark_test", default=None,
                      help="One benchmark tests to run. Tests are files in JSON "
                           "format describing a set of test, input data requirements and what "
                           "metrics to record.")
    parser.add_option("-s", "--benchmarks-dir",
                      action="store", type="string", dest="benchmarks_dir", default=None,
                      help="A path to a directory with benchmark tests to run. All files "
                           "with .json extension will be loaded and run."
                           "metrics to record.")
    parser.add_option("-i", "--benchmark-id",
                      action="store", type="string", dest="benchmark_id", default=None,
                      help="(mandatory).What name/id to use for this recording.")
    parser.add_option("-d", "--download-dir",
                      action="store", type="string", dest="download_dir", default='download',
                      help=" (Default: ./download).The directory path where all input "
                           "data used by the benchmarks should be stored.")
    parser.add_option("-b", "--bowtie-dir",
                      action="store", type="string", dest="bowtie_dir", default='bowtie2',
                      help="(Default: ./bowtie2).Bowtie directory.")
    parser.add_option("-o", "--output-dir",
                      action="store", type="string", dest="out_dir", default='benchmark_rezults',
                      help="(Default: ./out). Directory where tests results are going to "
                           "be written.")
    parser.add_option("-l", "--log-file",
                      action="store", type="string", dest="log_fname", default=None,
                      help="(Default: stderr). Log file name if desired.")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="Print more info about each step.")
    parser.add_option("-u", "--debug",
                      action="store_true", dest="debug", default=False,
                      help="Print extra debug info. Used only for diagnosing purpose.")

    (options, args) = parser.parse_args()

    if options.benchmark_id is None:
        print("Mandatory option is missing (--benchmark-id)")
        parser.print_help()
        exit(-1)

    return options


if __name__ == "__main__":
    options = parse_args()

    if options.log_fname is not None:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            filename=options.log_fname,
                            level=logging.ERROR)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.ERROR)

    if options.verbose:
        logging.getLogger().setLevel(level=logging.INFO)

    if options.debug:
        logging.getLogger().setLevel(level=logging.DEBUG)

    curr_path = os.getcwd()

    if not os.path.exists(options.out_dir):
        logging.debug("Creating output dir: %s" % options.out_dir)
        os.mkdir(options.out_dir)

    if not os.path.exists(options.download_dir):
        logging.debug("Creating download dir %s" % options.download_dir)
        os.mkdir(options.download_dir)

    if options.benchmarks_dir is not None:
        if not os.path.exists(options.benchmarks_dir):
            logging.error("Cannot find benchmark directory %s" % options.benchmarks_dir)

    batch_benchmarks = bm.Benchmarks(benchmarks_dir=options.benchmarks_dir,
                                     benchmark_test=options.benchmark_test,
                                     data_dir=options.download_dir,
                                     output_dir=options.out_dir,
                                     bin_dir=options.bowtie_dir,
                                     benchmark_id=options.benchmark_id)

    for test_set in batch_benchmarks:
        if not test_set.input_data_loaded:
            test_set.load()
        test_set.run()






