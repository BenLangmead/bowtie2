#!/usr/bin/env python3
"""
A few items to deal with sets of benchmarks.
"""

import os
import re
import csv
import glob
import json
import logging
import subprocess
import samreader as Sr


class Benchmarks(object):
    """ 
    Iterable for all benchmarks found in our test directory. 
    """

    def __init__(self,
                 benchmarks_dir=None,
                 benchmark_test=None,
                 data_dir=None,
                 output_dir=None,
                 bin_dir=None,
                 benchmark_id=None):
        self.set_idx = 0
        self.values = list()
        self.benchmark_id = benchmark_id
        self.data_dir = data_dir
        self.bin_dir = bin_dir
        self.benchmarks_dir = benchmarks_dir
        self.output_dir = os.path.join(output_dir, self.benchmark_id)
        # 
        if os.path.exists(self.output_dir):
            logging.error("A benchmark with the same name already exists (%s)" % self.output_dir)
            raise OSError("Directory already exists: %s" % self.output_dir)
        #
        logging.debug("Creating test output directory: %s" % self.output_dir)
        os.mkdir(self.output_dir)

        if benchmark_test is not None:
            self.values.append(self._load_benchmark(benchmark_test))

        if self.benchmarks_dir is not None and os.path.isdir(self.benchmarks_dir):
            logging.debug("Parse all json files from: %s" % self.benchmarks_dir)
            all_tests = glob.glob(os.path.join(self.benchmarks_dir, "*.json"))
            if all_tests:
                for b_set in all_tests:
                    logging.debug("Adding test: %s" % b_set)
                    self.values.append(self._load_benchmark(b_set))
            else:
                logging.warn("No json files found!")

    def __iter__(self):
        self.set_idx = 0
        return self

    def __next__(self):
        if self.set_idx == len(self.values):
            raise StopIteration

        value = self.values[self.set_idx]
        self.set_idx += 1
        return value

    def _load_benchmark(self, bench_fname):
        logging.debug("Loading test: %s" % bench_fname)
        with open(bench_fname) as fp:
            set_data = json.load(fp)
            test_list = set_data["tests"]

            for test_desc in test_list:
                data_load_list = test_desc["input_data"]["loading"]
                for i, load_cmd in enumerate(data_load_list):
                    load_cmd = re.sub(r'##BT2DIR##', self.bin_dir, load_cmd)
                    load_cmd = re.sub(r'##DATADIR##', self.data_dir, load_cmd)
                    data_load_list[i] = load_cmd

                runable_opt = test_desc["runable"]
                for key in ["options", "parameters", "outfiles"]:
                    try:
                        runable_prop = runable_opt[key]
                    except KeyError:
                        continue

                    for i, item in enumerate(runable_prop):
                        item = re.sub(r'##BT2DIR##', self.bin_dir, item)
                        item = re.sub(r'##DATADIR##', self.data_dir, item)
                        item = re.sub(r'##OUTDIR##', self.output_dir, item)
                        runable_prop[i] = item

            return BenchmarkSet(set_data, self.data_dir, self.output_dir, self.bin_dir)


class BenchmarkSet(object):
    """ A Benchmark item
    """

    def __init__(self, data, data_dir, out_dir, bin_dir):
        self.data = data
        self.data_dir = data_dir
        self.out_dir = out_dir
        self.bin_dir = bin_dir
        self.input_data_loaded = False

    def run(self):
        logging.info("Running benchmark set: %s" % self.data["description"])
        all_tests = self.data["tests"]
        for test in all_tests:
            bench_cmd = globals()[test["metric"]](self, test)
            bench_cmd.launch()

    def load(self):
        if self.input_data_loaded:
            return
        logging.info("Loading data for %s" % self.data["name"])
        # load data
        test_list = self.data["tests"]

        for test in test_list:
            data_is_here = True
            data_set = test["input_data"]["files"]
            for data_file in data_set:
                if not os.path.isfile(os.path.join(self.data_dir, data_file)):
                    data_is_here = False

            if not data_is_here:
                logging.info("Generate data for %s" % test["name"])
                for cmd in test["input_data"]["loading"]:
                    logging.info("running: %s" % cmd)
                    subprocess.check_call(cmd, shell=True)

        self.input_data_loaded = True


class Runable(object):
    """ cmd line and run helper """

    def __init__(self, main_set, test):
        """ start """
        self.benchmark_set = main_set
        self.test = test
        self.prologue = ''
        self.err_log = os.path.join(self.benchmark_set.out_dir, test["name"]) + ".metric"

    def launch(self):
        """ builds cmd and launch"""
        logging.debug("Building command.")
        cmd = self._build_cmd()
        logging.debug("Running command: %s" % cmd)
        self._run_cmd(cmd)
        logging.debug("Calling report formating.")
        self._format_report()

    def _build_cmd(self):
        """ build command line"""
        space = " "
        test = self.test
        prg = os.path.join(self.benchmark_set.bin_dir, test["runable"]["program"])
        cmd = self.prologue + space + prg

        for opt in test["runable"]["options"]:
            cmd = cmd + space + opt

        for parm in test["runable"]["parameters"]:
            cmd = cmd + space + parm

        return cmd

    def _run_cmd(self, cmd):
        """ running """
        logging.info("Start Benchmark %s" % self.test["name"])
        logging.info("Running: %s" % cmd)
        with open(self.err_log, 'w') as errlog:
            subprocess.check_call(cmd, shell=True, stderr=errlog)

    def _format_report(self, cmd):
        """ not always required """
        pass


class TestTime(Runable):
    """ Time benchmarks """

    def __init__(self, main_set, test):
        super(TestTime, self).__init__(main_set, test)
        self.prologue = "/usr/bin/time -f %U,%S,%E "

    def _format_report(self):
        """ formats data in csv format"""
        csv_file = os.path.join(self.benchmark_set.out_dir, self.test["name"]) + ".csv"
        p1 = subprocess.Popen(["tail -1 %s" % self.err_log], shell=True, stdout=subprocess.PIPE)
        line = p1.communicate()[0].rstrip()

        with open(csv_file, 'w') as csvf:
            writer = csv.writer(csvf)
            writer.writerow(["Test", "User Time", "System Time", "Wall Time"])
            row = [self.test["name"]]
            row.extend(line.split(","))
            writer.writerow(row)


class TestAccuracy(Runable):
    """ accuracy """

    def __init__(self, main_set, test):
        super(TestAccuracy, self).__init__(main_set, test)

    def _format_report(self):
        """ collect data in csv format"""
        in_sam_file = self._get_first_sam_input_file()
        in_sam_file = os.path.join(self.benchmark_set.data_dir, in_sam_file)
        initial_data = dict()
        mapq_summary = dict()
        with open(in_sam_file, "r") as fh_sin:
            in_reader = Sr.SamReader(fh_sin)
            for rec in in_reader:
                fl = rec.flag & 255
                if fl % 2 == 0:
                    logging.error("Initial data does not look like a paired search.")
                if fl > 128:
                    q_name = rec.qname + "_2"
                elif fl > 64:
                    q_name = rec.qname + "_1"
                else:
                    logging.error("Again, initial data does not look like a paired search.")
                    q_name = rec.qname
                initial_data[q_name] = (rec.rname, rec.pos)

        out_sam_file = self.test["runable"]["outfiles"][0]
        with open(out_sam_file, "r") as fh_out:
            in_reader = Sr.SamReader(fh_out)
            for rec in in_reader:
                fl = rec.flag & 255
                if fl % 2 == 0:
                    logging.error("This does not look like a paired search.")
                if fl > 128:
                    q_name = rec.qname + "_2"
                elif fl > 64:
                    q_name = rec.qname + "_1"
                else:
                    logging.error("Again, initial data does not look like a paired search.")
                    q_name = rec.qname
                orig = initial_data[q_name]
                if orig[0] == rec.rname and orig[1] - rec.pos < 3:
                    delta = [1, 0]
                else:
                    delta = [0, 1]
                    logging.debug("%s: missed (pos:%d vs %d)" % (q_name, orig[1], rec.pos))
                try:
                    mapq_summary[rec.mapq] = list(map(sum, list(zip(delta, mapq_summary[rec.mapq]))))
                except KeyError:
                    mapq_summary[rec.mapq] = delta

        csv_file = os.path.join(self.benchmark_set.out_dir, self.test["name"]) + ".csv"
        with open(csv_file, 'w') as csvf:
            writer = csv.writer(csvf)
            writer.writerow(["Test name", "MAPQ", "No. Correct", "No. Misses"])
            for k in mapq_summary:
                row = [self.test["name"], k]
                row.extend(mapq_summary[k])
                writer.writerow(row)

    def _get_first_sam_input_file(self):
        all_input_files = self.test["input_data"]["files"]
        for fname in all_input_files:
            if fname[-4:] == ".sam":
                logging.debug("Compare with origin SAM file: %s" % fname)
                return fname

        raise LookupError("No SAM data input file defined for this test!")


