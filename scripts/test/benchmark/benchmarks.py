#!/usr/bin/env python
"""
A few items to deal with sets of benchmarks.
"""

import os
import re
import glob
import json
import logging
import subprocess


class Benchmarks(object):
    """ 
    Iterable for all benchmarks found in our test directory. 
    """    
    
    def __init__(self, benchmarks_dir = None, 
               data_dir = None,
               output_dir = None,
               bin_dir = None,
               commit = None,
               benchmark_id = None):
        self.set_idx = 0
        self.values  = list()
        self.commit = commit
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
        os.mkdir(self.output_dir)
        for b_set in glob.glob(os.path.join(self.benchmarks_dir,"*.bench")):
            self.values.append(self._load_benchmark(b_set))
        
        
        
    def __iter__(self):
        return self
    
    
    def next(self):
        if self.set_idx == len(self.values):
            raise StopIteration
        
        value = self.values[self.set_idx]
        self.set_idx += 1
        return value
    
    
    def _load_benchmark(self, bench_fname):
        with open(bench_fname) as fp:
            set_data = json.load(fp)
            test_list = set_data["tests"]

            for test_desc in test_list:
                data_load_list = test_desc["input_data"]["loading"]
                for i,load_cmd in enumerate(data_load_list):
                    load_cmd = re.sub(r'##BT2DIR##',self.bin_dir,load_cmd)
                    load_cmd = re.sub(r'##DATADIR##',self.data_dir,load_cmd)
                    data_load_list[i] = load_cmd

                runable_opt = test_desc["runable"]["options"]
                for i,run_opt in enumerate(runable_opt):
                    run_opt = re.sub(r'##BT2DIR##',self.bin_dir,run_opt)
                    run_opt = re.sub(r'##DATADIR##',self.data_dir,run_opt)
                    runable_opt[i] = run_opt
                
                runable_param  = test_desc["runable"]["parameters"]
                for i,run_parm in enumerate(runable_param):
                    run_parm = re.sub(r'##BT2DIR##',self.bin_dir,run_parm)
                    run_parm = re.sub(r'##DATADIR##',self.data_dir,run_parm)
                    runable_param[i] = run_parm

            return BenchmarkSet(set_data,self.data_dir,self.output_dir,self.bin_dir)
    
 
    
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
            logging.info("Starting: %s" % test["name"])
            logging.info("Metric: %s" % test["metric"])
            prologue = ''
            space    = " "
            err_out  = None
            
            if test["metric"] == "time":
                prologue = "/usr/bin/time -f %U,%S,%E " 
                err_out  = open(os.path.join(self.out_dir,test["name"])+".metric",'w')
            
            cmd = prologue + space + os.path.join(self.bin_dir,test["runable"]["program"])
            
            for opt in test["runable"]["options"]:
                cmd = cmd + space + opt   
            
            for parm in test["runable"]["parameters"]:
                cmd = cmd + space + parm  
            
            logging.info("running: %s" % cmd)
            subprocess.check_call(cmd,shell=True,stderr=err_out)
                 
            if err_out:
                err_out.close()    
            
            
    
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
                if not os.path.isfile(os.path.join(self.data_dir,data_file)):
                    data_is_here = False
            
            if not data_is_here:
                logging.info("Generate data for %s" % test["name"])
                for cmd in test["input_data"]["loading"]:
                    logging.info("running: %s" % cmd)
                    subprocess.check_call(cmd,shell=True)
                
        self.input_data_loaded = True
    
        
    
    
    
    
        