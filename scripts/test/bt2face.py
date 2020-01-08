#!/usr/bin/env python3

import os
import logging
import subprocess


class BowtieSuite(object):

    def __init__(self,bt2_path=''):
        curr_path           = os.path.realpath(bt2_path)
        self.bowtie_bin     = os.path.join(curr_path,'bowtie2')
        self.bowtie_build   = os.path.join(curr_path,'bowtie2-build')
        self.bowtie_inspect = os.path.join(curr_path,'bowtie2-inspect')
        logging.info('bowtie2 path: ' + self.bowtie_bin)
        logging.info('bowtie2-build path: ' + self.bowtie_build)
        logging.info('bowtie2-inspect path: ' + self.bowtie_inspect)


    def run(self, *args):
        cmd = self.bowtie_bin + " " + " ".join([i for i in args])
        logging.debug('align cmd: ' + cmd)
        return(subprocess.call(cmd,shell=True))


    def silent_run(self, *args):
        cmd = self.bowtie_bin + " " + " ".join([i for i in args])
        logging.debug('align cmd (silent): ' + cmd)
        return(subprocess.call(cmd,shell=True,stderr=open(os.devnull, 'w')))


    def build(self, *args):
        cmd = self.bowtie_build + " " + " ".join([i for i in args])
        curr_dir = os.getcwd()
        os.chdir(os.path.dirname(cmd.split()[-2]))
        logging.debug('build cmd: ' + cmd)
        ret = subprocess.check_call(cmd,shell=True)
        os.chdir(curr_dir)
        return(ret)
