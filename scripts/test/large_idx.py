#!/usr/bin/env python

import os
import gzip
import urllib2
import inspect
import unittest
import logging
import bt2face
import btdata
import subprocess
from optparse import OptionParser

              

class TestLargeIndex(unittest.TestCase):
    """
    Main fixture for large index testing.
    """
    
    def test_human(self):
        wdir = g_bdata.data_dir_path
        wdir = os.path.join(wdir,'human')
        rdir = g_bdata.reads_dir_path
        genome = g_bdata.genomes['human']
        genome_fasta = os.path.join(wdir,genome['ref_name'])
        genome_index = os.path.join(wdir,'human')
        reads        = os.path.join(rdir,'human_reads.fa')
        ret = g_bt.build("%s human" % genome_fasta)
        self.assertEqual(ret,0)
        args = "-x %s -f -U %s" % (genome_index,reads)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
    
    
    def test_mouse(self):
        wdir = g_bdata.data_dir_path
        wdir = os.path.join(wdir,'mouse')
        rdir = g_bdata.reads_dir_path
        genome = g_bdata.genomes['mouse']
        genome_fasta = os.path.join(wdir,genome['ref_name'])
        genome_index = os.path.join(wdir,'mouse')
        reads        = os.path.join(rdir,'mouse_reads.fa')
        ret = g_bt.build("%s mouse" % genome_fasta)
        self.assertEqual(ret,0)
        args = "-x %s -f -U %s" % (genome_index,reads)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
    
    
    def test_large_index(self):
        wdir = g_bdata.data_dir_path
        wdir = os.path.join(wdir,'ms_hum')
        rdir = g_bdata.reads_dir_path
        genome = g_bdata.joint_genomes['ms_hum']
        genome_fasta = os.path.join(wdir,genome['ref_name'])
        genome_index = os.path.join(wdir,'ms_hum')
        reads_human  = os.path.join(rdir,'human_reads.fa')
        reads_mouse  = os.path.join(rdir,'mouse_reads.fa')
        ret = g_bt.build("%s ms_hum" % genome_fasta)
        self.assertEqual(ret,0)
        args = "-x %s -f -U %s" % (genome_index,reads_human)
        ret = g_bt.run(args)
        args = "-x %s -f -U %s" % (genome_index,reads_mouse)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)


   
def get_suite():
    tests = ['test_human','test_mouse','test_large_index']
    return unittest.TestSuite(map(TestLargeIndex,tests))

    
            
def parse_args():
    usage = " %prog [options] \n\n"
    usage += "Warning, this test runs some VERY resource consuming tests.\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", 
                    action="store_true",dest="verbose", default=False,
                    help="Print more info about each test.")

    (options, args) = parser.parse_args()
    return options
    
    
g_bdata = None
g_bt    = None

if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.ERROR)    
    options = parse_args()

    runner = unittest.TextTestRunner()
    if options.verbose:
        logging.getLogger().setLevel(level=logging.INFO)
        runner = unittest.TextTestRunner(verbosity=2)
   
    src_file_path  = os.path.realpath(inspect.getsourcefile(parse_args))
    curr_path      = os.path.dirname(src_file_path)
    bw2_subdir     = 'bowtie2'

    i = curr_path.find(bw2_subdir)
    bt2_path = curr_path[:i+len(bw2_subdir)] 
    
    g_bdata = btdata.LargeTestsData(bt2_path)
    g_bt    = bt2face.BowtieSuite(bt2_path)
    runner.run(get_suite())

