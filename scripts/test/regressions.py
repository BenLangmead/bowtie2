#!/usr/bin/env python

import os
import inspect
import unittest
import logging
import bt2face
import dataface
import btdata
from optparse import OptionParser


class TestRegressions(unittest.TestCase):
    """
    Main regression fixture.
    """
    
    def test_288(self):
        """ Check if --un option works correctly when used with --no-unal
        """
        # run --un and record unaligned file size and orig file size
        # run --no-unal and record file size
        # check if --no-unali + --un files size equal the orig file size
        out_sam     = 'test288_out.sam'
        not_algn    = 'test288_nal.fastq'

        ref_index = os.path.join(g_bdata.index_dir_path,'lambda_virus')
        reads     = os.path.join(g_bdata.reads_dir_path,'longreads.fq')

        args = "--quiet -x %s -U %s -a --un %s > %s" % (ref_index,reads,not_algn,out_sam)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
        sam_size     = dataface.SamFile(out_sam).size()
        no_algn_size = dataface.FastaQFile(not_algn).size()
        args = "--quiet -x %s -U %s -a --un %s --no-unal > %s" % (ref_index,reads,not_algn,out_sam)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
        no_al_sam_size = dataface.SamFile(out_sam).size()
        self.assertEqual(no_al_sam_size + no_algn_size, sam_size)
        os.remove(out_sam)
        os.remove(not_algn)
    
    
   
def get_suite():
    tests = ['test_288']
    return unittest.TestSuite(map(TestRegressions,tests))

            
def parse_args():
    usage = " %prog [options] \n\n"
    usage += "Some regression tests.\n"
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
    
    g_bdata = btdata.ExampleData(bt2_path)
    g_bt    = bt2face.BowtieSuite(bt2_path)
    runner.run(get_suite())

