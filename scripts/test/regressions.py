#!/usr/bin/env python

import os
import inspect
import unittest
import logging
import shutil
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

        args = "--quiet -x %s -U %s -a --un %s -S %s" % (ref_index,reads,not_algn,out_sam)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
        sam_size     = dataface.SamFile(out_sam).size()
        no_algn_size = dataface.FastaQFile(not_algn).size()
        args = "--quiet -x %s -U %s -a --un %s --no-unal -S %s" % (ref_index,reads,not_algn,out_sam)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
        no_al_sam_size = dataface.SamFile(out_sam).size()
        self.assertEqual(no_al_sam_size + no_algn_size, sam_size)
        os.remove(out_sam)
        os.remove(not_algn)


    def test_279(self):
        """ Check if --un-conc option works correctly when used with --no-unal
        """
        out_no_conc = 'test279_nconc.fq'
        out_conc    = 'test279_conc.fq'
        out_sam     = 'test279.sam'
        out_p1_nc   = 'test279_nconc.1.fq'
        out_p2_nc   = 'test279_nconc.2.fq'
        out_p1_conc = 'test279_conc.1.fq'
        out_p2_conc = 'test279_conc.2.fq'
        ref_index   = os.path.join(g_bdata.index_dir_path,'lambda_virus')
        pairs_1     = os.path.join(g_bdata.reads_dir_path,'reads_1.fq')
        pairs_2     = os.path.join(g_bdata.reads_dir_path,'reads_2.fq')

        args = "--quiet -x %s -1 %s -2 %s --un-conc %s --dovetail --al-conc %s -S %s" % (ref_index,pairs_1,pairs_2,out_no_conc,out_conc,out_sam)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
        conc_p1_size_step1 = dataface.FastaQFile(out_p1_conc).size()
        conc_p2_size_step1 = dataface.FastaQFile(out_p2_conc).size()
        nc_p1_size_step1   = dataface.FastaQFile(out_p1_nc).size()
        nc_p2_size_step1   = dataface.FastaQFile(out_p2_nc).size()
        self.assertEqual(conc_p1_size_step1, conc_p2_size_step1, "Number of concordant sequences in pass 1 should match.")
        self.assertEqual(nc_p1_size_step1, nc_p2_size_step1, "Number of non-concordant sequences in pass 1 should match.")

        args = "--quiet -x %s -1 %s -2 %s --no-unal --un-conc %s --dovetail --al-conc %s -S %s" % (ref_index,pairs_1,pairs_2,out_no_conc,out_conc,out_sam)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)
        conc_p1_size_step2 = dataface.FastaQFile(out_p1_conc).size()
        conc_p2_size_step2 = dataface.FastaQFile(out_p2_conc).size()
        nc_p1_size_step2   = dataface.FastaQFile(out_p1_nc).size()
        nc_p2_size_step2   = dataface.FastaQFile(out_p2_nc).size()
        self.assertEqual(conc_p1_size_step2, conc_p2_size_step2, "Number of concordant sequences in pass 2 should match.")
        self.assertEqual(nc_p1_size_step2, nc_p2_size_step2, "Number of non-concordant sequences in pass 2 should match.")
        self.assertEqual(conc_p1_size_step1, conc_p1_size_step2, "Number of concordant sequences in both steps should match.")
        self.assertEqual(conc_p2_size_step1, conc_p2_size_step2, "Number of concordant sequences in both steps should match.")
        self.assertEqual(nc_p1_size_step1, nc_p1_size_step2, "Number of non-concordant sequences in both steps should match.")
        self.assertEqual(nc_p2_size_step1, nc_p2_size_step2, "Number of non-concordant sequences in both steps should match.")

        self.assertTrue(1,"temp")
        os.remove(out_p1_conc)
        os.remove(out_p2_conc)
        os.remove(out_p1_nc)
        os.remove(out_p2_nc)
        os.remove(out_sam)
         
    
    def test_version(self):
        """ Check if --version option works correctly. 
        """
        args = "--version"
        ret  = g_bt.silent_run(args)
        self.assertEqual(ret,0)
        
       
    def test_x_option(self):
        """ Check if -x is indeed mandatory. 
        """
        ref_index = os.path.join(g_bdata.index_dir_path,'lambda_virus')
        reads     = os.path.join(g_bdata.reads_dir_path,'longreads.fq')
        out_sam   = 'test_x_option.sam'

        args = "%s -U %s -S %s " % (ref_index,reads,out_sam)
        ret = g_bt.silent_run(args)
        self.assertNotEqual(ret,0)
        args = "-x %s -U %s -S %s " % (ref_index,reads,out_sam)
        ret = g_bt.silent_run(args)
        self.assertEqual(ret, 0)
        os.remove(out_sam)


    def test_313(self):
        """ Check if --un-conc, --al-conc, --un, --al do not fail in case of a dot 
            character present in their path name.
        """
        ref_index   = os.path.join(g_bdata.index_dir_path,'lambda_virus')
        pairs_1     = os.path.join(g_bdata.reads_dir_path,'reads_1.fq')
        pairs_2     = os.path.join(g_bdata.reads_dir_path,'reads_2.fq')
        no_dot_dir  = 'nodotdir'
        dot_dir     = 'v.12_.test'
        un_basename = 'bname'
        out_spec    = os.path.join(no_dot_dir,un_basename)
        out_spec_dot= os.path.join(dot_dir,un_basename)
        args        = "--quiet -x %s -1 %s -2 %s --un-conc %s -S %s" % (ref_index,pairs_1,pairs_2,out_spec,os.devnull)
        
        try:
            os.makedirs(no_dot_dir)
            os.makedirs(dot_dir)
        except:
            pass

        ret = g_bt.silent_run(args)
        self.assertEqual(ret, 0)
        mate1_count = dataface.FastaQFile(out_spec + '.1').size()
        args        = "--quiet -x %s -1 %s -2 %s --un-conc %s -S %s" % (ref_index,pairs_1,pairs_2,out_spec_dot,os.devnull)
        ret = g_bt.silent_run(args)
        self.assertEqual(ret, 0, "--un-conc should work with dots in its path.")
        mate1_count_dot = dataface.FastaQFile(out_spec_dot + '.1').size()
        self.assertEqual(mate1_count, mate1_count_dot, "Number of seqs should be the same.")

        args        = "--quiet -x %s -1 %s -2 %s --un-conc %s -S %s" % (ref_index,pairs_1,pairs_2,dot_dir,os.devnull)
        ret = g_bt.silent_run(args)
        self.assertEqual(ret, 0, "--un-conc should work with a path only.")
        mate1_fs    = os.path.join(dot_dir,"un-conc-mate.1")
        self.assertTrue(os.path.isfile(mate1_fs), "--un-conc output is missing for path only case.")
        mate1_count_dot_only = dataface.FastaQFile(mate1_fs).size()
        self.assertEqual(mate1_count, mate1_count_dot_only, "Number of seqs should be the same.")

        args= "--quiet -x %s -1 %s -2 %s --un-conc %s --al-conc %s --un %s -S %s" % (ref_index,pairs_1,pairs_2,dot_dir,dot_dir,dot_dir,os.devnull)
        ret = g_bt.run(args)
        self.assertEqual(ret, 0)
        mate1_al    = os.path.join(dot_dir,"al-conc-mate.1")
        mate1_un    = os.path.join(dot_dir,"un-seqs")
        self.assertTrue(
            os.path.isfile(mate1_fs) and 
            os.path.isfile(mate1_al) and 
            os.path.isfile(mate1_un),
            "--un-conc, --al-conc and --un should not clobber"
        )
        shutil.rmtree(no_dot_dir)
        shutil.rmtree(dot_dir)
        

   
def get_suite():
    loader = unittest.TestLoader()
    return loader.loadTestsFromTestCase(TestRegressions)
    
            
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

