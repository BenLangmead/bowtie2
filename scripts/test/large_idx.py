#!/usr/bin/env python

import os
import gzip
import urllib2
import inspect
import unittest
import logging
import subprocess
from optparse import OptionParser


class LargeTestsData(object):
    """
    Large tests use quite big datasets. Make sure these
    are present before starting the time consuming tests.
    """
    
    def __init__(self):
        self.data_dir       = 'big_data'
        src_file_path       = os.path.realpath(inspect.getsourcefile(parse_args))
        curr_path           = os.path.dirname(src_file_path)
        self.data_dir_path  = os.path.join(curr_path,self.data_dir)
        self.reads_dir_path = os.path.join(curr_path,'reads')
        
        try:
            os.stat(self.data_dir_path)
        except:
            logging.error("Cannot find the working datadir %s!" % self.data_dir_path)
            raise 
        
        self.genomes = dict()
        self.genomes['human'] = dict()
        hm = self.genomes['human']
        hm['link'] = "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
        hm['ref_name'] = 'human.fa'
        hm['chromosomes'] = []
        chromosomes = hm['chromosomes']
        for i in range(1,22):
            chromosomes.append('chr%d' % i)
        chromosomes.extend(['chrX', 'chrY', 'chrM'])
        
        self.genomes['mouse'] = dict()
        ms = self.genomes['mouse']
        ms['link'] = "ftp://hgdownload.cse.ucsc.edu/goldenPath/rn4/chromosomes"
        ms['ref_name'] = 'mouse.fa'
        ms['chromosomes'] = []
        chromosomes = ms['chromosomes']
        for i in range(1,21):
            chromosomes.append('chr%d' % i)
        chromosomes.extend(['chrX', 'chrM'])
           
        self.joint_genomes = dict()
        self.joint_genomes['ms_hum'] = dict()
        mh = self.joint_genomes['ms_hum']
        mh['link'] = None
        mh['ref_name'] = 'ms_hum.fa'
        mh['genomes'] = ['human','mouse']

        self.init_data()
        


    def init_data(self):
        """ Try and init the data we need.
        """
        for genome,gdata in self.genomes.iteritems():
            gn_path  = os.path.join(self.data_dir_path,genome)
            gn_fasta = os.path.join(gn_path,gdata['ref_name'])
            if not os.path.exists(gn_fasta):
                self._get_genome(genome)
                self._build_genome(genome)
        
        for genome,gdata in self.joint_genomes.iteritems():
            gn_path  = os.path.join(self.data_dir_path,genome)
            gn_fasta = os.path.join(gn_path,gdata['ref_name'])
            if not os.path.exists(gn_fasta):
                self._build_joint_genome(genome)

            
        
    def _get_genome(self,genome):
        g        = self.genomes[genome]
        gn_path  = os.path.join(self.data_dir_path,genome)
        
        if not os.path.exists(gn_path):
            os.mkdir(gn_path)
            
        logging.info("Downloading genome: %s " % genome)
        
        for chrs in g['chromosomes']:
            chr_file = chrs + ".fa.gz"
            fname = os.path.join(gn_path,chr_file)
            
            if os.path.exists(fname):
                logging.info("Skip %s (already present)" % chr_file)
                continue
            
            uri = g['link'] + r"/" + chr_file 
            logging.info("file: %s" % chr_file)
            
            try:
                f = open(fname,'wb')
                u = urllib2.urlopen(uri)
                f.write(u.read())
            except:
                f.close()
                os.remove(fname)
                os.close(u.fileno())
                raise
            else:
                os.close(u.fileno())
                u.close()
                f.close()
                
        
        
    def _build_genome(self,genome):
        g        = self.genomes[genome]
        gn_path  = os.path.join(self.data_dir_path,genome)
        gn_fasta = os.path.join(gn_path,g['ref_name'])

        logging.info("Building fasta file for genome: %s" % genome)

        f_gn = open(gn_fasta,'wb')
        
        for chrs in g['chromosomes']:
            chr_file = chrs + ".fa.gz"
            fname = os.path.join(gn_path,chr_file)
        
            try:
                f_chr = gzip.open(fname,'rb')
                f_gn.write(f_chr.read())
            except:
                f_chr_close()
                f_gn.close()
                os.remove(gn_fasta)
                raise
            else:
                f_chr.close()

        f_gn.close()

        

    def _build_joint_genome(self,genome):
        jg        = self.joint_genomes[genome]
        jgn_path  = os.path.join(self.data_dir_path,genome)
        jgn_fasta = os.path.join(jgn_path,jg['ref_name'])

        if not os.path.exists(jgn_path):
            os.mkdir(jgn_path)
         
        logging.info("Building fasta file for genome: %s" % genome)

        f_jg = open(jgn_fasta,'wb')
        for g in jg['genomes']:
            gn_path    = os.path.join(self.data_dir_path,g)
            fasta_file = os.path.join(gn_path,self.genomes[g]['ref_name'])
            try:
                fin = open(fasta_file,'rb')
                f_jg.write(fin.read())
            except:
                fin.close()
                f_jg.close()
                os.remove(jgn_fasta)
                raise
            else:
                fin.close()
                
        f_jg.close()
                
              

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
        reads_human  = os.path.join(rdir,'human_reads.fa'
        reads_mouse  = os.path.join(rdir,'mouse_reads.fa'
        ret = g_bt.build("%s ms_hum" % genome_fasta)
        self.assertEqual(ret,0)
        args = "-x %s -f -U %s" % (genome_index,reads_human)
        ret = g_bt.run(args)
        args = "-x %s -f -U %s" % (genome_index,reads_mouse)
        ret = g_bt.run(args)
        self.assertEqual(ret,0)



class BowtieSuite(object):
    
    def __init__(self):
        src_file_path       = os.path.realpath(inspect.getsourcefile(parse_args))
        curr_path           = os.path.dirname(src_file_path)
        bw2_subdir = 'bowtie2'
        i = curr_path.find(bw2_subdir)
        bw2_path = curr_path[:i+len(bw2_subdir)]
        self.bowtie_bin     = os.path.join(bw2_path,'bowtie2')
        self.bowtie_build   = os.path.join(bw2_path,'bowtie2-build')


    def run(self, *args):
        cmd = self.bowtie_bin + " " + " ".join([i for i in args])
        return(subprocess.check_call(cmd,shell=True))        
        
        
    def build(self, *args):
        cmd = self.bowtie_build + " " + " ".join([i for i in args])
        curr_dir = os.getcwd()
        os.chdir(os.path.dirname(cmd.split()[-2]))
        ret = subprocess.check_call(cmd,shell=True)
        os.chdir(curr_dir)
        return(ret)
        
    
   
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
    g_bdata = LargeTestsData()
    g_bt    = BowtieSuite()
    runner.run(get_suite())

