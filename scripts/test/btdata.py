#!/usr/bin/env python
"""
Note: This would look so much better replaced by XML or at least JSON. But 
      is not worth to do it for now.
"""

import os
import gzip
import urllib2
import logging


class LargeTestsData(object):
    """
    Large index tests use quite big datasets. Make sure these
    are present before starting the time consuming tests.
    """
    
    def __init__(self,bt2_path=''):
        self.data_dir       = 'big_data'
        curr_path           = os.path.realpath(bt2_path)

        curr_path           = os.path.join(curr_path,'scripts')
        curr_path           = os.path.join(curr_path,'test')
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
                
  
  
class ExampleData(object):
    """ The example data.
    """
    
    def __init__(self,bt2_path=''):
        curr_path           = os.path.realpath(bt2_path)
        curr_path           = os.path.join(curr_path,'example')
        self.index_dir_path = os.path.join(curr_path,'index')
        self.reads_dir_path = os.path.join(curr_path,'reads')
        self.ref_dir_path   = os.path.join(curr_path,'reference')
        
        try:
            os.stat(curr_path)
        except:
            logging.error("Cannot find the example datadir %s!" % curr_path)
            raise 
    
    
    
    
    
                

