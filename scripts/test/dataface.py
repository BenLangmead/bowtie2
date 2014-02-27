#!/usr/bin/env python

import os
import logging


class DataFace(object):
    """ Some data IFace bowtie can work with.
    """
    
    def size(self):
        raise NotImplementedError("size() needs to be implemented!")
 
 
 
class SamFile(DataFace):
    
    def __init__(self,sam_desc):
        self.file_desc = sam_desc
        
        
    def size(self):
        if hasattr(self,'no_frags'):
            return self.no_frags
        
        count = 0
        try:
            fh = open(self.file_desc,"r")
            for line in fh:
                if line[0] != "@":
                    count += 1
        except:
            logging.error("Exception reading sam file!")
            fh.close()
            raise
        else:
            fh.close()
    
        self.no_frags = count
        return count
            


class FastaQFile(DataFace):
    
    def __init__(self,fastq_desc):
        self.file_desc = fastq_desc
        
        
    def size(self):
        if hasattr(self,'no_seqs'):
            return self.no_seqs
        
        count = 0
        try:
            fh = open(self.file_desc,"r")
            lno = 0
            for line in fh:
                if lno == 0:
                    if line[0] != '@':
                        raise TypeError("This does not look like a fastq file!")
                    count += 1
                if lno == 2:
                    if line[0] != '+':
                        raise TypeError("This does not look like a fastq file!")
                lno = (lno + 1) % 4 
        except:
            logging.error("Exception reading fastq file!")
            fh.close()
            raise
        else:
            fh.close()
    
        self.no_seqs = count
        return self.no_seqs



    
    
    