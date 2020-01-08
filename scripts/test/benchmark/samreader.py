#!/usr/bin/env python3
"""
A reader of SAM format.

    import samreader
    
    # ...
    # Attach it to stdin and print the header.
    sam_reader = SamReader(sys.stdin)
    print sam_reader.header
    
    # ...
    # Open a SAM file and store each read starting position.
    with fopen(myfile) as fp:
        sr = SamReader(fp)
        for rec in sr:
            start_coord[rec.qname] = rec.pos
            

"""


class CigarString(object):
    """ 
    A basic Cigar. 
    """    
    
    def __init__(self, cigar):
        self.cigar = cigar
        self.set_idx = 0
        
        
    def __str__(self):
        return self.cigar
        

class SamHeader(object):
    """ 
    Sam Header. 
    """    
    
    def __init__(self, file_handler):
        self.header_lines = list()
        self.curr_idx = 0
        self._source_fh = file_handler   
        self.end_header_pointer = self._load_header()

    def __iter__(self):
        self.curr_idx = 0
        return self

    def __next__(self):
        if self.curr_idx == len(self.header_lines):
            raise StopIteration
            
        value = self.header_lines[self.curr_idx]
        self.curr_idx += 1
        return value

    def __str__(self):
        return "\n".join(self.header_lines)

    def _load_header(self):
        fh = self._source_fh
        fh.seek(0)
        
        last_pos = fh.tell()
        line = fh.readline().rstrip()
        while line[0] == '@':
            self.header_lines.append(line)
            last_pos = fh.tell()
            line = fh.readline().rstrip()
            
        fh.seek(last_pos)
        return last_pos
    

class SamRecord(object):
    """ 
    Record Item for SAM. 
    """    
    
    def __init__(self, sam_line):
        all_tokens = sam_line.rstrip().split('\t')
        # NOTE: not care about optional fields for now.
        self.rec_keys = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
        rec = {k: all_tokens[i] for i, k in enumerate(self.rec_keys)}
        rec['POS'] = int(rec['POS'])
        rec['MAPQ'] = int(rec['MAPQ'])
        rec['CIGAR'] = CigarString(rec['CIGAR'])
        rec['PNEXT'] = int(rec['PNEXT'])
        rec['TLEN'] = int(rec['TLEN']) 
        self.record = rec
        self.pos = rec['POS']
        self.qname = rec['QNAME']
        self.mapq = rec['MAPQ']
        self.rname = rec['RNAME']
        self.flag = int(rec['FLAG'])

    def __str__(self):
        return "\t".join([str(self.record[self.rec_keys[i]]) for i in range(len(self.rec_keys))])
    

class SamReader(object):
    """ 
    Iterable for all SAM records. 
    """    
    
    def __init__(self, file_handle):
        self._source_fh = file_handle
        self.header = SamHeader(file_handle)

    def __iter__(self):
        self._source_fh.seek(self.header.end_header_pointer)
        return self

    def __next__(self):
        line = self._source_fh.readline()
        if not line:
            raise StopIteration
        
        return SamRecord(line)


    
    


