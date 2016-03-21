#!/usr/bin/env python

__author__ = 'vantones'

from optparse import OptionParser


class FqRecord(object):
    """
    Record Item for FQ.
    """

    def __init__(self, seq_id, seq, qual):
        self.seq_id = seq_id
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return "\t".join([self.seq_id, self.seq, self.qual])


class FqReader(object):
    """
    Iterable for all FQ records.
    """

    def __init__(self, file_handle):
        self._source_fh = file_handle

    def __iter__(self):
        self._source_fh.seek(0)
        return self

    def next(self):
        seq_id = self._source_fh.readline().rstrip()[1:]
        seq = self._source_fh.readline().rstrip()
        useless_line = self._source_fh.readline()
        qual = self._source_fh.readline().rstrip()
        if not (seq_id or seq or useless_line or qual):
            raise StopIteration

        return FqRecord(seq_id, seq, qual)


def parse_args():
    usage = " %prog [options] \n\n"
    usage += "Extracts a specific number of FastQ records and writes them to stdout\n"
    usage += "as C records.\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--fastq-file",
                      action="store", type="string", dest="fastq_file", default=None,
                      help="(mandatory).FastQ file to read from.")
    parser.add_option("-c", "--count",
                      action="store", type="int", dest="seqs_count", default=100,
                      help="(Default: 100). Number of sequences to pull.")

    (options, args) = parser.parse_args()

    if options.fastq_file is None:
        print("Mandatory option is missing (--fastq-file)")
        parser.print_help()
        exit(-1)

    return options


if __name__ == "__main__":
    options = parse_args()

    num_seqs = 1
    with open(options.fastq_file, "r") as fin:
        fq = FqReader(fin)
        for rec in fq:
            if num_seqs > options.seqs_count:
                break
            num_seqs += 1
            print "  {"
            print "   \"" + rec.seq_id + "\","
            print "   \"" + rec.seq + "\","
            print "   \"" + rec.qual.replace('\\', '\\\\').replace('"', '\\"').replace('?', '\\?') + "\""
            print "  },"
