#!/usr/bin/env python

"""
sa.py

Parse and possibly sanity-check a .sa file output by bowtie2-build in --sa
mode.  These files have a very simple format: first is a uint32_t containing
the length of the suffix array, the rest is an array of that many uint32_ts
containing the suffix array.
"""

import sys
import struct

def loadBowtieSa(fh):
	""" Load a .sa file from handle into an array of ints """
	nsa = struct.unpack('I', fh.read(4))[0]
	return [ struct.unpack('I', fh.read(4))[0] for i in xrange(0, nsa) ]

def loadBowtieSaFilename(fn):
	""" Load a .sa file from filename into an array of ints """
	with open(fn, 'rb') as fh:
		return loadBowtieSa(fh)

def loadFasta(fns):
	""" Load the concatenation of all the A/C/G/T characters """
	falist = []
	dna = set(['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'])
	for fn in fns:
		with open(fn, 'r') as fh:
			for line in fh:
				if line[0] == '>':
					continue
				for c in line:
					if c in dna:
						falist.append(c)
	return ''.join(falist)

if __name__ == "__main__":
	import argparse
	
	parser = argparse.ArgumentParser(\
		description='Parse suffix array built from bowtie2-build')
	
	parser.add_argument(\
		'--sa', metavar='string', required=True, type=str,
		help='Suffix array file')
	parser.add_argument(\
		'--fa', metavar='string', type=str, nargs='+', help='FASTA file')

	args = parser.parse_args()
	
	def go():
		ref = None
		if args.fa is not None:
			ref = loadFasta(args.fa)
		sas = loadBowtieSaFilename(args.sa)
		# Suffix array is in sas; note that $ is considered greater than all
		# other characters
		if ref is not None:
			for i in xrange(1, len(sas)):
				sa1, sa2 = sas[i-1], sas[i]
				assert sa1 != sa2
				# Sanity check that suffixes are really in order
				while sa1 < len(ref) and sa2 < len(ref):
					if ref[sa1] < ref[sa2]:
						break
					assert ref[sa1] == ref[sa2]
					sa1 += 1
					sa2 += 1
				else:
					# Note: Bowtie treats $ as greater than all other
					# characters; so if these strings are tied up to the end of
					# one or the other, the longer string is prior
					assert sa1 < sa2, "%d, %d" % (sas[i-1], sas[i])
			assert sas[-1] == len(ref)
	
	go()
	