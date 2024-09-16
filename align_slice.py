#!/usr/bin/env python3

###############################################
# Author: B.M. Anderson
# Date: Sep 2024
# Description: extract a portion of a fasta alignment as a new file in the current directory with prefix "mod_"
###############################################


import argparse
from Bio import AlignIO
import os
import sys


# instantiate the parser and add args
parser = argparse.ArgumentParser(description = 'A script to slice a fasta alignment from start to end positions provided. ' +
	'Please provide start and/or end.')
parser.add_argument('align', type = str, help = 'The fasta alignment file', nargs = 1)
parser.add_argument('-s', type = int, dest = 'start', help = 'The start position (1-based); default = 1')
parser.add_argument('-e', type = int, dest = 'end', help = 'The end position (1-based, inclusive); default = end of alignment')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
align = args.align[0]
start = args.start
end = args.end

if not align:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not start:
	if not end:
		print('Neither start nor end specified! Exiting')
		sys.exit(1)
	start = 1


# read in the alignment and slice it, then write output
with open(align, 'r') as infile, open('mod_' + os.path.basename(align), 'w') as outfile:
	alignment = AlignIO.read(infile, 'fasta')
	if not end:
		end = alignment.get_alignment_length()
	elif end > alignment.get_alignment_length():
		print('End specified longer than alignment! Exiting')
		sys.exit(1)
	if start < 0:
		print('Start specified as negative! Exiting')
		sys.exit(1)
	elif start >= alignment.get_alignment_length() - 1:
		print('Start specified as longer than alignment! Exiting')
		sys.exit(1)

	slice = alignment[:, start - 1: end]
	AlignIO.write(slice, outfile, 'fasta')
