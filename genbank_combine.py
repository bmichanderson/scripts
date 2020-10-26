#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 17 June 2019
# Modified: Oct 2020
# Description: combine multiple genbank files into a single multigenbank file
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files

def help():
	print('A script to combine multiple genbank files into a single multigenbank file.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' file1 file2 ...')
	print('')


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# capture list of files
file_list = []

for arg in sys.argv[1:]:
	file_list.append(arg)


# read in each file and combine into a single multigenbank file
gbfiles = []

for infile in file_list:
	gbfile = SeqIO.read(open(infile, 'r'), 'gb')
	gbfiles.append(gbfile)

filename = '_'.join(file_list[0].split('_')[:-1]) + '_multi.gb'		# assuming naming convention of 'genus_specificep_accession.gb'

with open(filename, 'w') as outfile:
	for gbfile in gbfiles:
		SeqIO.write(gbfile, outfile, 'gb')
