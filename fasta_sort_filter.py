#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 3 Feb 2020
# Modified: Oct 2020
# Description: read a multifasta file (arg1), sort the entries by length, keep those >= or <= a specified length (arg2)
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files


def help():
	print('A script to filter a fasta file of contigs based on length, given a user specified length and whether min or max filtering.')
	print('The filtered entries are output as a mutlifasta file in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' fasta_file length_cutoff min/max')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 2:
	print('Arguments are fasta_file:  ' + str(sys.argv[1]) + ', length cutoff: ' + str(sys.argv[2]) + ', min/max: ' + str(sys.argv[3]))
else:
	sys.exit(help())


len_cut = int(sys.argv[2])

if sys.argv[3].lower() == 'max':
	mode = 'max'
elif sys.argv[3].lower() == 'min':
	mode = 'min'
else:
	sys.exit(help())


contig_list = []
len_exclude = 0
len_include = 0


with open(sys.argv[1], 'r') as contig_file, open('len_' + mode + str(len_cut) + '.fasta', 'w') as output_file:
	contig_iter = SeqIO.parse(contig_file, 'fasta')

	for contig in contig_iter:
		if any([all([mode == 'min', len(contig) >= len_cut]), all([mode == 'max', len(contig) <= len_cut])]):
			contig_list.append(contig)
			len_include = len_include + 1
		else:
			len_exclude = len_exclude + 1

	for contig in sorted(contig_list, key=lambda x: len(x), reverse=True):
		SeqIO.write(contig, output_file, 'fasta')

	# print a summary of reported lengths
	print('Found ' + str(len_exclude + len_include) + ' contigs, of which ' + str(len_include) + ' passed the filter.')
