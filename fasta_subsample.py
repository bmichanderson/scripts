#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 12 Feb 2019
# Updated: Oct 2020
# Description: subsample from a multifasta (arg1) so that no more than a set number (arg2) are present
##########################

import sys			# allows access to command line arguments
import random			# allows random sampling
from Bio import SeqIO		# for reading and writing Seq Record objects

def help():
	print('A script to subsample from a multifasta file to have no more than a given number.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' multi_fasta max')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
	print('Arguments are multi_fasta =  ' + str(sys.argv[1]) + ' and max = ' + str(sys.argv[2]))
	max_reads = int(sys.argv[2])
else:
	sys.exit(help())


# Subsample and output a new file
with open(sys.argv[1], 'r') as reads_file, open('subsampled_reads.fasta', 'w') as out_file:

	out_list = []
	reads_list = []

	fastas = SeqIO.parse(reads_file, 'fasta')

	for index, fasta in enumerate(fastas, start=0):
		if index == 0:		# first entry is the reference
			out_list.append(fasta)
		else:
			reads_list.append(fasta)

	if len(reads_list) > max_reads:
		rand_list = random.sample(reads_list, max_reads)
	else:
		rand_list = reads_list

	for random_read in rand_list:
		out_list.append(random_read)

	for out_fasta in out_list:
		SeqIO.write(out_fasta, out_file, 'fasta')
