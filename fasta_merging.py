#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 18 May 2020
# Modified: Nov 2020
# Description: read a multifasta file and output a fasta file with the sequences merged together into one
##########################

import sys				# allows access to command line arguments
import re				# for use of regular expressions
from Bio import SeqIO			# for reading and writing sequence records

def help():
	print('A script to combine entries of a multifasta into a single fasta.')
	print('The fasta is output as a single file called merged.fasta in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' file')
	print('')


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# read in and merge sequences
seqs = []
with open(sys.argv[1], 'r') as fasta_file, open('merged.fasta', 'w') as out_file:
	fasta_iter = SeqIO.parse(fasta_file, 'fasta')
	index = 1
	for fasta in fasta_iter:
		if index == 1:
			index = index + 1
			fasta1 = fasta
			continue
		index = index + 1
		sequence = fasta.seq
		#sequence = str(sequence).replace('*', 'N')	# substitute for the strange * characters in some NOVOPlasty outputs (?)
		#sequence = re.sub('\d', 'N', str(sequence))	# substitute any digits in some NOVOPlasty outputs (?)
		if len(sequence) > 0:
			seqs.append(sequence)
		else:
			continue

	concatenated = fasta1.seq
	for sequence in seqs:
		concatenated = concatenated + sequence

	fasta1.seq = concatenated
	fasta1.description = 'merged ' + fasta1.description
	fasta1.id = 'merged'
	fasta1.name = 'merged'

	SeqIO.write(fasta1, out_file, 'fasta')
