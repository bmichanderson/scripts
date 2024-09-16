#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: 19 Oct 2018
# Modified: Oct 2020, Sep 2024
# Description: read a multifasta file (arg1) and output a fasta file for each entry
##########################


from Bio import SeqIO
import sys


# check that a multifasta was provided
if len(sys.argv[1:]) == 0:
	print('This script needs a single argument (the multifasta file to split)')
	sys.exit(1)


# open the input file and write a fasta file for each entry
with open(sys.argv[1], 'r') as fasta_file:
	for fasta in SeqIO.parse(fasta_file, 'fasta'):
		outname = '_'.join(fasta.description.strip().split()) + '.fasta'
		with open(outname, 'w') as outfile:
			SeqIO.write(fasta, outfile, 'fasta')
