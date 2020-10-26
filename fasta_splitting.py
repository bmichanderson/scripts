#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 19 Oct 2018
# Modified: Oct 2020
# Description: read a multifasta file (arg1) and output a fasta file for each entry
##########################


import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files


with open(sys.argv[1], 'r') as fasta_file:
	fasta_iter = SeqIO.parse(fasta_file, 'fasta')
	for entry in fasta_iter:
		with open('_'.join(entry.description.rstrip().split()) + '.fasta', 'w') as output_file:
			SeqIO.write(entry, output_file, 'fasta')
