#!/usr/bin/env python
##########################
# Author: B. Anderson
# Date: 28 May 2020
# Description: read a genbank file (arg1) and write a fasta
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files

def help():
	print('This script will convert a (multi)genbank file (arg1) into a single (multi)fasta file')
	print('')

if len(sys.argv) != 2:		# print help if incorrect number of arguments provided
	sys.exit(help())


# check the extension on the genbank file
filepath = sys.argv[1]
if filepath.split('.')[-1] != '.gb':
	filepath = ''.join(filepath.split('.')[:-1]) + '.gb'

# open the genbank, then write the fasta
with open(sys.argv[1], 'r') as gbfile, open(filepath.replace('.gb', '.fasta'), 'w') as out_file:
	records = SeqIO.parse(gbfile, 'genbank')
	for record in records:
		SeqIO.write(record, out_file, 'fasta')
