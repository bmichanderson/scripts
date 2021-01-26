#!/usr/bin/env python
##########################
# Author: B. Anderson
# Date: 26 Jan 2021
# Description: read a genbank or embl file (arg1) and write an embl or genbank (convert either way)
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files

def help():
	print('This script will convert a genbank file or embl file (arg1) into embl or genbank, respectively')
	print('The file should have the extension \".gb\" or \".embl\"')
	print('')

if len(sys.argv) != 2:		# print help if incorrect number of arguments provided
	sys.exit(help())


# check the extension on the file
genbank_present = False
embl_present = False
filepath = sys.argv[1]
if filepath.split('.')[-1] == 'gb':	#
	genbank_present = True
elif filepath.split('.')[-1] == 'embl':
	embl_present = True
else:
	sys.exti(help())


# determine which way to convert, then write the file
if genbank_present:
	with open(sys.argv[1], 'r') as gbfile, open(filepath.replace('.gb', '.embl'), 'w') as out_file:
		records = SeqIO.parse(gbfile, 'genbank')
		for record in records:
			SeqIO.write(record, out_file, 'embl')
else:
	with open(sys.argv[1], 'r') as emblfile, open(filepath.replace('.embl', '.gb'), 'w') as out_file:
		records = SeqIO.parse(emblfile, 'embl')
		for record in records:
			SeqIO.write(record, out_file, 'genbank')
