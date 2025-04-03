#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: 28 May 2020
# Modified: Apr 2025 (cleaned up and made more resilient to file naming)
# Description: read a (multi)genbank file (arg1) and write a (multi)fasta
##########################


import sys			# allows access to command line arguments
import argparse		# for parsing arguments
import os			# for interacting with the operating system and paths
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files


# instantiate the parser
parser = argparse.ArgumentParser(description = 'Convert (multi)genbank file(s) into (multi)fasta file(s). ' +
	'The output file(s) will share the same name as the input but with extension \".fasta\"')


# add arguments to parse
parser.add_argument('genbanks', type=str, help='The (multi)genbank files', nargs='*')


# parse the command line
if len(sys.argv[1:]) == 0:			# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
genbanks = args.genbanks


# read in the genbank(s) and write to fasta(s)
for genbank in genbanks:
	prefix, extension = os.path.splitext(genbank)
	with open(genbank, 'r') as gbfile, open(prefix + '.fasta', 'w') as outfile:
		records = SeqIO.parse(gbfile, 'genbank')
		for record in records:
			SeqIO.write(record, outfile, 'fasta')
