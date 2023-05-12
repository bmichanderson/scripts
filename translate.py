#!/usr/bin/env python3

##########################
# Author: Benjamin Anderson
# Date: March 2023
# Updated: May 2023
# Description: translate nucleotide (multi)fasta files to protein
# Note: it will create a file(s) with the same name, adding a "_prot" before the extension
##########################


import sys
import argparse
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to translate fasta files from nucleotide to protein')


# add arguments to parse
parser.add_argument('fastas', type=str, help='The (multi)fasta file(s)', nargs='*')
parser.add_argument('-c', type=int, dest = 'code', help='The NCBI genetic translation table to use [default = 1]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
fastas = args.fastas
code = args.code

if not code:
	code = 1


# read in the (multi)fasta(s) and translate
counter = 0
for fasta_file in fastas:
	counter = counter + 1
	count_entries = 0
	with open(fasta_file, 'r') as infile, \
		open(fasta_file.replace('.fasta', '_prot.fasta'), 'w') as outfile:
		entries = SeqIO.parse(infile, 'fasta')
		for entry in entries:
			entry.seq = entry.seq.translate(table = code)
			SeqIO.write(entry, outfile, 'fasta')
			count_entries = count_entries + 1
	print('Translated ' + str(count_entries) + ' samples from ' + fasta_file)

# record completion
print('Finished after translating ' + str(counter) + ' (multi)fasta files')
