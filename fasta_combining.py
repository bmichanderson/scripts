#!/usr/bin/env python3

##########################
# Author: B. Anderson
# Date: 10 Mar 2020
# Modified: April 2023; Oct 2020
# Description: combine fasta files, removing duplicate sequences
##########################


import sys
import argparse
from Bio import SeqIO		# for reading and writing sequence records


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to combine fasta entries in multiple files into a single file.')


# add arguments to parse
parser.add_argument('fasta_files', type=str, help='The fasta files', nargs='*')
parser.add_argument('-o', type = str, dest = 'out_file', help = 'The output file to create [default: combined.fasta]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
fasta_files = args.fasta_files
out_file = args.out_file

if not out_file:
	out_file = 'combined.fasta'


# extract the fastas from the files
fastas = []
for f_file in fasta_files:
	with open(f_file, 'r') as fasta_file:
		for entry in SeqIO.parse(fasta_file, 'fasta'):
			fastas.append(entry)


# remove exact duplicate entries
dlist = []
seqlist = []
to_remove = []
for num, fasta in enumerate(fastas):
	if fasta.description in dlist:
		d_indices = [i for i in range(len(dlist)) if dlist[i] == fasta.description]
		found = 'False'
		for index in d_indices:
			if fasta.seq == seqlist[index]:
				to_remove.append(num)
				found = 'True'
				break

		if found == 'False':
			dlist.append(fasta.description)
			seqlist.append(fasta.seq)
	else:
		dlist.append(fasta.description)
		seqlist.append(fasta.seq)

for index in sorted(to_remove, reverse = True):
	del fastas[index]


# output the fastas
with open(out_file, 'w') as output_file:
	for entry in fastas:
		SeqIO.write(entry, output_file, 'fasta')
