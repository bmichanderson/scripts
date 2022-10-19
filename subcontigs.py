#!/usr/bin/env python3

##########################
# Author: B. Anderson
# Date: 19 Oct 2022
# Description: break up a single fasta sequence into multiple contigs using a user-supplied text file
# 		with start and end coordinates, one per line (1-based and end-inclusive)
# Note:	the contigs will be saved in the current directory with the specified output prefix and a number
##########################


import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to break a single fasta sequence into multiple contigs \
	using a user-supplied text file with start and end points, tab separated and one per line (1-based and end-inclusive)')


# add arguments to parse
parser.add_argument('fasta', type=str, help='A fasta file to break up (single entry)')
parser.add_argument('-c', type=str, dest='coords', help='Coordinates text file with starts and ends, \
	tab separated and one start and end per line (1-based and end-inclusive)')
parser.add_argument('-o', type=str, dest='outpre', help='The output prefix for saving files [default output]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
fasta_file = args.fasta
coords = args.coords
outpre = args.outpre

if not fasta_file:
	sys.exit('Please provide a fasta file (no option)')

if not coords:
	sys.exit('Please specify a coordinates file with option -c')

if not outpre:
	outprefix = 'output'
else:
	outprefix = outpre


# read in the coordinates
contig_list = []
with open(coords, 'r') as coord_file:
		for line in coord_file:
			pieces = line.rstrip().split()
			if len(pieces) != 2:
				sys.exit('Coordinates specified incorrectly (should have two numbers per line)')
			else:
				contig_list.append(pieces)


# locate and extract the contigs, saving each to a new file
with open(fasta_file, 'r') as f_file:
	fasta = SeqIO.read(f_file, 'fasta')
	index = 1
	for contig in contig_list:
		start = int(contig[0])
		end = int(contig[1])

		if any([start < 0, start > len(fasta), end < 0, end > len(fasta)]):
			sys.exit('Coordinates out of range')

		if start > end:
			sys.exit('Please specify coordinates with start < end')

		new_contig = SeqRecord(
			fasta.seq[(start - 1): end],
			id = 'contig' + str(index),
			name = 'contig' + str(index),
			description = 'contig' + str(index) + ' from '  + fasta.id
		)

		with open(str(outprefix) + '_' + str(index) + '.fasta', 'w') as outfile:
			SeqIO.write(new_contig, outfile, 'fasta')

		index = index + 1
