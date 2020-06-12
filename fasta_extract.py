#!/usr/bin/env python
##########################
# Author: B. Anderson
# Date: 12 June 2020
# Description: extract a sequence from a fasta given start and end coordinates (1-based and end-inclusive)
##########################

import argparse
import sys			# allows access to command line arguments
from Bio import SeqIO		# for reading and writing Seq Record objects


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract a sequence from a fasta given start and end coords (1-based and end-inclusive)')


# add arguments to parse
parser.add_argument('fasta', type=str, help='A fasta file to extract from (single entry)')
parser.add_argument('-c', type=str, dest='coords', help='Coordinates in the form start..end, with start > end for compliment (1-based and end-inclusive)')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

fasta_file = args.fasta
coords = args.coords

if not coords:
	sys.exit('Please specify coordinates with option -c')


# locate and extract the region
with open(fasta_file, 'r') as f_file, open('extract.fasta', 'w') as out_file:
	if '..' in coords:
		start = int(coords.split('..')[0])
		end = int(coords.split('..')[1])
	else:
		sys.exit('Coordinates specified incorrectly')

	fasta = SeqIO.read(f_file, 'fasta')

	if any([start > len(fasta), end > len(fasta)]):
		sys.exit('Coordinates out of range')

	if start < end:		# forward strand
		out_file.write('>%s from %s\n%s\n' % ('sequence_' + coords, fasta.id, fasta.seq[(start - 1): end]))

	elif start > end:	# reverse strand
		out_file.write('>%s from %s\n%s\n' % ('sequence_' + coords, fasta.id, fasta.seq[end-1: start].reverse_complement()))

	else:			# if they are equal
		sys.exit('Please specify different start and end coordinates')
