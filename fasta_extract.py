#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: 12 June 2020
# Modified: Oct 2020, Mar 2021, Aug 2024, Apr 2025 (adjust output fasta headers)
# Description: extract a sequence from a fasta given start and end coordinates (1-based and end-inclusive) or
#	extract an entry from a multifasta based on string matching
##########################


import argparse
import sys
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract a sequence from a fasta ' +
	'given start and end coords (1-based and end-inclusive). Alternatively, specify a string to extract an entry ' +
	'from a multifasta')


# add arguments to parse
parser.add_argument('fasta', type=str, help='A fasta file to extract from (single entry) using coordinates, ' +
	'or a multifasta using a search string')
parser.add_argument('-c', type=str, dest='coords', help='Coordinates in the form start..end, ' +
	'with start > end for compliment (1-based and end-inclusive)')
parser.add_argument('-s', type=str, dest='search', help='A string to search against fasta descriptions ' +
	'in a multifasta file')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
fasta_file = args.fasta
coords = args.coords
search = args.search

if not fasta_file:
	print('Please specify a fasta file to parse\n')
	parser.print_help(sys.stderr)
	sys.exit(1)

if all([not coords, not search]):
	print('Please specify coordinates with option -c, or a string to search with option -s\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# locate and extract the region or the fastas
with open(fasta_file, 'r') as f_file, open('extract.fasta', 'w') as out_file:
	if coords:
		if '..' in coords:
			start = int(coords.split('..')[0])
			end = int(coords.split('..')[1])
		else:
			print('Coordinates specified incorrectly\n')
			parser.print_help(sys.stderr)
			sys.exit(1)

		fasta = SeqIO.read(f_file, 'fasta')

		if any([start > len(fasta), end > len(fasta)]):
			print('Coordinates out of range\n')
			parser.print_help(sys.stderr)
			sys.exit(1)

		if start < end:		# forward strand
			fasta.seq = fasta.seq[(start - 1): end]
			fasta.id = fasta.id + '_sequence_' + coords
			fasta.description = fasta.description + '_sequence_' + coords
			SeqIO.write(fasta, out_file, 'fasta')
		elif start > end:	# reverse strand
			fasta.seq = fasta.seq[end-1: start].reverse_complement()
			fasta.id = fasta.id + '_sequence_' + coords
			fasta.description = fasta.description + '_sequence_' + coords
			SeqIO.write(fasta, out_file, 'fasta')
		else:			# if they are equal
			print('Please specify different start and end coordinates\n')
			parser.print_help(sys.stderr)
			sys.exit(1)

	elif search:
		fastas = SeqIO.parse(f_file, 'fasta')
		found = 0
		for fasta in fastas:
			if str(search) in fasta.description:
				found = found + 1
				SeqIO.write(fasta, out_file, 'fasta')
		print('Found ' + str(found) + ' entries matching ' + str(search))

	else:
		print('This shouldn\'t happen...')
		sys.exit(1)
