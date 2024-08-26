#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Aug 2024
# Description: extend a linear representation of a circular contig (fasta) by a specified number of basepairs
# Note: This might be used to map to a longer representation of a plastome (to extend the overlaps)  
#	Will output a longer fasta in the current directory named `new_contig.fasta`
##########################


import argparse
from Bio import SeqIO
import sys



# instantiate the parser
parser = argparse.ArgumentParser(
	description = 'A script to extend a linear representation of a circular sequence (fasta file) by a specified amount')


# add arguments to parse
parser.add_argument('contig', type=str, help='The circular contig to extend')
parser.add_argument('-b', type=int, dest='bps', help='How much to extend each end by (default 150)')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

contig = args.contig
bps = args.bps


# assign variables
if not bps:
	bps = 150


# create new extended fasta
with open(contig, 'r') as contig_file, open('new_contig.fasta', 'w') as out_file:
	fasta = SeqIO.read(contig_file, 'fasta')
	add_front = fasta.seq[len(fasta.seq) - bps: ]		# add the end of the circle before the start
	add_back = fasta.seq[0: bps]		# add the start of the circle to the end
	fasta.seq = add_front + fasta.seq + add_back
	SeqIO.write(fasta, out_file, 'fasta')
