#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 18 May 2020
# Modified: Nov 2020
# Description: read a multifasta file and output a fasta file with the sequences merged together into one
##########################


import argparse
import sys				# allows access to command line arguments
#import re				# for use of regular expressions
from Bio import SeqIO			# for reading and writing sequence records
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to combine entries of a multifasta into a single fasta, with the option of a spacer')


# add arguments to parse
parser.add_argument('contigs', type=str, help='The (multi)fasta file to combine')
parser.add_argument('-s', type=str, dest='spacer', help='Specify a text string of bases to be used as a spacer between contigs, e.g. N')
parser.add_argument('-r', type=int, dest='repeats', help='Specify the number of times the string is repeated in the spacer, e.g. 100 [default = 1]')


# parse the command line
if len(sys.argv[1:]) == 0:              # if there are no arguments
        parser.print_help(sys.stderr)
        sys.exit(1)

args = parser.parse_args()

con = args.contigs
spacer = args.spacer
repeats = args.repeats

if not repeats:
	repeats = 1

if not spacer:
        spacer_present = False
else:
	spacer_present = True
	spacer = Seq(spacer*repeats)


# read in and merge sequences
fastas = []
with open(con, 'r') as fasta_file, open('merged.fasta', 'w') as out_file:
	for fasta in SeqIO.parse(fasta_file, 'fasta'):
		#sequence = str(fasta.seq).replace('*', 'N')	# substitute for the strange * characters in some NOVOPlasty outputs (?)
		#sequence = re.sub('\d', 'N', str(fasta.seq))	# substitute any digits in some NOVOPlasty outputs (?)
		#fasta.seq = Seq(sequence)
		if len(fasta.seq) > 0:
			fastas.append(fasta)
		else:
			continue

	concatenated = fastas[0]

	if spacer_present:
		concatenated.seq = concatenated.seq + spacer + spacer.join([x.seq for x in fastas[1:]])
	else:
		for fasta in fastas[1:]:
			concatenated.seq = concatenated.seq + fasta.seq

	concatenated.description = 'merged ' + concatenated.description
	concatenated.id = 'merged'
	concatenated.name = 'merged'

	SeqIO.write(concatenated, out_file, 'fasta')
