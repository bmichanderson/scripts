#!/usr/bin/env python3

###############################################
# Author: B.M. Anderson
# Date: Sept 2024
# Description: a script to remove Ns from fasta sequences
# Note: fastas are output in the current directory with new prefix "mod_"
# Warning: if there are internal Ns, this will create false sequences
###############################################


import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to remove Ns from sequences in (multi)fasta files. ' +
	'The new fasta files are output in the current directory with prefix \"mod_\"')


# add arguments to parse
parser.add_argument('fasta', type = str, help = 'The (multi)fasta file(s)', nargs = '*')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
fasta = args.fasta


# read in the each (multi)fasta file and remove Ns
for fasta_file in fasta:
	with open(fasta_file, 'r') as input_file, \
		open('mod_' + os.path.basename(fasta_file), 'w') as output_file:
		for fasta in SeqIO.parse(input_file, 'fasta'):
			fasta.seq = Seq(''.join(str(fasta.seq).split('N')))
			SeqIO.write(fasta, output_file, 'fasta')
