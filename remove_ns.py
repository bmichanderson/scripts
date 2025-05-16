#!/usr/bin/env python3

###############################################
# Author: B.M. Anderson
# Date: Sept 2024
# Modified: May 2025 (add gaps and argument to specify)
# Description: a script to remove all Ns and/or gaps ("-") from fasta sequences
# Note: fastas are output in the current directory with new prefix "mod_"
# Warning: all will be removed, so ensure this is desired behaviour (may create incorrect sequences)
###############################################


import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to remove Ns and/or gaps from sequences in (multi)fasta files. ' +
	'The new fasta files are output in the current directory with prefix \"mod_\"')


# add arguments to parse
parser.add_argument('fasta', type = str, help = 'The (multi)fasta file(s)', nargs = '*')
parser.add_argument('-r', type = str, dest = 'remove', help = 'Which to remove: \"Ns\" [default], \"gaps\", or \"both\"')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
fasta = args.fasta
remove = args.remove

if remove:
	if remove == 'Ns':
		remove_index = [0]
	elif remove == 'gaps':
		remove_index = [1]
	elif remove == 'both':
		remove_index = [0, 1]
	else:
		parser.print_help(sys.stderr)
		sys.exit(1)
else:
	remove_index = [0]


# read in the each (multi)fasta file and remove elements
remove_list = ['N', '-']
for fasta_file in fasta:
	with open(fasta_file, 'r') as input_file, \
		open('mod_' + os.path.basename(fasta_file), 'w') as output_file:
		if len(remove_index) > 1:		# two to remove
			for fasta in SeqIO.parse(input_file, 'fasta'):
				string1 = ''.join(str(fasta.seq).split(remove_list[remove_index[0]]))
				fasta.seq = Seq(''.join(string1.split(remove_list[remove_index[1]])))
				SeqIO.write(fasta, output_file, 'fasta')
		else:
			for fasta in SeqIO.parse(input_file, 'fasta'):
				fasta.seq = Seq(''.join(str(fasta.seq).split(remove_list[remove_index[0]])))
				SeqIO.write(fasta, output_file, 'fasta')
