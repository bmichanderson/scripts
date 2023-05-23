#!/usr/bin/env python3

#################
# Author: B. Anderson
# Date: May 2023
# Description: count the number and percentage of ambiguities in (multi)fastas
#################


import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to count ambiguities in fasta sequences')


# add arguments to parse
parser.add_argument('fastas', type=str, help='The (multi)fasta file(s)', nargs='*')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
fastas = args.fastas


# set up a list of IUPAC ambiguities  
amb_list = ['K', 'k', 'M', 'm', 'R', 'r', 'S', 's', 'W', 'w', 'Y', 'y']


# for each fasta (and entry within) report the counts  
counter = 0
print('\t'.join(['File', 'ID', 'Length', 'Length_nogapN', 'Length_ambig', '%_ambig']))
for fasta_file in fastas:
	filename = os.path.basename(fasta_file)
	seqs = []
	counter = counter + 1
	count_entries = 0
	with open(fasta_file, 'r') as infile:
		entries = SeqIO.parse(infile, 'fasta')
		for entry in entries:
			unknown_length = entry.seq.count('N') + entry.seq.count('n') + entry.seq.count('-')
			compare_length = len(entry.seq) - unknown_length
			ambig_length = sum([entry.seq.count(ambig) for ambig in amb_list])
			count_entries = count_entries + 1
			if (compare_length > 0):
				print('\t'.join([filename, entry.id, str(len(entry.seq)),
					str(compare_length), str(ambig_length),
					'%.2f' % (float(100*ambig_length)/compare_length)]))
			else:
				print('\t'.join([filename, entry.id, str(len(entry.seq)),
					str(compare_length), str(ambig_length), 'N/A']))
			seqs.append(entry.seq)
		if len(seqs) > 1:
			cum_seq = Seq('').join(seqs)
			unknown_length = cum_seq.count('N') + cum_seq.count('n') + cum_seq.count('-')
			compare_length = len(cum_seq) - unknown_length
			ambig_length = sum([cum_seq.count(ambig) for ambig in amb_list])
			if (compare_length > 0):
				print('\t'.join([filename, 'Total(' + str(count_entries) + ')', str(len(cum_seq)),
					str(compare_length), str(ambig_length),
					'%.2f' % (float(100*ambig_length)/compare_length)]))
			else:
				print('\t'.join([filename, 'Total(' + str(count_entries) + ')', str(len(cum_seq)),
					str(compare_length), str(ambig_length), 'N/A']))
