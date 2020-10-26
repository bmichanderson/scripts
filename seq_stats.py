#!/usr/bin/env python

#################
# Author: B. Anderson
# Date: 27 Apr 2020
# Modified: Oct 2020
# Description: output basic sequence statistics for a set of input fasta files
#################

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq

def help():
	print('A script to output basic sequence statistics for a set of input fasta files')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' fasta1 fasta2 fasta3 ...')
	print('')


# print help if the script is called without arguments
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# for each of the input files, read the sequence(s) then calculate and output stats
cum_seqs = []
for fasta in sys.argv[1:]:
	print(fasta)

	fsa = SeqIO.parse(open(fasta, 'r'), 'fasta')
	entries = []
	for entry in fsa:
		entries.append(entry)

	seqs = []
	for entry in entries:
		seqs.append(entry.seq)
		print('ID: ' + entry.id)
		n_count = entry.seq.count('N') + entry.seq.count('n')
		print('Length: ' + str(len(entry.seq)) + ', GC content: ' + ('%.3f' % GC(entry.seq)) + ', proportion Ns: ' + ('%.2f' % (float(100*n_count)/len(entry.seq))))


	cum_seq = Seq('').join(seqs)
	cum_seqs.append(cum_seq)
	n_count = cum_seq.count('N') + cum_seq.count('n')
	if len(seqs) > 1:
		print('Cumulative length: ' + str(len(cum_seq)) + ', GC content: ' + ('%.3f' % GC(cum_seq)) + ', proportion Ns: ' + ('%.2f' % (float(100*n_count)/len(cum_seq))))
	print('')

overall_seq = Seq('').join(cum_seqs)
n_count = overall_seq.count('N') + overall_seq.count('n')
if len(cum_seqs) > 1:
	print('Overall length: ' + str(len(overall_seq)) + ', GC content: ' + ('%.3f' % GC(overall_seq)) + ', proportion Ns: ' + ('%.2f' % (float(100*n_count)/len(overall_seq))))
