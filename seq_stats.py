#!/usr/bin/env python

#################
# Author: B. Anderson
# Date: 27 Apr 2020
# Modified: Mar 2023 (add gaps), Nov 2020 (changed from proportion to percentage Ns), Oct 2022 (changed screen output format)
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
index = 1
print('\t'.join(['File', 'ID', 'Length', 'GC', '% N', '% Gaps']))
for fasta in sys.argv[1:]:
	fsa = SeqIO.parse(open(fasta, 'r'), 'fasta')
	entries = []
	for entry in fsa:
		entries.append(entry)

	seqs = []
	for entry in entries:
		seqs.append(entry.seq)
		n_count = entry.seq.count('N') + entry.seq.count('n')
		gap_count = entry.seq.count('-')
		print('\t'.join([fasta, entry.id, str(len(entry.seq)), '%.2f' % GC(entry.seq),
			'%.2f' % (float(100*n_count)/len(entry.seq)),
			'%.2f' % (float(100*gap_count)/len(entry.seq))]))

	cum_seq = Seq('').join(seqs)
	cum_seqs.append(cum_seq)
	n_count = cum_seq.count('N') + cum_seq.count('n')
	gap_count = cum_seq.count('-')
	if len(seqs) > 1:
		print('\t'.join([fasta, 'Total (' + str(len(seqs)) + ')', str(len(cum_seq)),
			'%.2f' % GC(cum_seq), '%.2f' % (float(100*n_count)/len(cum_seq)),
			'%.2f' % (float(100*gap_count)/len(cum_seq))]))
	index = index + 1

overall_seq = Seq('').join(cum_seqs)
n_count = overall_seq.count('N') + overall_seq.count('n')
gap_count = overall_seq.count('-')
if len(cum_seqs) > 1:
	print('')
	print('\t'.join(['Summary', 'Length', 'GC', '% N', '% Gaps']))
	print('\t'.join([str(index - 1) + ' files', str(len(overall_seq)),
		'%.2f' % GC(overall_seq), '%.2f' % (float(100*n_count)/len(overall_seq)),
		'%.2f' % (float(100*gap_count)/len(overall_seq))]))
