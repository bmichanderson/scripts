#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 8-9 Jan 2019
# Modified: Oct 2020
# Description: generate a consensus sequence from a fasta file containing a multiple sequence alignment (e.g. from MAFFT) given a threshold
##########################

import sys				# allows access to command line arguments
from Bio import SeqIO			# for reading and writing sequence records
from Bio import AlignIO			# for reading and writing alignment records
from Bio.Align import AlignInfo		# for calculating summary information from alignments
from Bio import Alphabet		# for generating alphabet types
from Bio.Alphabet import IUPAC		# for interpreting DNA sequence characters
from Bio.SeqRecord import SeqRecord	# for generating a sequence record from a sequence


# a function for when the script is called incorrectly or without arguments
def help():
        print('A script to generate a consensus sequence from a fasta file (arg1) using a threshold (arg2)')
        print('')
        print('Usage: ' + str(sys.argv[0]) + ' fasta_file threshold')
        print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
        print('Threshold: ' + float(sys.argv[2]))
else:
        sys.exit(help())


# import the first command argument as the alignment file and parse, and write the output
thresh = float(sys.argv[2])			# threshold for calculating consensus position

with open(sys.argv[1], 'r') as align_file, open('consensus.fasta', 'w') as out_file:
	align = AlignIO.read(align_file, 'fasta', alphabet = Alphabet.Gapped(IUPAC.ambiguous_dna))
	sum_align = AlignInfo.SummaryInfo(align)
	consensus = sum_align.dumb_consensus(thresh, 'n', IUPAC.ambiguous_dna)
	record = SeqRecord(consensus.upper(), id = sys.argv[1][0:sys.argv[1].find('_')], description = 'consensus using threshold of ' + str(thresh))
	SeqIO.write(record, out_file, 'fasta')
