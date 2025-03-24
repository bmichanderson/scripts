#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Jan 2019
# Modified: Oct 2020, Jun 2023 (major), Jul 2023, Mar 2025 (score Ns less frequently if base data above threshold)
# Description: generate a simple consensus sequence from a fasta multiple sequence alignment given a threshold
##########################


import sys
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to generate a consensus sequence (\"consensus.fasta\") ' +
    'from a fasta multiple sequence alignment using a threshold for counting bases')


# add arguments to parse
parser.add_argument('fasta', type = str, help = 'The multiple sequence alignment fasta file')
parser.add_argument('-t', type = str, dest = 'thresh', help = 'The minimum (>) frequency threshold ' +
	' (of the total non-missing characters) for counting a base at a position [default 0.2]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
fasta_file = args.fasta
threshold = args.thresh

if threshold:
	threshold = float(threshold)
else:
	threshold = 0.2


# create an ambiguity dictionary
amb_dict = {
	'AC': 'M',
	'AG': 'R',
	'AT': 'W',
	'CG': 'S',
	'CT': 'Y',
	'GT': 'K',
	'ACG': 'V',
	'ACT': 'H',
	'AGT': 'D',
	'CGT': 'B',
	'ACGT': 'N'
}


# read in the alignment and output the consensus
with open(fasta_file, 'r') as align_file, open('consensus.fasta', 'w') as out_file:
	alignment = AlignIO.read(align_file, 'fasta')
	consensus = []

	# convert to upper case
	for entry in alignment:
		entry.seq = entry.seq.upper()

	# for each position, use the bases that are present at > threshold to determine the consensus
	gaps = 0
	for pos in range(alignment.get_alignment_length()):
		align_slice = alignment[:, pos]
		nbases = [nbase for nbase in align_slice if nbase in ['-', 'N']]
		bases = [base for base in align_slice if base not in ['-', 'N']]
		base_list = []

		if len(bases) != 0:		# a residue is there
			if float(bases.count('A') / len(bases)) > threshold:
				base_list.append('A')

			if float(bases.count('C') / len(bases)) > threshold:
				base_list.append('C')

			if float(bases.count('G') / len(bases)) > threshold:
				base_list.append('G')

			if float(bases.count('T') / len(bases)) > threshold:
				base_list.append('T')

			if  any([len(base_list) == 0, len(bases) / (len(nbases) + len(bases)) < threshold]):
				# no bases above threshold or present data compared to total data less than threshold
				base = 'N'
			elif len(base_list) == 1:	# a single base above threshold
				base = base_list[0]
			else:
				base = amb_dict[''.join(base_list)]

			consensus.append(base)

		elif 'N' in nbases:
			consensus.append('N')

		else:
			gaps = gaps + 1

	# output the consensus sequence
	new_record = SeqRecord(Seq(''.join(consensus)), id = 'consensus', name = 'consensus', description = 'consensus')
	SeqIO.write(new_record, out_file, 'fasta')

	# report if there were gaps removed
	if gaps > 0:
		print('Removed ' + str(gaps) + ' gap positions')
