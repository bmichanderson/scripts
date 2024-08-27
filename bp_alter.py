#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: 20 Jul 2020
# Version update: Aug 2024
# Description: modify sequences in a (mulit)fasta file given an input text file specifying
#	contigs, locations and changes
#	The modified sequences are output in the current directory as a file "output_alter.fasta"
##########################


import argparse
from Bio import SeqIO					# for reading and writing sequence records
from Bio.SeqRecord import SeqRecord		# for generating a sequence record from a sequence
import sys


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to modify positions in sequences in a (multi)fasta file ' +
	'using a user-generated text file specifying changes. ' +
	'The text file should be in the form (tab-separated columns without a header): '
	'contig_id	position(1-based)	existing_base	alt_base(use \"-\" to delete it). '
	'If wanting to make an insertion, specify two bases, e.g. C -> CT would be: contig position C CT')


# add arguments to parse
parser.add_argument('-f', type=str, dest='fasta', help='The (multi)fasta file with contigs')
parser.add_argument('-p', type=str, dest='positions', help='The user-generated text file with positions to change')
parser.add_argument('-a', type=str, dest='ambi', help='Substitute ambiguity codes instead of the alternative base pairs (yes or no [default])')
parser.add_argument('-b', type=str, dest='both', help='Generate both possible contigs as copy1 and copy2 (yes or no [default])')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
fasta_file = args.fasta
positions_file = args.positions
ambiguities = args.ambi
both = args.both

use_ambi = False
do_both = False

if not fasta_file:
	parser.print_help(sys.stderr)
	sys.exit(0)
if not positions_file:
	parser.print_help(sys.stderr)
	sys.exit(0)

if ambiguities:
	if ambiguities.lower() == 'yes':
		use_ambi = True
	else:
		use_ambi = False

if both:
	if both.lower() == 'yes':
		do_both = True
	else:
		do_both = False


# Define an ambiguity dictionary
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


# Read the contig file and store the fastas
fasta_list = []
fastas = SeqIO.parse(open(fasta_file, 'r'), 'fasta')
for fasta in fastas:
	fasta_list.append(fasta)


# Read the positions file and store as a list
change_list = []
with open(positions_file, 'r') as pos_file:
	for line in pos_file:
		elements = line.strip().split('\t')
		con_id = elements[0]
		position = elements[1]
		orig_bp = elements[2]
		alt_bp = elements[3]
		if use_ambi:
			new_bp = amb_dict[''.join(sorted(list(orig_bp) + list(alt_bp)))]
		else:
			new_bp = alt_bp
		change_list.append([con_id, position, orig_bp, new_bp])


# Open an output file and cycle through the stored fastas, making specified changes
with open('output_alter.fasta', 'w') as outfile:
	nochanges = 0
	changed = 0
	for fasta in fasta_list:
		these_changes = []
		old_fasta = SeqRecord(fasta.seq, id = fasta.id, name = fasta.name, description = fasta.description)
		for change in change_list:
			if change[0] == fasta.id:
				these_changes.append(change)

		if len(these_changes) < 1:		# no changes for that contig
			nochanges = nochanges + 1
			SeqIO.write(fasta, outfile, 'fasta')
			continue
		else:
			changed = changed + 1

		if do_both:
			old_fasta.id = fasta.id + '_copy1'
			old_fasta.description = fasta.description.replace(fasta.id, fasta.id + '_copy1')
			SeqIO.write(old_fasta, outfile, 'fasta')

		# sort the changes in reverse so that deletions won't affect indexing
		sorted_changes = sorted(these_changes, key = lambda x: int(x[1]), reverse = True)
		for change in sorted_changes:
			position = int(change[1])
			old_bp = change[2]
			new_bp = change[3]
			if str(fasta.seq[position - 1]) == old_bp:		# the base matches user input
				if new_bp == '-':		# a deletion
					fasta.seq = fasta.seq[0: position - 1] + fasta.seq[position: ]
				else:
					fasta.seq = fasta.seq[0: position - 1] + new_bp + fasta.seq[position: ]
			else:
				print('Warning: incorrectly specified existing base for ' + fasta.id + '!')
				print('Exiting')
				sys.exit(1)

		# output the new fasta
		if do_both:
			old_id = fasta.id
			fasta.id = fasta.id + '_copy2'
			fasta.description = fasta.description.replace(old_id, fasta.id)

		SeqIO.write(fasta, outfile, 'fasta')


# Report completion
print('Finished altering sequence(s): ' + str(changed) + ' altered, ' + str(nochanges) + ' unaltered')
