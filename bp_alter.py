#!/usr/bin/env python
##########################
# Author: B. Anderson
# Date: 20 Jul 2020
# Description: modify base pairs in contigs in a (mulit)fasta file given an input text file specifying locations and changes
##########################

import sys				# allows access to command line arguments
import argparse
from Bio import SeqIO			# for reading and writing sequence records
from Bio.Seq import Seq			# for creating sequence objects
from Bio.SeqRecord import SeqRecord	# for generating a sequence record from a sequence


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to modify base pairs in contigs in a (multi)fasta file using a user-generated text file specifying changes')


# add arguments to parse
parser.add_argument('-f', type=str, dest='fasta', help='The (multi)fasta file with contigs')
parser.add_argument('-p', type=str, dest='positions', help='The user-generated text file with positions to change')
parser.add_argument('-a', type=str, dest='ambi', help='Substitute an ambiguity code instead of the alternative base pair (yes or no [default])')
parser.add_argument('-b', type=str, dest='both', help='Generate both possible contigs as copy1 and copy2 (yes or no [default])')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	print('Format of the tab-delimited positions file:')
	print('contig_id	position(1-based)	existing_base/cp1	alt_base(- to delete)/cp2')
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
	print('Please specify a fasta file of contigs')
	sys.exit(0)

if not positions_file:
	parser.print_help(sys.stderr)
	print('Please specify a positions file')
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
amb_dict = {'AC':'M', 'AG':'R', 'AT':'W', 'CG':'S', 'CT':'Y', 'GT':'K'}


# Read the contig file and store the fastas
fasta_list = []
fastas = SeqIO.parse(open(fasta_file, 'r'), 'fasta')
for f in fastas:
	fasta_list.append(f)


# Read the positions file and store a list of changes
change_list = []
with open(positions_file, 'r') as pos_file:
	for line in pos_file:
		elements = line.strip().split('\t')
		con_id = elements[0]
		position = elements[1]
		orig_bp = elements[2]
		alt_bp = elements[3]

		if use_ambi:
			new_bp = amb_dict[''.join(sorted([orig_bp, alt_bp]))]
			change_list.append([con_id, position, orig_bp, new_bp])
		else:
			change_list.append([con_id, position, orig_bp, alt_bp])


# Collect the changes for each contig together
changes_list = []
unchanged_list = []
for f in fasta_list:
	gather_list = []
	for change in change_list:
		if change[0] == f.id:
			gather_list.append(change)
	if gather_list:
		changes_list.append([f, gather_list])
	else:
		unchanged_list.append(f)


# Alter each contig and output
with open('out_mod.fasta', 'w') as outfile:
	for combo in changes_list:
		f_chars = []
		gaps = []
		for char in str(combo[0].seq):
			f_chars.append(char)

		if do_both:
			f2_chars = f_chars.copy()
			for change in combo[1]:
				p = int(change[1])
				bp1 = change[2]
				bp2 = change[3]
				f_chars[p-1] = bp1
				f2_chars[p-1] = bp2

			SeqIO.write(SeqRecord(Seq(''.join(f_chars)), id=combo[0].id, description=combo[0].description + '_mod'), outfile, 'fasta')
			SeqIO.write(SeqRecord(Seq(''.join(f2_chars)), id=combo[0].id, description=combo[0].description + '_mod2'), outfile, 'fasta')

		else:
			for change in combo[1]:
				p = int(change[1])
				bp1 = change[2]
				bp2 = change[3]
				if use_ambi:
					f_chars[p-1] = bp2
				elif f_chars[p-1] == bp1:
					if bp2 == '-':
						gaps.append(p)
					else:
						f_chars[p-1] = bp2
				else:
					print('Warning: incorrectly specified existing bp for ' + combo[0].id + '!')

			if gaps:
				for index in sorted(gaps, reverse=True):
					del f_chars[index-1]

			SeqIO.write(SeqRecord(Seq(''.join(f_chars)), id=combo[0].id, description=combo[0].description + '_mod'), outfile, 'fasta')

	for f in unchanged_list:
		SeqIO.write(f, outfile, 'fasta')
