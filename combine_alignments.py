#!/usr/bin/env python

################
# Author: B. Anderson
# Date: 15 Aug 2020
# Modified: 30 Jan 2021
# Description: combine fasta alignments for a single concatenated version with all taxa
# Change: added a parser and made the script able to combine alignments with limited header information, e.g. Genus_species only
################


import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to concatenate multiple sequence alignments, inserting dashes for missing taxa')


# add arguments to parse
parser.add_argument('alignments', type=str, help='The multiple sequence alignments', nargs='*')
parser.add_argument('-f', type=str, dest='format',
			help='Specify name format: original [default] = space delimited, with taxa as fourth and fifth field; ' +
			'simple = underscore delimited, with taxa as first and second field')


# parse the command line
if len(sys.argv[1:]) == 0:              # if there are no arguments
        parser.print_help(sys.stderr)
        sys.exit(1)

args = parser.parse_args()

alignments = args.alignments		# a list
format = args.format
if format:
	if format.lower() == 'simple':
		format_original = False
	else:
		format_original = True
else:
	format_original = True


files = []
for file in alignments:
	files.append(file)


# read each file and capture data
entry_list = []
align_num = 1
for multifasta in files:
	with open(multifasta, 'r') as align:
		fastas = SeqIO.parse(align, 'fasta')
		for fasta in fastas:
			entry_list.append([align_num, fasta])
	align_num = align_num + 1

print('Read in ' + str(align_num - 1) + ' alignments')


# determine all the taxa present
taxa_list = []
for entry in entry_list:
	if format_original:
		head = entry[1].description.split()
		acc = head[2]
		if (head[3] + ' ' + head[4]) in (x[0] for x in taxa_list):
			continue
		else:
			taxa_list.append([head[3] + ' ' + head[4], acc])
	else:
		head = entry[1].description.strip()
		if head in taxa_list:
			continue
		else:
			taxa_list.append(head)

taxa_list.sort()
print('Detected ' + str(len(taxa_list)) + ' taxa present across alignments')


# for each region, create fasta files for missing taxa that are "-" and as long as the others
for num in range(1, align_num):
	taxa_minilist = []
	length = 0
	count = 0
	for entry in entry_list:
		if entry[0] == num:
			if length == 0:
				length = len(entry[1])
			if len(entry[1]) > length:
				sys.exit('Alignment entries are not equal length!')
			if format_original:
				head = entry[1].description.split()
				taxon = head[3] + ' ' + head[4]
			else:
				taxon = entry[1].description.strip()
			if taxon in taxa_minilist:
				print('Taxon represented more than once for alignment ' + str(num))
				continue
			else:
				taxa_minilist.append(taxon)
	for tentry in taxa_list:
		if format_original:
			if tentry[0] in taxa_minilist:
				continue
			else:
				new_fasta = SeqRecord(Seq('-' * length), id = 'gaps', name = 'gaps', description = 'gaps from ' + tentry[1] + ' ' + tentry[0])
				entry_list.append([num, new_fasta])
				count = count + 1
		else:
			if tentry in taxa_minilist:
				continue
			else:
				new_fasta = SeqRecord(Seq('-' * length), id = tentry, name = tentry, description = tentry)
				entry_list.append([num, new_fasta])
				count = count + 1

	print('Added ' + str(count) + ' gap entries for missing taxa for alignment ' + str(num))


# sort the entry list and start concatenating
entry_list = sorted(entry_list, key=lambda x: x[0])
out_list = []
if format_original:
	for tentry in taxa_list:
		new_fasta = SeqRecord(Seq(""), id = 'concat', name = 'concat', description = 'concat from ' + tentry[1] + ' ' + tentry[0])
		for entry in entry_list:
			head = entry[1].description.split()
			if (head[3] + ' ' + head[4]) == tentry[0]:
				new_fasta.seq = new_fasta.seq + entry[1].seq
		out_list.append([tentry[0], new_fasta])
else:
	for tentry in taxa_list:
		new_fasta = SeqRecord(Seq(""), id = tentry, name = tentry, description = tentry + ' concat')
		for entry in entry_list:
			head = entry[1].description.strip()
			if head == tentry:
				new_fasta.seq = new_fasta.seq + entry[1].seq
		out_list.append([tentry, new_fasta])


# write the final concatenated fastas into a single file
with open('combine_out.fasta', 'w') as outfile:
	for outf in sorted(out_list, key=lambda x: x[0]):
		SeqIO.write(outf[1], outfile, 'fasta')

