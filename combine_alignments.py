#!/usr/bin/env python

################
# Author: B. Anderson
# Date: 15 Aug 2002
# Description: combine fasta alignments for a single concatenated version with all taxa
################

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# check that there are arguments
if len(sys.argv[1:]) == 0:
	sys.exit('Please indicate fasta alignments to combine')


# parse the command line for the input fasta alignments
files = []
for file in sys.argv[1:]:
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
	head = entry[1].description.split()
	acc = head[2]
	if (head[3] + ' ' + head[4]) in (x[0] for x in taxa_list):
		continue
	else:
		taxa_list.append([head[3] + ' ' + head[4], acc])

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
			head = entry[1].description.split()
			if (head[3] + ' ' + head[4]) in taxa_minilist:
				print('Taxon represented more than once for alignment ' + str(num))
				continue
			else:
				taxa_minilist.append(head[3] + ' ' + head[4])
	for tentry in taxa_list:
		if tentry[0] in taxa_minilist:
			continue
		else:
			new_fasta = SeqRecord(Seq('-' * length), id = 'gaps', name = 'gaps', description = 'gaps from ' + tentry[1] + ' ' + tentry[0])
			entry_list.append([num, new_fasta])
			count = count + 1

	print('Added ' + str(count) + ' gap entries for missing taxa for alignment ' + str(num))


# sort the entry list and start concatenating
entry_list = sorted(entry_list, key=lambda x: x[0])
out_list = []
for tentry in taxa_list:
	new_fasta = SeqRecord(Seq(""), id = 'concat', name = 'concat', description = 'concat from ' + tentry[1] + ' ' + tentry[0])
	for entry in entry_list:
		head = entry[1].description.split()
		if (head[3] + ' ' + head[4]) == tentry[0]:
			new_fasta.seq = new_fasta.seq + entry[1].seq
	out_list.append([tentry[0], new_fasta])


# write the final concatenated fastas into a single file
with open('combine_out.fasta', 'w') as outfile:
	for outf in sorted(out_list, key=lambda x: x[0]):
		SeqIO.write(outf[1], outfile, 'fasta')

