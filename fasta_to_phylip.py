#!/usr/bin/env python

###############################################
# Author: B. Anderson
# Date: 27 Jun 2019
# Modified: Oct 2020, Jun 2021 (to make more flexible)
# Description: read in a fasta alignment (with a particular naming convention) and convert it to phylip for running in RAxML
###############################################

import sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def help():
	print('A script to convert a fasta alignment to a phylip alignment.')
	print('The converted alignment is output as a phylip file in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' fasta_alignment')
	print('')
	print('Note that the fasta naming convention should follow:')
	print('>region from ID Genus specific_ep ...')


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# read in the alignment
fasta_list = []

fastas = SeqIO.parse(open(sys.argv[1], 'r'), 'fasta')

for fasta in fastas:
	fasta_list.append(fasta)


# change the names to relaxed phylip format
label_list = []
for fasta in fasta_list:
	if len(fasta.description.split()) > 4:
		label = fasta.description.split()[3] + '_' + fasta.description.split()[4]
	else:
		label = fasta.description.split()[0]

	index = 2
	orig_label = label
	while label in label_list:
		label = orig_label + str(index)
		index += 1

	label_list.append(label)

	fasta.id = label


# output the new alignment
with open(sys.argv[1].replace('.fasta', '.phy'), 'w') as out_file:
	alignment = MultipleSeqAlignment(fasta_list)
	AlignIO.write(alignment, out_file, 'phylip-relaxed')
