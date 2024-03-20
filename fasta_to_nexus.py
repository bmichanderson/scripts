#!/usr/bin/env python3

###############################################
# Author: B.M. Anderson
# Date: 20 Mar 2024
# Description: read in a fasta alignment (*.fasta) and convert it to NEXUS for use with MrBayes
###############################################


import sys
import os
from Bio import SeqIO


def help():
	print('A script to convert a fasta alignment (named \'*.fasta\') to a NEXUS alignment')
	print('The converted alignment is output as a \'*.nex\' file in the current directory')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' fasta_alignment')
	print('')


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())
else:
	fasta_file = sys.argv[1]


# read in the alignment
fasta_list = []
fastas = SeqIO.parse(open(fasta_file, 'r'), 'fasta')
for fasta in fastas:
	fasta_list.append(fasta)


# construct the data block information
data_block = ('Begin data;\n\tDimensions ntax=' + str(len(fasta_list)) +
	' nchar=' + str(len(fasta_list[0].seq)) + ';\n\t' +
	'Format datatype=DNA missing=-;\n\tMatrix\n\n')


# write the output
with open(os.path.basename(fasta_file).replace('.fasta', '.nex'), 'w') as outfile:
	outfile.write('#NEXUS\n\n')
	outfile.write(data_block)
	for fasta in fasta_list:
		outfile.write(fasta.id + '\t' + str(fasta.seq) + '\n')
	outfile.write('\t;\nEnd;\n')
