#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 10 Mar 2020
# Modified: Oct 2020
# Description: combine fasta files, removing duplicate sequences
##########################

import sys
from Bio import SeqIO		# for reading and writing sequence records

def help():
	print('A script to combine fasta entries in multiple files into a single file.')
	print('The fastas are output as a single multifasta file called combined.fasta in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' file1 file2 ...')
	print('')


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# capture list of files
file_list = []
for arg in sys.argv[1:]:
	file_list.append(arg)


# extract the fastas from the files
fastas = []
for f_file in file_list:
	with open(f_file, 'r') as fasta_file:
		fasta_iter = SeqIO.parse(fasta_file, 'fasta')
		for entry in fasta_iter:
			fastas.append(entry)


# remove duplicate entries
dlist = []
seqlist = []
to_remove = []

for num, fasta in enumerate(fastas):
	if fasta.description in dlist:
		d_indices = [i for i in range(len(dlist)) if dlist[i] == fasta.description]
		found = 'False'
		for index in d_indices:
			if fasta.seq == seqlist[index]:
				to_remove.append(num)
				found = 'True'
				break

		if found == 'False':
			dlist.append(fasta.description)
			seqlist.append(fasta.seq)

	else:
		dlist.append(fasta.description)
		seqlist.append(fasta.seq)

for index in sorted(to_remove, reverse = True):
	del fastas[index]


# sort the fastas, then output to a new file
#fastas = sorted(fastas, key=lambda x: len(x.seq), reverse=True)

with open('combined.fasta', 'w') as output_file:
	for entry in fastas:
		SeqIO.write(entry, output_file, 'fasta')
