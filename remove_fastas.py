#!/usr/bin/env python
###############################################
# Author: B. Anderson
# Date: 28 Jun 2019
# Description: a script to simply remove entries from multifastas based on a string matches (arg1)
###############################################

import sys
from Bio import SeqIO

def help():
	print('A script to remove fasta entries based on string matching.')
	print('The new fasta files are output in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' string1,string2,... multifasta1 multifasta2 ...')
	print('')

# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# parse the command line
strings = sys.argv[1]
string_list = []
for part in strings.split(','):
	string_list.append(part)

file_list = []
for arg in sys.argv[2:]:
	file_list.append(arg)


# read in the each multifasta file and remove entries
for multifasta in file_list:

	fasta_list = []

	fastas = SeqIO.parse(open(multifasta, 'r'), 'fasta')

	for fasta in fastas:
		fasta_list.append(fasta)

	remove_indices = []

	for i,fasta in enumerate(fasta_list):
		for string in string_list:
			if string in fasta.description:
				remove_indices.append(i)

	for i in sorted(remove_indices, reverse = True):
		del fasta_list[i]

	print('Removed ' + str(len(remove_indices)) + ' entries from ' + multifasta)

	with open('mod_' + multifasta, 'w') as outfile:
		for fasta in fasta_list:
			SeqIO.write(fasta, outfile, 'fasta')
