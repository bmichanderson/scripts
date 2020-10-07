#!/usr/bin/env python

###############################################
# Author: B. Anderson
# Date: 12 Mar 2019
# Modified: Oct 2020
# Description: read in an alignment and remove positions with more than a specified number or percentage of gaps
###############################################

import sys
from Bio import SeqIO

def help():
	print('A script to remove missing positions in an alignment.')
	print('The cleaned alignment is output as a fasta file in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' options(-... -...) fasta_alignment')
	print('')
	print('Options:')
	print('	-g	Maximum number of gaps in the alignment at any position')
	print('')
	print('	-p	Maximum percentage (%) of gap positions in the alignment at any position')

opts = ['-g', '-p']
opt_valdic = {}


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# assign arguments to variables and evaluate options
fasta_alignment = sys.argv[-1:][0]	# retrieve the last argument slice, then entry as string

for index, arg in enumerate(sys.argv[1:-1]):	# don't include last entry
	if arg[0] == '-':
		if arg in opts:
			if index < len(sys.argv[1:-1]) - 1:
				if sys.argv[1:-1][index + 1][0] == '-':
					sys.exit('No argument provided for ' + str(arg))
				opt_valdic[arg] = sys.argv[1:-1][index + 1]
			else:
				sys.exit(help())
		else:
			sys.exit('The argument ' + str(arg) + ' is not a valid option')
	else:
		continue

if '-g' in opt_valdic:
	gap_number = True
	num_gaps = int(opt_valdic['-g'])
else:
	gap_number = False

if '-p' in opt_valdic:
	gap_percent = True
	perc_gaps = int(opt_valdic['-p'])
else:
	gap_percent = False


# read in the alignment
fasta_list = []
seq_list = []

fastas = SeqIO.parse(open(fasta_alignment, 'r'), 'fasta')

for fasta in fastas:
	fasta_list.append(fasta)
	seq_list.append(str(fasta.seq))


# now that the data has been read, we need to identify positions with too many gaps
del_cols = []
index = 0
for align_col in zip(*seq_list):		# e.g. [[a, b, c], [a, b, c]] ---> [(a, a), (b, b) ...
	gap_count = align_col.count('-')

	if gap_number:
		if gap_count > num_gaps:
			del_cols.append(index)
			index += 1
		else:
			index += 1

	elif gap_percent:
		gap_count_perc = round(100*float(gap_count)/len(align_col), 0)
		if gap_count_perc > perc_gaps:
			del_cols.append(index)
			index += 1
		else:
			index += 1

	else:
		sys.exit('Specify gap threshold')


# delete those positions in all the fastas
for fasta in fasta_list:
	new_seq = fasta.seq
	for col in sorted(del_cols, reverse=True):		# have to start from the highest index!
		new_seq = new_seq[:col] + new_seq[col+1:]
	fasta.seq = new_seq


# output the new alignment
with open(fasta_alignment.replace('.fasta', '_clean.fasta'), 'w') as out_file:
	for fasta in fasta_list:
		SeqIO.write(fasta, out_file, 'fasta')
