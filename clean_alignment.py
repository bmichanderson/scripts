#!/usr/bin/env python

###############################################
# Author: B. Anderson
# Date: 12 Mar 2019
# Modified: Jun 2021 (updated for phylip and using AlignIO), Oct 2020
# Description: read in an alignment and remove positions with more than a specified number or percentage of gaps or Ns
###############################################


import sys
from Bio import AlignIO


def help():
	print('A script to remove missing positions in an alignment.')
	print('The cleaned alignment is output as a fasta/phylip file in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' options(-... -...) alignment')
	print('')
	print('Options:')
	print('	-f	File type (fasta [default] or phylip)')
	print('')
	print('	-g	Maximum number of gaps/Ns in the alignment at any position')
	print('')
	print('	-p	Maximum percentage (%) of gap/N positions in the alignment at any position')
	print('')


opts = ['-f', '-g', '-p']
opt_valdic = {}


# print help if the script is called without args

if len(sys.argv[1:]) == 0:
	sys.exit(help())


# assign arguments to variables and evaluate options

alignment = sys.argv[-1:][0]		# retrieve the last argument slice, then entry as string

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


if '-f' in opt_valdic:
	align_specified = True
	align_type = opt_valdic['-f']
else:
	align_specified = False
	align_type = 'fasta'


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

if align_type == 'fasta':

	align = AlignIO.read(open(alignment, 'r'), 'fasta')

elif align_type == 'phylip':

	align = AlignIO.read(open(alignment, 'r'), 'phylip-relaxed')

else:
	sys.exit('Please specify alignment type as "fasta" or "phylip"')


# now that the data has been read, we need to identify positions with too many gaps or Ns

keep_cols = []
counter = 1
align_length = len(align)		# this is number of rows; columns is get_alignment_length()

for index in range(0, align.get_alignment_length()):
	gap_count = align[:, index].count('-')
	n_count = align[:, index].upper().count('N')

	if gap_number:
		if gap_count + n_count > num_gaps:
			continue
		else:
			keep_cols.append(index)

	elif gap_percent:
		gap_count_perc = round(100*float(gap_count + n_count)/align_length, 0)
		if gap_count_perc > perc_gaps:
			continue
		else:
			keep_cols.append(index)

	else:
		sys.exit('Specify gap/N threshold')

	# insert a counter to enable tracking of large alignment files
	if index // 100000 == counter:
		print(str(index))
		counter = counter + 1


# slice the alignments to keep the correct columns

print('Keeping ' + str(len(keep_cols)) + ' columns')

new_alignment = align[:, keep_cols[0]: keep_cols[0] + 1]	# make a single column alignment

index = 1
counter = 1
for column in keep_cols[1:]:
	new_alignment = new_alignment + align[:, column: column + 1]		# add single column alignments
	index = index + 1
	if index // 10000 == counter:
		print(str(index))
		counter = counter + 1


# output the new alignment

if align_type == 'fasta':

	with open(alignment.replace('.f', '_clean.f'), 'w') as out_file:
		AlignIO.write(new_alignment, out_file, 'fasta')

elif align_type == 'phylip':

	with open(alignment.replace('.p', '_clean.p'), 'w') as out_file:
		AlignIO.write(new_alignment, out_file, 'phylip-relaxed')
