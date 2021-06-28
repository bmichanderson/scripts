#!/usr/bin/env python

###############################################
# Author: B. Anderson
# Date: 12 Mar 2019
# Modified: Jun 2021 (updated for phylip and using AlignIO, numpy and pandas), Oct 2020
# Description: read in an alignment and remove positions with more than a specified number or percentage of gaps or Ns
###############################################


import sys
import pandas
import numpy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


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


# capture the sequence ids

ids = []

for sequence in align:
	ids.append(sequence.id)

#print('Sequence ids are: ')
#for id in ids:
#	print(id)

len_align = align.get_alignment_length()
print('Alignment length is: ' + str(len_align))


# now that the data has been read, we need to convert the sequences to an array, then to a dataframe
# due to memory issues, I need to do this in chunks when the file is large

# define a function to convert to upper case
def upperit(x):
	return(str(x).upper())


# if the alignment is super large
if len_align > 500000:
	align_list = []
	increment = 200000

	for start in range(0, len_align, increment):
		end = min(len_align, start + increment)

		print('Filtering from positions: ' + str(start) + ' to ' + str(end) + ' of ' + str(len_align))

		n_array = numpy.array([list(sample) for sample in align[:, start: end]])

		p_df = pandas.DataFrame(n_array)
		rows = float(len(p_df))

		# convert to upper case for search

		p_df = p_df.applymap(upperit)

		# now we need to filter the dataframe by columns, to keep columns which pass filters
		# from https://stackoverflow.com/questions/31614804/how-to-delete-a-column-in-pandas-dataframe-based-on-a-condition/31618820

		if gap_number:
			filtp_df = p_df.loc[:, (p_df.eq('-').sum() + p_df.eq('N').sum() <= num_gaps)]

		elif gap_percent:
			filtp_df = p_df.loc[:, (round(100*(p_df.eq('-').sum() + p_df.eq('N').sum())/rows, 0) <= perc_gaps)]

		else:
			sys.exit('Specify gap/N threshold')

		# now that the dataframe is filtered, we need to convert it back into an alignment

		new_records = []

		for index, sequence in enumerate(filtp_df.values.tolist()):
			new_sequence = Seq(''.join(sequence))
			new_record = SeqRecord(new_sequence, id = ids[index])
			new_records.append(new_record)

		# append it to a list
		align_list.append(MultipleSeqAlignment(new_records))


	# stick the alignments together and output the new alignment

	new_alignment = align_list[0]
	for piece in align_list[1:]:
		new_alignment = new_alignment + piece

	print('Filtered alignment is ' + str(new_alignment.get_alignment_length()) + ' bp')

	if align_type == 'fasta':

		with open(alignment.replace('.f', '_clean.f'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'fasta')

	elif align_type == 'phylip':

		with open(alignment.replace('.p', '_clean.p'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'phylip-relaxed')


else:
	n_array = numpy.array([list(sample) for sample in align])
	print('Dimensions of array: ' + str(n_array.shape))

	p_df = pandas.DataFrame(n_array)
	rows = float(len(p_df))

	# convert to upper case for search

	p_df = p_df.applymap(upperit)


	# now we need to filter the dataframe by columns, to keep columns which pass filters
	# from https://stackoverflow.com/questions/31614804/how-to-delete-a-column-in-pandas-dataframe-based-on-a-condition/31618820

	if gap_number:
		filtp_df = p_df.loc[:, (p_df.eq('-').sum() + p_df.eq('N').sum() <= num_gaps)]

	elif gap_percent:
		filtp_df = p_df.loc[:, (round(100*(p_df.eq('-').sum() + p_df.eq('N').sum())/rows, 0) <= perc_gaps)]

	else:
		sys.exit('Specify gap/N threshold')

	#print('New data frame has ' + str(rows) + ' rows and ' + len(filtp_df.columns) + ' columns')


	# now that the dataframe is filtered, we need to convert it back into an alignment

	new_records = []

	for index, sequence in enumerate(filtp_df.values.tolist()):
		new_sequence = Seq(''.join(sequence))
		new_record = SeqRecord(new_sequence, id = ids[index])
		new_records.append(new_record)

	new_alignment = MultipleSeqAlignment(new_records)

	print('Filtered alignment is ' + str(new_alignment.get_alignment_length()) + ' bp')

	# output the new alignment

	if align_type == 'fasta':

		with open(alignment.replace('.f', '_clean.f'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'fasta')

	elif align_type == 'phylip':

		with open(alignment.replace('.p', '_clean.p'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'phylip-relaxed')
