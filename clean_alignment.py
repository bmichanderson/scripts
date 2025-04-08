#!/usr/bin/env python3

###############################################
# Author: B.M. Anderson
# Date: 12 Mar 2019
# Modified: Oct 2020; Jun 2021 (updated for phylip and using AlignIO, numpy and pandas);
#	Mar 2023 (add sequence removal option); Sep 2024 (made default no filtering; updated for pandas function deprecation);
#	Apr 2025 (added better arg parsing; updated for pandas version differences; changed output reporting slightly)
# Description: read in an alignment and remove positions/sequences with more than a specified number or percentage of gaps or Ns
###############################################


import argparse
import sys
import pandas
import numpy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


# define a function to convert to upper case
# updated: x is now a dataframe
def upperit(x):
	if int(''.join(pandas.__version__.split('.')[0:2])) > 20:
		return(x.map(lambda y: str(y).upper()))
	else:
		return(x.applymap(lambda y: str(y).upper()))


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to remove positions/sequences in an alignment ' +
	'based on number/percentage of gaps and Ns. Choosing both a number and a percentage overrides the percentage. ' +
	'The cleaned alignment is output as a fasta/phylip file in the current directory.')


# add arguments to parse
parser.add_argument('fasta', type = str, help = 'The multiple sequence alignment fasta file')
parser.add_argument('-f', type = str, dest = 'format', help = 'The format of the alignment: fasta [default] or phylip')
parser.add_argument('-g', type = int, dest = 'gaps', help = 'The maximum number of gaps and Ns allowed (<=) in ' +
	'the alignment at any position [default: no maximum]')
parser.add_argument('-p', type = float, dest = 'percent', help = 'The maximum percentage (%) of gaps and Ns allowed (<=) in ' +
	'the alignment at any position [default: 100]')
parser.add_argument('-s', type = float, dest = 'seq_percent', help = 'The maximum percentage (%) of gaps and Ns allowed (<=) in ' +
	'a sequence to keep it (after position filtering) [default: 100]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
alignment = args.fasta
format = args.format
gaps = args.gaps
percent = args.percent
seq_percent = args.seq_percent

if format:
	align_type = str(format)
else:
	align_type = 'fasta'

if gaps:
	num_gaps = gaps

if percent:
	perc_gaps = percent

if seq_percent:
	seq_perc_gaps = seq_percent


# read in the alignment
if align_type == 'fasta':
	align = AlignIO.read(open(alignment, 'r'), 'fasta')
elif align_type == 'phylip':
	align = AlignIO.read(open(alignment, 'r'), 'phylip-relaxed')
else:
	sys.exit('Please specify alignment type as "fasta" or "phylip"')


# capture the sequence ids and descriptions
ids = []
descriptions = []
for sequence in align:
	ids.append(sequence.id)
	descriptions.append(sequence.description)

len_align = align.get_alignment_length()
print('Unfiltered alignment: ' + str(len(align)) + ' sequences of length ' + str(len_align) + ' bp')


# now that the data has been read, we need to convert the sequences to an array, then to a dataframe
# due to memory issues, I need to do this in chunks when the file is large
if len_align > 500000:	# alignment is super large
	align_list = []
	increment = 200000

	for start in range(0, len_align, increment):
		end = min(len_align, start + increment)

		print('Filtering from positions: ' + str(start) + ' to ' + str(end) + ' of ' + str(len_align))

		n_array = numpy.array([list(sample) for sample in align[:, start: end]])
		p_df = pandas.DataFrame(n_array)
		rows = float(len(p_df))

		# convert to upper case for search
		p_df = upperit(p_df)

		# now we need to filter the dataframe by columns, to keep columns which pass filters
		# from https://stackoverflow.com/questions/31614804/how-to-delete-a-column-in-pandas-dataframe-based-on-a-condition/31618820
		if gaps:
			filtp_df = p_df.loc[:, (p_df.eq('-').sum() + p_df.eq('N').sum() <= num_gaps)]
		elif percent:
			filtp_df = p_df.loc[:, (round(100*(p_df.eq('-').sum() + p_df.eq('N').sum())/rows, 0) <= perc_gaps)]
		else:
			# No gap/N threshold specified, so keeping all positions
			filtp_df = p_df

		# now that the dataframe is filtered, we need to convert it back into an alignment
		new_records = []
		for index, sequence in enumerate(filtp_df.values.tolist()):
			new_sequence = Seq(''.join(sequence))
			new_record = SeqRecord(new_sequence, id = ids[index], description = descriptions[index])
			new_records.append(new_record)

		# append it to a list
		align_list.append(MultipleSeqAlignment(new_records))


	# stick the alignments together and output the new alignment
	new_alignment = align_list[0]
	for piece in align_list[1:]:
		new_alignment = new_alignment + piece

	# if there is a flag to check gap composition of sequences, do so
	if seq_percent:
		temp_record_list = []
		for align_rec in new_alignment:
			this_count = align_rec.seq.count('-') + align_rec.seq.count('N')
			if 100 * float(this_count) / len(align_rec.seq) > seq_perc_gaps:
				continue
			
			temp_record_list.append(align_rec)

		new_alignment = MultipleSeqAlignment(temp_record_list)


	print('Filtered alignment: ' + str(len(new_alignment)) + ' sequences of length ' +
		str(new_alignment.get_alignment_length()) + ' bp')

	if align_type == 'fasta':
		with open(alignment.replace('.f', '_clean.f'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'fasta')
	elif align_type == 'phylip':
		with open(alignment.replace('.p', '_clean.p'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'phylip-relaxed')

else:		# alignment is reasonable size
	n_array = numpy.array([list(sample) for sample in align])
	#print('Dimensions of array: ' + str(n_array.shape))
	p_df = pandas.DataFrame(n_array)
	rows = float(len(p_df))

	# convert to upper case for search
	p_df = upperit(p_df)

	# now we need to filter the dataframe by columns, to keep columns which pass filters
	# from https://stackoverflow.com/questions/31614804/how-to-delete-a-column-in-pandas-dataframe-based-on-a-condition/31618820
	if gaps:
		filtp_df = p_df.loc[:, (p_df.eq('-').sum() + p_df.eq('N').sum() <= num_gaps)]
	elif percent:
		filtp_df = p_df.loc[:, (round(100*(p_df.eq('-').sum() + p_df.eq('N').sum())/rows, 0) <= perc_gaps)]
	else:
		print('No gap/N threshold specified, so keeping all positions')
		filtp_df = p_df

	# now that the dataframe is filtered, we need to convert it back into an alignment
	new_records = []
	for index, sequence in enumerate(filtp_df.values.tolist()):
		new_sequence = Seq(''.join(sequence))
		# if there is a flag to check gap composition of sequences, do so
		if seq_percent:
			this_count = new_sequence.count('-') + new_sequence.count('N')
			if 100 * float(this_count) / len(new_sequence) > seq_perc_gaps:
				continue

		new_record = SeqRecord(new_sequence, id = ids[index], description = descriptions[index])
		new_records.append(new_record)

	new_alignment = MultipleSeqAlignment(new_records)
	print('Filtered alignment: ' + str(len(new_alignment)) + ' sequences of length ' +
		str(new_alignment.get_alignment_length()) + ' bp')

	# output the new alignment
	if align_type == 'fasta':
		with open(alignment.replace('.f', '_clean.f'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'fasta')
	elif align_type == 'phylip':
		with open(alignment.replace('.p', '_clean.p'), 'w') as out_file:
			AlignIO.write(new_alignment, out_file, 'phylip-relaxed')
