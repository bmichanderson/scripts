#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 11 June 2020
# Description: parse tabular (-outfmt 6) BLASTN/X results
##########################

import argparse
import sys			# allows access to command line arguments
import re			# allows use of regular expressions
from Bio import SeqIO		# for reading and writing sequence record objects



# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to parse BLASTN/X results in tabular format (-outfmt 6); use \'6 std qseq sseq qlen slen stitle\' if you wish to have alignments with gaps')


# add arguments to parse
parser.add_argument('blast_out', type=str, help='The blast output to parse')
parser.add_argument('-b', type=str, dest='plot', help='Create (yes or no [default]) a BED-like format file for plotting with e.g. R package Sushi; only for BLASTN')
parser.add_argument('-e', type=int, dest='extend', help='Extension (length of extra alignment to print; requires -q and -s) [default 50 bp]')
parser.add_argument('-l', type=int, dest='len_thresh', help='Length threshold (bp) (only HSPs longer than this will be kept) [default 0]')
parser.add_argument('-p', type=int, dest='pid_thresh', help='Percent identity, e.g. 75; only HSPs with higher identity than this will be kept [default 0]')
parser.add_argument('-q', type=str, dest='query', help='Query (multi)fasta file (needed to show longer alignments)')
parser.add_argument('-s', type=str, dest='sbjct', help='Subject (multi)fasta file (needed to show longer alignments)')
parser.add_argument('-t', type=str, dest='type', help='Type of BLAST search: x or n [default]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

blast_out = args.blast_out
create_plot = args.plot
extension = args.extend
len_thresh = args.len_thresh
pid_thresh = args.pid_thresh
query_file = args.query
sbjct_file = args.sbjct
type = args.type.lower()


# assign variables depending upon whether present or not

if create_plot:
	if create_plot.lower() == 'yes':
		create_plot = 'yes'
	else:
		create_plot = 'no'
else:
	create_plot = 'no'

if not extension:
	extension = 50

if not len_thresh:
	len_thresh = 0

if not pid_thresh:
	pid_thresh = 0

if not type:
	type = 'n'
elif any([type == 'x', type == 'n']):
	print('BLAST' + type.upper() + ' specified')
else:
	parser.print_help(sys.stderr)
	sys.exit(1)


print('Parsing blast_out: ' + str(blast_out) + ', length_threshold: ' + str(len_thresh) + ', pid_threshold: ' + str(pid_thresh)
		+ ', extension: ' + str(extension) + ', create_plot: ' + str(create_plot))


# Parse the BLAST results
hsp_list = []
with open(blast_out, 'r') as b_file:
	# check if there are any hits in the file
	if len(b_file.read(1)) == 0:
		print('No BLAST hits found for ' + blast_out)
		sys.exit()
	else:
		# filter the hits
		b_file.seek(0)
		for line in b_file:
		# actual line is 'query', 'sbjct', 'pident', 'length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'sbjct_start', 'sbjct_end', 'evalue', 'bitsore', 'qseq', 'sseq', 'qlen', 'slen', 'stitle'
			data = line.strip().split('\t')

			pid = float(data[2])
			align_length = int(data[3])

			if align_length < len_thresh:
				continue
			elif pid < pid_thresh:
				continue
			else:
				hsp_list.append(data)

	# check if there are any hits in the list
	if len(hsp_list) == 0:		# no hits passed the filters
		print('No BLAST hits passed filters for ' + blast_out)
		sys.exit()


# summarize and output the filtered blast results
print('')
print('After filtering, there are ' + str(len(hsp_list)) + ' blast high scoring segment pairs (hsps) retained')
print('')


# if a plotting file was requested, create it (only for blastn)
if all([create_plot == 'yes', type == 'n']):
	too_short = 0
	print('Creating BED-like plotting file plot_data.tab')
	plot_file = open('plot_data.tab', 'w')

	# write the header to accommodate standard and longer outfmt
	plot_file.write('\t'.join(['sbjct', 'sbjct_start', 'sbjct_end', 'query', 'query_start', 'query_end', 'score', 'sstrand', 'qlen', 'slen']) + '\n')

	# write the output for each line, after adjusting coordinates to 0-based and end-inclusive (BED)
	# actual line is 'query', 'sbjct', 'pident', 'length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'sbjct_start', 'sbjct_end', 'evalue', 'bitsore', 'qseq', 'sseq', 'qlen', 'slen', 'stitle'
	for hsp in hsp_list:
		# determine subject strand
		if int(hsp[8]) < int(hsp[9]):		# positive strand
			sstrand = '+'
			sbjct_start = str(int(hsp[8]) - 1)		# in BED format this is 0-based (1-based in BLAST)
			sbjct_end = hsp[9]				# in BED format this is 1-based (1-based in BLAST); length of feature is end - start

		else:
			sstrand = '-'
			sbjct_start = str(int(hsp[9]) - 1)		# we still need the smaller number first in BED
			sbjct_end = hsp[8]


		# determine if the hit is longer than min length
		if int(sbjct_end) - int(sbjct_start) < len_thresh:
			too_short = too_short + 1
			continue

		# determine query strand (NOT NECESSARY SINCE ALWAYS POSITIVE)
#		if int(hsp[6]) < int(hsp[7]):		# positive strand
#			qstrand = '+'
#			query_start = str(int(hsp[6]) - 1)
#			query_end = hsp[7]
#		else:
#			qstrand = '-'
#			query_start = str(int(hsp[7]) - 1)
#			query_end = hsp[6]

		# assign other variables
		sbjct = hsp[1]
		query = hsp[0]
		query_start = str(int(hsp[6]) - 1)
		query_end = hsp[7]
		score = hsp[2]			# calling the percent identical the score

		if len(hsp) > 12:		# assume qseq sseq qlen slen stitle
			qlen = hsp[14]
			slen = hsp[15]
		else:
			qlen = '?'
			slen = '?'

		# print the line to the file
		plot_file.write('\t'.join([sbjct, sbjct_start, sbjct_end, query, query_start, query_end, score, sstrand, qlen, slen]) + '\n')


	# report how many shorter hits were removed (passed previously based on gaps)
	if too_short > 0:
		print('Did not copy ' + str(too_short) + ' hits based on too short subject coverage')

	plot_file.close()


# Output alignments
with open('hits.txt', 'w') as hits_file:
	hit_list = []

	# hsp is 'query', 'sbjct', 'pident', 'length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'sbjct_start', 'sbjct_end', 'evalue', 'bitsore', 'qseq', 'sseq', 'qlen', 'slen', 'stitle'
	for hsp in hsp_list:
		if query_file:
			# locate the query fasta in the (multi)fasta file
			found_query = False
			for q_fasta in SeqIO.parse(open(query_file, 'r'), 'fasta'):
				if q_fasta.description.split()[0] == hsp[0]:
					found_query = True
					if int(hsp[6]) < int(hsp[7]):		# positive strand
						seq_str = str(q_fasta.seq)
					else:
						seq_str = str(q_fasta.seq.reverse_complement())
					query_len = len(seq_str)
					break

			# extract the sequence, adding necessary spaces or extending beyond the hit
			if found_query:
				if int(hsp[6]) < int(hsp[7]):		# positive strand
					new_start = int(hsp[6])
					new_end = int(hsp[7])
				else:
					new_start = query_len - int(hsp[6]) + 1
					new_end = query_len - int(hsp[7]) + 1



				if new_start < extension + 1:
					spaces = extension + 1 - new_start
					nuc_query = (' '*spaces + seq_str[:(new_start - 1)] + ' ' + seq_str[(new_start - 1):new_end] +
							' ' + seq_str[new_end:(new_end + extension)])
				else:
					nuc_query = (seq_str[(new_start - extension - 1):(new_start - 1)] + ' ' + seq_str[(new_start - 1):new_end] +
							' ' + seq_str[new_end:(new_end + extension)])

				# insert gaps if we have that information
				if len(hsp) > 12:		# assume qseq and sseq present
					if type == 'n':
						gap_index = [g.start() for g in re.finditer('-', hsp[12])]		# find gaps in the query sequence
						for gap_start in gap_index:
							nuc_query = nuc_query[0:(gap_start + extension)] + '-' + nuc_query[(gap_start + extension):]
					else:		# blastx
						gap_index = [g.start() for g in re.finditer('-', '  '.join(list(hsp[12])))]		# find gaps in the query sequence
						for gap_start in gap_index:
							nuc_query = nuc_query[0:(gap_start + extension + 1)] + ' '*3 + nuc_query[(gap_start + extension + 1):]


				if all([type == 'x', len(hsp) > 12]):
					prot_query = ' '*(extension + 2) + '  '.join(list(hsp[12]))
				else:
					prot_query = ' '

			else:
				sys.exit('Problem with query fasta file!')
		else:
			if len(hsp) > 12:		# assume qseq and sseq present
				nuc_query = ' '*(extension + 1) + hsp[12]
				query_len = int(hsp[14])
				if type == 'x':
					prot_query = ' '*(extension + 2) + '  '.join(list(hsp[12]))
					nuc_query = ' '
			else:
				nuc_query = ' '
				query_len = '?'
				prot_query = ' '


		if sbjct_file:
			# locate the sbjct fasta in the (multi)fasta file
			found_sbjct = False
			for s_fasta in SeqIO.parse(open(sbjct_file, 'r'), 'fasta'):
				if s_fasta.description.split()[0] == hsp[1]:
					found_sbjct = True
					stitle = s_fasta.description
					if int(hsp[8]) < int(hsp[9]):		# positive strand
						seq_str = str(s_fasta.seq)
					else:
						seq_str = str(s_fasta.seq.reverse_complement())
					sbjct_len = len(seq_str)
					break

			# extract the sequence, adding necessary spaces or extending beyond the hit
			if found_sbjct:
				if int(hsp[8]) < int(hsp[9]):		# positive strand
					new_start = int(hsp[8])
					new_end = int(hsp[9])
				else:
					new_start = sbjct_len - int(hsp[8]) + 1
					new_end = sbjct_len - int(hsp[9]) + 1

				if new_start < extension + 1:
					spaces = extension + 1 - new_start
					nuc_sbjct = (' '*spaces + seq_str[:(new_start - 1)] + ' ' + seq_str[(new_start - 1):new_end] +
							' ' + seq_str[new_end:(new_end + extension)])
				else:
					nuc_sbjct = (seq_str[(new_start - extension - 1):(new_start - 1)] + ' ' + seq_str[(new_start - 1):new_end] +
							' ' + seq_str[new_end:(new_end + extension)])

				# insert gaps if we have that information
				if len(hsp) > 12:		# assume qseq and sseq present
					gap_index = [g.start() for g in re.finditer('-', hsp[13])]		# find the gaps in the subject sequence
					for gap_start in gap_index:
						nuc_sbjct = nuc_sbjct[0:(gap_start + extension)] + '-' + nuc_sbjct[(gap_start + extension):]

				if type == 'x':
					prot_sbjct = ' '

			else:
				sys.exit('Problem with subject fasta file!')
		else:
			if len(hsp) > 12:			# assume qseq and sseq present
				nuc_sbjct = ' '*(extension + 1) + hsp[13]
				sbjct_len = int(hsp[15])
				stitle = hsp[16]
				if type == 'x':
					prot_sbjct = ' '*(extension + 2) + '  '.join(list(hsp[13]))
			else:
				prot_sbjct = ' '
				nuc_sbjct = ' '
				sbjct_len = '?'
				stitle = hsp[1]


		# append the sequences and information to the hit list depending on the type of blast search
		# hsp is 'query', 'sbjct', 'pident', 'length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'sbjct_start', 'sbjct_end', 'evalue', 'bitsore', 'qseq', 'sseq', 'qlen', 'slen', 'stitle'
		if type == 'n':
			hit_list.append((	str(hsp[3]) + ' bp, ' + str(hsp[2]) + ' % identical',
						'Query: ' + hsp[0],
						'Query loc: ' + str(hsp[6]) + '..' + str(hsp[7]) + ' of ' + str(query_len),
						'Sbjct: ' + stitle,
						'Sbjct loc: ' + str(hsp[8]) + '..' + str(hsp[9]) + ' of ' + str(sbjct_len),
						nuc_query,
						nuc_sbjct + '\n',
						))
		else:
			hit_list.append((	str(hsp[3]) + ' residues, ' + str(hsp[2]) + ' % identical',
						'Query: ' + hsp[0],
						'Query loc: ' + str(hsp[6]) + '..' + str(hsp[7]) + ' of ' + str(query_len),
						'Sbjct: ' + stitle,
						'Sbjct loc: ' + str(hsp[8]) + '..' + str(hsp[9]) + ' of ' + str(sbjct_len),
						prot_query,
						nuc_query,
						prot_sbjct + '\n',
						))


	# write the output
	hits_file.write('\n'.join('\n'.join(i) for i in hit_list))
