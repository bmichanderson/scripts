#!/usr/bin/env python

#####################
# Author: B. Anderson
# Date: 15 Apr 2020, 21 May 2020 (updated to now handle the more BED-like output of the new blastn_parse script)
# Description: parse a plot_data.tab file from blastn_parse to calculate the percentage cover
#####################

import sys


def help():
	print('A script to parse a plot_data.tab file from parse_blastn to calculate the percentage cover')
	print('')
	print('Usage: ' + sys.argv[0] + ' plot_data.tab')
	print('')


# print help if the script is called without arguments
if (len(sys.argv[1:])) == 0:
	sys.exit(help())


# read in the plot_data.tab file line by line and capture the hit information
hit_dict = {}
sbjct_len_dict = {}
cover_dict = {}
sbjct_list = []
query_list = []
with open(sys.argv[1], 'r') as input:
	for lineno, line in enumerate(input, start=1):
		if lineno == 1:
			continue

		sbjct = line.split()[0]
		sbjct_start = int(line.split()[1])
		sbjct_end = int(line.split()[2])
		query = line.split()[3]
		query_start = int(line.split()[4])
		query_end = int(line.split()[5])
		#score = float(line.split()[6])
		#strand = line.split()[7]
		if len(line.split()) > 8:
			query_len = int(line.split()[8])
			sbjct_len = int(line.split()[9])
		else:
			query_len = 0
			sbjct_len = 0


		# ignore lines with full length hits with same query and sbjct
		if all([query_start == 0, sbjct_start == 0, sbjct == query]):
			continue


		# capture the hit information by subject
		if sbjct not in hit_dict:
			#hit_dict[sbjct] = [sorted((sbjct_start, sbjct_end))]
			# NOTE: since these are now proper BED coordinates, the actual base positions are automatically lesser to greater and 0-based and non-inclusive of end position
			hit_dict[sbjct] = [(sbjct_start, sbjct_end)]
		else:
			#hit_dict[sbjct].append(sorted((sbjct_start, sbjct_end)))
			hit_dict[sbjct].append((sbjct_start, sbjct_end))


		# capture the subject lengths and sbjct/query presence
		if sbjct not in sbjct_len_dict:
			sbjct_len_dict[sbjct] = sbjct_len

		if sbjct not in sbjct_list:
			sbjct_list.append(sbjct)

		if query not in query_list:
			query_list.append(query)


# First convert the hits into non-overlapping ranges, then calculate the length of the ranges
for sbjct in hit_dict:
	hit_ranges = []
	for begin, end in sorted(hit_dict[sbjct]):		# from https://stackoverflow.com/questions/15273693/union-of-multiple-ranges
		#if hit_ranges and hit_ranges[-1][1] >= begin - 1:
		if hit_ranges and hit_ranges[-1][1] >= begin:		# had to change this since BED ranges don't include the last position
			hit_ranges[-1][1] = max(hit_ranges[-1][1], end)
		else:
			hit_ranges.append([begin, end])

	total_length = 0
	for begin, end in hit_ranges:
		#total_length = total_length + len(range(begin, end + 1))
		total_length = total_length + len(range(begin, end))		# another change since BED lengths are the same as range lengths (i.e. end - begin)

	if sbjct not in cover_dict:
		cover_dict[sbjct] = total_length
	else:
		print('Duplicate subject?')


# report the findings
cover_len = 0
total_len = 0
for sbjct in cover_dict:
	cover_len = cover_len + cover_dict[sbjct]
for sbjct in sbjct_len_dict:
	total_len = total_len + sbjct_len_dict[sbjct]

print('Summary for BLAST results with ' + str(len(sbjct_len_dict)) + ' subjects and ' + str(len(query_list)) + ' queries.')
print('')
print('Hit coverage:')
for sbjct in sorted(sbjct_list):
	if sbjct in cover_dict:
		print(sbjct + ': ' + str(cover_dict[sbjct]) + ' of ' + str(sbjct_len_dict[sbjct]))
print('')
print('Total covered: ' + str(cover_len) + ' of ' + str(total_len))
