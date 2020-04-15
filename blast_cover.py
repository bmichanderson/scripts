#!/usr/bin/env python

#####################
# Script: blast_cover.py
# Author: B. Anderson
# Date: 15 Apr 2020
# Description: parse a plot.tab file from parse_blastn to calculate the percentage cover
#####################

import sys

def help():
	print('A script to parse a plot.tab file from parse_blastn to calculate the percentage cover')
	print('')
	print('Usage: ' + sys.argv[0] + ' plot.tab')
	print('')


# print help if the script is called without arguments
if (len(sys.argv[1:])) == 0:
	sys.exit(help())


# read in the plot.tab file line by line and capture the hit information
hit_dict = {}
sbjct_len_dict = {}
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
			hit_dict[sbjct] = [sorted((sbjct_start, sbjct_end))]
		else:
			hit_dict[sbjct].append(sorted((sbjct_start, sbjct_end)))


		# capture the subject lengths and query presence
		if sbjct not in sbjct_len_dict:
			sbjct_len_dict[sbjct] = sbjct_len

		if query not in query_list:
			query_list.append(query)



# Now process the capture data to avoid overlapping hits
keep_dict = {}
for sbjct in hit_dict:
	keep_hits = []
	hit_dict[sbjct].sort(key=lambda x: x[1]-x[0], reverse=True)
	keep_hits.append(hit_dict[sbjct][0])
	for hit in hit_dict[sbjct][1:]:
		skip = False
		for keeper in keep_hits:
			if all([hit[0] > keeper[0], hit[1] < keeper[1]]):	# if the hit is within a previous hit
				skip = True
		if skip:
			continue
		else:
			keep_hits.append(hit)

	keep_dict[sbjct] = keep_hits


# Now calculate the coverage of the remaining hits, one subject at a time
tally = 0
cover_dict = {}
for sbjct in keep_dict:
	hit_cover = 0
	for hit in keep_dict[sbjct]:
		hit_cover = hit_cover + (hit[1] - hit[0])
		tally = tally + (hit[1] - hit[0])

	cover_dict[sbjct] = hit_cover


# report the findings
total_len = 0
for key, value in sbjct_len_dict.iteritems():
	total_len = total_len + value

print('Summary for BLAST results with ' + str(len(sbjct_len_dict)) + ' subjects and ' + str(len(query_list)) + ' queries.')
print('')
print('Hits covered:')
for key, value in sorted(cover_dict.iteritems()):
	print(key + ': ' + str(value) + ' of ' + str(sbjct_len_dict[key]))
print('')
print('Total covered: ' + str(tally) + ' of ' + str(total_len))
