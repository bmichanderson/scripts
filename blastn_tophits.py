#!/usr/bin/env python

#####################
# Author: B. Anderson
# Date: 1, 3, 4 June 2020
# Modified: May 2021 to prioritize length and keep more hits
# Description: parse a plot_data.tab file from blastn_parse with multiple queries to keep the top hits
#####################

import argparse
import sys


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to parse plot_data.tab from blastn_parse.py with multiple queries to keep the top hits')


# add arguments to parse
parser.add_argument('plot_data', type=str, help='The plot_data.tab output to parse')
parser.add_argument('-d', type=int, dest='id_diff', help='The maximum difference in percent identity to keep a hit ecompassed by another [default 0]')
parser.add_argument('-l', type=int, dest='min_length', help='The minimum length of a hit or overlap to consider complete and equivalent [default 10]')
parser.add_argument('-p', type=int, dest='min_id', help='The minimum percent identity to keep a hit [default 50]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

plot_data = args.plot_data

if args.id_diff:
	id_diff  = args.id_diff
else:
	id_diff = 0

if args.min_length:
	min_length = args.min_length
else:
	min_length = 10

if args.min_id:
	min_id = args.min_id
else:
	min_id = 50



# read in the plot_data.tab file line by line and capture the hit information
hit_dict = {}
with open(plot_data, 'r') as input:
	for lineno, line in enumerate(input, start=1):
		if lineno == 1:
			continue

		line_list = line.split()

		sbjct = line_list[0]
		sbjct_start = int(line_list[1])
		sbjct_end = int(line_list[2])
		query = line_list[3]
		query_start = int(line_list[4])
		query_end = int(line_list[5])
		score = float(line_list[6])
		sstrand = line_list[7]
		qlen = int(line_list[8])
		slen = int(line_list[9])

		# capture the hit information by subject
		if sbjct not in hit_dict:
			hit_dict[sbjct] = [(sbjct, sbjct_start, sbjct_end, query, query_start, query_end, score, sstrand, qlen, slen)]
		else:
			hit_dict[sbjct].append((sbjct, sbjct_start, sbjct_end, query, query_start, query_end, score, sstrand, qlen, slen))


# For each subject, determine the top hits and keep them
top_hits_dict = {}
for sbjct in hit_dict:
	counter = 2
	pass_no = 1
	sbjct_list = hit_dict[sbjct][:]
	identicals = 0
	too_short = 0

	while counter > 0:
		keep_hits = []
		ranges_covered = []

		# sort the hits by score, then by length (prioritizing length)
		sbjct_list.sort(key=lambda x: x[6], reverse=True)
		sbjct_list.sort(key=lambda x: (x[2] - x[1]), reverse=True)

		# for each hit, calculate whether it overlaps another hit and the comparative score
		keep_hits.append(sbjct_list[0])
		ranges_covered.append([sbjct_list[0][1], sbjct_list[0][2], sbjct_list[0][6], sbjct_list[0][3]])		# start, end, score, query

		for hit in sbjct_list[1:]:
			keep = True
			query = hit[3]
			start = hit[1]
			end = hit[2]
			score = hit[6]
			for index, range in enumerate(ranges_covered):
				if score < min_id:		# hit doesn't fit filter
					keep = False
					break

				if end - start < min_length:	# hit is too short
					too_short = too_short + 1
					keep = False
					break

				if any([start > range[1], end < range[0]]):		# hit is entirely outside range, so keep it
					continue

				elif all([start >= range[0], end <= range[1]]):		# hit is completely within range
					if score == range[2]:
						if query == range[3]:			# same query, drop it
							keep = False
							break
						else:
							identicals = identicals + 1
							break
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False
						break

				elif any([start < range[0], end > range[1]]):		# hit overlaps one or both sides
					if range[0] - start >= min_length:		# longer overlap, so keep as distinct hit
						continue
					elif end - range[1] >= min_length:		# longer overlap, so keep as distinct hit
						continue
					elif score == range[2]:
						identicals = identicals + 1
						break
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False
						break
				else:
					continue

			if keep:
				if hit not in keep_hits:
					keep_hits.append(hit)
				if [start, end, score, query] not in ranges_covered:
					ranges_covered.append([start, end, score, query])


		# re-set the sbjct list
		sbjct_list = keep_hits

		# possibly re-iterate
		print('Pass ' + str(pass_no) + ' complete for subject ' + sbjct)
		print('Found ' + str(identicals) + ' hits with identical score that entirely matched or overlapped with minimal difference')
		if too_short > 0:
			print('Chuck ' + str(too_short) + ' hits with length < min_length')
		counter = counter - 1
		pass_no = pass_no + 1
		identicals = 0
		too_short = 0

	# Now that the while loop is finished, output the final sbjct_list
	top_hits_dict[sbjct] = sbjct_list


# Now that top hits have been collected for each sbjct, we can output a new BED-like file with only the top hits
print('Plotting top hits to new_plot_data.tab')
with open('new_plot_data.tab', 'w') as outfile:
	outfile.write('\t'.join(['sbjct', 'sbjct_start', 'sbjct_end', 'query', 'query_start', 'query_end', 'score', 'sstrand', 'qlen', 'slen']) + '\n')
	for sbjct in top_hits_dict:
		for hit in top_hits_dict[sbjct]:
			outfile.write('\t'.join(str(x) for x in hit) + '\n')
