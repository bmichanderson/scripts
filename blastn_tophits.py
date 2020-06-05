#!/usr/bin/env python

#####################
# Author: B. Anderson
# Date: 1, 3, 4 June 2020
# Description: parse a plot_data.tab file from blastn_parse with multiple queries to keep the top hits
#####################

import argparse
import sys


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to parse plot_data.tab from blastn_parse.py with multiple queries to keep the top hits')


# add arguments to parse
parser.add_argument('plot_data', type=str, help='The plot_data.tab output to parse')
parser.add_argument('-d', type=int, dest='id_diff', help='The maximum difference in percent identity to keep a hit at the same location [default 0]')
parser.add_argument('-l', type=int, dest='min_length', help='The minimum length of a partial hit (trimmed) [default 50]')
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
	min_length = 50

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
	counter = 3
	pass_no = 1
	sbjct_list = hit_dict[sbjct][:]
	identicals = 0

	while counter > 0:
		keep_hits = []
		ranges_covered = []

		# sort the hits by length then by score (i.e. score is priority)
		sbjct_list.sort(key=lambda x: (x[2] - x[1]), reverse=True)
		sbjct_list.sort(key=lambda x: x[6], reverse=True)

		# for each hit, calculate whether it overlaps another hit and the comparative score
		keep_hits.append(sbjct_list[0])
		ranges_covered.append([sbjct_list[0][1], sbjct_list[0][2], sbjct_list[0][6]])		# start, end, score
		for hit in sbjct_list[1:]:
			keep = True
			keep_remain = True
			start = hit[1]
			end = hit[2]
			score = hit[6]
			frag1 = False
			len1 = 1000000
			frag1_end = 0
			frag2 = False
			len2 = 1000000
			frag2_start = 0
			for index, range in enumerate(ranges_covered):
				if score < min_id:		# hit doesn't fit filter
					keep = False
					break

				if score > range[2] + id_diff:			# hit has a higher score
					counter = counter + 1		# queue another pass to re-sort the list
					break

				if all([start == range[0], end == range[1]]):		# hit coordinates are identical
					if score == range[2]:				# hit score is also identical, so keep it
						identicals = identicals + 1
						break
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False

				elif any([start > range[1], end < range[0]]):		# hit is outside range, so keep it
					continue

				elif all([start > range[0], end < range[1]]):		# hit is completely within range
					if score == range[2]:
						identicals = identicals + 1
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False

				elif start == range[0]:			# end is incongruous
					if all([end > range[1], end - range[1] >= min_length]):		# if extends and longer than min_length
						frag2 = True
						if end - range[1] < len2:		# if this is the shortest fragment
							frag2_start = range[1]
							len2 = end - frag2_start
							if range[2] - score > id_diff:
								keep_remain = False
					elif score == range[2]:
						identicals = identicals + 1
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False

				elif end == range[1]:			# start is incongruous
					if all([start < range[0], range[0] - start >= min_length]):	# if extends and longer than min_length
						frag1 = True
						if range[0] - start < len1:		# if this is the shortest fragment
							frag1_end = range[0]
							len1 = frag1_end - start
							if range[2] - score > id_diff:
								keep_remain = False
					elif score == range[2]:
						identicals = identicals + 1
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False

				elif all([start < range[0], end < range[1]]):		# hit overlaps left side
					if range[0] - start >= min_length:		# if longer than min_length
						frag1 = True
						if range[0] - start < len1:		# if this is the shortest fragment
							frag1_end = range[0]
							len1 = frag1_end - start
							if range[2] - score > id_diff:
								keep_remain = False
					elif score == range[2]:
						identicals = identicals + 1
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False

				elif all([start > range[0], end > range[1]]):		# hit overlaps right side
					if end - range[1] >= min_length:			# if longer than min_length
						frag2 = True
						if end - range[1] < len2:		# if this is the shortest fragment
							frag2_start = range[1]
							len2 = end - frag2_start
							if range[2] - score > id_diff:
								keep_remain = False
					elif score == range[2]:
						identicals = identicals + 1
					elif range[2] - score > id_diff:		# hit score is less by more than id_diff
						keep = False

				elif all([start < range[0], end > range[1]]):		# hit fully encompases range
					if range[2] - score <= id_diff:
						continue
					else:
						keep_remain = False
						if range[0] - start >= min_length:		# if longer than min_length
							frag1 = True
							if range[0] - start < len1:		# if this is the shortest fragment
								frag1_end = range[0]
								len1 = frag1_end - start

						if end - range[1] >= min_length:			# if longer than min_length
							frag2 = True
							if end - range[1] < len2:		# if this is the shortest fragment
								frag2_start = range[1]
								len2 = end - frag2_start

						if all([range[0] - start < min_length, end - range[1] < min_length]):
							keep = False

				else:
					continue

			if keep:
				if any([frag1, frag2]):
					if all([frag1, frag2]):
						new_hit1 = list(hit)
						new_hit1[2] = frag1_end
						if new_hit1 not in keep_hits:
							keep_hits.append(new_hit1)
						if [start, frag1_end, score] not in ranges_covered:
							ranges_covered.append([start, frag1_end, score])

						new_hit2 = list(hit)
						new_hit2[1] = frag2_start
						if new_hit2 not in keep_hits:
							keep_hits.append(new_hit2)
						if [frag2_start, end, score] not in ranges_covered:
							ranges_covered.append([frag2_start, end, score])

						if keep_remain:
							remainder = list(hit)
							remainder[1] = frag1_end
							remainder[2] = frag2_start
							if remainder not in keep_hits:
								keep_hits.append(remainder)
							if [frag1_end, frag2_start, score] not in ranges_covered:
								ranges_covered.append([frag1_end, frag2_start, score])

					elif frag1:
						new_hit1 = list(hit)
						new_hit1[2] = frag1_end
						if new_hit1 not in keep_hits:
							keep_hits.append(new_hit1)
						if [start, frag1_end, score] not in ranges_covered:
							ranges_covered.append([start, frag1_end, score])

						if keep_remain:
							remainder = list(hit)
							remainder[1] = frag1_end
							if remainder not in keep_hits:
								keep_hits.append(remainder)
							if [frag1_end, end, score] not in ranges_covered:
								ranges_covered.append([frag1_end, end, score])

					elif frag2:
						new_hit2 = list(hit)
						new_hit2[1] = frag2_start
						if new_hit2 not in keep_hits:
							keep_hits.append(new_hit2)
						if [frag2_start, end, score] not in ranges_covered:
							ranges_covered.append([frag2_start, end, score])

						if keep_remain:
							remainder = list(hit)
							remainder[2] = frag2_start
							if remainder not in keep_hits:
								keep_hits.append(remainder)
							if [start, frag2_start, score] not in ranges_covered:
								ranges_covered.append([start, frag2_start, score])

				else:
					if hit not in keep_hits:
						keep_hits.append(hit)
					if [start, end, score] not in ranges_covered:
						ranges_covered.append([start, end, score])

		# re-set the sbjct list
		sbjct_list = keep_hits

		# possibly re-iterate
		print('Pass ' + str(pass_no) + ' complete for subject ' + sbjct)
		print('Found ' + str(identicals) + ' hits with identical score that entirely matched or overlapped another')
		counter = counter - 1
		pass_no = pass_no + 1


	# Now that the while loop is finished, output the final sbjct_list
	top_hits_dict[sbjct] = sbjct_list


# Now that top hits have been collected for each sbjct, we can output a new BED-like file with only the top hits
print('Plotting top hits to new_plot_data.tab')
with open('new_plot_data.tab', 'w') as outfile:
	outfile.write('\t'.join(['sbjct', 'sbjct_start', 'sbjct_end', 'query', 'query_start', 'query_end', 'score', 'sstrand', 'qlen', 'slen']) + '\n')
	for sbjct in top_hits_dict:
		for hit in top_hits_dict[sbjct]:
			outfile.write('\t'.join(str(x) for x in hit) + '\n')
