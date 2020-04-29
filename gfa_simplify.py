#!/usr/bin/env python

#####################
# Author: B. Anderson
# Date: 29 April 2020
# Description: attempt to simplify a complicated gfa graph structure from Bandage and output contigs
#####################

import sys
#from Bio import SeqIO
from Bio.Seq import Seq

def help():
	print('Simplify complicated gfa graph structures from Bandage and output contigs')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' gfa_file')
	print('')


# print help if no arguments provided
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# read in the gfa file and capture the sequences, their lengths, and the links
seq_list = []
link_list = []
length_tally = 0
with open(sys.argv[1], 'r') as gfa_file:
	for line in gfa_file:
		if line.split()[0] == 'S':		# a sequence line
			seq_num = int(line.split()[1])
			seq_str = line.split()[2]
			seq_len = int(line.split()[3].split(':')[2])
			if len(seq_str) != seq_len:
				print('Problem with sequence length for sequence ' + str(seq_num) + '!')
			length_tally = length_tally + seq_len
			if seq_num not in seq_list:
				seq_list.append([seq_num, seq_str])
			else:
				print('Problem with sequence labelling for sequence ' + str(seq_num) + '!')

		elif line.split()[0] == 'L':		# a link line
			seq_1 = int(line.split()[1])
			seq_1_dir = line.split()[2]
			seq_2 = int(line.split()[3])
			seq_2_dir = line.split()[4]
			overlap = line.split()[5]
			if overlap != '0M':
				print('Overlap detected between ' + str(seq_1) + ' and ' + str(seq_2) + '!')
			link_list.append([seq_1, seq_1_dir, seq_2, seq_2_dir])

			# also append the corresponding connection in the other direction
			if seq_2_dir == '+':
				seq_1b_dir = '-'
			else:
				seq_1b_dir = '+'
			if seq_1_dir == '+':
				seq_2b_dir = '-'
			else:
				seq_2b_dir = '+'
			link_list.append([seq_2, seq_1b_dir, seq_1, seq_2b_dir])
		else:
			continue

print('Read in ' + str(len(seq_list)) + ' nodes and ' + str(int(len(link_list)/2)) + ' links')
print('Total length of assembly graph: ' + str(length_tally))
print('')

# define a function to retrieve the links for a given sequence number
def find_links(seq_num, link_list):
	hit = False
	miss = False
	links = []
	for link in sorted(link_list):
		if seq_num != link[0]:
			continue

		if link[0] == link[2]:		# circles on itself
			return 'Circle'

		if all([hit, miss]):		# if we have hit a link previously and then missed on the last one (avoid searching all the sorted links)
			return links

		if seq_num == link[0]:		# a link including our sequence number
			hit = True
			links.append(link)

		if all([hit, seq_num != link[0]]):
			miss = True
	if not links:
		return 'Not'
	else:
		return links


# define a function for recursive search of paths (depends on find_links)
def recursive_paths(seq_num, seq_dir, link_list, path_list, visited_list):
	global all_paths		# needs to be defined each time prior to calling this function

	# determine links for the input seq_num
	next_links = find_links(seq_num, link_list)
	if next_links == 'Not':		# shouldn't happen
		sys.exit('Hit a \"Not\" return from links!')		# troubleshooting
	else:
		hit = False
		for link in next_links:
			if not link:		# shouldn't happen
				print('Hit a non-link!')
				continue

			if link[1] == seq_dir:		# a continuation from the previous contig
				hit = True
				#if (str(link[2]) + link[3]) in visited_list:		# circling back to a previously visited contig = end of this path
				if link[2] in visited_list:
#					print('Circled back to ' + str(link[2]) + link[3])
					all_paths.append(path_list)
				else:
					link_path = path_list[:]
					link_path.append([link[2], link[3]])
					#visited_list.append(str(link[2]) + link[3])
					visited_list.append(link[2])
					recursive_paths(link[2], link[3], link_list, link_path, visited_list)
			else:
				continue

		if not hit:
#			print('Last link search did not find a continuation')
			all_paths.append(path_list)


# starting with the longest node, follow links to determine all possible paths
# find longest paths and non-redundant representation of the possible paths
circular = []
linear = []
filt_paths = []
for seq_num, seq_str in seq_list[:20]:		# the seq_list is automatically sorted by length/line in the gfa file
	links = find_links(seq_num, link_list)

	print(str(seq_num) + ':')

	if links == 'Circle':
		print('Circle')
		print('')
		circular.append(seq_num)
		continue
	elif links == 'Not':
		print('Not')
		print('')
		continue
	else:
		linear.append(seq_num)

	for link in links:
		all_paths = []		# a global variable to be updated for each link in recursive_paths
		path_list = []
		path_list.append([link[0], link[1]])
		path_list.append([link[2], link[3]])
		visited_list = []
		#visited_list.append(str(link[0]) + link[1])
		visited_list.append(link[0])
		#visited_list.append(str(link[2]) + link[3])
		visited_list.append(link[2])

		# recursively determine possible paths, stopping when hitting a contig visited before
		recursive_paths(link[2], link[3], link_list, path_list, visited_list)

		# evaluate the paths captured
		keep_paths = []
		touched = []
		dropped = 0
		for path in sorted(all_paths, key=len, reverse=True):
			keep = False
			nums = [step[0] for step in path]
			for num in nums:
				if num in touched:
					continue
				else:
					touched.append(num)
					keep = True
			if keep:
				keep_paths.append(path)
			else:
				dropped = dropped + 1

#		print(set(range(1, len(seq_list))) - set(touched))
#		print('Dropped: ' + str(dropped) + ' of ' + str(len(all_paths)))

		for path in keep_paths:
			filt_paths.append(path)


# Now find longest paths and non-redundanct representation across all paths from every starting point
keep_paths = []
touched = []
dropped = 0
for path in sorted(filt_paths, key=len, reverse=True):
	keep = False
	nums = [step[0] for step in path]
	for num in nums:
		if num in touched:
			continue
		else:
			touched.append(num)
			keep = True

	if keep:
		keep_paths.append(path)
	else:
		dropped = dropped + 1

print('Overall, there are ' + str(len(keep_paths)) + ' paths retained')
print('Dropped: ' + str(dropped) + ' of ' + str(len(filt_paths)))
print('')
print('Difference between possible and included contigs:')
print(set(range(1, len(seq_list))) - set(touched))
print('')
print('Circular contigs:')
print(circular)
#print(set(range(1, len(seq_list))) - set(linear))
print('')


# Define a function for simplifying overlapping series
def remove_overlap(paths_list):
	num_series = [[step[0] for step in paths_list[0]]]
	new_paths = [paths_list[0]]
	for path in paths_list[1:]:
		new_1 = []
		new_2 = []
		match_length = 0
		nums = [step[0] for step in path]
		for series in num_series:
			matches = []
			i = 0
			j = 0
			j_start = 0
			j_end = 0
			for num in nums:
				if num in series:
					i = series.index(num)
					j_start = j
					break
				else:
					j = j + 1

			while i < len(series) and j < len(nums):
				if series[i] == nums[j]:
					matches.append(j)
					i = i + 1
					j = j + 1
				else:
					j_end = j
					break
			if all([len(matches) > match_length, len(matches) > 5]):
				print(matches)
				match_length = len(matches)
				new_1 = path[0:j_start]
				new_2 = path[j_end:]
#				print(new_1)
#				print(new_2)
#				print('')

		if any([len(new_1) > 1, len(new_2) > 1]):
			if len(new_1) > 1:
				new_paths.append(new_1)
				num_series.append([step[0] for step in new_1])
			if len(new_2) > 1:
				new_paths.append(new_2)
				num_series.append([step[0] for step in new_2])
		else:
			new_paths.append(path)

	return new_paths


# Further simplifying to remove sections of paths that are fully contained in longer paths (iteratively)
new_paths = keep_paths
i = 0
while i < 3:
	new_paths = remove_overlap(new_paths)
	print(len(new_paths))
	i = i + 1


# Now create output paths and contigs from the filtered paths
final_paths = new_paths

touched = []
for path in sorted(final_paths, key=len, reverse=True):
	nums = [step[0] for step in path]
	for num in nums:
		if num in touched:
			continue
		else:
			touched.append(num)

print('Final difference between possible and included contigs:')
print(set(range(1, len(seq_list))) - set(touched))
print('')


with open('paths_out.txt', 'w') as outfile:
	print('Outputting path sequences...')
	index = 1
	for path in final_paths:
		link_str = ''
		for step in path:
			link_str = link_str + str(step[0]) + step[1] + ','
		outfile.write('Path ' + str(index) + ': ' + link_str.rstrip(',') + '\n')
		index = index + 1

with open('contigs_out.fasta', 'w') as outfile:
	print('Outputting contigs...')
	index = 1
	for path in final_paths:
		contig_seq = ''
		for step in path:
			seq_num = step[0]
			seq_dir = step[1]
			if seq_dir == '+':
				if seq_num == seq_list[seq_num - 1][0]:
					step_seq = seq_list[seq_num - 1][1]
				else:
					print('Problem with sequence ordering!')
			else:
				if seq_num == seq_list[seq_num - 1][0]:
					step_seq = str(Seq(seq_list[seq_num - 1][1]).reverse_complement())
				else:
					print('Problem with sequence ordering!')
			contig_seq = contig_seq + step_seq
		outfile.write('>contig' + str(index) + '\n' + contig_seq + '\n')
		index = index + 1
