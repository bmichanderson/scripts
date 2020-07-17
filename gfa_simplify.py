#!/usr/bin/env python

#####################
# Author: B. Anderson
# Date: 16 July 2020
# Description: attempt to simplify a complicated gfa graph structure from Bandage and output contigs
#####################

import sys
import argparse
from Bio.Seq import Seq
from collections import Counter


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to simplify complicated gfa graph structures from Bandage and output contigs')


# add arguments to parse
parser.add_argument('gfa_file', type=str, help='Required input of a gfa file representing the assembly structure')
parser.add_argument('-p', type=str, dest='paths_file', help='An option to read in an existing modified paths file for only exporting contigs')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
gfa_file = args.gfa_file
paths_file = args.paths_file

if paths_file:
	pathsfile_present = True
else:
	pathsfile_present = False


# read in the gfa file and capture the sequences, their lengths, and the links
seq_list = []
link_list = []
length_tally = 0
with open(gfa_file, 'r') as gfafile:
	for line in gfafile:
		line_elements = line.split()
		if line_elements[0] == 'S':		# a sequence line
			seq_num = line_elements[1]
			seq_str = line_elements[2]
			seq_len = int(line_elements[3].split(':')[2])
			if len(seq_str) != seq_len:
				print('Problem with sequence length for sequence ' + str(seq_num) + '!')
			length_tally = length_tally + seq_len
			if seq_num not in seq_list:
				seq_list.append([seq_num, seq_str])
			else:
				print('Problem with sequence labelling for sequence ' + str(seq_num) + '!')

		elif line_elements[0] == 'L':		# a link line
			seq_1 = line_elements[1]
			seq_1_dir = line_elements[2]
			seq_2 = line_elements[3]
			seq_2_dir = line_elements[4]
			overlap = int(line_elements[5].strip('M'))

			link_list.append([seq_1, seq_1_dir, seq_2, seq_2_dir, overlap])

			# also append the corresponding connection in the other direction
			if seq_2_dir == '+':
				seq_1b_dir = '-'
			else:
				seq_1b_dir = '+'
			if seq_1_dir == '+':
				seq_2b_dir = '-'
			else:
				seq_2b_dir = '+'

			link_list.append([seq_2, seq_1b_dir, seq_1, seq_2b_dir, overlap])

		else:		# skipping possible path lines
			continue

print('Read in ' + str(len(seq_list)) + ' nodes and ' + str(int(len(link_list)/2)) + ' links')
print('Total length of assembly graph: ' + str(length_tally))
print('')



#################################
# Define functions
#################################


## find_links
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


## recursive_paths
# define a function for recursive search of paths (depends on find_links)
def recursive_paths(seq_num, seq_dir, link_list, path_list, visited_list):
	global all_paths		# needs to be defined each time prior to calling this function

	# determine links for the input seq_num
	next_links = find_links(seq_num, link_list)
	if any([next_links == 'Not', next_links == 'Circle']):		# shouldn't happen
		sys.exit('Hit a \"Not\" or \"Circle\" return from links for ' + seq_num + '!')		# troubleshooting
	else:
		hit = False
		for link in next_links:
			if not link:		# shouldn't happen
				print('Hit a non-link!')		# troubleshooting
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


## longest_match
# Define a function for finding the longest hit (indices of first argument) between two lists of numbers
def longest_match(series1, series2):
	longest_matches = []
	for index, num in enumerate(series1):		# iterate through the series to find matching numbers in the other series
		match_list = []
		starts = [j for j, x in enumerate(series2) if x == num]
		if len(starts) > 0:
			for start in starts:		# see how long we can match between the two series
				matches = []
				i = index
				j = start
				while i < len(series1) and j < len(series2):
					if series1[i] == series2[j]:
						matches.append(i)
						i = i + 1
						j = j + 1
					else:
						break
				match_list.append(matches)
		else:
			continue
		if len(match_list) > 0:
			longest_match = sorted(match_list, key=len, reverse=True)[0]
			longest_matches.append(longest_match)
		else:
			continue
	if len(longest_matches) > 0:
		return sorted(longest_matches, key=len, reverse=True)[0]
	else:
		return []


## remove_overlap
# Define a function for simplifying overlapping paths (depends on longest_match)
def remove_overlap(paths_list):
	paths_list = sorted(paths_list, key=len, reverse=True)		# sort the list by length
	new_paths = [paths_list[0]]	# automatically keep the first (longest) path
	num_series = [[step[0] for step in paths_list[0]]]		# turn the first path into a series of numbers (no orientation)
	singletons = []
	for path in paths_list[1:]:	# iterate through the remaining paths, comparing them to the growing set of number series
		if len(path) == 1:
			singletons.append(path)
			continue
		new_1 = []	# initialize potential sub-paths created by removing matching sequences
		new_2 = []
		best_match_length = 1			# set this to determine minimum match retained
		nums = [step[0] for step in path]	# turn the path into a series of numbers
		nums_rev = nums[::-1]			# reverse the list to check for reverse hits
		for series in num_series:		# iterate over the existing processed paths to look for overlaps
			match_length = 0
			match = longest_match(nums, series)
			match_rev = longest_match(nums_rev, series)

			if len(match) >= len(match_rev):
				match_length = len(match)
				chosen_match = match
			else:
				match_length = len(match_rev)
				corr_match_rev = []			# we need to switch the indices back to the original path indices (but now forward too)
				for index in match_rev:
					orig_index = len(nums_rev) - index - 1
					corr_match_rev.append(orig_index)
				corr_match_rev.reverse()
				chosen_match = corr_match_rev

			if match_length > best_match_length:
				best_match_length = match_length
				new_1 = path[0:chosen_match[0]]
				new_2 = path[chosen_match[-1] + 1:]
			else:
				continue

		if any([len(new_1) > 0, len(new_2) > 0]):		# was split into at least one new path (the best will have assigned new_1 and new_2)
			if len(new_1) > 0:
				if len(new_1) == 1:
					singletons.append(new_1)
				else:
					new_paths.append(new_1)
					num_series.append([step[0] for step in new_1])
			if len(new_2) > 0:
				if len(new_2) == 1:
					singletons.append(new_2)
				else:
					new_paths.append(new_2)
					num_series.append([step[0] for step in new_2])
		elif len(path) == best_match_length:		# full length match
			continue
		else:				# no sufficiently long hits
			new_paths.append(path)
			num_series.append([step[0] for step in path])


	# Deal with singletons (duplicates and those present in existing paths)
	singletons = [list(i) for i in set(tuple(x[0]) for x in singletons)]		# remove duplicates

	remove_indices = []
	for index, singleton in enumerate(singletons):
		single_num = singleton[0]
		found = False
		for series in num_series:
			if single_num in series:		# if hit, regardless of orientation
				found = True
				break
		if found:
			remove_indices.append(index)
	remove_indices = list(set(remove_indices))		# remove duplicate indices
	for index in sorted(remove_indices, reverse=True):
		del singletons[index]
	print('Deleted ' + str(len(remove_indices)) + ' singletons')

	for singleton in singletons:
		new_paths.append([singleton])

	# Deal with duplicate paths of any length
	remove_indices = []
	for path in new_paths:
		hits = []
		for index, other_path in enumerate(new_paths):
			if path == other_path:
				hits.append(index)
		if len(hits) > 1:
			for hit in hits[1:]:
				remove_indices.append(hit)
	remove_indices = list(set(remove_indices))		# remove duplicate indices
	for index in sorted(remove_indices, reverse=True):
		del new_paths[index]
	print('Deleted ' + str(len(remove_indices)) + ' exact duplicate paths')


	# Finish and return the new paths without overlaps
	return new_paths






##########################
# Continue execution
##########################
# if there is no paths file provided, proceed to generate all possible paths
if not pathsfile_present:

	# starting with the longest node, follow links to determine all possible paths
	# find longest paths and non-redundant representation of the possible paths
	circular = []
	disconnected = []
	filt_paths = []
	seq_list = sorted(seq_list, key=lambda x: len(x[1]), reverse=True)		# not always automatically sorted
	for seq_num, seq_str in seq_list:
		links = find_links(seq_num, link_list)

		#print(str(seq_num) + ':')

		if links == 'Circle':
			#print('Circle')
			#print('')
			circular.append(seq_num)
			continue
		elif links == 'Not':
			#print('No links found')
			#print('')
			disconnected.append(seq_num)
			continue

		for link in links:
			all_paths = []		# a global variable to be updated for each link in recursive_paths
			path_list = []
			path_list.append([link[0], link[1]])
			path_list.append([link[2], link[3]])
			visited_list = []
			visited_list.append(link[0])
			visited_list.append(link[2])

			# recursively determine possible paths, stopping when hitting a contig visited before
			recursive_paths(link[2], link[3], link_list, path_list, visited_list)

			# evaluate the paths captured to remove those with all contigs hit in previous paths
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

			for path in keep_paths:
				filt_paths.append(path)


	# Now find longest paths and non-redundanct representation across *all* paths from every starting point
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
	print(set([x[0] for x in seq_list]) - set(touched) - set(circular) - set(disconnected))
	print('')
	if circular:
		print('Circular contigs:')
		print(circular)
		print('')
	if disconnected:
		print('Disconnected contigs:')
		print(disconnected)
		print('')


	# Now we need to reduce path overlap to try to get the smallest representation of the mess of connections
	# Simplifying to remove sections of paths that are contained in longer paths (iteratively)
	new_paths = keep_paths[:]
	i = 0
	while i < 10:		# set number of iterations to 10
		new_paths = remove_overlap(new_paths)
		print(len(new_paths))
		i = i + 1

	print('')


	# Further simplifying by merging any paths that end with the start of another path
	new_paths = sorted(new_paths, key=len, reverse=True)[:]
	merge_count = 0
	hit_end = False
	iteration_counter = 1
	while not hit_end:
		for index, path in enumerate(new_paths):
			hit_merge = False
			to_merge = ()

			#print('Index: ' + str(index) + ', Length: ' + str(len(new_paths)))
			if index >= len(new_paths) - 1:
				hit_end = True

			if len(path) == 1:
				continue

			for other_index, other_path in enumerate(new_paths):
				if len(other_path) == 1:
					continue
				if index == other_index:		# the same path
					continue
				if path[-1] == other_path[0]:		# exact same end and start point
					to_merge = (index, other_index)
					hit_merge = True
					break
				else:
					continue

			if hit_merge:
				break
			else:
				continue

		if len(to_merge) > 0:
			new_paths[to_merge[0]] = new_paths[to_merge[0]][:-1] + new_paths[to_merge[1]]
			del new_paths[to_merge[1]]
			merge_count = merge_count + 1

		iteration_counter = iteration_counter + 1

	print('Merged ' + str(merge_count) + ' pairs of paths')
	print('')


	# Additionally simplify by removing start or end path links that are already present in other paths (trimming shorter contigs)
	new_paths = sorted(new_paths, key=len, reverse=True)[:]
	trim_start_count = 0
	trim_end_count = 0
	hit_end = False
	while not hit_end:
		for index, path in enumerate(new_paths):
			tstart = False
			tend = False
			trim_start = 0
			trim_end = 0
			if index == len(new_paths) - 1:
				hit_end = True
			if len(path) == 1:
				continue

			nums = [step[0] for step in path]	# turn the path into a series of numbers

			for other_index, other_path in enumerate(new_paths):
				if len(other_path) == 1:
					continue
				if index == other_index:		# the same path
					continue

				other_nums = [step[0] for step in other_path]	# turn the path into a series of numbers

				if all([not tstart, nums[0] in other_nums]):		# start present inside another path
					tstart = True
					trim_start = index
				if all([not tend, nums[-1] in other_nums]):		# end present in another path
					tend = True
					trim_end = index

				if all([tstart, tend]):
					break
				else:
					continue
			if any([tstart, tend]):
				break
			else:
				continue

		if tstart:
			if len(new_paths[trim_start]) > 1:
				new_paths[trim_start] = new_paths[trim_start][1:]
				trim_start_count = trim_start_count + 1
			else:
				del new_paths[trim_start]
		if tend:
			if len(new_paths[trim_end]) > 1:
				new_paths[trim_end] = new_paths[trim_end][:-1]
				trim_end_count = trim_end_count + 1
			else:
				del new_path[trim_end]


	print('Trimmed ' + str(trim_start_count) + ' starts and ' + str(trim_end_count) + ' ends from paths')
	print('')


	# Delete any exact duplicates possibly created (singletons), and any singletons created that are now present in paths
	remove_indices = []
	is_singleton = False
	for path in new_paths:
		is_singelton = False
		hits = []
		if len(path) == 1:
			is_singleton = True
		for index, other_path in enumerate(new_paths):
			if path == other_path:
				hits.append(index)
			elif is_singleton:
				nums = [step[0] for step in other_path]
				if path[0][0] in nums:
					hits.append(index)
		if len(hits) > 1:
			for hit in hits[1:]:
				remove_indices.append(hit)
	remove_indices = list(set(remove_indices))		# remove duplicate indices
	for index in sorted(remove_indices, reverse=True):
		del new_paths[index]

	if len(remove_indices) > 0:
		print('Deleted ' + str(len(remove_indices)) + ' exact duplicate paths and singletons in other paths')
		print('')


	# If possible, re-connect singletons that have been orphaned
	new_paths = sorted(new_paths, key=len, reverse=True)[:]
	singletons = []
	singleton_nums = []
	path_starts = []
	path_ends = []
	linked = 0
	unlinked = 0
	to_link = []
	to_delete = []
	for index, path in enumerate(new_paths):
		if len(path) == 1:
			singletons.append((index, path))
			singleton_nums.append(path[0][0])
		else:
			path_starts.append(path[0])
			path_ends.append(path[-1])

	for end in path_ends:
		links = find_links(end[0], link_list)
		for link in links:
			if link[1] == end[1]:		# same direction
				if link[2] in singleton_nums:
					print('Found a hit from a path to a singleton!')
					print(link)
					for i, singleton in sorted(singletons, key=lambda x: x[1][0]):		# by contig number
						if link[2] == singleton[0][0]:
							for index, path in enumerate(new_paths):
								if [link[0], link[1]] == path[-1]:
									to_link.append((index, link))
									to_delete.append(i)


	for i, singleton in sorted(singletons, key=lambda x: x[1][0]):		# sort by contig number, so from largest actual sequence
		links = find_links(singleton[0][0], link_list)
		for link in links:
			if [link[2], link[3]] in path_starts:		# a direct connection
				print('Found a hit from a singleton to a path!')
				print(link)
				for index, path in enumerate(new_paths):
					if [link[2], link[3]] == path[0]:
						to_link.append((index, link))
						to_delete.append(i)

			elif link[2] in singleton_nums:
				print('Found a possible hit between singletons!')
				print(link)


	for index, link in to_link:
		if new_paths[index][0][0] == link[2]:		# the link is showing singleton --> path_start
			new_paths[index] = [[link[0], link[1]]] + new_paths[index]		# concatenate
			linked = linked + 1
		elif new_paths[index][-1][0] == link[0]:		# the link is showing path_end --> singleton
			new_paths[index] = new_paths[index] + [[link[2], link[3]]]
			linked = linked + 1

	for index in sorted(to_delete, reverse=True):
		del new_paths[index]

	for path in new_paths:
		if len(path) == 1:
			unlinked = unlinked + 1

	print('Singletons without links to other singletons or paths: ' + str(unlinked))
	print('Singletons re-linked with existing paths: ' + str(linked))
	print('')



########################
# Now, read in the paths file if present; otherwise, proceed with the existing paths
if pathsfile_present:
	circular = []
	disconnected = []
	paths = []
	with open(paths_file, 'r') as pfile:
		for line in pfile:
			path_string = line.strip().split()[2]
			path_list = path_string.split(',')
			path = [[x[0:-1], x[-1]] for x in path_list]

			if line.strip().split()[0] == 'Circle':
				circular.append(path[0][0])
			elif line.strip().split()[0] == 'Disc':
				disconnected.append(path[0][0])
			else:
				paths.append(path)
	final_paths = sorted(paths, key=len, reverse=True)


# Count and summarize the final path set
if not pathsfile_present:
	final_paths = sorted(new_paths, key=len, reverse=True)[:]

touched = []
for path in sorted(final_paths, key=len, reverse=True):
	nums = [step[0] for step in path]
	for num in nums:
		if num in touched:
			continue
		else:
			touched.append(num)

print('Total number of final paths: ' + str(len(final_paths)))
print('Final difference between possible and included contigs:')
print(set([x[0] for x in seq_list]) - set(touched) - set(circular) - set(disconnected))
print('')


# Create output paths and contigs from the filtered paths
contig_counting = []
with open('paths_out.txt', 'w') as outfile:
	print('Outputting path sequences to paths_out.txt')
	index = 1
	for path in final_paths:
		link_str = ''
		for step in path:
			link_str = link_str + step[0] + step[1] + ','
			contig_counting.append(step[0])
		outfile.write('Path ' + str(index) + ': ' + link_str.rstrip(',') + '\n')
		index = index + 1

	# need to output paths for circles and disconnected pieces
	if circular:
		index = 1
		for circ in circular:
			outfile.write('Circle ' + str(index) + ': ' + circ + '+' + '\n')
			index = index + 1
	if disconnected:
		index = 1
		for disc in disconnected:
			outfile.write('Disc ' + str(index) + ': ' + disc + '+' + '\n')
			index = index + 1


	# check if there are over-represented contigs
	over_contigs = []
	counts = Counter(contig_counting)
	for entry in counts:
		if counts[entry] > 1:
			over_contigs.append((entry, counts[entry]))
	two_times = []
	lfive = []
	gfive = []
	gten = []
	for entry, count in over_contigs:
		if count == 2:
			two_times.append(entry)
		elif count <= 5:
			lfive.append(entry)
		elif any([count > 5, count > 10]):
			if count <= 10:
				gfive.append(entry)
			else:
				gten.append(entry)
		else:
			print('Problem!')
	if any([two_times, lfive, gfive, gten]):
		print('Over-represented contigs:')
		if two_times:
			print('Two times: ' + str(sorted(two_times)))
		if lfive:
			print('<=5: ' + str(sorted(lfive)))
		if gfive:
			print('>5: ' + str(sorted(gfive)))
		if gten:
			print('>10: ' + str(sorted(gten)))


# Generate contigs
with open('contigs_out.fasta', 'w') as outfile:
	print('Outputting contigs to contigs_out.fasta')
	index = 1
	for path in final_paths:
		contig_seq = ''
		for i, step in enumerate(path):
			seq_num = step[0]
			seq_dir = step[1]
			seq_ind = [x[0] for x in seq_list].index(seq_num)
			seq_str = seq_list[seq_ind][1]
			if i < len(path) - 1:		# if the step continues to another, there is a chance of overlap
				for link in link_list:
					if all([link[0] == seq_num, link[2] == path[i + 1][0]]):
						overlap = link[4]
						break
				if seq_dir == '+':
					step_seq = seq_str[:len(seq_str)-overlap]
				else:
					step_seq = str(Seq(seq_str).reverse_complement()[:len(seq_str)-overlap])
			else:		# if the step is a singleton or the last step
				if len(path) == 1:	# if the step is a singleton (don't want to duplicate overlaps from paths that have been axed)
					max_overlap1 = 0
					max_overlap2 = 0
					hit_1 = False
					hit_2 = False
					for link in link_list:
						if link[0] == seq_num:
							if link[1] == '+':
								hit_1 = True
								if link[4] > max_overlap1:
									max_overlap1 = link[4]
							elif link[1] == '-':
								hit_2 = True
								if link[4] > max_overlap2:
									max_overlap2 = link[4]
					if all([hit_1, hit_2]):
						step_seq = seq_str[max_overlap2:len(seq_str)-max_overlap1]
					elif hit_1:
						step_seq = seq_str[:len(seq_str)-max_overlap1]
					elif hit_2:
						step_seq = str(Seq(seq_str).reverse_complement()[:len(seq_str)-max_overlap2])
					else:
						print('Singleton has no path links!?!?')		# troubleshooting
				else:		# if the step is the last step in a longer path
					if seq_dir == '+':
						step_seq = seq_str
					else:
						step_seq = str(Seq(seq_str).reverse_complement())
			contig_seq = contig_seq + step_seq
		outfile.write('>contig_p' + str(index) + '\n' + contig_seq + '\n')
		index = index + 1

	if circular:
		index = 1
		for circ in circular:
			for link in link_list:
				if all([link[0] == circ, link[2] == circ]):
					overlap = link[4]
					break
			seq_ind = [x[0] for x in seq_list].index(circ)
			contig_seq = seq_list[seq_ind][1][:-overlap]
			outfile.write('>contig_c' + str(index) + '\n' + contig_seq + '\n')
			index = index + 1

	if disconnected:
		index = 1
		for disc in disconnected:
			seq_ind = [x[0] for x in seq_list].index(disc)
			contig_seq = seq_list[seq_ind][1]
			outfile.write('>contig_d' + str(index) + '\n' + contig_seq + '\n')
			index = index + 1
