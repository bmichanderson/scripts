#!/usr/bin/env python3

###############################################
# Author: B. Anderson
# Date: May 2023
# Modified: Jun-July 2023
# Description: Assess heterozygous positions in hybrid sequences (sample ID file supplied with -f) in
#	fasta multiple sequence alignment(s) to determine likely parent samples
# Note: Loci with fewer than 8 heterozygous positions for a given hybrid will not be counted
###############################################


import sys
import os
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to assess hybrid heterozygous positions in fasta multiple sequence' +
	'alignments to determine potential parent samples')


# add arguments to parse
parser.add_argument('alignments', type = str, nargs = '*', help = 'The multiple sequence alignment(s)')
parser.add_argument('-f', type = str, dest = 'hybrids_file', help = 'The file with sample IDs of hybrids in the alignments (one per line)')
parser.add_argument('-t', type = int, dest = 'thresh', help = 'The threshold (greater than or equal to) for minimum percent of ' +
	'the total for reporting a top parent count [default 20]; e.g., only report if the top count is at least 20 percent of ' +
	'the total number of ambiguities')
parser.add_argument('-r', type = int, dest = 'rep', help = 'The reporting level for number of parents; if there are more than this ' +
	'number of unique parents for a given count, only report \">x\" [default x = 10]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
alignments = args.alignments
hybrids_file = args.hybrids_file
threshold = args.thresh
rep_thresh = args.rep

if any([not alignments, not hybrids_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)

if not threshold:
	threshold = 20

if not rep_thresh:
	rep_thresh = 10


# define functions

## a function to separate specified hybrids from other samples in a multiple sequence alignment
## It returns a list of two lists: a list of hybrid sequences, and a list of other sequences (+ hybrids)
def process_alignment(alignment, hybrid_list):
	hybrid_seqs = []
	other_seqs = []
	for sequence in alignment:
		if sequence.id in hybrid_list:
			hybrid_seqs.append(sequence)
		else:
			other_seqs.append(sequence)
	other_seqs.extend(hybrid_seqs)		# in case a hybrid may be a parent
	return([hybrid_seqs, other_seqs])


## A function to find ambiguities in a list of hybrid sequences
## It returns a list of lists, with each sub-list consisting of a hybrid ID and
## a list of the location and base of ambiguities in the sequence
def get_ambiguities(hybrid_seqs):
	ambiguity_list = ['K', 'k', 'M', 'm', 'R', 'r', 'S', 's', 'W', 'w', 'Y', 'y']
	output_list = []
	for hybrid_seq in hybrid_seqs:
		ambiguities = []
		for index, base in enumerate(hybrid_seq.seq):
			if base in ambiguity_list:
				ambiguities.append([index, base])
		output_list.append([hybrid_seq.id, ambiguities])
	return(output_list)


## a function to determine groups of parent samples for a given set of ambiguities
## The provided input "ambiguities" is a list of [position, base] of ambiguities in an alignment
## The input "alignment" is a multiple sequence alignment of potential parents
## It returns a list of groups of parents
def get_groups(ambiguities, alignment):
	# define an ambiguity dictionary
	amb_dict = {
		'K' : ['G', 'T'],
		'k' : ['g', 't'],
		'M' : ['A', 'C'],
		'm' : ['a', 'c'],
		'R' : ['A', 'G'],
		'r' : ['a', 'g'],
		'S' : ['C', 'G'],
		's' : ['c', 'g'],
		'W' : ['A', 'T'],
		'w' : ['a', 't'],
		'Y' : ['C', 'T'],
		'y' : ['c', 't']
	}
	output_list = []
	for ambiguity in ambiguities:
		base1 = amb_dict[ambiguity[1]][0]
		base2 = amb_dict[ambiguity[1]][1]
		this_column = alignment[:, ambiguity[0]]		# slice all rows (samples), one position
		# determine if the column includes both possible bases for the ambiguity
		#if ambiguity[1] in set(this_column):		# found other identical ambiguities
		#	continue
		if all([base1 in set(this_column), base2 in set(this_column)]):
			# collect the two groups of parents
			parents1 = []
			parents2 = []
			for index, base in enumerate(this_column):
				if base == base1:
					parents1.append(alignment[index].id)
				elif base == base2:
					parents2.append(alignment[index].id)
				else:
					continue
			output_list.append([parents1, parents2])
		else:		# couldn't find both bases
			continue
	return(output_list)


## a function to get counts of parents occurring in opposite groups
## The input "parents" are a list of all possible parents to compare
## The input "groups" are the groupings to assess (a list of lists each with two lists of parents)
## It returns a list (sorted by count) of lists with parent combinations and their counts
def get_opposites(parents, groups):
	output_list = []
	for index, parent in enumerate(parents):
		if index < (len(parents) - 1):	# if at least two more parents to compare
			for parent2 in parents[index + 1: ]:
				count = 0
				for group in groups:
					if any([
						all([parent in group[0], parent2 in group[1]]),
						all([parent in group[1], parent2 in group[0]])
						]):		# different groups
						count = count + 1
					else:
						continue
				output_list.append([[parent, parent2], count])
	# now sort the parent combinations by how often they group separately
	output_list.sort(key = lambda x: x[1], reverse = True)
	return(output_list)


## a function to create a parent result list
## The six inputs are:
## The hybrid ID, the filename of the alignment, the total number of ambiguities,
## and lists for each of the top three scoring parent combos,
## which include a list of the count and the combos (list of two samples) with that count
## It returns a list with the following fields for reporting (12):
## 'Hybrid', 'Alignment', 'Ambiguities', 'count1', 'Parent1a', 'Parent1b',
## 'count2', 'Parent2a', 'Parent2b', 'count3', 'Parent3a', 'Parent3b'
def parent_report(hybrid, alignment_name, total, parent1_list, parent2_list, parent3_list):
	output_list = []
	output_list.append(str(hybrid))
	output_list.append(str(alignment_name))
	output_list.append(str(total))

	for parent_list in [parent1_list, parent2_list, parent3_list]:
		output_list.append(str(parent_list[0]))		# the count

		# determine if there are multiple combos
		if len(parent_list[1]) > 1:
			if len(set(sum(parent_list[1], []))) > rep_thresh:		# more than x unique parents
				output_list.append('>' + str(rep_thresh))
				output_list.append('>' + str(rep_thresh))
			else:
				parenta_list = []
				parentb_list = []
				for combo in parent_list[1]:
					if any([combo[0] in parenta_list, combo[0] in parentb_list]):
						if any([combo[1] in parenta_list, combo[1] in parentb_list]):
							continue
						elif combo[0] in parenta_list:
							parentb_list.append(combo[1])
						else:
							parenta_list.append(combo[1])
					elif combo[1] in parenta_list:
						parentb_list.append(combo[0])
					elif combo[1] in parentb_list:
						parenta_list.append(combo[0])
					else:
						parenta_list.append(combo[0])
						parentb_list.append(combo[1])

				parenta_string = '/'.join(list(set(parenta_list)))
				parentb_string = '/'.join(list(set(parentb_list)))
				output_list.append(parenta_string)
				output_list.append(parentb_string)
		else:
			output_list.append(str(parent_list[1][0][0]))
			output_list.append(str(parent_list[1][0][1]))

	return(output_list)




###########################
# Run
###########################


# get the hybrid sample IDs
hybrids = []
with open(hybrids_file, 'r') as infile:
	for line in infile:
		hybrids.append(line.strip())


# create a list to keep track of the results for each hybrid across loci
# The idea is to retain:
#	the hybrid ID
#	the alignment name
#	the number of ambiguities for that hybrid in the alignment
#	the top count and parent combos if there are enough ambiguities to evaluate them
#	the second and third counts and parent combos
result_list = []


# print the categories for the output to stdout
print('\t'.join(['Hybrid', 'Alignment', 'Ambiguities', 'count1', 'Parent1a', 'Parent1b',
	'count2', 'Parent2a', 'Parent2b', 'count3', 'Parent3a', 'Parent3b']))


# read in each alignment, extract hybrid sequences, retain non-hybrid samples in a new alignment,
# then process the hybrids and their ambiguities
for alignment_file in alignments:
	#print('\nProcessing alignment ' + os.path.basename(alignment_file))
	sequences = process_alignment(AlignIO.read(open(alignment_file, 'r'), 'fasta'), hybrids)

	# turn the non-hybrid sequences into an alignment
	clean_alignment = MultipleSeqAlignment(sequences[1])

	# get the locations of and bases for ambiguities in the hybrid sequences
	hybrids_ambiguities = get_ambiguities(sequences[0])

	# for each hybrid, determine the parent groups for its valid ambiguities
	# ambiguities are valid if there are samples carrying both bases in that column
	for hybrid_entry in hybrids_ambiguities:
		# go through the hybrid's ambiguities and determine groups of parents
		groups = get_groups(hybrid_entry[1], clean_alignment)

		# collate parent samples into a list
		parents = list(set(sum(sum(groups, []), [])))

		if any([len(parents) < 2, len(hybrid_entry[1]) < 8]):	# too few parents or sites
			result_list.append(parent_report(hybrid_entry[0], os.path.basename(alignment_file),
				len(hybrid_entry[1]), ['-', [['-', '-']]], ['-', [['-', '-']]], ['-', [['-', '-']]]))
			continue		# move on to the next hybrid and its ambiguities
		else:
			# get counts of parent combos when they occur in opposite groups
			opposite_count = get_opposites(parents, groups)

			# determine the top parent combos
			count_vals = sorted(list(set([item[1] for item in opposite_count])), reverse = True)
			count1 = count_vals[0]

			if (count1 / len(hybrid_entry[1])) * 100 < threshold:
				result_list.append(parent_report(hybrid_entry[0], os.path.basename(alignment_file),
					len(hybrid_entry[1]), ['-', [['-', '-']]], ['-', [['-', '-']]], ['-', [['-', '-']]]))
			else:
				if len(count_vals) > 1:
					count2 = count_vals[1]
				else:
					count2 = 0
				if len(count_vals) > 2:
					count3 = count_vals[2]
				else:
					count3 = 0

				parents1 = [item[0] for item in opposite_count if item[1] == count1]
				parents2 = [item[0] for item in opposite_count if item[1] == count2]
				parents3 = [item[0] for item in opposite_count if item[1] == count3]

				if len(parents3) > 0:
					result_list.append(parent_report(hybrid_entry[0], os.path.basename(alignment_file),
						len(hybrid_entry[1]), [count1, parents1], [count2, parents2], [count3, parents3]))
				elif len(parents2) > 0:		# only two sets of parents
					result_list.append(parent_report(hybrid_entry[0], os.path.basename(alignment_file),
						len(hybrid_entry[1]), [count1, parents1], [count2, parents2], [count3, [['-', '-']]]))
				else:	# only one set of parents
					result_list.append(parent_report(hybrid_entry[0], os.path.basename(alignment_file),
						len(hybrid_entry[1]), [count1, parents1], [count2, [['-', '-']]], [count3, [['-', '-']]]))


# output the result to stdout, first sorting by hybrid ID and then alignment name
outlist = sorted(result_list, key = lambda x: (x[0], x[1]))
for entry in outlist:
	print('\t'.join(entry))


# if there is more than one locus, filter the results to also report locus counts for top parents
# This only reports parents that were unambiguously top (ignores results wih >x)
if len(alignments) > 1:
	print('\n\n')
	print('\t'.join(['Hybrid', 'Loci', 'Loci_gt_thresh', 'count1', 'Parent1', 'count2',
	'Parent2', 'count3', 'Parent3', 'count4', 'Parent4']))
	for hybrid in hybrids:
		numloci = 0
		numloci_gt_thresh = 0
		top_parents = []
		counts = []
		for result in result_list:
			if result[0] == hybrid:
				numloci = numloci + 1
				if all(['-' not in result[3], int(result[2]) > 0]):
					if (int(result[3]) / int(result[2])) * 100 >= threshold:
						numloci_gt_thresh = numloci_gt_thresh + 1
						if '>' not in result[4]:
							if '/' in result[4]:
								top_parents.extend(result[4].split('/'))
							else:
								top_parents.append(result[4])
							if '/' in result[5]:
								top_parents.extend(result[5].split('/'))
							else:
								top_parents.append(result[5])

		# now count the occurrences of each parent in top parents
		for parent in set(top_parents):
			counts.append([parent, top_parents.count(parent)])

		counts.sort(key = lambda x: x[1], reverse = True)

		# report the top four across all loci
		if len(counts) > 0:
			count1 = counts[0][1]
			parent1 = counts[0][0]
			count2 = '-'
			parent2 = '-'
			count3 = '-'
			parent3 = '-'
			count4 = '-'
			parent4 = '-'

			if len(counts) > 3:
				count4 = counts[3][1]
				parent4 = counts[3][0]

			if len(counts) > 2:
				count3 = counts[2][1]
				parent3 = counts[2][0]

			if len(counts) > 1:
				count2 = counts[1][1]
				parent2 = counts[1][0]

			print('\t'.join([str(hybrid), str(numloci), str(numloci_gt_thresh), str(count1),
				str(parent1), str(count2), str(parent2), str(count3), str(parent3),
				str(count4), str(parent4)]))
