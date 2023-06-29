#!/usr/bin/env python3

###############################################
# Author: B. Anderson
# Date: May 2023
# Description: Assess heterozygous positions in hybrid sequences (sample ID file supplied with -f) in
#	fasta multiple sequence alignment(s) to determine likely parent samples
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


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
alignments = args.alignments
hybrids_file = args.hybrids_file

if any([not alignments, not hybrids_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)


# define functions

## a function to separate specified hybrids from other samples in a multiple sequence alignment
## It returns a list of two lists: a list of hybrid sequences, and a list of other sequences
def process_alignment(alignment, hybrid_list):
	hybrid_seqs = []
	other_seqs = []
	for sequence in alignment:
		if sequence.id in hybrid_list:
			hybrid_seqs.append(sequence)
		else:
			other_seqs.append(sequence)
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
		# determine if the column includes both possible bases for the ambiguity, and no other similar ambiguities
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
## The input "groups" are the various groupings to assess
## (a list of lists each with two lists of parents)
## It returns a list (sorted by count) of lists with parent combinations
## and their counts
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
	# now sort the parents by how often they group separately
	output_list.sort(key = lambda x: x[1], reverse = True)
	return(output_list)


## a function to create a report string
def parent_report(parent_list):
	if len(set(sum(parent_list, []))) > 10:		# more than 10 unique parents in top comparisons
		report = 'more than 10 (' + str(len(set(sum(parent_list, [])))) + ')'
	else:
		report = '/'.join([(item[0] + ';' + item[1]) for item in parent_list])
	return(report)




###########################
# Run
###########################


# get the hybrid sample IDs
hybrids = []
with open(hybrids_file, 'r') as infile:
	for line in infile:
		hybrids.append(line.strip())


# read in each alignment, extract hybrid sequences, retain non-hybrid samples in a new alignment,
# then process the hybrids and their ambiguities
for alignment_file in alignments:
	print('\nProcessing alignment ' + os.path.basename(alignment_file))
	sequences = process_alignment(AlignIO.read(open(alignment_file, 'r'), 'fasta'), hybrids)

	print('\n' + '\t'.join(['Hybrid', 'Ambiguities', 'Top_count', '2nd_count', '3rd_count',
		'Top_parents', '2nd_parents', '3rd_parents']))

	# turn the non-hybrid sequences into an alignment
	clean_alignment = MultipleSeqAlignment(sequences[1])

	# get the locations of and bases for ambiguities in the hybrid sequences
	hybrids_ambiguities = get_ambiguities(sequences[0])

	# for each hybrid, determine the parent groups for its valid ambiguities
	# ambiguities are valid if there are samples carrying both bases in that column
	# ambiguities are invalid if there is only one of the other bases (uninformative)
	for hybrid_entry in hybrids_ambiguities:
		# go through the hybrid's ambiguities and determine groups of parents
		groups = get_groups(hybrid_entry[1], clean_alignment)

		# determine how many parents in total
		tally = []
		for group in groups:
			for parent_list in group:
				tally = tally + parent_list

		parents = list(set(tally))

		if len(parents) < 2:	# not enough
			print('\t'.join([str(hybrid_entry[0]), str(len(hybrid_entry[1])), '-', '-', '-',
				'(none)', '(none)', '(none)']))
			continue		# move on to the next hybrid and its ambiguities
		else:
			# get counts of parents occurring in opposite groups
			opposite_count = get_opposites(parents, groups)

			# now evaluate how often the parents group separately
			count_vals = sorted(list(set([item[1] for item in opposite_count])), reverse = True)
			max_count = count_vals[0]
			top_parents = [item[0] for item in opposite_count if item[1] == max_count]
			top_report = parent_report(top_parents)

			if len(count_vals) > 1:
				second_count = count_vals[1]
				second_parents = [item[0] for item in opposite_count if item[1] == second_count]
				second_report = parent_report(second_parents)
			else:
				second_count = 0
				second_report = '(none)'

			if len(count_vals) > 2:
				third_count = count_vals[2]
				third_parents = [item[0] for item in opposite_count if item[1] == third_count]
				third_report = parent_report(third_parents)
			else:
				third_count = 0
				third_report = '(none)'

			# report
			if len(hybrid_entry[1]) < 4:		# if there are few heterozygous sites for the hybrid
				print('\t'.join([str(hybrid_entry[0]), str(len(hybrid_entry[1])), 'fewer than 4', '-', '-',
					'(none)', '(none)', '(none)']))
			else:
				if max_count > len(hybrid_entry[1])/2:		# more than half the ambiguities support the top comparisons
					top_count = str(max_count) + '*'
				else:
					top_count = str(max_count)

				print('\t'.join([str(hybrid_entry[0]), str(len(hybrid_entry[1])), top_count,
					str(second_count), str(third_count), top_report, second_report, third_report]))


# HOW TO EXTEND ACROSS ALIGNMENTS? PASS SCORES TO SETS OF PARENTS?
# TOP BY SIMPLE COUNT OF LOCI (TOP in EACH, WHEN LESS THAN 10)?
