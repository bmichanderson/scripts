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


# define an ambiguity dictionary and list
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
amb_list = ['K', 'k', 'M', 'm', 'R', 'r', 'S', 's', 'W', 'w', 'Y', 'y']


# get the hybrid sample IDs
hybrids = []
with open(hybrids_file, 'r') as infile:
	for line in infile:
		hybrids.append(line.strip())


# read in each alignment, extract hybrid sequences, retain non-hybrid samples in a new alignment
for alignment_file in alignments:
	print('\nProcessing alignment ' + os.path.basename(alignment_file))
	hybrid_seqs = []
	keep_seqs = []
	keep_samples = []
	for sequence in AlignIO.read(open(alignment_file, 'r'), 'fasta'):
		if sequence.id in hybrids:
			hybrid_seqs.append(sequence)
		else:
			keep_seqs.append(sequence)
			keep_samples.append(sequence.id)

	clean_alignment = MultipleSeqAlignment(keep_seqs)

	# for each hybrid in the alignment, determine where and what the ambiguities are
	eval_hybrids = []
	for hybrid_seq in hybrid_seqs:
		ambiguities = []
		for index, base in enumerate(hybrid_seq.seq):
			if base in amb_list:
				ambiguities.append([index, base])

		#print('Hybrid ' + str(hybrid_seq.id) + ' has ' + str(len(ambiguities)) + ' ambiguities')
		eval_hybrids.append([hybrid_seq.id, ambiguities])

	# for each hybrid to evaluate, create sample groups at informative heterozygous sites
	print('\n' + '\t'.join(['Hybrid', 'Ambiguities', 'Top_count', '2nd_count', '3rd_count',
		'Top_parents', '2nd_parents', '3rd_parents', 'Unique_parents']))
	for entry in eval_hybrids:
		parents = []
		groups_list = []
		for ambiguity in entry[1]:
			amb_base = ambiguity[1]
			# slice the position in the alignment
			this_column = clean_alignment[:, ambiguity[0]]		# slice all rows (samples), one position
			# determine if the column includes both possible bases for the ambiguity, and no other similar ambiguities
			if amb_base in set(this_column):
				#print('Found other instances of that ambiguous base')
				continue
			elif all([amb_dict[amb_base][0] in set(this_column), amb_dict[amb_base][1] in set(this_column)]):
				# collect the two groups of parents
				parents1 = []
				parents2 = []
				for index, base in enumerate(this_column):
					if base == amb_dict[amb_base][0]:
						parents1.append(keep_samples[index])
					elif base == amb_dict[amb_base][1]:
						parents2.append(keep_samples[index])
					else:
						#print('This position does not match one of the alternatives')
						continue
				
				groups_list.append([parents1, parents2])
				parents = list(set(parents + parents1 + parents2))
			else:
				#print('Could not find both alternative bases in the alignment')
				continue

		# determine when parents occur in opposite groups
		opposite_count = []
		if len(parents) < 2:	# not enough parents recovered to report
			print('\t'.join([str(entry[0]), str(len(entry[1])), '-', '-', '-',
				'(none)', '(none)', '(none)', '(none)']))
			continue
		for index, parent in enumerate(parents):	# cycle through each parent sample
			if index < (len(parents) - 1):	# if at least two more parents to compare
				for parent2 in parents[index + 1: ]:
					count = 0
					for group in groups_list:
						if any([
								all([parent in group[0], parent2 in group[1]]),
								all([parent in group[1], parent2 in group[0]])
							]):		# different groups
							count = count + 1
						else:
							continue
					
					opposite_count.append([[parent, parent2], count])

		# now sort the parents by how often they group separately
		opposite_count.sort(key = lambda x: x[1], reverse = True)
		count_vals = sorted(list(set([item[1] for item in opposite_count])), reverse = True)
		max_count = count_vals[0]

		## make a function to create a report string
		def parent_report(parent_list):
			if len(set(sum(parent_list, []))) > 10:		# more than 10 unique parents in top comparisons
				report = 'more than 10 (' + str(len(set(sum(parent_list, [])))) + ')'
			else:
				report = '/'.join([(item[0] + ';' + item[1]) for item in parent_list])
			return(report)

		# assess the sorted parent combos
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
		if len(entry[1]) < 4:		# if there are few heterozygous sites for the hybrid
			print('\t'.join([str(entry[0]), str(len(entry[1])), 'fewer than 4', '-', '-',
				'(none)', '(none)', '(none)', '(none)']))
		else:
			if max_count > len(entry[1])/2:		# more than half the ambiguities support the top comparisons
				top_count = str(max_count) + '*'
			else:
				top_count = str(max_count)

			print('\t'.join([str(entry[0]), str(len(entry[1])),	top_count,
				str(second_count), str(third_count), top_report, second_report, third_report]))



# HOW TO EXTEND ACROSS ALIGNMENTS? PASS SCORES TO SETS OF PARENTS?
# TOP BY SIMPLE COUNT OF LOCI (TOP in EACH, WHEN LESS THAN 10)?
