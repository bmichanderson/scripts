#!/usr/bin/env python3

###############################################
# Author: B. Anderson
# Date: Aug 2023
# Description: parse the results of a series of RELAX runs to select optimal values
#	This assumes the results have already been tabulated in csv files, with a column of samples,
#	and three columns per gene corresponding to K, p and LogL of the RELAX alternative model
#	Arguments are the csv files (one per run; at least two), which should have the same unique taxa and genes
###############################################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to parse tabulated (csv file) RELAX outputs')


# add arguments to parse
parser.add_argument('results', type = str, nargs = '*', help = 'The csv result files; these need to have ' +
	'the same samples and genes in each.')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
results = args.results

if len(results) < 2:
	print('Please provide at least two run results files\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# Read in the results
gene_list = []
taxa = []
master_list = []
for fileno, result_file in enumerate(results):
	with open(result_file, 'r') as infile:
		for lineno, line in enumerate(infile):
			parts = line.strip().split(',')
			if lineno == 0:
				gene_entries = parts[1: ]
				if fileno == 0:		# first file
					for gene in gene_entries:
						if gene not in gene_list:
							gene_list.append(gene)
						else:
							continue
				elif len(gene_list) != (len(gene_entries) / 3):
					print('Wrong number of gene entries in ' + str(result_file))
					parser.print_help(sys.stderr)
					sys.exit(1)
			else:
				taxon = parts[0]
				if fileno == 0:		# first file
					taxa.append(taxon)
				elif taxon not in taxa:
					print('Wrong taxon in ' + str(result_file))
					parser.print_help(sys.stderr)
					sys.exit(1)
				
				master_list.append([taxon, parts[1: ]])


# Go through the samples in the list and choose the best values for each sample
# Determine the p-value significance and indicate it where appropriate
# Round the K value
# Count how many significant runs had K < 1 and how many had K > 1 (to assess bimodality)
out_list = []
for taxon in taxa:
	top_list = []
	result_list = [item[1] for item in master_list if item[0] == taxon]		# grab the runs for that taxon
	for index, gene in enumerate(gene_list):
		if result_list[0][index * 3] == '-':		# gene missing
			top_list.append(['-', '-'])
		else:
			Kvals = [float(item[index * 3 + 0]) for item in result_list]
			pvals = [float(item[index * 3 + 1]) for item in result_list]
			Lvals = [float(item[index * 3 + 2]) for item in result_list]
			# determine top likelihood and corresponding K and p for reporting
			max_index = Lvals.index(max(Lvals))
			if pvals[max_index] < 0.05:
				if pvals[max_index] < 0.01:
					if pvals[max_index] < 0.001:
						signif = '< 0.001'
						top_K = str(round(Kvals[max_index], 1)) + '***'
					else:
						signif = '< 0.01'
						top_K = str(round(Kvals[max_index], 1)) + '**'
				else:
					signif = '< 0.05'
					top_K = str(round(Kvals[max_index], 1)) + '*'
			else:
				signif = 'n.s.'
				top_K = str(round(Kvals[max_index], 1))
			# determine how many significant runs had divergent K values (< 1 and > 1)
			sig_indices = [i for i, pval in enumerate(pvals) if pval < 0.05]
			if len(sig_indices) == 0:		# no significant runs
				count_blurb = 'n.s.'
			else:
				sig_Kvals = [Kvals[i] for i in sig_indices]
				count_greater = len([item for item in sig_Kvals if item >= 1])
				count_less = len([item for item in sig_Kvals if item < 1])
				count_blurb = str(count_less) + '/' + str(count_greater) + \
					' of ' + str(count_greater + count_less)
			# record
			top_list.append([top_K, count_blurb])
	out_list.append([taxon, top_list])


# Print the output as csv (can be directed to a file with ">")
print('Sample,' + ','.join([(item + ',' + item) for item in gene_list]))
for entry in out_list:
	print(entry[0] + ',' + ','.join([','.join(item) for item in entry[1]]))
