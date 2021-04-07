#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 10 Mar 2020
# Modified: Oct 2020, April 2021 (updated argument handling and made more flexible to name of gene with _ in it)
# Description: take multifasta files as input, extract genes, sort based on name, combine all similar genes into single files as output
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files
import argparse

# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract specified regions and combine them across multifasta files ' +
						'extracted from genbank files. The extracted regions are output as single multifasta ' +
						'files in the current directory. It is assumed that no species are duplicated ' +
						'(one accession number per species).')


# add arguments to parse
parser.add_argument('extract_files', type=str, nargs='+', help='The multifasta files to extract from')
parser.add_argument('-f', type=str, dest='regions_file', help='File (required) with regions for extracting (one per line)')
parser.add_argument('-d', type=str, dest='dict_file', help='Dictionary (optional) of alternative region names in the form: ' +
								'"alt_name	region_name", one per line')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

extract_files = args.extract_files
regions_file = args.regions_file
dict_file = args.dict_file


# capture list of files
file_list = []
for file in extract_files:
	file_list.append(file)


# read in a list of regions for keeping
if regions_file:
	gene_list = []
	with open(regions_file, 'r') as gl:
		for line in gl:
			gene_list.append(line.rstrip().lower())
else:
	parser.print_help(sys.stderr)
	sys.exit(1)


# read in a dictionary of region name alternatives
gene_name_dict = {}
if dict_file:
	with open(dict_file, 'r') as dl:
		for line in dl:
			entry = str(line.rstrip().split()[0].lower())
			equiv = str(line.rstrip().split()[1].lower())
			gene_name_dict[entry] = equiv


# read in each extract file and sort and filter into list for later output
cap_list = []

if len(gene_list) > 0:
	for multifile in file_list:
		with open(multifile, 'r') as mf:
			fastas = SeqIO.parse(mf, 'fasta')
			for fasta in fastas:
				first_element = fasta.description.split()[0].lower()
				if '_part' in first_element:
					gene_element = first_element.split('_part')[0]		# in case this is a _part gene listing
				else:
					gene_element = first_element
				gene_name = gene_name_dict.get(gene_element, gene_element)	# will simply return the gene_name if no dictionary hit
				if gene_name in gene_list:
					cap_list.append((gene_name, fasta))
else:
	sys.exit('There is something wrong with the regions file you specified')


# combine all the regions into single files, removing exact duplicates and sorting
misses = []

for gene in gene_list:
	hits = []
	for locus in cap_list:
		if str(locus[0]) == str(gene):
			hits.append(locus[1])		# a fasta file is appended for each time the gene_name matches


	# find exact duplicates and remove them
	org_list = []
	seq_list = []
	remove_indices = []

	for num, hit in enumerate(hits):
		# assuming naming of regions is from genbank_parse.py with >gene from accession Genus specific_ep
		org = hit.description.split()[3] + '_' + hit.description.split()[4]
		if org in org_list:
			org_indices = [i for i in range(len(org_list)) if org_list[i] == org]
			found = 'False'
			for index in org_indices:
				if hit.seq == seq_list[index]:
					remove_indices.append(num)
					found = 'True'
					break

			if found == 'False':
				org_list.append(org)
				seq_list.append(hit.seq)

		else:
			org_list.append(org)
			seq_list.append(hit.seq)

	for index in sorted(remove_indices, reverse = True):
		del hits[index]


	# output the regions
	hits = sorted(hits, key = lambda x: x.description.split()[3])		# sort by genus

	if len(hits) > 0:
		with open(gene + '.fasta', 'w') as out_file:
			for hit in hits:
				SeqIO.write(hit, out_file, 'fasta')
	else:
		misses.append(gene)

if len(misses) > 0:
	print('No sequences found for ' + ', '.join(misses))
