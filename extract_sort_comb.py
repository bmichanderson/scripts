#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: 10 Mar 2020
# Modified: Oct 2020, April 2021, Aug 2024
# Description: take multifasta files as input (e.g. output from genbank_parse.py), extract genes, sort based on sample name, 
#	remove exact duplicates, and combine all similar genes into single files as output
# Note: this script expects the fasta entries in the multifastas to have a format similar to the output from genbank_parse.py,
#	e.g. >genename from <accession> <genus> <species>
#	or to have an analogous format of at least three fields
#	e.g. >genename from <field>
#	If it detects that there are fewer than five fields, it will sort based on the third, otherwise based on genus
##########################


import argparse
from Bio import SeqIO
import sys


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract specified regions and combine them across multifasta files ' +
	'extracted from genbank files. The extracted regions are output as single multifasta files in the current directory.')


# add arguments to parse
parser.add_argument('extract_files', type=str, nargs='*', help='The multifasta files to extract regions from')
parser.add_argument('-f', type=str, dest='regions_file', help='File (optional) with regions for combining (one per line)')
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
gene_list = []
if regions_file:
	with open(regions_file, 'r') as gl:
		for line in gl:
			gene_list.append(line.rstrip().lower())
else:
	print('Regions file not specified, so will include all regions\n')


# read in a dictionary of region name alternatives
gene_name_dict = {}
if dict_file:
	with open(dict_file, 'r') as dl:
		for line in dl:
			pieces = line.strip().split()
			entry = pieces[0].lower()
			equiv = pieces[1].lower()
			gene_name_dict[entry] = equiv


# read in each extract file and sort and filter into list for later output
cap_list = []
species_present = 0
species_absent = 0
for file in file_list:
	with open(file, 'r') as mf:
		fastas = SeqIO.parse(mf, 'fasta')
		for fasta in fastas:
			description_parts = fasta.description.split()
			if len(description_parts) > 4:		# original naming convention: >gene from accession genus species
				species_present = species_present + 1
			elif len(description_parts) > 2:	# there's a third field to sort on (sample)
				species_absent = species_absent + 1
			else:	# wrong naming convention
				print('Incorrect fasta description formatting\n')
				parser.print_help(sys.stderr)
				sys.exit(1)
			first_element = description_parts[0].lower()
			if '_part' in first_element:
				gene_element = first_element.split('_part')[0]		# in case this is a _part gene listing
			else:
				gene_element = first_element
			gene_name = gene_name_dict.get(gene_element, gene_element)	# will simply return the gene_name if no dictionary hit
			if regions_file:
				if gene_name in gene_list:
					cap_list.append((gene_name, fasta))
			else:
				if gene_name in gene_list:
					cap_list.append((gene_name, fasta))
				else:
					gene_list.append(gene_name)
					cap_list.append((gene_name, fasta))

if species_present > 0:
	if species_absent > 0:
		print('Inconsistent fasta description lengths\n')
		parser.print_help(sys.stderr)
		sys.exit(1)
	else:
		use_genus = True
elif species_absent > 0:
	use_genus = False
else:
	print('Incorrect fasta description formatting\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# combine all the regions into single files, removing exact duplicates and sorting
misses = []
for gene in gene_list:
	hits = []
	for locus in cap_list:
		if str(locus[0]) == str(gene):
			hits.append(locus[1])		# a fasta file is appended for each time the gene_name matches

	# find exact duplicates (based on third field) and remove them
	sample_list = []
	seq_list = []
	remove_indices = []

	for num, hit in enumerate(hits):
		sample = hit.description.split()[2]
		if sample in sample_list:
			sample_indices = [i for i in range(len(sample_list)) if sample_list[i] == sample]
			found = 'False'
			for index in sample_indices:
				if hit.seq == seq_list[index]:
					remove_indices.append(num)
					found = 'True'
					break

			if found == 'False':
				sample_list.append(sample)
				seq_list.append(hit.seq)

		else:
			sample_list.append(sample)
			seq_list.append(hit.seq)

	if len(remove_indices) > 0:
		for index in sorted(remove_indices, reverse = True):
			del hits[index]
		print('Removed ' + str(len(remove_indices)) + ' entries for gene ' + str(gene))

	# output the regions
	if use_genus:
		hits = sorted(hits, key = lambda x: x.description.split()[3])
	else:
		hits = sorted(hits, key = lambda x: x.description.split()[2])

	if len(hits) > 0:
		with open(gene + '.fasta', 'w') as out_file:
			for hit in hits:
				SeqIO.write(hit, out_file, 'fasta')
		print('Output ' + str(len(hits)) + ' entries for gene ' + str(gene))
	else:
		misses.append(gene)


print('\nOutput generated for ' + str(len(gene_list) - len(misses)) + ' regions')
if len(misses) > 0:
	print('No sequences found for ' + ', '.join(misses))
