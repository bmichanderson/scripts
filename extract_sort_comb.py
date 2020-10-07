#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 10 Mar 2020
# Modified: Oct 2020
# Description: take multifasta files as input, extract genes, sort based on name, combine all similar genes into single files as output
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files

def help():
	print('A script to extract specified regions and combine them across multifasta files extracted from genbank files.')
	print('The extracted regions are output as single multifasta files in the current directory.')
	print('It is assumed that no species are duplicated (one accession number per species).')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' options(-... -...) file1 file2 ...')
	print('')
	print('Options:')
	print('	-f,	File (required) with master list of regions for extracting (lower case), one per line')
	print('')
	print('	-d,	Dictionary (optional) of alternative region names for regions in the first file')
	print('		This should be in the form: "alternative_name	master_name", one per line')

options = ['-f', '-d']
opt_valdic = {}


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# process command line
for index, arg in enumerate(sys.argv[1:]):
	if arg[0] == '-': 	# if an option
		if arg in options:
			if sys.argv[1:][index + 1][0] == '-':
				sys.exit('No argument provided for ' + str(arg))
			opt_valdic[arg] = sys.argv[1:][index + 1]
		else:
			sys.exit('The argument ' + str(arg) + ' is not a valid option')


# capture list of files
file_list = []
for arg in sys.argv[1:]:
	if arg not in opt_valdic:
		file_list.append(arg)


# read in a list of regions for keeping
if '-f' in opt_valdic:
	gene_list = []
	with open(opt_valdic['-f'], 'r') as gl:
		for line in gl:
			gene_list.append(line.rstrip())


# read in a dictionary of region name alternatives
gene_name_dict = {}
if '-d' in opt_valdic:
	with open(opt_valdic['-d'], 'r') as dl:
		for line in dl:
			entry = str(line.rstrip().split()[0])
			equiv = str(line.rstrip().split()[1])
			gene_name_dict[entry] = equiv


# read in each extract file and sort and filter into list for later output
cap_list = []

if len(gene_list) > 0:
	for multifile in file_list:
		with open(multifile, 'r') as mf:
			fastas = SeqIO.parse(mf, 'fasta')
			for fasta in fastas:
				first_element = fasta.description.split()[0].lower()
				gene_element = first_element.split('_')[0]		# in case this is a _part gene listing
				gene_name = gene_name_dict.get(gene_element, gene_element)	# will simply return the gene_name if no dictionary hit
				if gene_name in gene_list:
					cap_list.append((gene_name, fasta))
else:
	sys.exit('Specify file with list of regions using the -f option')


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
