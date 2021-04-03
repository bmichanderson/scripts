#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 27 Jan 2021
# Modified: April 2021 (make product addition more specific and based on input file)
# Description: add product to a features table tRNA or rRNA entry; remove fragment genes and optionally also remove specified entries
##########################

import sys			# allows access to command line arguments
import argparse			# allows argument parsing


# instantiate the parser
parser = argparse.ArgumentParser(description = 	'A script to add product information to a feature table for tRNA or rRNA entries, ' +
						'and optionally remove genes. Outputs a table as features_edited.tbl ' +
						'NOTE: this only works on feature tables with gene in each feature. ' +
						'NOTE: this also removes any features with \"fragment\" in their name.')


# add arguments to parse
parser.add_argument('feat_tbl', type=str, help='The feature table to add to')
parser.add_argument('-k', type=str, dest='key', help='The tab-delimited keys file with headings gene_name, formal_name, type, product')
parser.add_argument('-r', type=str, dest='remove', help='(Optional) A comma-delimited list of gene names to remove from the table')
parser.add_argument('-p', type=str, dest='prods', help='A file with genes to add product to, one gene name per line')


# parse the command line
if len(sys.argv[1:]) == 0:              # if there are no arguments
        parser.print_help(sys.stderr)
        sys.exit(1)

args = parser.parse_args()

feat_tbl = args.feat_tbl
key = args.key
remove = args.remove
prods = args.prods

if remove:
	remove_requested = True
	print('Will remove: ' + ' '.join(remove.split(',')))
else:
	remove_requested = False


# Read in the key file and create a dictionary for looking up name, type and product
if not key:
        parser.print_help(sys.stderr)
        sys.exit(1)

key_dict = {}
with open(key, 'r') as key_file:
	for line_no, line in enumerate(key_file, start=1):
		if line_no == 1:	# first line is the header row
			continue
		entry = line.rstrip().split('\t')
		key_dict[entry[0].lower()] = '\t'.join(entry[1:4])


# Read in the remove list if present
remove_list = []

if remove_requested:
	for item in remove.strip().split(','):
		remove_list.append(item.lower())


# Read in the prods list (genes which need a product added) if present
prods_list = []

if prods:
	with open(prods, 'r') as infile:
		for line in infile:
			gene_name = line.rstrip()
			prods_list.append(gene_name.lower())


# Read in the features table line by line, optionally remove lines, then read through remaining and update products
lines = []
remove_lines = []
to_remove = False
start = 0
add_product = False
with open(feat_tbl, 'r') as feat_file, open('features_edited.tbl', 'w') as out_file:
	for index, line in enumerate(feat_file):
		elements = line.rstrip().split('\t')
		if elements[0].startswith('>'):			# a new contig
			lines.append(line)
		elif len(elements) == 3:		# hit the next feature
			if to_remove:
				for num in range(start, index):
					remove_lines.append(num)
				to_remove = False
			start = index
			lines.append(line)
		elif len(elements) >= 4:			# properties of a feature
			if elements[3] == 'gene':
				if elements[4].lower() in remove_list:		# if this feature is supposed to be removed
					to_remove = True
				if 'fragment' in elements[4].lower():		# if this feature is a fragment (no product)
					print('Fragment: ' + elements[4].lower())
					to_remove = True
			lines.append(line)
		elif len(elements) == 2:			# coordinates
			lines.append(line)
		else:
			sys.exit('Problematic line! See ' + str(index+1))

	if len(remove_lines) > 0:
		print('Removing ' + str(len(remove_lines)) + ' lines that are specified or contain "fragment" in their name')
		for line in sorted(remove_lines, reverse = True):
			del lines[line]

	for line in lines:
		elements = line.rstrip().split('\t')
		if elements[0].startswith('>'):		# a new contig
			out_file.write(line)
		elif len(elements) == 3:				# start of a feature
			if any([elements[2] == 'tRNA', elements[2] == 'rRNA', elements[2] == 'CDS']):
#				add_product = True
				prod_feature = True
				out_file.write(line)
			else:
#				add_product = False
				prod_feature = False
				out_file.write(line)
		elif len(elements) >=4:				# properties of a feature
			if elements[3] == 'gene':
#				if add_product:
				if prod_feature:
					gene = elements[4].lower()
					if all([prods, gene in prods_list]):
						product = key_dict[gene].split('\t')[2]
						out_file.write(line)
						out_file.write('\t\t\tproduct\t' + product + '\n')
					else:
						out_file.write(line)
				else:
					out_file.write(line)
			else:
				out_file.write(line)
		elif len(elements) == 2:
			out_file.write(line)
		else:
			sys.exit('Data problem!!')


