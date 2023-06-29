#!/usr/bin/env python3

#####################
# Author: B. Anderson
# Date: May 2021
# Modified: June 2023
# Description: Retain a subset (-s list.txt) of fasta samples from multifastas
#####################


import sys
import os
import argparse
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to retain a subset (-s list.txt) of fasta samples from multifasta inputs.' +
	' Output files will have the same name but starting with \"subset_\".' +
	' NOTE: it is up to you to ensure the subset names only match once per multifasta.' +
	' By default, the taxon names should be the 4th and 5th fields in the descriptions.' +
	' To simply match the id, use \"-i yes\".')


# add arguments to parse
parser.add_argument('multifastas', type=str, nargs='+', help='The multifasta files to be filtered')
parser.add_argument('-s', type=str, dest='subset', help='File (required) with names of samples to match EXACTLY in the multifastas')
parser.add_argument('-i', type=str, dest='idmatch', help='Whether (\"yes\" or \"no\") to match the ID rather than the original design [default no]')


# parse the command line
if len(sys.argv[1:]) == 0:		# no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
multifastas = args.multifastas
subset = args.subset
idmatch = args.idmatch

if not idmatch:
	idmatch = 'no'


# capture list of files for filtering
file_list = []
for multifasta in multifastas:
	file_list.append(multifasta)


# read in the list of samples to keep
if subset:
	sample_list = []
	with open(subset, 'r') as sub_list:
		for line in sub_list:
			if idmatch == 'yes':
				sample_list.append(line.strip())
			else:
				sample_list.append(" ".join(line.rstrip().split()))
else:
	parser.print_help(sys.stderr)
	sys.exit(1)


# filter each multifasta file
for multifasta in file_list:
	omit_count = 0
	retained_count = 0
	with open(multifasta, 'r') as infile, open('subset_' + os.path.basename(multifasta), 'w') as outfile:
		fastas = SeqIO.parse(infile, 'fasta')
		for fasta in fastas:
			if idmatch == 'yes':
				name = fasta.id
			else:
				name = " ".join(fasta.description.split()[3:5])

			if name in sample_list:
				SeqIO.write(fasta, outfile, 'fasta')
				retained_count = retained_count + 1
			else:
				omit_count = omit_count + 1
	print('Retained ' + str(retained_count) + ' and omitted ' + str(omit_count) + ' entries from ' + multifasta)
