#!/usr/bin/env python
###############################################
# Author: B. Anderson
# Date: 28 Jun 2019
# Updated: Oct 2021 (can now take a text file with strings for string matching, one per line; also cleaned up and added parser)
# Updated: Jun 2022
# Description: a script to simply remove entries from multifastas based on string matches (arg1)
###############################################


import sys
import os
import argparse
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to remove fasta entries based on string matching. ' +
	'The new fasta files are output in the current directory with prefix \"mod_\".')


# add arguments to parse
parser.add_argument('strings', type = str, help = 'The comma delimited strings to search for ' +
	'when removing fastas from the multifastas (or nothing if file specified), ' +
	'and then the space-delimited multifastas', nargs = '*')
parser.add_argument('-f', type = str, dest = 'string_file', help = 'A text file with a list of strings, ' +
	'one per line (optional). If not specified, then the first space-delimited command line argument ' +
	'will be interpreted as the comma-delimited search strings.')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

strings = args.strings
string_file = args.string_file

if not strings:
	parser.print_help(sys.stderr)
	sys.exit(1)

string_list = []
file_list = []
if not string_file:			# if there is no file specified with '-f'
	search_strings = strings[0]
	for part in search_strings.split(','):
		string_list.append(part)
	for fasta in strings[1:]:
		file_list.append(fasta)
else:
	with open(string_file, 'r') as str_file:
		for line in str_file:
			string_list.append(line.rstrip())
	for fasta in strings:
		file_list.append(fasta)


# read in the each multifasta file and remove entries
for multifasta in file_list:
	fasta_list = []
	fastas = SeqIO.parse(open(multifasta, 'r'), 'fasta')
	for fasta in fastas:
		fasta_list.append(fasta)
	remove_indices = []
	for i, fasta in enumerate(fasta_list):
		for string in string_list:
			if string in fasta.description:
				remove_indices.append(i)
	for i in sorted(set(remove_indices), reverse = True):
		del fasta_list[i]
	print('Removed ' + str(len(set(remove_indices))) + ' entries from ' + os.path.basename(multifasta))
	with open('mod_' + os.path.basename(multifasta), 'w') as outfile:
		for fasta in fasta_list:
			SeqIO.write(fasta, outfile, 'fasta')
