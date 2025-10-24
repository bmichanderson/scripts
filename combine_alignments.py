#!/usr/bin/env python3

################
# Author: B.M. Anderson
# Date: 15 Aug 2020
# Modified: 30 Jan 2021 (added parser and flexibility for header info), April 2021 (adding support for generating a partitions file)
#	Sep 2021 (added additional naming convention for broader applicability and made region names more flexible and a report for many alignments)
#	Apr 2025 (added arg for search strings of taxa to drop; improved speed and efficiency for large numbers of alignments)
#	Oct 2025 (changed file naming for when files aren't in the starting directory; fixed bug for last locus missing)
# Description: combine fasta alignments for a single concatenated version with all desired taxa
################


import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to concatenate multiple sequence alignments, inserting dashes for missing taxa ' +
	'and optionally to create a nexus format partitions file. ' +
	'Files are assumed to be named: align_{region}_... in the original [default] naming convention. ' +
	'If named otherwise, file names (minus the suffix) are used for region names in any partition file created.')


# add arguments to parse
parser.add_argument('alignments', type=str, help='The multiple sequence alignments', nargs='*')
parser.add_argument('-f', type=str, dest='format',
	help='Specify entry name format in the alignment files: ' +
	'original [default] = space delimited, with taxa as fourth and fifth field; ' +
	'simple = underscore delimited, with taxa as first and second field; ' +
	'single = no delimitation, a single unique identifier per individual')
parser.add_argument('-p', type=str, dest='partitions', help='Create a partitions file, \"locus\" [default] or ' +
	'\"cds\" for by codon position assuming regions begin with position 1')
parser.add_argument('-s', type=str, dest='strings', help='A text file with strings to search against ' +
	'taxa to specify which to ignore when combining alignments')


# parse the command line
if len(sys.argv[1:]) == 0:			# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

alignments = args.alignments		# a list

format = args.format
format_original = False
format_simple = False
format_single = False
if format:
	if format.lower() == 'simple':
		format_simple = True
	elif format.lower() == 'single':
		format_single = True
	else:
		format_original = True
else:
	format_original = True

partitions = args.partitions
if partitions:
	if partitions.lower() == 'locus':
		cds = False
	elif partitions.lower() == 'cds':
		cds = True
	else:
		parser.print_help(sys.stderr)
		sys.exit(1)

strings = args.strings
if strings:
	str_list = []
	with open(strings, 'r') as str_file:
		for line in str_file:
			str_list.append(line.strip())

files = []
for file in alignments:
	files.append(file)


# read each file and capture data
entry_list = []
align_num = 1
regions = []
for multifasta in files:
	filename = os.path.basename(multifasta)
	with open(multifasta, 'r') as align:
		fastas = SeqIO.parse(align, 'fasta')
		len_list = []
		for fasta in fastas:
			entry_list.append([align_num, fasta])
			if len(fasta) not in len_list:
				len_list.append(len(fasta))
		if len(len_list) > 1:
			sys.exit('Problem with aligned sequence lengths for ' + multifasta)
		else:
			if format_original:
				regions.append([align_num, filename.split('_')[1], len_list[0]])	# assuming files are named: align_region_...
			else:
				regions.append([align_num, os.path.splitext(filename)[0], len_list[0]])	# strips out extension and uses file name for region
	align_num = align_num + 1


# determine all the taxa present
taxa_list = []
for entry in entry_list:
	head = entry[1].description.split()
	if format_simple or format_single:
		if head[0] in taxa_list:
			continue
		else:
			taxa_list.append(head[0])
	else:		# format_original
		acc = head[2]
		if (head[3] + ' ' + head[4]) in (x[0] for x in taxa_list):
			continue
		else:
			taxa_list.append([head[3] + ' ' + head[4], acc])

taxa_list.sort()
print('Detected ' + str(len(taxa_list)) + ' taxa present across ' + str(align_num - 1) + ' alignments')


# if a strings file was provided, drop taxa that match
if strings:
	keep_list = []
	drop_list = []
	for taxon in taxa_list:
		found = False
		for string in str_list:
			if str(string) == str(taxon):
				drop_list.append(taxon)
				found = True
				break
		if not found:
			keep_list.append(taxon)

	count = len(taxa_list) - len(keep_list)
	if count != len(drop_list):
		print('PROBLEM with taxa')
		sys.exit(1)
	if count > 0:
		print('Ignoring ' + str(count) + ' taxa based on provided list of strings')

	taxa_list = keep_list


# for each region, create fasta files for missing taxa that are "-" and as long as the others
align_count = 0
filtered_list = []
num = 1		# start at the first alignment
length = 0
count = 0
taxa_minilist = []
missing = []
for entry in entry_list:
	if entry[0] == num:		# continue evaluating the alignment for this next entry
		length = len(entry[1])
		if format_original:
			head = entry[1].description.split()
			taxon = head[3] + ' ' + head[4]
		else:
			taxon = entry[1].description.split()[0]		# just in case there is extra text
		if taxon not in taxa_list:		# because it is set to be ignored
			continue
		if taxon in taxa_minilist:
			print('\nNOTE!\n' + taxon + ' is represented more than once for alignment ' + str(num) + ': ' +
				regions[num-1][1] + '\n')
			continue
		else:
			taxa_minilist.append(taxon)
			filtered_list.append([taxon, entry])

	else:					# finished an alignment, now fill in missing sequences and start the next
		if len(taxa_minilist) != len(taxa_list):		# some missing taxa
			for tentry in taxa_list:
				if format_original:
					if tentry[0] in taxa_minilist:
						continue
					else:
						missing.append(tentry[0])
						new_fasta = SeqRecord(Seq('-' * length), id = 'gaps', name = 'gaps', description = 'gaps from ' +
							tentry[1] + ' ' + tentry[0])
						filtered_list.append([tentry[0], [num, new_fasta]])
						count = count + 1
				else:
					if tentry in taxa_minilist:
						continue
					else:
						missing.append(tentry)
						new_fasta = SeqRecord(Seq('-' * length), id = tentry, name = tentry, description = tentry)
						filtered_list.append([tentry, [num, new_fasta]])
						count = count + 1

		if len(missing) > 0 and not format_single:
			print('Added ' + str(count) + ' gap entries for missing taxa for alignment ' + str(num) + ': ' + regions[num-1][1])
			print('Missing taxa are: ' + '\t'.join(missing))

		align_count = align_count + 1
		if align_count % 200 == 0:
			print('Evaluated ' + str(align_count) + ' alignments of ' + str(align_num - 1) + ' alignments')

		# reset vars and start on the next
		num = entry[0]
		count = 0
		taxa_minilist = []
		missing = []

		# also catch the info from the current entry
		length = len(entry[1])
		if format_original:
			head = entry[1].description.split()
			taxon = head[3] + ' ' + head[4]
		else:
			taxon = entry[1].description.split()[0]
		if taxon not in taxa_list:		# because it is set to be ignored
			continue
		if taxon in taxa_minilist:
			print('\nNOTE!\n' + taxon + ' is represented more than once for alignment ' + str(num) + ': ' +
				regions[num-1][1] + '\n')
			continue
		else:
			taxa_minilist.append(taxon)
			filtered_list.append([taxon, entry])


# after completing the pass of the entry list, finish processing the last alignment
if len(taxa_minilist) != len(taxa_list):		# some missing taxa
	for tentry in taxa_list:
		if format_original:
			if tentry[0] in taxa_minilist:
				continue
			else:
				missing.append(tentry[0])
				new_fasta = SeqRecord(Seq('-' * length), id = 'gaps', name = 'gaps', description = 'gaps from ' +
					tentry[1] + ' ' + tentry[0])
				filtered_list.append([tentry[0], [num, new_fasta]])
				count = count + 1
		else:
			if tentry in taxa_minilist:
				continue
			else:
				missing.append(tentry)
				new_fasta = SeqRecord(Seq('-' * length), id = tentry, name = tentry, description = tentry)
				filtered_list.append([tentry, [num, new_fasta]])
				count = count + 1

if len(missing) > 0 and not format_single:
	print('Added ' + str(count) + ' gap entries for missing taxa for alignment ' + str(num) + ': ' + regions[num-1][1])
	print('Missing taxa are: ' + '\t'.join(missing))

align_count = align_count + 1
if align_count % 200 == 0:
	print('Evaluated ' + str(align_count) + ' alignments of ' + str(align_num - 1) + ' alignments')


# concatenate regions for each taxon
out_list = []
if format_original:
	for taxon in taxa_list:
		sublist = [entry[1] for entry in filtered_list if entry[0] == taxon[0]]
		sort_sublist = sorted(sublist, key = lambda x: x[0])
		new_fasta = SeqRecord(Seq(""), id = 'concat', name = 'concat', description = 'concat from ' + taxon[1] + ' ' + taxon[0])
		for entry in sort_sublist:
			new_fasta.seq = new_fasta.seq + entry[1].seq
		out_list.append([taxon[0], new_fasta])
else:
	for taxon in taxa_list:
		sublist = [entry[1] for entry in filtered_list if entry[0] == taxon]
		sort_sublist = sorted(sublist, key = lambda x: x[0])
		new_fasta = SeqRecord(Seq(""), id = taxon, name = taxon, description = taxon + ' concat')
		for entry in sort_sublist:
			new_fasta.seq = new_fasta.seq + entry[1].seq
		out_list.append([taxon, new_fasta])


# write the final concatenated fastas into a single file
with open('combine_out.fasta', 'w') as outfile:
	for outf in out_list:
		SeqIO.write(outf[1], outfile, 'fasta')


# if a partitions file was requested, create it
if partitions:
	running_pos = 1
	with open('combine_partitions.nex', 'w') as partfile:
		partfile.write('#nexus\nbegin sets;\n')
		if cds:
			for region in regions:
				partfile.write('\tcharset ' + region[1] + '_1 = ' + str(running_pos) + '-' +
						str(running_pos + region[2] - 1) + '\\3;\n')
				partfile.write('\tcharset ' + region[1] + '_2 = ' + str(running_pos + 1) + '-' +
						str(running_pos + region[2] - 1) + '\\3;\n')
				partfile.write('\tcharset ' + region[1] + '_3 = ' + str(running_pos + 2) + '-' +
						str(running_pos + region[2] - 1) + '\\3;\n')
				running_pos = running_pos + region[2]
		else:
			for region in regions:
				partfile.write('\tcharset ' + region[1] + ' = ' + str(running_pos) + '-' +
						str(running_pos + region[2] - 1) + ';\n')
				running_pos = running_pos + region[2]
		partfile.write('end;\n')




###################
# OLD VERSION
#for num in range(1, align_num):
#	align_count = align_count + 1
#	if align_count % 200 == 0:
#		print('Evaluated ' + str(align_count) + ' alignments of ' + str(align_num - 1) + ' alignments')
#	taxa_minilist = []
#	length = 0
#	count = 0
#	missing = []
#	found = False
#	for entry in entry_list:
#		if all([found, entry[0] != num]):
#			break
#		if entry[0] == num:
#			found = True
#		else:
#			continue
#		if length == 0:
#			length = len(entry[1])
#		if format_original:
#			head = entry[1].description.split()
#			taxon = head[3] + ' ' + head[4]
#		else:
#			taxon = entry[1].description.strip()
#		if taxon not in taxa_list:		# because it is set to be ignored
#			continue
#		if taxon in taxa_minilist:
#			print('\nNOTE!\n' + taxon + ' is represented more than once for alignment ' + str(num) + ': ' +
#				regions[num-1][1] + '\n')
#			continue
#		else:
#			taxa_minilist.append(taxon)
#			filtered_list.append([taxon, entry])
#	if len(taxa_minilist) == len(taxa_list):		# no missing
#		continue
#	else:
#		for tentry in taxa_list:
#			if format_original:
#				if tentry[0] in taxa_minilist:
#					continue
#				else:
#					missing.append(tentry[0])
#					new_fasta = SeqRecord(Seq('-' * length), id = 'gaps', name = 'gaps', description = 'gaps from ' +
#						tentry[1] + ' ' + tentry[0])
#					filtered_list.append([num, new_fasta])
#					count = count + 1
#			else:
#				if tentry in taxa_minilist:
#					continue
#				else:
#					missing.append(tentry)
#					new_fasta = SeqRecord(Seq('-' * length), id = tentry, name = tentry, description = tentry)
#					filtered_list.append([num, new_fasta])
#					count = count + 1
#
#		if len(missing) > 0 and not format_single:
#			print('Added ' + str(count) + ' gap entries for missing taxa for alignment ' + str(num) + ': ' + regions[num-1][1])
#			print('Missing taxa are: ' + '\t'.join(missing))
