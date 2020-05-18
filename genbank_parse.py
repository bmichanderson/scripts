#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 18 May 2020
# Description: Parse a genbank file to extract information
##########################


import sys			# allows access to command line arguments
import argparse
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract specified regions from a genbank/multigenbank file or list features present')


# add arguments to parse
parser.add_argument('gb_file', type=str, help='A genbank/multigenbank file to parse')
parser.add_argument('-f', type=str, dest='regions_file', help='File with regions for extracting (lower case, one per line)')
parser.add_argument('-t', type=str, dest='type', help='The type of extract output: nucl [default] or prot (for use with the -f option)')
parser.add_argument('-l', type=str, dest='list', help='Specify a specific type of feature to list: CDS, rRNA, or tRNA')
parser.add_argument('-g', type=str, dest='gene', help='Specify the name of a gene to extract the sequence(s)')
parser.add_argument('-c', type=str, dest='coords', help='Specify the coordinates of a sequence to extract: start..end, with start > end for compliment')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

gb_file = args.gb_file
regions_file = args.regions_file
type_extract = args.type
list_type = args.list
gene_name = args.gene
coords = args.coords


# if -l is present, run the listing and exit

if list_type:
	if list_type.lower() == 'cds':
		feature_type = 'CDS'
	elif list_type.lower() == 'rrna':
		feature_type = 'rRNA'
	elif list_type.lower() == 'trna':
		feature_type = 'tRNA'
	else:
		sys.exit('Specify a feature type to list as CDS, rRNA or tRNA')

	print_list = []
	gbks = SeqIO.parse(gb_file, 'genbank')
	for gbk in gbks:
		if any((feature_type == 'CDS', feature_type == 'rRNA')):
			for feature in gbk.features:
				if feature.type == feature_type:
					if 'gene' in feature.qualifiers:
						if 'pseudo' in feature.qualifiers:
							print_list.append(''.join(feature.qualifiers['gene'][0].lower().split()) + '-pseudo')
						else:
							print_list.append(''.join(feature.qualifiers['gene'][0].lower().split()))	# for removing spaces
					elif 'product' in feature.qualifiers:
						print_list.append(''.join(feature.qualifiers['product'][0].lower().split()))
		# tRNA
		elif feature_type == 'tRNA':
			for feature in gbk.features:
				if feature.type == feature_type:
					if 'gene' in feature.qualifiers:
						name = ''.join(feature.qualifiers['gene'][0].lower().split())
					elif 'product' in feature.qualifiers:
						name = ''.join(feature.qualifiers['product'][0].lower().split())
					else:
						continue

					if '(' in name:
							mod_name = name.replace('(', '-')
							mod_name = mod_name.replace(')', '')
					else:
							mod_name = name

					if 'pseudo' in feature.qualifiers:
						print_list.append(mod_name + '-pseudo')
						continue

					if 'anticodon' in feature.qualifiers:
						anti_list = list(feature.qualifiers['anticodon'][0].rstrip(')').split(','))
						anti_codon = anti_list[-1].split(':')[-1]

						bases = list(('a', 'c', 'g', 't', 'u'))
						if not all((mod_name[-2:][0] in bases, mod_name[-2:][1] in bases)):		# if the end of the name isn't already a codon
							if mod_name[-3:] == '-cp':
								mod_name = mod_name[:-3] + '-' + anti_codon + '-cp'
							else:
								mod_name = mod_name + '-' + anti_codon
					elif 'note' in feature.qualifiers:
						note = str(feature.qualifiers['note'])
						if 'anticodon' in note:
							start = note.find('anticodon:')
							anti_codon = note[start + 10: start + 13].lower()

							bases = list(('a', 'c', 'g', 't', 'u'))
							if not all((mod_name[-2:][0] in bases, mod_name[-2:][1] in bases)):		# if the end of the name isn't already a codon
								if mod_name[-3:] == '-cp':
									mod_name = mod_name[:-3] + '-' + anti_codon + '-cp'
								else:
									mod_name = mod_name + '-' + anti_codon

					if 'note' in feature.qualifiers:
						if all(('plast' in str(feature.qualifiers['note']), '-cp' not in mod_name)):
							mod_name = mod_name + '-cp'

					print_list.append(mod_name)

	for gene in sorted(print_list):
		print(gene)
	sys.exit()



# If -g or -c are present, extract the region and exit

if any([gene_name, coords]):
	gbks = []
	genbanks = SeqIO.parse(gb_file, 'genbank')
	for gb in genbanks:
		gbks.append(gb)
	with open(gbks[0].annotations['organism'].split()[0] + '_extract.fasta', 'w') as out_file:
		if gene_name:
			for gbk in gbks:
				for feature in gbk.features:
					if feature.type == 'gene':
						name = ''.join(feature.qualifiers['gene'][0].lower().split())
						if name == gene_name.lower():
							out_file.write(">%s from %s\n%s\n" % (name, gbk.name + ' ' + gbk.annotations['organism'],
												feature.location.extract(gbk).seq))

		elif coords:
			if len(gbks) > 1:		# if this is a multigenbank file, coords shouldn't work
				sys.exit('Cannot specify coordinates to extract from a multigenbank file')
			elif '..' in coords:
				start = int(coords.split('..')[0])
				end = int(coords.split('..')[1])
				if start < end:
					out_file.write(">%s from %s\n%s\n" % ('sequence_' + coords, gbks[0].name + ' ' + gbks[0].annotations['organism'],
										gbks[0].seq[start-1:end]))
				elif start > end:	# opposite strand format
					out_file.write(">%s from %s\n%s\n" % ('sequence_' + coords, gbks[0].name + ' ' + gbks[0].annotations['organism'],
										gbks[0].seq[end-1:start].reverse_complement()))
				else:
					sys.exit('Ensure the start and end coordinates are not the same')
			else:
				sys.exit('Coordinates specified incorrectly. Need to be in the form start..end')

	sys.exit()


# If a file is provided to extract named regions, proceed

if regions_file:
	genes = []
	with open(regions_file, 'r') as gene_list:
		for gene in gene_list:
			genes.append(gene.rstrip())
	genes = list(set(genes))		# remove duplicates

	if type_extract:
		if type_extract.lower() == 'prot':
			type_extract = 'prot'
		else:
			type_extract = 'nucl'
	else:
		type_extract = 'nucl'

	gbks = []
	genbanks = SeqIO.parse(gb_file, 'genbank')
	for gb in genbanks:
		gbks.append(gb)
	with open(gbks[0].annotations['organism'].split()[0] + '_' + str(type_extract) + '_extract.fasta', 'w') as out_file:
		copy_num = {}
		ref_needed = 'False'		# a flag for when a multi-record genbank has features referencing different sequences
		multi_parts = []
		for gbk in gbks:
			for feature in gbk.features:

				# CDS
				if all((feature.type == 'CDS', 'gene' in feature.qualifiers)):

					if 'pseudo' in feature.qualifiers:
						name = ''.join(feature.qualifiers['gene'][0].lower().split()) + '-pseudo'
					else:
						name = ''.join(feature.qualifiers['gene'][0].lower().split())

					if name in genes:
						if name in copy_num:		# if this region has multiple copies
							copy_num[name] = copy_num[name] + 1
							print('Duplicate feature: ' + name + ' detected')
						else:
							copy_num[name] = 1

						if copy_num[name] > 1:
							name = name + '_' + str(copy_num[name])

						if type_extract == 'prot':
							if 'translation' in feature.qualifiers:
								out_file.write(">%s from %s\n%s\n" % (name, gbk.name + ' ' + gbk.annotations['organism'],
													feature.qualifiers['translation'][0]))
						elif type_extract == 'nucl':
							ref_present = 'False'
							for part in feature.location.parts:
								if part.ref:		# if there is a reference to another sequence
									ref_present = 'True'
									ref_needed = 'True'		# this flag will now trigger another pass over the list

							if ref_present == 'True':
								parts = []
								for part in feature.location.parts:		# note that the order of parts is important esp. if trans
									if part.ref:
										parts.append(part)
									else:		# the part is found in this genbank, so ref = None; need to add
										part.ref = gbk.id
										parts.append(part)
								multi_parts.append((name, parts))
							else:
								out_file.write(">%s from %s\n%s\n" % (name, gbk.name + ' ' + gbk.annotations['organism'],
													feature.location.extract(gbk).seq))
						else:
							sys.exit('Type of extraction (nucl or prot) specified incorrectly')


				# rRNA
				elif (feature.type == 'rRNA'):

					if 'gene' in feature.qualifiers:
						name = ''.join(feature.qualifiers['gene'][0].lower().split())

						if name in genes:
							out_file.write(">%s from %s\n%s\n" % (name, gbk.name + ' ' + gbk.annotations['organism'],
												feature.location.extract(gbk).seq))

					elif 'product' in feature.qualifiers:
						name = ''.join(feature.qualifiers['product'][0].lower().split())

						if name in genes:
							out_file.write(">%s from %s\n%s\n" % (name, gbk.name + ' ' + gbk.annotations['organism'],
												feature.location.extract(gbk).seq))


				# tRNA
				elif feature.type == 'tRNA':

					if 'gene' in feature.qualifiers:
						name = ''.join(feature.qualifiers['gene'][0].lower().split())
					elif 'product' in feature.qualifiers:
						name = ''.join(feature.qualifiers['product'][0].lower().split())
					else:
						continue

					if '(' in name:
							mod_name = name.replace('(', '-')
							mod_name = mod_name.replace(')', '')
					else:
							mod_name = name

					if 'pseudo' in feature.qualifiers:
						mod_name = mod_name + '-pseudo'
						if mod_name in genes:
							out_file.write(">%s from %s\n%s\n" % (mod_name, gbk.name + ' ' + gbk.annotations['organism'],
												feature.location.extract(gbk).seq))

					if 'anticodon' in feature.qualifiers:
						anti_list = list(feature.qualifiers['anticodon'][0].rstrip(')').split(','))
						anti_codon = anti_list[-1].split(':')[-1]

						bases = list(('a', 'c', 'g', 't', 'u'))
						if not all((mod_name[-2:][0] in bases, mod_name[-2:][1] in bases)):		# if the end of the name isn't already a codon
							if mod_name[-3:] == '-cp':
								mod_name = mod_name[:-3] + '-' + anti_codon + '-cp'
							else:
								mod_name = mod_name + '-' + anti_codon
					elif 'note' in feature.qualifiers:
						note = str(feature.qualifiers['note'])
						if 'anticodon' in note:
							start = note.find('anticodon:')
							anti_codon = note[start + 10: start + 13].lower()

							bases = list(('a', 'c', 'g', 't', 'u'))
							if not all((mod_name[-2:][0] in bases, mod_name[-2:][1] in bases)):		# if the end of the name isn't already a codon
								if mod_name[-3:] == '-cp':
									mod_name = mod_name[:-3] + '-' + anti_codon + '-cp'
								else:
									mod_name = mod_name + '-' + anti_codon

					if 'note' in feature.qualifiers:
							if all(('plast' in str(feature.qualifiers['note']), '-cp' not in mod_name)):
								mod_name = mod_name + '-cp'

					if mod_name in genes:
						out_file.write(">%s from %s\n%s\n" % (mod_name, gbk.name + ' ' + gbk.annotations['organism'],
											feature.location.extract(gbk).seq))

		if all((ref_needed == 'True', type_extract == 'nucl')):			# if a multi-record genbank had a CDS feature referencing another sequence
			for entry in multi_parts:
				seq = ''
				name = entry[0]
				parts = entry[1]
				for part in parts:
					for gbk in gbks:
						if gbk.id == part.ref:
							part.ref = ''
							if seq:
								seq = seq + part.extract(gbk.seq)
							else:
								seq = part.extract(gbk.seq)		# for the first part, to initiate a sequence object

				out_file.write(">%s from %s\n%s\n" % (name, 'multi-record ' + gbks[0].annotations['organism'], seq))

else:
	parser.print_help(sys.stderr)
	sys.exit(1)

