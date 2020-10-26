#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 3 Apr 2020
# Modified: Oct 2020
# Description: generate a Sequin format table of annotations for use with a .fsa file and tbl2asn then asn2gb to ultimately create a GenBank flatfile
##########################

import sys			# allows access to command line arguments

def help():
	print('A script to generate a Sequin format table of annotations from a file with names of genes and locations, one per line.')
	print('')
	print('Format of the tab-delimited locations file (after headers):')
	print('gene_name	note(; delimited)	start1	end1	start2	end2	start3	end3	...')
	print('')
	print('Introns should be annotated similarly to coding sequences, with keyword \'intron\' in the note')
	print('')
	print('Also include a tab-delimited key file for naming conventions and gene products of four columns with headers, in the format:')
	print('gene_name	formal_name	type	product')
	print('')
	print('Gene type should be one of CDS, rRNA or tRNA')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' locations_file key_file')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
	print('Arguments are locations_file:  ' + str(sys.argv[1]) + ', and key_file: ' + str(sys.argv[2]))
else:
	sys.exit(help())


# Read in the key file and create a dictionary for looking up name, type and product
key_dict = {}
with open(sys.argv[2], 'r') as key_file:
	for line_no, line in enumerate(key_file, start=1):
		if line_no == 1:	# first line is the header row
			continue
		entry = line.rstrip().split('\t')
		key_dict[entry[0].lower()] = '\t'.join(entry[1:4])


# Read in the locations file and output the Sequin format table
with open(sys.argv[1], 'r') as loc_file, open('features.tbl', 'w') as out_file:

	out_file.write('>Features SeqID\n')		# write header

	locus_index = 1

	for line_no, line in enumerate(loc_file, start=1):

		if line_no == 1:	# first line is the header row
			continue

		old_entries = line.rstrip().split('\t')
		entries = [entry.rstrip() for entry in old_entries]		# to get rid of whitespace
		gene = entries[0].lower()

		pieces = []

		if len(entries) > 4:			# i.e. if the gene/feature is in more than one piece
			loc_iter = iter(entries[4:])
			for loc in loc_iter:
				pieces.append((loc, next(loc_iter)))

		note_entry = ''

		if len(str(entries[1])) > 0:		# parse notes if the second entry has text

			notes = str(entries[1]).split(';')

			for note in notes:

				note = note.lstrip().rstrip().lower()

				if 'pseudo' in note:
					note_entry = note_entry + '\t\t\tpseudo\n'

				elif 'trans-splicing' in note:
					note_entry = note_entry + '\t\t\texception\ttrans-splicing\n'

				elif 'codon_start' in note:
					note_entry = note_entry + '\t\t\tcodon_start\t' + note.split(' ')[1] + '\n'

				else:
					note_entry = note_entry + '\t\t\tnote\t' + note + '\n'


		start1 = entries[2]
		start2 = entries[3]

		gene_name = key_dict[gene].split('\t')[0]
		gene_type = key_dict[gene].split('\t')[1]
		gene_prod = key_dict[gene].split('\t')[2]


		if not any((gene_type == 'CDS', gene_type == 'rRNA', gene_type == 'tRNA')):
			 sys.exit('Incorrectly specified gene type')


		if gene_type == 'CDS':			# prepend a protein product identifier when annotating a CDS

			note_entry = '\t\t\tprotein_id\tgnl|DBNAME|LOCUSTAG_' + format(locus_index, '03d') + '\n' + note_entry


		if len(pieces) > 0:

			if 'trans-splicing' in note_entry:

				out_file.write(str(start1) + '\t' + str(start2) + '\tgene\n' +
						'\n'.join('\t'.join(x) for x in pieces) + '\n' +
						'\t\t\tgene\t' + gene_name + '\n' +
						'\t\t\tlocus_tag\tLOCUSTAG_' + format(locus_index, '03d') + '\n' +
						'\t\t\texception\ttrans-splicing\n' +
						str(start1) + '\t' + str(start2) + '\t' + gene_type + '\n' +
						'\n'.join('\t'.join(x) for x in pieces) + '\n' +
						'\t\t\tproduct\t' + gene_prod + '\n' +
						note_entry +
						str(start1) + '\t' + str(start2) + '\texon\n' +
						'\t\t\tnumber\t1\n' +
						'\n'.join(('\t'.join(x) + '\texon\n\t\t\tnumber\t' + str(index)) for index, x in enumerate(pieces, start=2)) + '\n')

			elif 'intron' in note_entry:

				out_file.write(str(start1) + '\t' + str(start2) + '\tintron\n' +
						'\t\t\tnumber\t1\n' +
						'\n'.join(('\t'.join(x) + '\tintron\n\t\t\tnumber\t' + str(index)) for index, x in enumerate(pieces, start=2)) + '\n')

			else:
				out_file.write(str(start1) + '\t' + str(pieces[-1][1]) + '\tgene\n' +
						'\t\t\tgene\t' + gene_name + '\n' +
						'\t\t\tlocus_tag\tLOCUSTAG_' + format(locus_index, '03d') + '\n' +
						str(start1) + '\t' + str(start2) + '\t' + gene_type + '\n' +
						'\n'.join('\t'.join(x) for x in pieces) + '\n' +
						'\t\t\tproduct\t' + gene_prod + '\n' +
						note_entry +
						str(start1) + '\t' + str(start2) + '\texon\n' +
						'\t\t\tnumber\t1\n' +
						'\n'.join(('\t'.join(x) + '\texon\n\t\t\tnumber\t' + str(index)) for index, x in enumerate(pieces, start=2)) + '\n')

		else:
			if 'intron' in note_entry:

				out_file.write(str(start1) + '\t' + str(start2) + '\tintron\n' +
						'\t\t\tnumber\t1\n')


			else:
				out_file.write(str(start1) + '\t' + str(start2) + '\tgene\n' +
						'\t\t\tgene\t' + gene_name + '\n' +
						'\t\t\tlocus_tag\tLOCUSTAG_' + format(locus_index, '03d') + '\n' +
						str(start1) + '\t' + str(start2) + '\t' + gene_type + '\n' +
						'\t\t\tproduct\t' + gene_prod + '\n' +
						note_entry)

		if 'intron' not in note_entry:
			locus_index = locus_index + 1

