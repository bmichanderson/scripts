#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 6 Feb 2020
# Modified: Oct 2020
# Description: read a multifasta file (arg1), keep those >= or <= a specified coverage (arg2)
##########################

import sys			# allows access to command line arguments
import re			# allows use of regular expression searching
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files

def help():
	print('A script to filter a fasta file of contigs based on reported coverage, given a user specified coverage and whether min or max filtering.')
	print('The filtered entries are output as a mutlifasta file in the current directory.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' fasta_file cover_cutoff min/max')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 2:
	print('Arguments are fasta_file:  ' + str(sys.argv[1]) + ', cover_cutoff: ' + str(sys.argv[2]) + ', min/max: ' + str(sys.argv[3]))
else:
	sys.exit(help())



cov_cut = int(sys.argv[2])


if sys.argv[3].lower() == 'max':
	mode = 'max'
elif sys.argv[3].lower() == 'min':
	mode = 'min'
else:
	sys.exit(help())


fasta_list = []
cover_list = []

with open(sys.argv[1], 'r') as contig_file:
	fastas = SeqIO.parse(contig_file, 'fasta')
	for fasta in fastas:
		if 'cov=' in fasta.description:				# format is either direct from Tadpole (>contig_884,length=62594,cov=11.6,gc=0.506) or 
			descrip = fasta.description.replace(',', '')	# via Bandage (>NODE_contig_10013length=829cov=4.8gc=0.302+_length_829_cov_1)
			cov_search = re.search('.*(cov=.*)gc.*', descrip)
			if cov_search:
				cover_str = cov_search.group(1)
				cover = float(cover_str.split('=')[1])
				cover_list.append(cover)
				if any([all([mode == 'min', cover >= cov_cut]), all([mode == 'max', cover <= cov_cut])]):
					fasta_list.append(fasta)


		elif fasta.description[0] == 'N':			# format is >NODE_42_length_506_cov_15.8063
			cover = float(fasta.description.split('_')[-1])
			cover_list.append(cover)
			if any([all([mode == 'min', cover >= cov_cut]), all([mode == 'max', cover <= cov_cut])]):
				fasta_list.append(fasta)

		elif fasta.description[0].isdigit:			# format is >42 506 7998
			cover = int(fasta.description.split()[2])/int(fasta.description.split()[1])
			cover_list.append(cover)
			if any([all([mode == 'min', cover >= cov_cut]), all([mode == 'max', cover <= cov_cut])]):
				fasta_list.append(fasta)

		else:
			print('Format not recognized')
			exit()

with open('cover_' + mode + str(cov_cut) + '.fasta', 'w') as output_file:
	for fasta in fasta_list:
		SeqIO.write(fasta, output_file, 'fasta')

# print the found cover values
	for cover in sorted(cover_list):
		print(cover)

