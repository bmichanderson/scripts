#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 29 Jan 2019
# Modified: Oct 2020
# Description: filter a bam file for reads that map for greater than a given % (arg2) of their length
##########################

import sys			# allows access to command line arguments
import re			# allows the use of regular expressions
import pysam			# allows the reading/writing and querying of sam/bam files
from Bio import SeqIO		# allows the manipulation of Seq Record objects

def help():
	print('A script to filter an input sorted and indexed bam file for reads that hit a reference for a given % of their length.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' bam_file hit_percent(1--100) output_bam(yes or [no])')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
	if 'y' in sys.argv[3]:
		output_bam = 'yes'
	else:
		output_bam = 'no'
	hit_percent = float(sys.argv[2])
	print('Arguments are bam_file =  ' + str(sys.argv[1]) + ', hit_percent = ' + str(hit_percent) + ', and output_bam = ' + output_bam)
else:
	sys.exit(help())


# Filter the bam file
bamfile = pysam.AlignmentFile(sys.argv[1], 'rb')
rec_keep = []
for rec in bamfile.fetch():
	if rec.reference_length >= ((hit_percent/100) * rec.query_length):
		rec_keep.append(rec)


# if an output bam is desired, write it
if output_bam == 'yes':
	out_file = pysam.AlignmentFile('reads_out.bam', 'w', template = bamfile)
	for rec in rec_keep:
		out_file.write(rec)
	out_file.close()


# write the output reads as a fasta
with open('reads_out.fasta', 'w') as out_file:
	for rec in rec_keep:
		out_file.write('>' + rec.query_name + '\n' + rec.query_sequence + '\n')
