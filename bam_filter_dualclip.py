#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 28-29 Jan 2019
# Modified: Oct 2020
# Description: filter a bam file for reads that map across an entire reference extract (arg2) with soft clipping on both sides
##########################

import sys			# allows access to command line arguments
import re			# allows the use of regular expressions
import pysam			# allows the reading/writing and querying of sam/bam files
from Bio import SeqIO		# allows the manipulation of Seq Record objects

def help():
	print('A script to filter an input sorted and indexed bam file for reads that completely span a reference extract with soft clipping on both sides.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' bam_file extract_fasta output_bam(yes or [no])')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
	if 'y' in sys.argv[3]:
		output_bam = 'yes'
	else:
		output_bam = 'no'
	print('Arguments are bam_file =  ' + str(sys.argv[1]) + ', extract_fasta = ' + str(sys.argv[2]) + ', and output_bam = ' + output_bam)
else:
	sys.exit(help())


# Import the reference
ref = SeqIO.read(open(sys.argv[2], 'r'), 'fasta')
len_extract = len(ref.seq)


# Filter the bam file
bamfile = pysam.AlignmentFile(sys.argv[1], 'rb')
rec_list = []
rec_keep = []
for rec in bamfile.fetch():
	if rec.reference_length > (0.9 * len_extract):
		rec_list.append(rec)

for rec in rec_list:
	if re.match('^(\d+)S(.+)S$', rec.cigarstring):	# for soft clipping both ends needs to start with digits + S, and end with S
		rec_keep.append(rec)


# if an output bam is desired, write it
if output_bam == 'yes':
	out_file = pysam.AlignmentFile('reads_out.bam', 'w', template = bamfile)
	for rec in rec_keep:
		out_file.write(rec)
	out_file.close()


# write the output reads as a fasta with the reference
with open('reads_out.fasta', 'w') as out_file:
	out_file.write('>' + ref.name + '\n' + str(ref.seq) + '\n')
	for rec in rec_keep:
		out_file.write('>' + rec.query_name + '\n' + rec.query_sequence + '\n')

