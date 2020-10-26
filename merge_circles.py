#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 20 Sep 2019
# Modified: Oct 2020
# Description: merge two contigs (args 1 and 3) given start coordinates of an exact repeat shared between them (args 2 4)
##########################

import sys				# allows access to command line arguments
from Bio import SeqIO			# for reading and writing Seq Record objects
from Bio.SeqRecord import SeqRecord	# for creating SeqRecord objects

def help():
	print('A script to merge two circular contigs given start coordinates for a repeat region shared between them.')
	print('')
	print('Ensure the start and end coordinates are for hits that are on the same strand (+).')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' contig1 start1 contig2 start2')
	print('')

# Ensure there are command args
if len(sys.argv[1:]) > 3:
	print('Arguments are contig1 = ' + str(sys.argv[1]) + ', start1 = ' + str(sys.argv[2]) + ', contig2 = ' + str(sys.argv[3]) + ', start2 = ' + str(sys.argv[4]))
	start1 = int(sys.argv[2])
	start2 = int(sys.argv[4])
else:
	sys.exit(help())

# Read in the files, combine and output the new merged contig
with open(sys.argv[1], 'r') as contig1_file, open(sys.argv[3], 'r') as contig2_file, open('new_contig.fasta', 'w') as out_file:
	fasta1 = SeqIO.read(contig1_file, 'fasta')
	fasta2 = SeqIO.read(contig2_file, 'fasta')
	new_fasta = SeqRecord(fasta1.seq[0: start1-1] + fasta2.seq[start2-1: ] + fasta2.seq[0: start2-1] + fasta1.seq[start1-1: ], id= 'new_contig' , description= 'len=' + str(len(fasta1.seq) + len(fasta2.seq)))
	SeqIO.write(new_fasta, out_file, 'fasta')
