#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 28-29 Jan 2019
# Modified: Oct 2020
# Description: set the starting base position for a circular contig (arg 1) to a set value (arg 2) and optionally reverse complement the sequence
##########################

import sys			# allows access to command line arguments
from Bio import SeqIO		# for reading and writing Seq Record objects

def help():
	print('A script to first set the starting position for a circular contig (arg 1) and then optionally reverse complement the sequence.')
	print('To simply reverse complement, set the start_position argument to 1 and enter yes for the reverse_comp argument.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' contig_fasta start_position(1 is first) reverse_comp(yes or [no])' )
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
	if len(sys.argv[1:]) > 2:
		if 'y' in sys.argv[3]:
			reverse_comp = 'yes'
		else:
			reverse_comp = 'no'
	else:
		reverse_comp = 'no'
	print('Arguments are contig_fasta =  ' + str(sys.argv[1]) + ', start_position = ' + str(sys.argv[2]) + ', and reverse_comp = ' + reverse_comp)
	start = int(sys.argv[2])
else:
	sys.exit(help())


# Set a new starting position and output the file
with open(sys.argv[1], 'r') as contig_file, open('new_contig.fasta', 'w') as out_file:
	fasta = SeqIO.read(contig_file, 'fasta')
	if start > 1:
		fasta.seq = fasta.seq[start-1:] + fasta.seq[0: start-1]
	if reverse_comp == 'yes':
		fasta.seq = fasta.seq.reverse_complement()
	SeqIO.write(fasta, out_file, 'fasta')
