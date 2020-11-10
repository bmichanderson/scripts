#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 3 Dec 2018
# Modified: Nov 2020
# Description: read a (multi)fasta file and BLASTN contigs against a reference; keep contigs if they hit the reference
##########################

import sys
import argparse
import subprocess		# a module to allow calls to programs outside python

from Bio import SeqIO		# SeqIO is part of Biopython for parsing files
from Bio.Blast import NCBIXML	# the module for parsing results in xml format (default)


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to BLASTN a (multi)fasta file and keep entries that hit a reference. Requires NCBI command line tools, e.g. makeblastdb and blastn')


# add arguments to parse
parser.add_argument('contigs', type=str, help='The (multi)fasta file to BLAST')
parser.add_argument('-r', type=str, dest='reference', help='Specify a reference fasta file (required)')
parser.add_argument('-p', type=int, dest='percent', help='The percentage of the contig length that must be hit to retain [default 20]')


# parse the command line
if len(sys.argv[1:]) == 0:              # if there are no arguments
        parser.print_help(sys.stderr)
        sys.exit(1)

args = parser.parse_args()

con = args.contigs
ref = args.reference
per = args.percent

if not per:
	per = 20


# make a blast nucleotide database from the fasta reference
subprocess.call('makeblastdb -in ' + ref + ' -out tempdb -dbtype nucl', shell=True)


# blastn the contigs against the reference database
subprocess.call('blastn -task blastn -query ' + con + ' -out temp.xml -db tempdb -outfmt 5', shell=True)


# capture the fasta contigs into a list for evaluating length and choosing a subset
con_list = []
with open(con, 'r') as fasta_file:
	for con in SeqIO.parse(fasta_file, 'fasta'):
		con_list.append(con)

print('')
print('Read in ' + str(len(con_list)) + ' contigs')


# parse the resulting xml file to record what fasta entries hit for at least the set percentage of their length
# note that this only keeps hsps with identity > 70% and length > 100 bp
fasta_hits = []
no_hits = 0
no_len = 0
with open('temp.xml', 'r') as result_handle:
	for blast_record in NCBIXML.parse(result_handle):
		if blast_record.alignments:
			lengths = []
			for align in blast_record.alignments:
				for hsp in align.hsps:
					if all([float(hsp.identities)/float(hsp.align_length) * 100 > 70, hsp.align_length > 100]):
						lengths.append(hsp.align_length)

			for con in con_list:
				if con.description == blast_record.query:
					if sum(lengths) >= (len(con) * per / 100):
						fasta_hits.append(con.description)
						break
					else:
						no_len = no_len + 1
						break
		else:
			no_hits = no_hits + 1

print('A total of ' + str(no_hits) + ' contigs did not hit the reference and ' + str(no_len) + ' hit with insufficient length')


# create a new list of contigs that hit, then write them out
new_con_list = []
for con in con_list:
	if con.description in fasta_hits:
		new_con_list.append(con)

print('Found ' + str(len(new_con_list)) + ' contigs that hit the reference for more than ' + str(per) + ' percent of their length')

with open('hits.fasta', 'w') as outfile:
	for con in new_con_list:
		SeqIO.write(con, outfile, 'fasta')


# remove the temporary files created
subprocess.call('rm temp*', shell=True)
