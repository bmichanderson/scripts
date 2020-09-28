#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 28 Sep 2020
# Description: read a text file (arg1) of accession numbers (ideally with taxon names) and download genbank records
##########################

import sys			# allows access to command line arguments
import os			# allows interacting with the operating system directory structure

from Bio import Entrez		# for accessing records online
from Bio import SeqIO		# SeqIO is part of Biopython for parsing files

#from urllib.error import HTTPError 	# for Python3
#from urllib2 import HTTPError	# needed for the batch downloading
#import time			# for delaying execution?


def help():
	print('A script to download Genbank records given a file with a list of accession numbers and ideally taxon info.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' accessions_file')
	print('')
	print('If including taxon info, the file should be in tab-delimited format and one entry per line: accession_no, genus, specific_ep.')
	print('If there are only accession numbers, they should appear one per line.')
	print('')


if len(sys.argv) != 2:
	sys.exit(help())
else:
	acc_file = sys.argv[1]


Entrez.email = 'banderson2914@gmail.com'	# always tell NCBI who is requesting data


# read in and parse the file then decide how to proceed
acc_list_no = []
acc_list_taxa = []

with open(acc_file, 'r') as infile:
	for line in infile:
		if len(line.rstrip().split()) == 1:	# only accession number
			acc_num = line.rstrip().split()[0]
			acc_list_no.append(acc_num)
		elif len(line.rstrip().split()) > 2:		# including taxon info
			acc_num = line.rstrip().split()[0]
			taxon = line.rstrip().split()[1] + '_' + line.rstrip().split()[2]
			acc_list_taxa.append((acc_num, taxon))
		else:
			sys.exit(help())

if all((len(acc_list_no) < 25, len(acc_list_taxa) < 25)):		# if dealing with a small number of records
	if len(acc_list_no) > 0:		# if there are entries without taxon info
		print('Looking for ' + str(len(acc_list_no)) + ' Genbank files based on only accession numbers.')
		for acc in acc_list_no:
			print('Downloading Genbank file for accession number ' + str(acc))
			net_handle = Entrez.efetch(db = 'nucleotide', id = acc, rettype = 'gb', retmode = 'text')
			out_handle = open('temp.gb', 'w')
			out_handle.write(net_handle.read())
			out_handle.close()
			net_handle.close()

			gb_handle = open('temp.gb', 'r')
			gb = SeqIO.read(gb_handle, 'genbank')
			gb_handle.close()

			taxon = '_'.join(gb.annotations['organism'].split()[0:2])
			filename = taxon + '_' + acc.split('.')[0] + '.gb'

			if not os.path.isfile(filename):		# if the file doesn't already exist
				os.rename('temp.gb', filename)
			else:
				os.remove('temp.gb')

	if len(acc_list_taxa) > 0:		# if there are entries with taxon info
		print('Looking for ' + str(len(acc_list_taxa)) + ' Genbank files with taxon info.')
		for acc in acc_list_taxa:
			filename = acc[1] + '_' + acc[0].split('.')[0] + '.gb'		# include the split('.') in case format is NC_12345.1
			if not os.path.isfile(filename):
				print('Downloading Genbank file ' + filename)
				net_handle = Entrez.efetch(db = 'nucleotide', id = acc[0], rettype = 'gb', retmode = 'text')
				out_handle = open(filename, 'w')
				out_handle.write(net_handle.read())
				out_handle.close()
				net_handle.close()

else:		# if there are numerous records to download -- inspired by http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter:entrez
	if len(acc_list_no) > 0:		# if there are entries without taxon info
		posting = Entrez.epost(db = 'nucleotide', id = ','.join(acc_list_no))
		search_results = Entrez.read(posting)
		webenv = search_results['WebEnv']
		query_key = search_results['QueryKey']

		count = len(acc_list_no)
		if count > 10:
			batch_size = 10
		else:
			batch_size = count

		for start in range(0, count, batch_size):
			end = min(count, start + batch_size)
			print('Downloading Genbank files ' + str(start + 1) + ' to ' + str(end) + ' of ' + str(count) + ' based on only accession numbers.')

			# comment this block out if using the block below for HTTP errors
			net_handle = Entrez.efetch(db = 'nucleotide', retstart = start, retmax = batch_size, webenv = webenv,
							query_key = query_key, idtype = 'acc', rettype = 'gb', retmode = 'text')

			# a block to accommodate HTTP errors (apparently)
#			attempt = 0
#			while attempt < 3:
#				attempt += 1
#				try:
#					net_handle = Entrez.efetch(db = 'nucleotide', retstart = start, retmax = batch_size, webenv = webenv,
#									query_key = query_key, idtype = 'acc', rettype = 'gb', retmode = 'text')
#
#				except HTTPError as err:
#					if 500 <= err.code <= 599:
#						print('Received error from server ' + str(err))
#						print('Attempt ' + str(attempt) + ' of 3')
#						time.sleep(15)
#					else:
#						raise

			# write the data (batch_size genbank files)
			genbanks = SeqIO.parse(net_handle, 'gb')
			for position in range(start, end):
				out_handle = open('temp' + str(position) + '.gb', 'w')
				SeqIO.write(next(genbanks), out_handle, 'gb')
				out_handle.close()
			net_handle.close()

		# rename the genbank files
		for index in range(0, count):
			gb_handle = open('temp' + str(index) + '.gb', 'r')
			gb = SeqIO.read(gb_handle, 'genbank')
			gb_handle.close()

			taxon = '_'.join(gb.annotations['organism'].split()[0:2])
			filename = taxon + '_' + gb.annotations['accessions'][0].split('.')[0] + '.gb'

			if not os.path.isfile(filename):
				os.rename('temp' + str(index) + '.gb', filename)
			else:
				os.remove('temp' + str(index) + '.gb')


	if len(acc_list_taxa) > 0:		# if there are entries with taxon info
		acc_nums_taxa = []
		for acc in acc_list_taxa:
			acc_nums_taxa.append(acc[0])

		posting = Entrez.epost(db = 'nucleotide', id = ','.join(acc_nums_taxa))
		search_results = Entrez.read(posting)
		webenv = search_results['WebEnv']
		query_key = search_results['QueryKey']

		count = len(acc_list_taxa)
		if count > 10:
			batch_size = 10
		else:
			batch_size = count

		for start in range(0, count, batch_size):
			end = min(count, start + batch_size)
			print('Downloading Genbank files ' + str(start + 1) + ' to ' + str(end) + ' of ' + str(count) + ' with taxon info.')

			# comment this block out if using the block below for HTTP errors
			net_handle = Entrez.efetch(db = 'nucleotide', retstart = start, retmax = batch_size, webenv = webenv,
							query_key = query_key, idtype = 'acc', rettype = 'gb', retmode = 'text')

			# a block to accommodate HTTP errors (apparently)
#			attempt = 0
#			while attempt < 3:
#				attempt += 1
#				try:
#					net_handle = Entrez.efetch(db = 'nucleotide', retstart = start, retmax = batch_size, webenv = webenv,
#									query_key = query_key, idtype = 'acc', rettype = 'gb', retmode = 'text')
#
#				except HTTPError as err:
#					if 500 <= err.code <= 599:
#						print('Received error from server ' + str(err))
#						print('Attempt ' + str(attempt) + ' of 3')
#						time.sleep(15)
#					else:
#						raise

			# write the data (batch_size genbank files)
			genbanks = SeqIO.parse(net_handle, 'gb')
			for position in range(start, end):
				filename = acc_list_taxa[position][1] + '_' + acc_list_taxa[position][0].split('.')[0] + '.gb'
				if not os.path.isfile(filename):
					out_handle = open(filename, 'w')
					SeqIO.write(next(genbanks), out_handle, 'gb')
					out_handle.close()
				else:
					next(genbanks)

			net_handle.close()
