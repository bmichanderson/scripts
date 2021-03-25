#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 25 Feb 2020
# Modified: Oct 2020, Feb 2021, Mar 2021 (deal with UNVERIFIED)
# Description: search NCBI Genbank databases and report a summary for later use in downloading
##########################

import argparse
import sys			# allows access to command line arguments
from Bio import Entrez		# for interacting with the online Entrez databases


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to search Genbank and provide a summary with accession numbers for later downloading')


# add arguments to parse
parser.add_argument('search_str', type=str, help='The search string to submit to NCBI; put in quotations: \" \"')
parser.add_argument('-f', dest='full', action='store_true', help='Run the full summary and output to file [default: do not]')
parser.add_argument('-m', type=int, dest='max_recs', help='Specify the maximum number of records to return [default 50]')
parser.add_argument('-t', dest='tax', action='store_true', help='Retrieve higher taxonomic info for the records on a full summary [default: do not]')
parser.set_defaults(full=False)
parser.set_defaults(tax=False)


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

search_str = args.search_str
full = args.full
max_recs = args.max_recs
tax = args.tax

if not max_recs:
	max_recs = 50

if full:
	print('Running a full summary and outputting result to search_summary.tab')
elif tax:
	sys.exit('The taxonomy -t option only works with a full summary -f option')


# Perform the search and retrieve accession numbers
Entrez.email = 'benjamin.anderson@su.se'		# so they can notify me if my usage is improper; this is required

handle = Entrez.esearch(db = 'nucleotide', term = search_str, idtype = 'acc', retmax = max_recs)
record = Entrez.read(handle)
handle.close()

if int(record['Count']) > max_recs:		# this may be different from the number of accession numbers returned (note: default retmax = 20)
	print('Note that ' + record['Count'] + ' records were found, but only ' + str(max_recs) + ' are set to be retrieved')

acc_list = record['IdList']

print('Retrieved ' + str(len(acc_list)) + ' accession numbers for the input search string')


# Report result and optionally retrieve and output additional information to file
if not full:
	print('Accession numbers:')
	for acc_num in acc_list:
		print(acc_num)
else:
	rec_list = []
	genera_list = []
	genera_search = []
	count = len(acc_list)

	if count > 50:
		batch_size = 50
	else:
		batch_size = count

	for start in range(0, count, batch_size):
		end = min(count, start + batch_size)
		print('Searching Genbank accession details for ' + str(start + 1) + ' to ' + str(end) + ' of ' + str(count))

		sub_acc_list = acc_list[start: end]
		id_list = ','.join(sub_acc_list)	# comma-delimited list of accession numbers
		handle = Entrez.esummary(db = 'nucleotide', id = id_list)
		records = Entrez.read(handle)
		handle.close()

		for index, record in enumerate(records):
			title = record['Title']
			if title.split()[0] == 'UNVERIFIED:':
				genus = title.split()[1]
				specific_ep = title.split()[2]
			else:
				genus = title.split()[0]
				specific_ep = title.split()[1]
			length = str(record['Length'])
			tax_id = record['TaxId']

			if genus not in genera_list:
				genera_list.append(genus)
				genera_search.append((genus, tax_id))

			rec_list.append((genus, specific_ep, sub_acc_list[index], length, title))

	# Use the taxonomy database to extract order and family for each genus, if requested
	if tax:
		genera_dic = {}
		count = len(genera_search)

		if count > 50:
			batch_size = 50
		else:
			batch_size = count

		for start in range(0, count, batch_size):
			end = min(count, start + batch_size)
			print('Searching Genbank taxonomy details for ' + str(start + 1) + ' to ' + str(end) + ' of ' + str(count) + ' genera')

			sub_gen_list = genera_search[start: end]
			id_list = ','.join([str(x[1]) for x in sub_gen_list])
			handle = Entrez.efetch(db = 'Taxonomy', id = id_list)
			records = Entrez.read(handle)
			handle.close()

			for index, record in enumerate(records):
				for level in record['LineageEx']:
					if level['Rank'] == 'order':
						order = level['ScientificName']
					elif level['Rank'] == 'family':
						family = level['ScientificName']
					else:
						continue
				genera_dic[sub_gen_list[index][0]] = (order, family)

	# Output a summary
	out_list = []
	if tax:
		for rec in rec_list:
			out_list.append('\t'.join(genera_dic[rec[0]]) + '\t' + '\t'.join(rec))
	else:
		for rec in rec_list:
			out_list.append('\t'.join(rec))

	with open('search_summary.tab', 'w') as outfile:
		if tax:
			outfile.write('\t'.join(('order', 'family', 'genus', 'specific_ep', 'accession', 'length', 'full_title')) + '\n' + '\n'.join(out_list) + '\n')
		else:
			outfile.write('\t'.join(('genus', 'specific_ep', 'accession', 'length', 'full_title')) + '\n' + '\n'.join(out_list) + '\n')
