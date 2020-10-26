#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 25 Feb 2020
# Modified: Oct 2020
# Description: search NCBI Genbank databases and report a summary for later use in downloading
##########################

import sys			# allows access to command line arguments
from Bio import Entrez		# for interacting with the online Entrez databases

def help():
	print('A script to search Genbank and provide a summary with accession numbers for later downloading.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' options(- ... - ...) search_string')
	print('')
	print('Options:')
	print('	-f	Whether to run the full summary, yes or no [default]')
	print('')
	print('	-m	Maximum number of records to return/search [default 50]')
	print('')
	print('	-t	Whether to search for and report taxonomy, yes or no [default]')
	print('')
	print('The search string should be in parentheses, e.g. \"Cuscuta[Organism] AND mitochondrial[Title]\"')
	print('')


options = ['-f', '-m', '-t']
opt_valdic = {}


# Parse the command line
if len(sys.argv[1:]) > 0:

	search_string = sys.argv[-1:][0]		# retrieve the last argument slice, then entry as string

	for index, arg in enumerate(sys.argv[1:]):
		if arg[0] == '-': 			# if an option
			if arg in options:
				if sys.argv[1:][index + 1][0] == '-':
					sys.exit('No argument provided for ' + str(arg))
				opt_valdic[arg] = sys.argv[1:][index + 1]
			else:
				sys.exit('The argument ' + str(arg) + ' is not a valid option')

	if '-f' in opt_valdic:
		if opt_valdic['-f'].lower() == 'yes':
			full_summ = True
	else:
		full_summ = False

	if '-m' in opt_valdic:
		max_records = int(opt_valdic['-m'])
	else:
		max_records = 50

	if '-t' in opt_valdic:
		if opt_valdic['-t'].lower() == 'yes':
			search_tax = True
	else:
		search_tax = False

else:
	sys.exit(help())


# perform the search and retrieve accession numbers
Entrez.email = 'benjamin.anderson@su.se'		# so they can notify me if my usage is improper; this is required

handle = Entrez.esearch(db = 'nucleotide', term = search_string, idtype = 'acc', retmax = max_records)
record = Entrez.read(handle)
handle.close()

if int(record['Count']) > max_records:		# this may be different from the number of accession numbers returned (note: default retmax = 20)
	print('Note that ' + record['Count'] + ' records were found, but only ' + str(max_records) + ' are set to be retrieved.')

acc_list = record['IdList']

print('Retrieved ' + str(len(acc_list)) + ' accession numbers for the input search')

if not full_summ:
	print('Accession numbers:')
	for acc_num in acc_list:
		print(acc_num)


else:		# For each accession number, retrieve the record details and organism info to prepare the summary
	rec_list = []
	genera_list = []
	count = len(acc_list)

	if count > 50:
		batch_size = 50
	else:
		batch_size = count

	for start in range(0, count, batch_size):
		end = min(count, start + batch_size)
		print('Searching Genbank accession details for ' + str(start + 1) + ' to ' + str(end) + ' of ' + str(count))
		for position in range(start, end):
			handle = Entrez.esummary(db = 'nucleotide', id = acc_list[position])
			records = Entrez.read(handle)
			handle.close()

			title = records[0]['Title']
			genus = title.split()[0]
			specific_ep = title.split()[1]
			length = str(records[0]['Length'])
			tax_id = records[0]['TaxId']

			if genus not in genera_list:
				genera_list.append((genus, tax_id))

			rec_list.append((genus, specific_ep, acc_list[position], length, title))

	# Use the taxonomy database to extract order and family for each genus, if requested
	if search_tax:
		genera_dic = {}
		count = len(genera_list)

		if count > 50:
			batch_size = 50
		else:
			batch_size = count

		for start in range(0, count, batch_size):
			end = min(count, start + batch_size)
			print('Searching Genbank taxonomy details for ' + str(start + 1) + ' to ' + str(end) + ' of ' + str(count) + ' genera')
			for position in range(start, end):
				tax_id = genera_list[position][1]
				handle = Entrez.efetch(db = 'Taxonomy', id = tax_id)
				records = Entrez.read(handle)
				handle.close()

				for level in records[0]['LineageEx']:
					if level['Rank'] == 'order':
						order = level['ScientificName']
					elif level['Rank'] == 'family':
						family = level['ScientificName']
					else:
						continue

				genera_dic[genera_list[position][0]] = (order, family)


	# Output a summary
	out_list = []
	if search_tax:
		for rec in rec_list:
			out_list.append('\t'.join(genera_dic[rec[0]]) + '\t' + '\t'.join(rec))
	else:
		for rec in rec_list:
			out_list.append('\t'.join(rec))

	with open('search_summary.tab', 'w') as outfile:
		if search_tax:
			outfile.write('\t'.join(('order', 'family', 'genus', 'specific_ep', 'accession', 'length', 'full_title')) + '\n' + '\n'.join(out_list) + '\n')
		else:
			outfile.write('\t'.join(('genus', 'specific_ep', 'accession', 'length', 'full_title')) + '\n' + '\n'.join(out_list) + '\n')

