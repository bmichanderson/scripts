#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 18 May 2020
# Modified: June 2023
# Description: Parse the output of a HMMER search, and optionally filter for top/multiple hits
##########################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to parse the results file of a hmmsearch of a HMMER profile and ' +
	'a set of protein sequences; e.g., translated output from ORFfinder CDS fasta. If desired, the hits can be filtered ' +
	'using the \"-t\" argument.')


# add arguments to parse
parser.add_argument('hmmer', type = str, help = 'The HMMER output file')
parser.add_argument('-t', type = float, dest = 'thresh', help = 'The E-value threshold for filtering, e.g. \"1e-5\"')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
hmmer = args.hmmer
thresh = args.thresh

if thresh:
	filter_hits = True
	threshold = thresh
else:
	filter_hits = False
	threshold = 1


# parse the hmmer output line by line and capture info into lists
with open(hmmer, 'r') as hmmer_file, open(hmmer + '.parsed', 'w') as out_file:
	# first, capture all the hit information
	capture_hits = []
	hit_num = 0
	domain_num = 0
	index = 0

	for line_no, line in enumerate(hmmer_file, start = 1):
		if not line.strip():
			continue
		pieces = line.strip().split()
		if pieces[0] == 'Query:':	# indication of hit in hmmer database
			query = line.split()[1]
			mlen = line.split()[2].strip('[').strip(']').split('=')[1]
			hit_num = 0
		elif pieces[0] == '>>':		# indication of domain hit
			hit = pieces[1]
			index = line_no
			hit_num = hit_num + 1
			domain_num = 1
		elif all([pieces[0] == str(domain_num), line_no == index + 2 + domain_num]):	# hits are sequential starting 3 lines below >>
			# capture hit coordinate information of importance
			evalue = pieces[5]		# independent domain evalue, not the conditional evalue
			hmm_coord = pieces[6:9]		# two numbers and a string
			env_coord = pieces[12:15]	# two numbers and a string
			capture_hits.append([query, mlen, hit, evalue, hmm_coord, env_coord])
			domain_num = domain_num + 1
		elif pieces[0] == '//':		# indication of end of hmmer hit
			print('Processed query: ' + query)
		else:
			continue

	# second, format and output the captured info
	query_list = []
	out_file.write('Query\tMlen\tContig\tC_len\tE-value\t' +
		'H_coord_start\tH_coord_end\tH_indicator\tE_coord_start\t' +
		'E_coord_end\tE_indicator\tC_coord_start\tC_coord_end\t\n')
	print('Captured ' + str(len(capture_hits)) + ' hits')

	for capture in capture_hits:
		# re-assign variables for clarity
		query = capture[0]
		mlen = capture[1]
		hit = capture[2]
		evalue = capture[3]
		hmm_coord = capture[4]
		env_coord = capture[5]

		# ORFfinder output may contain hit location information
		locus_parts = hit.split(':')
		contig = locus_parts[0]
		if len(locus_parts) > 1:
			coord_start = locus_parts[1].split('-')[0]
			coord_end = int(locus_parts[1].split('-')[1])
			# convert the coordinates to match the envelope coords in terms of contig coords
			eq_coord = []
			start_adjust = (int(env_coord[0]) - 1) * 3
			end_adjust = (int(env_coord[1]) - int(env_coord[0])) * 3
			if str(coord_start).startswith('c'):		# reverse
				eq_coord.append(str(int(coord_start.split('c')[1]) - start_adjust))
				eq_coord.append(str(int(coord_start.split('c')[1]) - start_adjust - end_adjust - 2))
				orflen = (int(coord_start.split('c')[1]) - coord_end + 1) / 3
			else:	# normal
				eq_coord.append(str(int(coord_start) + start_adjust))
				eq_coord.append(str(int(coord_start) + start_adjust + end_adjust + 2))	# need to add two bp for rest of codon
				orflen = (coord_end - int(coord_start) + 1) / 3
		else:		# not ORFfinder output
			orflen = 'unknown'
			eq_coord = ['-', '-']

		# write the output
		out_file.write('\t'.join([query, str(mlen), contig, str(orflen), str(evalue)]) +
			'\t' + '\t'.join(hmm_coord) + '\t' + '\t'.join(env_coord) +
			'\t' + '\t'.join(eq_coord) + '\n')

		# record query
		if query not in query_list:
			query_list.append(query)

	# third, sort and filter the query hits, if requested
	if filter_hits:
		with open(hmmer + '.parsed.filtered', 'w') as out_file2:
			out_file2.write('Query\tORFs\n')
			for query in query_list:
				sublist = []
				found = False
				for capture in capture_hits:
					if capture[0] == query:
						found = True
						sublist.append(capture)
					elif found:		# assumes the query hits are all consecutive
						break
					else:
						continue

				# filter out hits with poor scores (> threshold)
				filtered_sublist = [item for item in sublist if float(item[3]) <= threshold]
				sublist = filtered_sublist
				if len(sublist) < 1:
					print('No hits passed filters for locus ' + str(query))
					continue

				# parse the hits for a query and determine most common contig, lowest E-value if ties
				contigs = [item[2].split(':')[0] for item in sublist]
				evalues = [float(item[3]) for item in sublist]

				counts = []
				contig_list = []
				evalue_list = []
				for contig in set(contigs):
					counts.append(contigs.count(contig))
					contig_list.append(contig)
				if len([item for item in counts if item == max(counts)]) > 1:		# ties are present
					indices = [index for index, item in enumerate(counts) if item == max(counts)]
					top_contigs = [contig_list[index] for index in indices]
					for contig in top_contigs:
						evalue_list.append(evalues[contigs.index(contig)])	# just the first evalue
					most_common = top_contigs[evalue_list.index(min(evalue_list))]
				else:
					most_common = max(set(contigs), key = contigs.count)

				# keep hits for the most common contig (lowest E-value in ties)
				indices = [index for index, item in enumerate(contigs) if item == most_common]
				keep_hits = [sublist[index] for index in indices]

				# order hits by HMMER profile order (first part of hmm_coord)
				keep_hits.sort(key = lambda x: int(x[4][0]))

				# write the query and comma-separated names of the ORFs to output, excluding dups
				final_orfs = []
				for hit in keep_hits:
					if hit[2] not in final_orfs:
						final_orfs.append(hit[2])
				out_file2.write(query + '\t' + ','.join(final_orfs) + '\n')
