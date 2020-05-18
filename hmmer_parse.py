#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 18 May 2020
# Description: Parse the output of a HMMER search
##########################


import sys			# allows access to command line arguments


def help():
	print('A script to parse the results of a hmmsearch of a HMMER profile and a protein sequence(s) file from ORFfinder.')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' hmmer_file.out')
	print('')


# print help if the script is called without args
if len(sys.argv[1:]) == 0:
	sys.exit(help())


# parse the hmmer output line by line and capture info into lists
with open(sys.argv[1], 'r') as hmmer_file, open(sys.argv[1].replace('.out', '') + '_parsed.out', 'w') as out_file:

	# first, capture all the hit information
	capture_hits = []
	hit_num = 0
	domain_num = 0
	index = 0

	for line_no, line in enumerate(hmmer_file, start=1):
		if line.startswith('#'):
			continue
		elif not line.strip():
			continue
		elif line.split()[0] == 'Query:':	# indication of first hit in hmmer database
			query = line.split()[1]
			mlen = line.split()[2].strip('[').strip(']').split('=')[1]
			hit_num = 0
		elif line.startswith('>>'):		# indication of domain hit
			hit = line.split()[1]
			index = line_no
			hit_num = hit_num + 1
			domain_num = 1
		elif all([line.lstrip().split()[0] == str(domain_num), line_no == index + 2 + domain_num]):	# hits are sequential starting 3 lines below >>
			# capture hit coordinate information of importance
			hmm_coord = line.split()[6:9]
			env_coord = line.split()[12:15]
			capture_hits.append([query, hit, hit_num, domain_num, hmm_coord, env_coord, mlen])
			domain_num = domain_num + 1
		elif line.startswith('//'):		# indication of end of hmmer hit
			print('Processed query: ' + query)
		else:
			continue


	# second, format the captured info for output
	query_list = []

	out_file.write('Contig\tHit_num\tDomain_num\tHmm_coord\tMlen\tEnv_coord\tORFlen\tEquiv_coord\t\n')

	print(len(capture_hits))

	for capture in capture_hits:

		# re-assign variables for clarity
		query = capture[0]
		hit = capture[1]
		hit_num = capture[2]
		domain_num = capture[3]
		hmm_coord = capture[4]
		env_coord = capture[5]
		mlen = capture[6]

		# this is specific to the formatting of ORFfinder output
		locus = hit.split('_')[1]
		contig = locus.split(':')[0]
		coord_start = locus.split(':')[1]
		coord_end = locus.split(':')[2]

		# convert the coordinates to match the envelope coords in terms of contig coords
		eq_coord = []
		actual_start = int(coord_start) + 1		# to correct for ORFfinder bug
		actual_end = int(coord_end) + 1		# to correct for ORFfinder bug
		start_adjust = (int(env_coord[0]) - 1) * 3
		end_adjust = (int(env_coord[1]) - int(env_coord[0])) * 3
		if actual_start < actual_end:		# positive strand
			eq_coord.append(str(actual_start + start_adjust))
			eq_coord.append(str(actual_start + start_adjust + end_adjust + 2))	# need to add two bp for rest of codon
			orflen = (actual_end - actual_start + 1) / 3
		else:					# negative strand
			eq_coord.append(str(actual_start - start_adjust))
			eq_coord.append(str(actual_start - start_adjust - end_adjust - 2))
			orflen = (actual_start - actual_end + 1) / 3

		if query not in query_list:
			query_list.append(query)
			out_file.write(query + '\n')
			out_file.write('\t'.join([contig, str(hit_num), str(domain_num)]) + '\t' + ' '.join(hmm_coord) + '\t' + str(mlen) + '\t' + ' '.join(env_coord) + '\t' + str(orflen) + '\t' + '\t'.join(eq_coord) + '\n')
		else:
			out_file.write('\t'.join([contig, str(hit_num), str(domain_num)]) + '\t' + ' '.join(hmm_coord) + '\t' + str(mlen) + '\t' + ' '.join(env_coord) + '\t' + str(orflen) + '\t' + '\t'.join(eq_coord) + '\n')
