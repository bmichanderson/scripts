#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: 15 Feb 2019
# Modified: Oct 2020
# Description: adjust an annotations locations file (arg1) by resetting the origin of the contig to a set value (arg2)
##########################

import sys			# allows access to command line arguments

def help():
	print('A script to update an annotations locations file after the contig it is based on has been reformed to a new origin.')
	print('')
	print('Format of the tab-delimited locations file (after header line):')
	print('gene_name	note	start1	end1	start2	end2	start3	end3	...')
	print('')
	print('Note that if reverse complement is done, it happens AFTER setting the new start (so new start should be relative to original contig).')
	print('')
	print('Usage: ' + str(sys.argv[0]) + ' locations_file contig_length new_start reverse_comp(yes or [no])')
	print('')


# Ensure there are command args
if len(sys.argv[1:]) > 1:
	if len(sys.argv[1:]) > 3:
		if 'y' in sys.argv[4]:
			reverse_comp = 'yes'
		else:
			reverse_comp = 'no'
	else:
		reverse_comp = 'no'

	print('Arguments are locations_file =  ' + str(sys.argv[1]) + ', contig_length = ' + str(sys.argv[2]) + ', new_start = ' + str(sys.argv[3]) + ', and reverse_comp = ' + reverse_comp)

	length = int(sys.argv[2])
	new_start = int(sys.argv[3])

else:
	sys.exit(help())


# Read in the locations file line by line and adjust appropriately, then output new file
out_list = []

with open(sys.argv[1], 'r') as locations_file, open('new_annotations.txt', 'w') as out_file:
	for lineno, line in enumerate(locations_file, start=1):
		if lineno == 1:
			out_list.append('\t'.join(line.rstrip().split('\t')))
		else:
			entries = line.rstrip().split('\t')
			for index, entry in enumerate(entries, start=0):
				if index < 2:		# first two columns are gene name and note
					continue
				elif int(entry) < new_start:
					entries[index] = str(int(length) + 1 - new_start + int(entry))
				elif int(entry) >= new_start:
					entries[index] = str(1 - new_start + int(entry))
				else:
					sys.exit("Error: incorrectly formatted locations file or specified length/start value")

			if reverse_comp == 'yes':
				for index, entry in enumerate(entries, start=0):
					if index < 2:		# first two columns are gene name and note
						continue
					else:
						entries[index] = str(int(length) + 1 - int(entry))
			out_list.append('\t'.join(entries))

	out_file.write('\n'.join(out_list))
