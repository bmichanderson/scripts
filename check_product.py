#!/usr/bin/env python

# check if product is present for each CDS in a genbank file (arg1)


import sys
from Bio import SeqIO


CDS_missing = []
RNA_missing = []


with open(sys.argv[1], 'r') as gbfile:
	genbanks = SeqIO.parse(gbfile, 'genbank')
	for gb in genbanks:
		for feature in gb.features:
			if feature.type == 'CDS':
				if 'product' not in feature.qualifiers:
					CDS_missing.append(feature.qualifiers['gene'][0])
			elif 'RNA' in feature.type:
				if 'product' not in feature.qualifiers:
					RNA_missing.append(feature.qualifiers['gene'][0])

#print('Features missing product qualifier:')
#print('\tCDS: ' + '\t'.join(CDS_missing))
#print('\tRNA: ' + '\t'.join(RNA_missing))
for item in CDS_missing:
	print(item)
for item in RNA_missing:
	print(item)
