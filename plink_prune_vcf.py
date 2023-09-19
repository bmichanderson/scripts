#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: September 2023
# Description: use a PLINK *.prune.out file to remove SNPs from a VCF
##########################


import argparse
import os
import sys


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to use the output of PLINK "--indep" to remove SNPs from a VCF')


# add arguments to parse
parser.add_argument('-p', type = str, dest = 'plink_file', help = 'The PLINK *.prune.out file with SNP IDs in the format CHROM:POS, one per line')
parser.add_argument('-v', type = str, dest = 'vcf_file', help = 'The VCF file to remove SNPs from; output will prefix the name with "mod_"')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
plink_file = args.plink_file
vcf_file = args.vcf_file
if any([not plink_file, not vcf_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)


# read in the SNPs to remove
snp_list = []
with open(plink_file, 'r') as pfile:
	for line in pfile:
		chrom = line.rstrip().split(':')[0]
		pos = line.rstrip().split(':')[1]
		snp_list.append(str(chrom) + '_' + str(pos))


# filter the VCF
dropped = 0
kept = 0
with open(vcf_file, 'r') as vcf, open('mod_' + os.path.basename(vcf_file), 'w') as outvcf:
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			outvcf.write(line)
		else:
			line_fields = line.rstrip().split()
			chrom = line_fields[0]
			pos = line_fields[1]
			search_str = str(chrom) + '_' + str(pos)
			if search_str in snp_list:
				dropped = dropped + 1
				continue
			else:
				kept = kept + 1
				outvcf.write(line)


print('Removed ' + str(dropped) + ' SNPs and kept ' + str(kept) + ' SNPs in file mod_' + os.path.basename(vcf_file))
