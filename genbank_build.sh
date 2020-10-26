#!/bin/bash

##########################
# Author: B. Anderson
# Date: 3 April 2020
# Modified: Oct 2020
# Description: using a fasta file (arg1), locations file (arg2) and a key file (arg3), and with a submission template (.sbt) in the directory, build a genbank file
##########################

# python script to call
gen_annot_script=~/scripts/annotations_tbl.py


# a function for when the script is called incorrectly
usage ()
{
	echo "This script will build a genbank file from a fasta, a table of feature locations, and a specified key file (for products of proteins)."
	echo
	echo "It requires a template (template.sbt) from NCBI in the current directory, and the NCBI programs tbl2asn and asn2gb in the path."
	echo
	echo "Note: the fasta description can contain modifiers for the created genbank"
	echo ">SeqID [organism=Cuscuta australis] [topology=circular] [gcode=11] chloroplast, complete genome"
	echo
	echo "Usage: $(basename $0) fasta_file locations_file key_file"
	echo
	exit 1
}

if [ $# -eq 0 ]; then
	usage
elif [ -z "$1" ]; then
	usage
elif [ -z "$2" ]; then
	usage
elif [ -z "$3" ]; then
	usage
fi


# assign command line arguments
fsa_file="$1"
loc_file="$2"
key_file="$3"


# generate the features table by calling the python script
"$gen_annot_script" "$loc_file" "$key_file"


# modify the features table to match the name of the contig
contig_name=$(grep -m 1 ">" "$fsa_file" | awk '{gsub(">",""); print $1}')
sed -i "s/SeqID/$contig_name/g" features.tbl


# build the Sequin format table file in the current directory
tbl2asn -i "$fsa_file" -f features.tbl -t template.sbt


# convert the resulting ANS.1 file into a GenBank file
asn2gb -i "${fsa_file/.*/.sqn}" -o "${fsa_file/.*/.gb}"
