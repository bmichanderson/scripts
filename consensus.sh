#!/bin/bash

##########################
# Author: B. Anderson
# Date: 8 Jan 2019
# Modified: Oct 2020
# Description: build a consensus of a fasta file containing sequences for a region
##########################

# python script to call
consensus_script=~/scripts/consensus.py

# a function for when the script is called incorrectly
usage ()
{
	echo "This script will build a consensus file from an input fasta file with sequences for a given region and using a threshold (arg2)."
	echo "It requires MAFFT to be in the path and Biopython."
	echo
	echo "Usage: $(basename $0) input threshold"
	echo
	exit 1
}

if [ $# -eq 0 ]; then
	usage
elif [ -z "$1" ]; then
	usage
elif [ -z "$2" ]; then
	usage
fi


# assign command line arguments
genes_file="$1"
threshold="$2"

gfile="$( basename $genes_file )"


# construct an alignment using MAFFT
mafft --auto "$genes_file" > "${gfile/.fasta/-alignment.fasta}"


# generate a consensus by calling a python script
"$consensus_script" "${gfile/.fasta/-alignment.fasta}" "$threshold"


# remove the intermediate file
rm "${gfile/.fasta/-alignment.fasta}"
