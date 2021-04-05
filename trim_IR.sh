#!/bin/bash

#####
# Author: B. Anderson
# Date: March 2021
# Modified: Apr 2021 (to make more consistent with updated reform_cp.sh)
# Description: using an input fasta file, determine a likely IR location at the end of the file and output a fasta with a single IR
#####


# python script for fasta manipulation
extract_script="/home/banderson/scripts/fasta_extract.py"


# Define a help function for an incorrect or empty call
help()
{
(       echo -e '\n''\t'"This script will attempt to trim the last IR from a chloroplast fasta file (arg1)"
        echo -e '\t'"It requires NCBI tools (makeblastdb, blastn)"
        echo -e '\t'"If successful, it will create a new copy of the chloroplast fasta as \"new_cp.fasta\""'\n'
        echo -e '\t'"Usage: trim_IR.sh cp.fasta "'\n'
	echo -e '\t'"It may not succeed if it cannot find a reciprocal hit easily or at the end of the chloroplast"'\n'
) 1>&2
        exit 1
}


# Check if there are the correct number of command arguments
if [ $# -ne 1 ]; then
        help
fi


# Set full paths for file argument
chloroplast=$(readlink -f "$1")


# Blast the file against itself
makeblastdb -in "$chloroplast" -out db -dbtype nucl -logfile temp.log && blastn -query "$chloroplast" -db db -outfmt "6 std slen" > temp_blast.out


# Evaluate the blast results to find a putative IR
length=$(head -n 1 temp_blast.out | cut -f 13)

hit2=$(awk ' NR==2 {print $4} ' temp_blast.out)
hit3=$(awk ' NR==3 {print $4} ' temp_blast.out)
hit4=$(awk ' NR==4 {print $4} ' temp_blast.out)
hit5=$(awk ' NR==5 {print $4} ' temp_blast.out)

IR_end=0

if [ $(echo "$hit3" - "$hit2" | bc | sed 's/-//') -gt 10 ]; then
	if [ $(echo "$hit4" - "$hit3" | bc | sed 's/-//') -gt 10 ]; then
		if [ $(echo "$hit5" - "$hit4" | bc | sed 's/-//') -gt 10 ]; then
			echo "No reciprocal hits in first five blast hits; not processing further"
		else
			IR_start=$(awk ' NR==4 {print $7} ' temp_blast.out)
			IR_end=$(awk ' NR==4 {print $8} ' temp_blast.out)
		fi
	else
		IR_start=$(awk ' NR==3 {print $7} ' temp_blast.out)
		IR_end=$(awk ' NR==3 {print $8} ' temp_blast.out)
	fi
else
	IR_start=$(awk ' NR==2 {print $7} ' temp_blast.out)
	IR_end=$(awk ' NR==2 {print $8} ' temp_blast.out)
fi

rm temp_blast.out temp.log temp.perf db.nhr db.nsq db.nin

if [ "$IR_end" -gt 0 ]; then
	if [ "$IR_end" -ne "$length" ]; then
		echo "Putative IR end point not at the end of the chloroplast; not processing further"
	else
		echo "Removing putative IR from end of chloroplast. Length removed = $(echo $IR_end - $IR_start | bc)"
		coord=$(echo "$IR_start" - 1 | bc)
		"$extract_script" -c 1.."$coord" "$chloroplast"
		cat <(head -n 1 "$chloroplast") <(tail -n +2 extract.fasta) > new_cp.fasta
		rm extract.fasta
	fi
fi
