#!/bin/bash

#################
# Author: B.M. Anderson
# Date: Sept 2024
# Modified: Apr 2025 (made pre-padding unnecessary, but there will always be Ns in output)
# Description: Reform a nuclear ribosomal (18S–ITS–26S) fasta assembly contig using supplied reference sequences
# Note: the script will keep up to 500 bp before 18S and beyond 26S (can be changed below with variable "retain")
#	Incomplete assemblies (no portions of either 18S or 26S) will be discarded without changes)
#	It requires NCBI tools (makeblastdb, blastn) and three Python scripts for fasta manipulation
#################


# python scripts used for fasta manipulation
extract_script=~/scripts/fasta_extract.py
reform_script=~/scripts/reform_contig.py
pad_script=~/scripts/extend_circle.py


# set how many bases to keep (max) in front of 18S and beyond 26S
retain=500


# Define a help function for an incorrect or empty call
help()
{
(	echo -e "\nThis script will reform a single nuclear ribosomal DNA assembly fasta file to set" \
	"consistent start position and orientation of the repeat unit\n"
	echo -e "It requires NCBI tools (makeblastdb, blastn) and two fasta references:"
	echo -e "\t1 - 18S in the orientation expected before ITS"
	echo -e "\t2 - 26S in the orientation expected after ITS"
	echo -e "\nIt will create a new copy of the ribosomal fasta as \"new_ribo.fasta\"\n"
	echo -e "Usage: reform_ribo.sh assembly.fasta ref1 ref2\n"
	echo -e "WARNING: will remove temp*.out temp*.fasta db* among other files"
) 1>&2
	exit 1
}


# Check if there are enough command arguments
if [ "$#" -lt 3 ]; then
	help
fi


# Set full paths for file arguments
assembly=$(readlink -f "$1")
ref1=$(readlink -f "$2")
ref2=$(readlink -f "$3")


# Step 1: generate a BLAST database for the assembly in the working directory and blast against itself
# Assuming it is unknown whether the sequence is circular, add 50 Ns to each side to
# represent an unknown gap after any adjustments
python3 "$pad_script" -b 50 -n yes "$assembly" && mv new_contig.fasta temp_ribo.fasta

makeblastdb -in temp_ribo.fasta -out db -dbtype nucl -logfile temp.log
blastn -query temp_ribo.fasta -db db -outfmt "6 std slen" | head > temp_blast.out
len=$(awk ' NR==1 {print $13} ' temp_blast.out)


# Step 2: determine the orientation of 18S and reverse complement the assembly if needed
#	This will take the top hit if 18S is broken up, and raise a warning if there are more hits
blastn -query "$ref1" -db db -outfmt 6 > temp_blast_ref1.out
if [ ! -s temp_blast_ref1.out ]; then
	echo "No hit to 18S in $(basename $assembly)! Exiting..."
	rm temp*.out temp.perf temp.log db*
	exit 1
fi
if [ $(cat temp_blast_ref1.out | wc -l) -gt 1 ]; then
	echo "More than one hit to 18S in $(basename $assembly)! Ordering by reference..."
	sort -k 7 -n temp_blast_ref1.out > temp && mv temp temp_blast_ref1.out
fi
ref1_start=$(awk ' NR==1 {print $9} ' temp_blast_ref1.out)
ref1_end=$(awk ' NR==1 {print $10} ' temp_blast_ref1.out)
if [ "$ref1_start" -gt "$ref1_end" ]; then
	python3 "$reform_script" temp_ribo.fasta 1 yes 
	mv new_contig.fasta temp_ribo.fasta
	rm db*
	makeblastdb -in temp_ribo.fasta -out db -dbtype nucl -logfile temp.log
	blastn -query "$ref1" -db db -outfmt 6 > temp_blast_ref1.out
	ref1_start=$(awk ' NR==1 {print $9} ' temp_blast_ref1.out)
	ref1_end=$(awk ' NR==1 {print $10} ' temp_blast_ref1.out)
fi


# Step 3: reset the start of the assembly to $retain bases before the start of 18S
#	If 18S doesn't have enough bases in front of it, the 100 Ns will represent an unknown gap
if [ "$ref1_start" -lt $((retain + 1)) ]; then
	new_start=$(echo "$len - $retain - 1 + $ref1_start" | bc)
else
	new_start=$(echo "$ref1_start - $retain" | bc)
fi
python3 "$reform_script" temp_ribo.fasta "$new_start" && mv new_contig.fasta temp_ribo.fasta


# Step 4: regenerate the database and check for the location of 26S
#	This will again use the top (first) hit
rm db*
makeblastdb -in temp_ribo.fasta -out db -dbtype nucl -logfile temp.log
blastn -query "$ref2" -db db -outfmt 6 > temp_blast_ref2.out
if [ ! -s temp_blast_ref2.out ]; then
	echo "No hit to 26S in $(basename $assembly)! Exiting..."
	rm temp*.out temp.perf temp.log db*
	exit
fi
ref2_start=$(awk ' NR==1 {print $9} ' temp_blast_ref2.out)
ref2_end=$(awk ' NR==1 {print $10} ' temp_blast_ref2.out)


# Step 5: extract just the desired portion of the assembly
#	If the end is less than the desired retained bases beyond the end of 26S, don't alter
if [ $(echo "$len - $ref2_end" | bc) -lt "$retain" ]; then
	echo -e "No need to trim end"
else
	new_end=$(echo "$ref2_end + $retain" | bc)
	python3 "$extract_script" temp_ribo.fasta -c 1.."$new_end"
	mv extract.fasta temp_ribo.fasta
fi

mv temp_ribo.fasta new_ribo.fasta
rm temp*.out temp.perf temp.log db*
