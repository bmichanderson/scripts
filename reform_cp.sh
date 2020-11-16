#!/bin/bash

#################
# Author: B. Anderson
# Date: Nov 2020
# Modified:
# Description: Reform a chloroplast fasta assembly (single circle) for consistent start position and orientation of LSC and SSC
#################


# python scripts used for fasta manipulation
extract_script="/home/banderson/scripts/fasta_extract.py"
merge_script="/home/banderson/scripts/fasta_merging.py"
reform_script="/home/banderson/scripts/reform_contig.py"


# Define a help function for an incorrect or empty call
help()
{
(	echo -e '\n''\t'"This script will reform a single chloroplast fasta file to set consistent start position and orientation of LSC and SSC"
	echo -e '\t'"It requires NCBI tools (makeblastdb, blastn) and two references:"
	echo -e '\t\t'"1 - a gene (preferably psbA) near the start of the LSC that normally has reverse orientation"
	echo -e '\t\t'"2 - a gene (preferably ccsA) in the SSC that normally has forward orientation"
	echo -e '\t'"It will create a new copy of the chloroplast fasta as \"new_cp.fasta\""'\n'
	echo -e '\t'"Usage: reform_cp.sh cp.fasta ref1 ref2"'\n'
) 1>&2
	exit 1
}


# Check if there are the correct number of command arguments
if [ $# -ne 3 ]; then
	help
fi


# Set full paths for file arguments
chloroplast=$(readlink -f "$1")
ref1=$(readlink -f "$2")
ref2=$(readlink -f "$3")


## Step 1: determine the IR locations and set a start at the beginning of the LSC

cp "$chloroplast" temp_cp.fasta
makeblastdb -in temp_cp.fasta -out db -dbtype nucl -logfile temp.log
blastn -query temp_cp.fasta -db db -outfmt 6 | head > temp_blast.out

len=$(awk ' NR==1 {print $4} ' temp_blast.out)
irlen=$(awk ' NR==2 {print $4} ' temp_blast.out)


# a: find largest possible IR (start may be in an IR) by setting a new contig at half length location and comparing blast results

"$reform_script" temp_cp.fasta $(echo "$len / 2" | bc) && mv new_contig.fasta temp_cp2.fasta
makeblastdb -in temp_cp2.fasta -out db2 -dbtype nucl -logfile temp.log
blastn -query temp_cp2.fasta -db db2 -outfmt 6 | head > temp_blast2.out

ir2len=$(awk ' NR==2 {print $4} ' temp_blast2.out)
if [ $ir2len -gt $irlen ]; then echo "found longer IR" && mv temp_cp2.fasta temp_cp.fasta && mv temp_blast2.out temp_blast.out; fi
rm db*


# b: find the locations of the IR boundaries and set start based on LSC

ira_start=$(awk ' NR==3 {print $7} ' temp_blast.out)
ira_end=$(awk ' NR==3 {print $8} ' temp_blast.out)
irb_start=$(awk ' NR==3 {print $10} ' temp_blast.out)
irb_end=$(awk ' NR==3 {print $9} ' temp_blast.out)
if [ $(echo "$irb_start - $ira_end" | bc) -gt $(echo "$len / 3" | bc) ]; then new_start=$(echo "$ira_end + 1" | bc); \
elif [ $irb_end -eq $len ]; then new_start=1; else new_start=$(echo "$irb_end + 1" | bc); fi

if [ $new_start -ne 1 ]; then "$reform_script" temp_cp.fasta "$new_start" && mv new_contig.fasta temp_cp.fasta; fi


# c: regenerate the database and set locations of IR boundaries

makeblastdb -in temp_cp.fasta -out db -dbtype nucl -logfile temp.log
blastn -query temp_cp.fasta -db db -outfmt 6 | head > temp_blast.out
ira_start=$(awk ' NR==3 {print $7} ' temp_blast.out)
ira_end=$(awk ' NR==3 {print $8} ' temp_blast.out)
irb_start=$(awk ' NR==3 {print $10} ' temp_blast.out)
irb_end=$(awk ' NR==3 {print $9} ' temp_blast.out)


## Step 2: extract the LSC; check orientation and location of ref1 in LSC; modify LSC as needed

"$extract_script" -c 1..$(echo "$ira_start - 1" | bc) temp_cp.fasta && mv extract.fasta LSC.fasta

blastn -query "$ref1" -db db -outfmt 6 > temp_blast_ref1.out
ref1_start=$(awk ' NR==1 {print $9} ' temp_blast_ref1.out)
ref1_end=$(awk ' NR==1 {print $10} ' temp_blast_ref1.out)
if [ $ref1_end -gt $ref1_start ]; then "$reform_script" LSC.fasta 1 yes && mv new_contig.fasta LSC.fasta; \
elif [ $ref1_end -gt 1000 ]; then echo "Uh oh! Script isn't working as intended!"; fi


## Step 3: extract the SSC; check orientation of ref2 in SSC; modify SSC as needed

"$extract_script" -c $(echo "$ira_end + 1" | bc)..$(echo "$irb_start - 1" | bc) temp_cp.fasta && mv extract.fasta SSC.fasta

blastn -query "$ref2" -db db -outfmt 6 > temp_blast_ref2.out
ref2_start=$(awk ' NR==1 {print $9} ' temp_blast_ref2.out)
ref2_end=$(awk ' NR==1 {print $10} ' temp_blast_ref2.out)		# this may only be a partial hit, but all I need is orientation
if [ $ref2_end -lt $ref2_start ]; then "$reform_script" SSC.fasta 1 yes && mv new_contig.fasta SSC.fasta; fi


## Step 4: extract the IRs, stitch all back together and remove intermediate files

"$extract_script" -c $ira_start..$ira_end temp_cp.fasta && mv extract.fasta IRA.fasta
"$extract_script" -c $irb_start..$irb_end temp_cp.fasta && mv extract.fasta IRB.fasta

cat LSC.fasta IRA.fasta SSC.fasta IRB.fasta > new_cp.fasta

"$merge_script" new_cp.fasta && mv merged.fasta new_cp.fasta

rm temp* db* IR* LSC.fasta SSC.fasta
