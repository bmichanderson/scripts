#!/bin/bash

#################
# Author: B. Anderson
# Date: Nov 2020
# Modified: Mar 2021, April 2021 (for more checking and flexibility when sequence is lower quality and hits are broken up)
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
	echo -e '\t'"WARNING: will remove temp*.out temp*.fasta db* among other files"'\n'
) 1>&2
	exit 1
}


# Check if there are the correct number of command arguments
if [ "$#" -ne 3 ]; then
	help
fi


# Set full paths for file arguments
chloroplast=$(readlink -f "$1")
ref1=$(readlink -f "$2")
ref2=$(readlink -f "$3")


## Step 1: determine the IR locations and set a start at the beginning of the LSC

cp "$chloroplast" temp_cp.fasta
makeblastdb -in temp_cp.fasta -out db -dbtype nucl -logfile temp.log
blastn -query temp_cp.fasta -db db -outfmt "6 std slen" | head > temp_blast.out

len=$(awk ' NR==1 {print $13} ' temp_blast.out)

# typically, the second largest hits are the IR, but it is good to check given the sequence may have errors or NNNs
hit2=$(awk ' NR==2 {print $4} ' temp_blast.out)
hit3=$(awk ' NR==3 {print $4} ' temp_blast.out)
hit4=$(awk ' NR==4 {print $4} ' temp_blast.out)
hit5=$(awk ' NR==5 {print $4} ' temp_blast.out)
hit6=$(awk ' NR==6 {print $4} ' temp_blast.out)

if [ $(echo "$hit3" - "$hit2" | bc | sed 's/-//') -gt 10 ]; then
	echo "Potential problem with file and blast hits" && cp temp_blast.out REVIEW_ME_blast.out
	if [ $(echo "$hit4" - "$hit3" | bc | sed 's/-//') -gt 10 ]; then
		if [ $(echo "$hit5" - "$hit4" | bc | sed 's/-//') -gt 10 ]; then
			if [ $(echo "$hit6" - "$hit5" | bc | sed 's/-//') -gt 10 ]; then
				echo "No reciprocal hits in first six blast hits; not processing further"
				exit 1
			else
				ir_firstrow=5
			fi
		else
			ir_firstrow=4
		fi
	else
		ir_firstrow=3
	fi
else
	ir_firstrow=2
fi

irlen=$(awk -v "row=$ir_firstrow" ' NR==row {print $4} ' temp_blast.out)
ir_identity=$(awk -v "row=$ir_firstrow" ' NR==row {print $3} ' temp_blast.out)


# a: find largest possible IR (start may be in an IR) by setting a new contig at half length location and comparing blast results

"$reform_script" temp_cp.fasta $(echo "$len / 2" | bc) && mv new_contig.fasta temp_cp2.fasta
makeblastdb -in temp_cp2.fasta -out db2 -dbtype nucl -logfile temp.log
blastn -query temp_cp2.fasta -db db2 -outfmt 6 | head > temp_blast2.out

hit2=$(awk ' NR==2 {print $4} ' temp_blast2.out)
hit3=$(awk ' NR==3 {print $4} ' temp_blast2.out)
hit4=$(awk ' NR==4 {print $4} ' temp_blast2.out)
hit5=$(awk ' NR==5 {print $4} ' temp_blast2.out)
hit6=$(awk ' NR==6 {print $4} ' temp_blast2.out)

if [ $(echo "$hit3" - "$hit2" | bc | sed 's/-//') -gt 10 ]; then
	echo "Potential problem with file and blast hits" && cp temp_blast2.out REVIEW_ME_blast2.out
	if [ $(echo "$hit4" - "$hit3" | bc | sed 's/-//') -gt 10 ]; then
		if [ $(echo "$hit5" - "$hit4" | bc | sed 's/-//') -gt 10 ]; then
			if [ $(echo "$hit6" - "$hit5" | bc | sed 's/-//') -gt 10 ]; then
				echo "No reciprocal hits in first six blast hits; not processing further"
				exit 1
			else
			        ir2_firstrow=5
			fi
		else
		        ir2_firstrow=4
		fi
	else
	        ir2_firstrow=3
	fi
else
        ir2_firstrow=2
fi

ir2len=$(awk -v "row=$ir2_firstrow" ' NR==row {print $4} ' temp_blast2.out)
ir2_identity=$(awk -v "row=$ir2_firstrow" ' NR==row {print $3} ' temp_blast2.out)

if [ "$ir2len" -gt $(echo "$irlen" + 20 | bc) ]; then
	echo "found longer IR: $ir2len ident: $ir2_identity vs $irlen ident: $ir_identity"
	mv temp_cp2.fasta temp_cp.fasta && mv temp_blast2.out temp_blast.out
	ir_firstrow="$ir2_firstrow"
fi

rm db*


# b: find the locations of the IR boundaries and set start based on LSC

irb_start=$(awk -v "row=$ir_firstrow" ' NR==row {print $7} ' temp_blast.out)
irb_end=$(awk -v "row=$ir_firstrow" ' NR==row {print $8} ' temp_blast.out)
ira_start=$(awk -v "row=$ir_firstrow" ' NR==row {print $10} ' temp_blast.out)
ira_end=$(awk -v "row=$ir_firstrow" ' NR==row {print $9} ' temp_blast.out)


# check whether the hits are in the right order
if [ "$ira_start" -gt "$irb_start" ]; then
	echo "Hits in the wrong order! Check the output" && cp temp_blast.out REVIEW_ME_blast3.out
	exit 1
fi
#ira_start=$(awk ' NR==4 {print $7} ' temp_blast.out) && \
#ira_end=$(awk ' NR==4 {print $8} ' temp_blast.out) && \
#irb_start=$(awk ' NR==4 {print $10} ' temp_blast.out) && \
#irb_end=$(awk ' NR==4 {print $9} ' temp_blast.out); fi


# check if the piece between the IRs in this orientation is the LSC or SSC
if [ $(echo "$irb_start - $ira_end" | bc) -gt $(echo "$len / 3" | bc) ]; then
	new_start=$(echo "$ira_end + 1" | bc)
elif [ "$irb_end" -eq "$len" ]; then
	new_start=1
else
	new_start=$(echo "$irb_end + 1" | bc)
fi


if [ "$new_start" -ne 1 ]; then
	"$reform_script" temp_cp.fasta "$new_start" && mv new_contig.fasta temp_cp.fasta
fi



# c: regenerate the database and set locations of IR boundaries

makeblastdb -in temp_cp.fasta -out db -dbtype nucl -logfile temp.log
blastn -query temp_cp.fasta -db db -outfmt 6 | head > temp_blast.out

hit2=$(awk ' NR==2 {print $4} ' temp_blast.out)
hit3=$(awk ' NR==3 {print $4} ' temp_blast.out)
hit4=$(awk ' NR==4 {print $4} ' temp_blast.out)
hit5=$(awk ' NR==5 {print $4} ' temp_blast.out)
hit6=$(awk ' NR==6 {print $4} ' temp_blast.out)

if [ $(echo "$hit3" - "$hit2" | bc | sed 's/-//') -gt 10 ]; then
	echo "Potential problem with file and blast hits" && cp temp_blast.out REVIEW_ME_blast4.out
	if [ $(echo "$hit4" - "$hit3" | bc | sed 's/-//') -gt 10 ]; then
		if [ $(echo "$hit5" - "$hit4" | bc | sed 's/-//') -gt 10 ]; then
			if [ $(echo "$hit6" - "$hit5" | bc | sed 's/-//') -gt 10 ]; then
				echo "No reciprocal hits in first six blast hits; not processing further"
				exit 1
			else
				ir_firstrow=5
			fi
		else
			ir_firstrow=4
		fi
	else
		ir_firstrow=3
	fi
else
	ir_firstrow=2
fi

irlen=$(awk -v "row=$ir_firstrow" ' NR==row {print $4} ' temp_blast.out)
ir_identity=$(awk -v "row=$ir_firstrow" ' NR==row {print $3} ' temp_blast.out)

irb_start=$(awk -v "row=$ir_firstrow" ' NR==row {print $7} ' temp_blast.out)
irb_end=$(awk -v "row=$ir_firstrow" ' NR==row {print $8} ' temp_blast.out)
ira_start=$(awk -v "row=$ir_firstrow" ' NR==row {print $10} ' temp_blast.out)
ira_end=$(awk -v "row=$ir_firstrow" ' NR==row {print $9} ' temp_blast.out)

# check whether the hits are in the right order
if [ "$ira_start" -gt "$irb_start" ]; then
	echo "Hits in the wrong order! Check the output" && cp temp_blast.out REVIEW_ME_blast5.out
	exit 1
fi


## Step 2: extract the LSC; check orientation and location of ref1 in LSC; modify LSC as needed

"$extract_script" -c 1..$(echo "$ira_start - 1" | bc) temp_cp.fasta && mv extract.fasta LSC.fasta

blastn -query "$ref1" -db db -outfmt 6 > temp_blast_ref1.out
ref1_start=$(awk ' NR==1 {print $9} ' temp_blast_ref1.out)
ref1_end=$(awk ' NR==1 {print $10} ' temp_blast_ref1.out)
if [ "$ref1_end" -gt "$ref1_start" ]; then
	"$reform_script" LSC.fasta 1 yes && mv new_contig.fasta LSC.fasta
elif [ "$ref1_end" -gt 1000 ]; then
	echo "Reference hit >1000 bp from start; take care if that is not expected!"
fi


## Step 3: extract the SSC; check orientation of ref2 in SSC; modify SSC as needed

"$extract_script" -c $(echo "$ira_end + 1" | bc)..$(echo "$irb_start - 1" | bc) temp_cp.fasta && mv extract.fasta SSC.fasta

blastn -query "$ref2" -db db -outfmt 6 > temp_blast_ref2.out
ref2_start=$(awk ' NR==1 {print $9} ' temp_blast_ref2.out)
ref2_end=$(awk ' NR==1 {print $10} ' temp_blast_ref2.out)		# this may only be a partial hit, but all I need is orientation
if [ "$ref2_end" -lt "$ref2_start" ]; then
	"$reform_script" SSC.fasta 1 yes && mv new_contig.fasta SSC.fasta
fi


## Step 4: extract the IRs, stitch all back together and remove intermediate files

"$extract_script" -c "$ira_start".."$ira_end" temp_cp.fasta && mv extract.fasta IRA.fasta
"$extract_script" -c "$irb_start".."$irb_end" temp_cp.fasta && mv extract.fasta IRB.fasta

cat LSC.fasta IRA.fasta SSC.fasta IRB.fasta > new_cp.fasta

"$merge_script" new_cp.fasta && mv merged.fasta new_cp.fasta

rm temp*.out temp*.fasta temp.log temp.perf db* IR* LSC.fasta SSC.fasta
