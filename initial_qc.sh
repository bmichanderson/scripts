#!/bin/bash

##############
# Author: Ben Anderson
# Date: Nov 2020
# Modified:
# Description: perform a series of initial QC steps on Illumina paired reads
##############


# Define a help function for an incorrect or empty call
help()
{
(	echo -e '\n''\t'"This script runs a series of QC steps on Illumina paired reads"
	echo -e '\t'"It requires BBMap scripts in the PATH, as well as fastqc"
	echo -e '\t'"It will create two corrected read files, as well as a merged and a pair of unmerged files"
	echo -e '\t'"It will also summarize read lengths in the merged and unmerged files, and create a log qc.log"'\n'
	echo -e '\t'"Usage: initial_qc.sh file1 file2"'\n'
) 1>&2
	exit 1
}


# Check if there are the correct number of command arguments
if [ $# -ne 2 ]; then
	help
fi


# Create a log file
echo -e "Quality control log file"'\n' > qc.log


# Report starting time
echo "Starting QC run at $(date) with the following input:"
echo "Working directory: $(pwd)"
echo "file1 = $1"
echo "file2 = $2"

start=$(date +%s)


# Step 1: Remove duplicates
echo -e '\n'"****"'\n'"Removing duplicates with clumpify"'\n'"****"'\n' | tee -a qc.log

clumpify.sh in1="$1" in2="$2" out1=dedup_1.fastq.gz out2=dedup_2.fastq.gz dedupe |& tee -a qc.log


# Step 2: Trim adapters, for quality and PhiX contamination
echo -e '\n'"****"'\n'"Trimming adapters and for quality and PhiX"'\n'"****"'\n' | tee -a qc.log

bbduk.sh in1=dedup_1.fastq.gz in2=dedup_2.fastq.gz out1=temp1_1.fastq.gz out2=temp1_2.fastq.gz ktrim=r \
literal="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA","AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" k=21 hdist=2 hdist2=1 mink=15 \
minlength=50 tbo tpe |& tee -a qc.log && \
bbduk.sh in1=temp1_1.fastq.gz in2=temp1_2.fastq.gz out1=temp2_1.fastq.gz out2=temp2_2.fastq.gz ktrim=r \
literal="AGATCGGAAGAGCAC","AGATCGGAAGAGCGT" k=8 restrictright=15 minlength=50 |& tee -a qc.log && \
bbduk.sh in1=temp2_1.fastq.gz in2=temp2_2.fastq.gz out1=temp3_1.fastq.gz out2=temp3_2.fastq.gz ktrim=r \
literal="AGATCGGA","AGATCGGA" k=6 restrictright=8 minlength=50 maxns=0 qtrim=r trimq=20 |& tee -a qc.log && \
bbduk.sh in1=temp3_1.fastq.gz in2=temp3_2.fastq.gz out1=trim_dedup_1.fastq.gz out2=trim_dedup_2.fastq.gz ref=phix k=31 hdist=2 |& tee -a qc.log && \
rm temp* dedup*


# Step 3: Interleave reads (for convenience) and correct sequencing errors
echo -e '\n'"****"'\n'"Correcting sequencing errors"'\n'"****"'\n' | tee -a qc.log

reformat.sh in1=trim_dedup_1.fastq.gz in2=trim_dedup_2.fastq.gz out=temp.fastq.gz |& tee -a qc.log && \
bbmerge.sh in=temp.fastq.gz out=ecco.fastq.gz ecco mix vstrict ordered prefilter=2 prealloc=t minlength=50 && \
mv ecco.fastq.gz temp.fastq.gz |& tee -a qc.log && \
clumpify.sh in=temp.fastq.gz out=eccc.fastq.gz ecc passes=4 && \
mv eccc.fastq.gz temp.fastq.gz |& tee -a qc.log && \
tadpole.sh in=temp.fastq.gz out=ecct.fastq.gz ecc k=62 ordered prefilter=2 prealloc=t errormult1=64 && \
mv ecct.fastq.gz temp.fastq.gz |& tee -a qc.log


# Step 4: Merge overlapping reads then limit read length and rename
echo -e '\n'"****"'\n'"Merging overlapping reads"'\n'"****"'\n' | tee -a qc.log

bbmerge-auto.sh in=temp.fastq.gz out=temp_merged.fastq.gz outu=temp_unmerged.fastq.gz strict \
k=93 extend2=80 rem ordered prefilter=2 prealloc=t |& tee -a qc.log && \
bbduk.sh in=temp_merged.fastq.gz out=merged.fastq.gz minlength=50 |& tee -a qc.log && \
bbduk.sh in=temp_unmerged.fastq.gz out1=unmerged_1.fastq.gz out2=unmerged_2.fastq.gz minlength=50 |& tee -a qc.log && \
bbduk.sh in=temp.fastq.gz out1=corrected_1.fastq.gz out2=corrected_2.fastq.gz minlength=50 |& tee -a qc.log && \
rm temp* trim_dedup*


# Step 5: Run a fastqc on the output
echo -e '\n'"****"'\n'"Running a fastqc on the output read files"'\n'"****"'\n' | tee -a qc.log

if [ ! -d "fastqc" ]; then
	mkdir "fastqc";
fi

fastqc -t 5 -o fastqc merged.fastq.gz unmerged_1.fastq.gz unmerged_2.fastq.gz corrected_1.fastq.gz corrected_2.fastq.gz


# Step 6: Assess read length and number of reads for the output
echo -e '\n'"****"'\n'"Assessing read length and number of reads for output"'\n'"****"'\n' | tee -a qc.log

readlength.sh in1=unmerged_1.fastq.gz in2=unmerged_2.fastq.gz out=readlength_qc_unmerged.txt |& tee -a qc.log && \
readlength.sh in=merged.fastq.gz out=readlength_qc_merged.txt |& tee -a qc.log && \
readlength.sh in1=corrected_1.fastq.gz in2=corrected_2.fastq.gz out=readlength_qc_corrected.txt |& tee -a qc.log


# Report time taken for the QC steps
end=$(date +%s)
duration=$(( $end - $start ))
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e '\n'"Finished at $(date), after running for $duration seconds, or $duration_mins minutes, or $duration_hours hours" | tee -a qc.log
