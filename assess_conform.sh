#!/bin/bash

##########################
# Author: B. Anderson
# Date: 12 Feb 2019
# Modified: Oct 2020
# Description: given an input contig file (arg1) and coordinates (args 2 3) for a repeat shared elsewhere, extract the region, map reads (arg 4) and filter for those that span repeat, then align
##########################

# paths to python scripts called
extract_script=~/scripts/fasta_extract.py
filt_bam_script=~/scripts/bam_filter_dualclip.py
subsamp_script=~/scripts/fasta_subsample.py


# a function for when the script is called incorrectly or without arguments
usage ()
{
	echo "This script will provide an alignment of long reads for a repeat region in an input contig."
	echo "It requires MAFFT, minimap2 and Samtools to be in the path, and Biopython and Pysam."
	echo
	echo "Usage: $(basename $0) contig_fasta start_position end_position long_reads"
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
elif [ -z "$4" ]; then
	usage
fi


# assign command line arguments
contig_fasta="$1"
start="$2"
end="$3"
long_reads="$4"


# length of the region
length=$(($end - $start + 1))


# extract the region
"$extract_script" "$contig_fasta" -c "$start".."$end"


# map long reads to the extracted region
minimap2 -ax map-pb -t 6 extract.fasta "$long_reads" > extract_minimap.sam


# filter the resulting sam file and produce a bam file with accompanying index (2308 = 4:mapped + 256:secondary + 2048:chimeric)
samtools view -bS -F 2308 extract_minimap.sam | samtools sort -o extract_minimap_sorted.bam
samtools index extract_minimap_sorted.bam


# query the resulting bam file in Python to only keep reads which span the contig and are soft clipped
"$filt_bam_script" extract_minimap_sorted.bam extract.fasta no


# (optional) sort and index the resulting output for evaluation in IGV
#samtools sort reads_out.bam -o reads_out_sorted.bam
#samtools index reads_out_sorted.bam


# subsample the reads for alignment so that there are no more than 50 (can adjust as you prefer)
"$subsamp_script" reads_out.fasta 50


# align the subsampled reads and the reference extract
# can set a different number of threads if desired
mafft --thread 6 --auto subsampled_reads.fasta > alignment.fasta


# remove the intermediate files
rm extract* reads_out.fasta
