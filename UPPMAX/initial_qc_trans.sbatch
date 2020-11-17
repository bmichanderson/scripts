#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 9 Jul 2019
# Modified: Oct 2020
# Description: run a series of qc steps on raw Illumina reads from transcriptome sequencing
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J trans_initial_qc


# set default input type
input_type="paired"


# define a function for when the script is called incorrectly
usage()
{
	echo "Usage: script.sbatch infile_1 [infile_2 for paired]" 1>&2
	exit 1
}

if [ $# -eq 0 ]; then
	usage
elif [ -z "$1" ]; then
	usage
elif [ -z "$2" ]; then
	input_type="single"
fi


# report command line arguments
if [ "$input_type" = "single" ]; then

	echo "Starting at $(date) with the following input parameters:"
	echo "input_type = $input_type infile_1 = $1"

else

	echo "Starting at $(date) with the following input parameters:"
	echo "input_type = $input_type infile_1 = $1 infile_2 = $2"

fi


# load modules (first the module for bioinformatic tools)
echo "Loading modules"

module load bioinfo-tools
module load bbmap/38.08 FastQC/0.11.5


# assign working directory to scratch
work_dir="$SNIC_TMP"


# copy input files to the temporary directory
echo "Copying files to working directory"

cp "$1" "$work_dir"

if [ "$input_type" = "paired" ]; then
	cp "$2" "$work_dir"
fi


# assign file names from (possibly) path arguments
file1="$( basename $1 )"

if [ "$input_type" = "paired" ]; then
	file2="$( basename $2 )"
fi


# change to the working directory and start qc
cd "$work_dir"


# run FastQC on the raw input
echo "Running initial FastQC"

if [ ! -d "fastqc" ]; then
	mkdir -p "fastqc"
fi

if [ "$input_type" = "paired" ]; then

	fastqc -t 2 -o fastqc "$file1" "$file2"

else

	fastqc -o fastqc "$file1"

fi


# remove duplicates
echo "Removing duplicates with clumpify"

if [ "$input_type" = "paired" ]; then

	clumpify.sh in1="$file1" in2="$file2" out1=dedup_"$file1" out2=dedup_"$file2" dedupe

else

	clumpify.sh in="$file1" out=dedup_"$file1" dedupe

fi


# trim adapters, quality, first 9 bp, and PhiX contamination
echo "Trimming adapters, for quality, and for PhiX contamination with bbduk"

if [ "$input_type" = "paired" ]; then

	bbduk.sh in1=dedup_"$file1" in2=dedup_"$file2" out1=temp1_1.fastq.gz out2=temp1_2.fastq.gz ktrim=r literal="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA","AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" k=21 hdist=2 hdist2=1 mink=15 minlength=50 tbo tpe

	bbduk.sh in1=temp1_1.fastq.gz in2=temp1_2.fastq.gz out1=temp2_1.fastq.gz out2=temp2_2.fastq.gz literal="AGATCGGAAGAGCAC","AGATCGGAAGAGCGT" ktrim=r k=8 restrictright=15 minlength=50

	bbduk.sh in1=temp2_1.fastq.gz in2=temp2_2.fastq.gz out1=temp3_1.fastq.gz out2=temp3_2.fastq.gz literal="AGATCGGA","AGATCGGA" ktrim=r k=6 restrictright=8 minlength=50 maxns=0 qtrim=r trimq=20

	bbduk.sh in1=temp3_1.fastq.gz in2=temp3_2.fastq.gz out1=trim_dedup_"$file1" out2=trim_dedup_"$file2" ref=phix k=31 hdist=2

else

	bbduk.sh in=dedup_"$file1" out=temp1.fastq.gz literal="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ktrim=r k=21 hdist=2 hdist2=1 mink=15 minlength=50 

	bbduk.sh in=temp1.fastq.gz out=temp2.fastq.gz literal="AGATCGGAAGAGCAC" ktrim=r k=8 restrictright=15 minlength=50

	bbduk.sh in=temp2.fastq.gz out=temp3.fastq.gz literal="AGATCGGA" ktrim=r k=6 restrictright=8 minlength=50 maxns=0 qtrim=r trimq=20

	bbduk.sh in=temp3.fastq.gz out=trim_dedup_"$file1" ref=phix k=31 hdist=2

fi


# normalize
# bbnorm uses kmer depth, with default k=31
# so, for 150 bp reads and default kmer depth of 100x, target read depth is Dr = 100 * 150 / (150-31+1) = 125
echo "Normalizing data for 100x kmer coverage, min 5"

if [ "$input_type" = "paired" ]; then

	bbnorm.sh in1=trim_dedup_"$file1" in2=trim_dedup_"$file2" out1=norm_trim_dedup_"$file1" out2=norm_trim_dedup_"$file2" target=100 mindepth=5

else

	bbnorm.sh in=trim_dedup_"$file1" out=norm_trim_dedup_"$file1" target=100 mindepth=5

fi


# run FastQC on the two main outputs
echo "Running final FastQC"

if [ ! -d "fastqc" ]; then
	mkdir -p "fastqc"
fi

if [ "$input_type" = "paired" ]; then

	fastqc -t 4 -o fastqc trim_dedup_"$file1" trim_dedup_"$file2" norm_trim_dedup_"$file1" norm_trim_dedup_"$file2"

else

	fastqc -t 2 -o fastqc trim_dedup_"$file1" norm_trim_dedup_"$file1"

fi


# remove raw files and intermediate dedup files and copy remainder to the storage directory
rm "$file1"

if [ "$input_type" = "paired" ]; then
	rm "$file2"
fi

rm dedup* temp*

echo "Copying output to the storage directory under ID: $SLURM_JOB_ID"

out_dir=/proj/snic2020-16-240/runs/"$SLURM_JOB_ID"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

cp -a . "$out_dir"


# copy CPU/RAM usage of the job to the storage directory
echo "Copying usage file to the storage directory"

node=$( hostname -s )
if [[ "$node" == r* ]]; then
	system=rackham
elif [[ "$node" == s* ]]; then
	system=snowy
else
	system=rackham
fi

cp /sw/share/slurm/"$system"/uppmax_jobstats/"$node"/"$SLURM_JOB_ID" "$out_dir"

echo "Finished at $(date)"
