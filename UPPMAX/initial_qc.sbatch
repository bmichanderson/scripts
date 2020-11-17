#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 28 Jan 2020
# Modified: Nov 2020
# Description: run a series of qc steps on raw Illumina reads
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 1-00:00:00
#SBATCH -J initial_qc


# determine max memory to specify
# if you have problems with memory, you can request 256G nodes or also 512G nodes on Snowy and a couple 1 TB nodes on Rackham
# use -C mem256GB

num_cpus="$SLURM_NPROCS"

if [ -n "$SLURM_MEM_PER_NODE" ]; then

	memory="$SLURM_MEM_PER_NODE"

elif [ -n "$SLURM_MEM_PER_CPU" ]; then

	mem_per_cpu="$SLURM_MEM_PER_CPU"
	memory=$((mem_per_cpu * num_cpus))

fi

max_mem=$((memory * 85 / 100000))g


# report starting values
echo
echo "Number of cpus: $num_cpus"
echo "Memory available: $memory"
echo "Max mem to assign: $max_mem"
echo

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
	echo

else

	echo "Starting at $(date) with the following input parameters:"
	echo "input_type = $input_type infile_1 = $1 infile_2 = $2"
	echo

fi


# load modules (first the module for bioinformatic tools)
echo
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


# Create a log file
echo -e "Quality control log file"'\n' > qc.log


# run FastQC on the raw input
echo -e '\n'"Running initial FastQC"'\n' | tee -a qc.log

if [ ! -d "fastqc" ]; then
	mkdir "fastqc"
fi

if [ "$input_type" = "paired" ]; then

	fastqc -t 2 -o fastqc "$file1" "$file2"

else

	fastqc -o fastqc "$file1"

fi


# remove duplicates
echo -e '\n'"Removing duplicates with clumpify"'\n' | tee -a qc.log

if [ "$input_type" = "paired" ]; then

	clumpify.sh in1="$file1" in2="$file2" out1=dedup_"$file1" out2=dedup_"$file2" -Xmx"$max_mem" dedupe |& tee -a qc.log

else

	clumpify.sh in="$file1" out=dedup_"$file1" -Xmx"$max_mem" dedupe |& tee -a qc.log

fi


# trim adapters, quality, and possible PhiX contamination
echo -e '\n'"Trimming adapters, for quality, and for PhiX contamination with bbduk"'\n' | tee -a qc.log

if [ "$input_type" = "paired" ]; then

	bbduk.sh in1=dedup_"$file1" in2=dedup_"$file2" out1=temp1_1.fastq.gz out2=temp1_2.fastq.gz ktrim=r \
	literal="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA","AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" k=21 hdist=2 hdist2=1 mink=15 \
	minlength=50 tbo tpe |& tee -a qc.log && \
	bbduk.sh in1=temp1_1.fastq.gz in2=temp1_2.fastq.gz out1=temp2_1.fastq.gz out2=temp2_2.fastq.gz ktrim=r \
	literal="AGATCGGAAGAGCAC","AGATCGGAAGAGCGT" k=8 restrictright=15 minlength=50 |& tee -a qc.log && \
	bbduk.sh in1=temp2_1.fastq.gz in2=temp2_2.fastq.gz out1=temp3_1.fastq.gz out2=temp3_2.fastq.gz ktrim=r \
	literal="AGATCGGA","AGATCGGA" k=6 restrictright=8 minlength=50 maxns=0 qtrim=r trimq=20 |& tee -a qc.log && \
	bbduk.sh in1=temp3_1.fastq.gz in2=temp3_2.fastq.gz out1=trim_dedup_"$file1" out2=trim_dedup_"$file2" ref=phix k=31 hdist=2 |& tee -a qc.log

else
	bbduk.sh in=dedup_"$file1" out=temp_trim1.fastq.gz ktrim=r \
	literal="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" k=21 hdist=2 hdist2=1 mink=15 minlength=50 |& tee -a qc.log && \
	bbduk.sh in=temp_trim1.fastq.gz out=temp_trim2.fastq.gz ktrim=r \
	literal="AGATCGGAAGAGCAC" k=8 restrictright=15 minlength=50 |& tee -a qc.log && \
	bbduk.sh in=temp_trim2.fastq.gz out=temp_trim3.fastq.gz ktrim=r \
	literal="AGATCGGA" k=6 restrictright=8 minlength=50 maxns=0 qtrim=r trimq=20 |& tee -a qc.log && \
	bbduk.sh in=temp_trim3.fastq.gz out=trim_dedup_"$file1" ref=phix k=31 hdist=2 |& tee -a qc.log

fi

rm temp* dedup*


# format reads as interleaved if paired
if [ "$input_type" = "paired" ]; then
	echo
	reformat.sh in1=trim_dedup_"$file1" in2=trim_dedup_"$file2" out=temp.fastq.gz |& tee -a qc.log
else
	ln -s trim_dedup_"$file1" temp.fastq.gz
fi


# error correction
echo -e '\n'"Correcting sequencing errors with bbmerge (for paired), clumpify and tadpole"'\n' | tee -a qc.log

if [ "$input_type" = "paired" ]; then

	bbmerge.sh in=temp.fastq.gz out=ecco.fastq.gz ecco mix vstrict ordered prefilter=2 prealloc=t minlength=50 -Xmx"$max_mem" |& tee -a qc.log && \
	mv ecco.fastq.gz temp.fastq.gz
	echo

fi

clumpify.sh in=temp.fastq.gz out=eccc.fastq.gz ecc passes=4 -Xmx"$max_mem" |& tee -a qc.log && \
mv eccc.fastq.gz temp.fastq.gz && echo && \
tadpole.sh in=temp.fastq.gz out=ecct.fastq.gz ecc k=62 ordered prefilter=2 prealloc=t -Xmx"$max_mem" errormult1=64 |& tee -a qc.log && \
mv ecct.fastq.gz temp.fastq.gz


if [ "$input_type" = "paired" ]; then

	# merging overlapping reads
	echo -e '\n'"Merging overlapping reads"'\n' | tee -a qc.log

	bbmerge-auto.sh in=temp.fastq.gz out=temp_merged.fastq.gz outu=temp_unmerged.fastq.gz strict k=93 \
	extend2=80 rem ordered prefilter=2 prealloc=t -Xmx"$max_mem" |& tee -a qc.log

fi


# limit the read length of output from merging and from correction, and rename corrected reads split into two files
echo -e '\n'"Limiting read length"'\n'

if [ "$input_type" = "paired" ]; then

	bbduk.sh in=temp_merged.fastq.gz out=merged.fastq.gz minlength=50 |& tee -a qc.log && \
	bbduk.sh in=temp_unmerged.fastq.gz out1=unmerged_1.fastq.gz out2=unmerged_2.fastq.gz minlength=50 |& tee -a qc.log && \
	bbduk.sh in=temp.fastq.gz out1=corrected_1.fastq.gz out2=corrected_2.fastq.gz minlength=50 |& tee -a qc.log

else

	bbduk.sh in=temp.fastq.gz out=corrected.fastq.gz minlength=50 |& tee -a qc.log

fi

rm temp* trim_dedup*

# run FastQC on the output
echo -e '\n'"Running a second FastQC"'\n'

if [ "$input_type" = "paired" ]; then

	fastqc -t 5 -o fastqc merged.fastq.gz unmerged_1.fastq.gz unmerged_2.fastq.gz corrected_1.fastq.gz corrected_2.fastq.gz

else

	fastqc -o fastqc corrected.fastq.gz

fi


# Assess read length and number of reads for the output
echo -e '\n'"Assessing read length and number of reads for output"'\n' | tee -a qc.log

if [ "$input_type" = "paired" ]; then

	readlength.sh in1=unmerged_1.fastq.gz in2=unmerged_2.fastq.gz out=readlength_qc_unmerged.txt |& tee -a qc.log && \
	readlength.sh in=merged.fastq.gz out=readlength_qc_merged.txt |& tee -a qc.log && \
	readlength.sh in1=corrected_1.fastq.gz in2=corrected_2.fastq.gz out=readlength_qc_corrected.txt |& tee -a qc.log

else

	readlength.sh in=corrected.fastq.gz out=readlength_qc_corrected.txt |& tee -a qc.log

fi


# remove raw files and copy remainder to the storage directory
rm "$file1"

if [ "$input_type" = "paired" ]; then
	rm "$file2"
fi


echo
echo "Copying output to the storage directory under ID: $SLURM_JOB_ID"

out_dir=/proj/snic2020-16-240/runs/"$SLURM_JOB_ID"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

cp -a . "$out_dir"


# copy CPU/RAM usage of the job to the storage directory
echo
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

echo
echo "Finished at $(date)"
