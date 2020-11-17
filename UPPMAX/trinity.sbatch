#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 15 Jul 2019
# Modified: Oct 2020
# Description: assemble a transcriptome from Illumina reads (single or paired-end); note: the single assumes reverse strand specific
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J trinity_run


# set default input type
input_type="paired"


# set maximum memory available for Trinity
max_mem=$(($SLURM_NPROCS * 6))


# assign command line arguments
if [ -z "$1" ]; then

	exit 1

elif [ -z "$2" ]; then

	input_type="single"
	file1="$( basename $1 )"
	echo "Starting at $(date) with the following input parameters:"
	echo "input_type = $input_type infile_1 = $file1"

else

	file1="$( basename $1 )"
	file2="$( basename $2 )"
	echo "Starting at $(date) with the following input parameters:"
	echo "input_type = $input_type infile_1 = $file1 infile_2 = $file2"

fi


# assign working directory to scratch
work_dir="$SNIC_TMP"


# copy input files to the temporary directory
echo "Copying files to working directory"

cp "$1" "$work_dir"

if [ "$input_type" = "paired" ]; then
	cp "$2" "$work_dir"
fi

cp /proj/snic2020-16-240/software/singularity-containers/trinity.img "$work_dir"


# change to the working directory
cd "$work_dir"


# start the Trinity run
echo "Starting Trinity run"

if [ "$input_type" = "paired" ]; then

	singularity run trinity.img --seqType fq --max_memory "$max_mem"G --left "$file1" --right "$file2" \
		--CPU "$SLURM_NPROCS" --output trinity_output --full_cleanup --no_version_check

else

	singularity run trininty.img --seqType fq --max_memory "$max_mem"G --single "$file1" --SS_lib_type R \
		--CPU "$SLURM_NPROCS" --output trinity_output --full_cleanup --no_version_check

fi


# remove raw files and copy remainder to the storage directory
rm "$file1"

if [ "$input_type" = "paired" ]; then
	rm "$file2"
fi

rm trinity.img

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
