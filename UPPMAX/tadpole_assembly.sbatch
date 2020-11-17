#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 22 January 2020
# Modified: Oct 2020
# Description: run an initial Tadpole assembly
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J tadpole


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


# set default input type (note that single merged reads or interleaved reads can also be used, since only extracting kmers)
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
echo "Loading modules"
echo

module load bioinfo-tools
module load bbmap/38.08


# assign working directory to scratch
work_dir="$SNIC_TMP"


# copy input files to the temporary directory
echo
echo "Copying files to working directory"
echo

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


# assemble for multiple values of K
echo
echo "Assembling for different values of K to find the optimum"
echo

if [ "$input_type" = "paired" ]; then

	tadwrapper.sh in="$file1","$file2" out=tadwrapper_contigs_%.fasta outfinal=tadwrapper_contigs k=40,81,102,124,217 bisect mincountretain=5 prefilter=2 prealloc=t -Xmx"$max_mem"

else

	tadwrapper.sh in="$file1" out=tadwrapper_contigs_%.fasta outfinal=tadwrapper_contigs k=40,81,102,124,217 bisect mincountretain=5 prefilter=2 prealloc=t -Xmx"$max_mem"

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
