#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 24 Feb 2020
# Modified: Oct 2020
# Description: assemble PacBio long reads using CANU
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J canu_assembly


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


# a function for when the script is called incorrectly or without arguments
usage()
{
	echo "This script will assemble long reads using the full CANU pipeline."
	echo
	echo "Usage: $(basename $0) options(- ... - ...) additional command line args"
	echo
	echo "Options:"
	echo
	echo "	-l	File with long reads (concatenate if you have multiple)"
	echo
	echo "	-g	Genome size estimate (including a k, m or g for kilobases, megabases or gigabases, without a space e.g. 550m)"
	echo
	exit 1
}


# parse the command line
if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi

reads="absent"
genome_est="absent"



FILES=()

while [[ $# -gt 0 ]]		# while the number of args is greater than 0
do
key="$1"

case $key in
	-l)
	reads="$2"
	shift
	shift
	;;
	-g)
	genome_est="$2"
	shift
	shift
	;;
	*)
	FILES+=("$1")		# capture arguments not connected to options in an array
	shift
	;;
esac
done

set -- "${FILES[@]}"		# restore the non-option args to command line args


if [ "$reads" = "absent" ]; then
	usage
fi

if [ "$genome_est" = "absent" ]; then
	usage
fi


# report command line arguments
echo "Starting at $(date) with the following input parameters:"
echo

echo "reads = $reads"
echo "genome estimate = $genome_est"
echo "additional arguments = $@"


# load modules
echo
echo "Loading modules"
echo

module load bioinfo-tools canu/1.8


# assign working directory to scratch
work_dir="$SNIC_TMP"


# copy input files to the temporary directory
echo
echo "Copying files to working directory"
echo

cp "$reads" "$work_dir"


# change to the working directory and start assembly
cd "$work_dir"


# run the full CANU pipeline, including error correction, read trimming and assembly
canu -p long_ass -d assembly genomeSize="$genome_est" -pacbio-raw "$( basename $reads)" "$@"


# remove raw files and copy remainder to the storage directory
echo
echo "Removing raw files and copying results to the storage directory under ID: $SLURM_JOB_ID"
echo

rm "$( basename $reads)"

out_dir=/proj/uppstore2018195/runs/"$SLURM_JOB_ID"

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
