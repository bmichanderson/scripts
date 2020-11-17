#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 29 Nov 2018
# Modified: Oct 2020
# Description: run ABySS assemblies on UPPMAX (execute this script with >sbatch [script] [args])
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J abyss_assembly


# set the OMP_NUM_THREADS variable to the number available
export OMP_NUM_THREADS="$SLURM_NPROCS"


# function for when the script is called incorrectly
usage()
{
	echo "Provide arguments: infile_1 infile_2 k (k_max) (k_step, default 2)" 1>&2
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
	kmax="$3"
	kstep=2
elif [ -z "$5" ]; then
	kstep=2
fi


# assign command line arguments
in1="$1"
in2="$2"
kmer="$3"
if [ -n "$4" ]; then
	kmax="$4"
	if [ "$kmax" -lt "$kmer" ]; then
		echo "Kmax must be larger than or equal to k"
		exit 1
	fi
fi
if [ -n "$5" ]; then
	kstep="$5"
	if [ "$kstep" -le 0 ]; then
		echo "Kstep must be greater than 0"
		exit 1
	fi
fi

echo "Starting at $(date) with the following input parameters:"
echo "infile_1 = $in1 infile_2 = $in2, kmer = $kmer, kmax = $kmax, kstep = $kstep"
echo "Loading modules"


# load modules (first the module for bioinformatic tools)
module load bioinfo-tools
module load bbmap/38.08
module load abyss/2.0.2-k128


# assign working directory to scratch
work_dir="$SNIC_TMP"


# merge the two input files into interleaved format, and put in the temporary directory
echo "Merging files and moving to working directory"
reformat.sh in1="$in1" in2="$in2" out="$work_dir"/input_interleaved.fastq.gz


# change to the working directory and start assemblies
cd "$work_dir"

for (( index=kmer; index<=kmax; index+=kstep )); do
	echo "Running ABySS for k = $index"
	mkdir k"$index"
	abyss-pe -C k"$index" name=$( basename -s _1.fastq.gz "$in1" ) k="$index" in="${work_dir}/input_interleaved.fastq.gz" np="$SLURM_NPROCS" q=15
done


# compute comparative assembly statistics
echo "Computing assembly statistics"
abyss-fac k*/$( basename -s _1.fastq.gz "$in1" )-contigs.fa


# delete input file and copy output to the storage directory
rm input_interleaved.fastq.gz
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
