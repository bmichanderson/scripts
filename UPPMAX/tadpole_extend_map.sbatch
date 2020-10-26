#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 6 February 2020
# Modified: Oct 2020
# Description: extend contigs using reads with Tadpole, then map the reads to the extended
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J extend_tadpole


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
usage ()
{
	echo "This script will extend contigs with given reads using Tadpole, then map those reads to the new contigs."
	echo
	echo "Usage: $(basename $0) options(- ... - ...)"
	echo
	echo "Options:"
	echo
	echo "	-c	File of contigs for extension"
	echo
	echo "	-1	File of forward reads (or interleaved paired, but then don't specify reverse)"
	echo
	echo "	-2	File of reverse reads"
	echo
	echo "	-s	File with single reads"
	echo
	exit 1
}

# function to join elements
function join_by { local IFS="$1"; shift; echo "$*"; }


# parse the command line
if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi

contigs="absent"
forward="absent"
reverse="absent"
single="absent"


FILES=()

while [[ $# -gt 0 ]]		# while the number of args is greater than 0
do
key="$1"

case $key in
	-c)
	contigs="$2"
	shift
	shift
	;;
	-1)
	forward="$2"
	shift
	shift
	;;
	-2)
	reverse="$2"
	shift
	shift
	;;
	-s)
	single="$2"
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


if [ "$contigs" = "absent" ]; then
	usage
fi


# report command line arguments
echo "Starting at $(date) with the following input parameters:"
echo

echo "contigs = $contigs"
echo "forward = $forward"
echo "reverse = $reverse"
echo "single = $single"


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

if [ "$contigs" != "absent" ]; then
	cp "$contigs" "$work_dir"
else
	usage
fi

if [ "$forward" != "absent" ]; then
	cp "$forward" "$work_dir"
fi

if [ "$reverse" != "absent" ]; then
	cp "$reverse" "$work_dir"
fi

if [ "$single" != "absent" ]; then
	cp "$single" "$work_dir"
fi


# change to the working directory and start extension
cd "$work_dir"


# create a tadpole command line
command_line="in=$( basename $contigs ) "
extra=""

if [ "$forward" != "absent" ]; then
	extra+="$( basename $forward ) "
fi

if [ "$reverse" != "absent" ]; then
	extra+="$( basename $reverse ) "
fi

if [ "$single" != "absent" ]; then
	extra+="$( basename $single ) "
fi

extra_joined=$(join_by , $extra)

command_line+="extra=$extra_joined "

command_line+="out=extended.fasta el=100000 er=100000 mode=extend k=31 -Xmx$max_mem"


# extend contigs with provided kmers from reads
echo
echo "Extending contigs"
echo

tadpole.sh "$command_line"


# map reads to the new contigs
echo
echo "Mapping reads to the new contigs"
echo

if [ "$forward" != "absent" ]; then

	if [ "$reverse" != "absent" ]; then

		bbmap.sh ref=extended.fasta in1="$( basename $forward )" in2="$( basename $reverse )" outm=mapu#.fastq.gz nodisk overwrite=t

	else	# perhaps interleaved paired reads

		bbmap.sh ref=extended.fasta in="$( basename $forward )" outm=mapu.fastq.gz nodisk overwrite=t

	fi

fi

if [ "$single" != "absent" ]; then

	bbmap.sh ref=extended.fasta in="$( basename $single )" outm=mapm.fastq.gz nodisk overwrite=t

fi


# remove raw files and copy remainder to the storage directory
if [ "$forward" != "absent" ]; then
	rm "$( basename $forward )"
fi

if [ "$reverse" != "absent" ]; then
	rm "$( basename $reverse )"
fi

if [ "$single" != "absent" ]; then
	rm "$( basename $single )"
fi

rm "$( basename $contigs )"


echo
echo "Copying output to the storage directory under ID: $SLURM_JOB_ID"

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
