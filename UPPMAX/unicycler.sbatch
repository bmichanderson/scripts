#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 6 February 2020
# Modified: Oct 2020
# Description: run a Unicycler analysis on UPPMAX
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J unicycler


# a function for when the script is called incorrectly or without arguments
usage()
{
	echo "This script will run a Unicycler assembly."
	echo
	echo "Usage: $(basename $0) options(- ... - ...)"
	echo
	echo "Options:"
	echo
	echo "	-1	File of forward reads"
	echo
	echo "	-2	File of reverse reads"
	echo
	echo "	-s	File with single reads"
	echo
	echo "	-l	File with long reads"
	echo
	exit 1
}

if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi

forward="absent"
reverse="absent"
single="absent"
long="absent"


FILES=()

while [[ $# -gt 0 ]]		# while the number of args is greater than 0
do
key="$1"

case $key in
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
	-l)
	long="$2"
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


# report command line arguments
echo "Starting at $(date) with the following input parameters:"
echo

echo "forward = $forward"
echo "reverse = $reverse"
echo "single = $single"
echo "long = $long"


# assign working directory to scratch
work_dir="$SNIC_TMP"


# copy input files and unicycler singularity container to the working directory
echo
echo "Copying files to working directory"
echo

if [ "$forward" != "absent" ]; then
	cp "$forward" "$work_dir"
fi

if [ "$reverse" != "absent" ]; then
	cp "$reverse" "$work_dir"
fi

if [ "$single" != "absent" ]; then
	cp "$single" "$work_dir"
fi

if [ "$long" != "absent" ]; then
	cp "$long" "$work_dir"
fi

cp /proj/uppstore2018195/software/singularity-containers/unicycler.img "$work_dir"


# change to the working directory and start unicycler
cd "$work_dir"


# create a unicycler command line
command_line=""

if [ "$forward" != "absent" ]; then
	command_line+="-1 $( basename $forward ) "
fi

if [ "$reverse" != "absent" ]; then
	command_line+="-2 $( basename $reverse ) "
fi

if [ "$single" != "absent" ]; then
	command_line+="-s $( basename $single ) "
fi

if [ "$long" != "absent" ]; then
	command_line+="-l $( basename $long ) "
fi

command_line+="-o . --threads $SLURM_NPROCS --no_pilon --no_rotate --keep 0 --verbosity 2 --no_correct --min_component_size 250 --min_dead_end_size 250"


echo
echo "Command line is:"
echo
echo "$command_line"
echo


# run unicycler
singularity exec unicycler.img unicycler $command_line


# remove input files and copy output to the storage directory
rm *.gz unicycler.img

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
