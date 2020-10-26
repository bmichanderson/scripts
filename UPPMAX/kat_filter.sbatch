#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 20 Feb 2020
# Modified: Oct 2020
# Description: filter reads based on kmer depth on UPPMAX
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J kmer_filter


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
	echo "This script will filter reads based on a kmer depth range."
	echo
	echo "Usage: $(basename $0) options(- ... - ...)"
	echo
	echo "Options:"
	echo
	echo "	-l	lower kmer depth"
	echo
	echo "	-h	higher kmer depth"
	echo
	echo "	-1	File of forward reads (or interleaved paired, but then don't specify reverse)"
	echo
	echo "	-2	File of reverse reads"
	echo
	echo "	-s	File with single reads"
	echo
	exit 1
}


# parse the command line
if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi
lowk="absent"
highk="absent"
forward="absent"
reverse="absent"
single="absent"


FILES=()

while [[ $# -gt 0 ]]		# while the number of args is greater than 0
do
key="$1"

case $key in
	-l)
	lowk="$2"
	shift
	shift
	;;
	-h)
	highk="$2"
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


if [ "$lowk" = "absent" ]; then
	if [ "$highk" = "absent" ]; then
		usage
	fi
fi


# report command line arguments
echo "Starting at $(date) with the following input parameters:"
echo

echo "low k = $lowk"
echo "high k = $highk"
echo "forward = $forward"
echo "reverse = $reverse"
echo "single = $single"


# assign working directory to scratch
work_dir="$SNIC_TMP"


# copy input files to the temporary directory
echo
echo "Copying files to working directory"
echo

cp /proj/uppstore2018195/software/singularity-containers/kat.img "$work_dir"


if [ "$forward" != "absent" ]; then
	cp "$forward" "$work_dir"
fi

if [ "$reverse" != "absent" ]; then
	cp "$reverse" "$work_dir"
fi

if [ "$single" != "absent" ]; then
	cp "$single" "$work_dir"
fi


# change to the working directory and start filtering
cd "$work_dir"


# create temporary files for reading into kat
temp="$(mktemp -u --suffix .fastq)"

input=""

if [ "$forward" != "absent" ]; then
	input+="$( basename $forward ) "
fi

if [ "$reverse" != "absent" ]; then
	input+="$( basename $reverse ) "
fi

if [ "$single" != "absent" ]; then
	input+="$( basename $single ) "
fi

mkfifo "$temp" && zcat $input > "$temp" &


# filter kmers
if [ "$lowk" = "absent" ]; then

	range=-"$highk"

	singularity exec kat.img kat filter kmer --high_count "$highk" -t "$num_cpus" -o filt"$range".kmer "$temp" && rm "$temp"

elif [ "$highk" = "absent" ]; then

	range="$lowk"-

	singularity exec kat.img kat filter kmer --low_count "$lowk" -t "$num_cpus" -o filt"$range".kmer "$temp" && rm "$temp"

else

	range="$lowk"-"$highk"

	singularity exec kat.img kat filter kmer --low_count "$lowk" --high_count "$highk" -t "$num_cpus" -o filt"$range".kmer "$temp" && rm "$temp"

fi


# filter the reads
if [ "$forward" != "absent" ]; then
	temp2="$(mktemp -u --suffix .fastq)"
	mkfifo "$temp2" && zcat "$( basename $forward )" > "$temp2" &
fi

if [ "$reverse" != "absent" ]; then
	temp3="$(mktemp -u --suffix .fastq)"
	mkfifo "$temp3" && zcat "$( basename $reverse )" > "$temp3" &
fi

if [ "$single" != "absent" ]; then
	temp1="$(mktemp -u --suffix .fastq)"
	mkfifo "$temp1" && zcat "$( basename $single )" > "$temp1" &
fi


if [ "$forward" != "absent" ]; then

	if [ "$reverse" != "absent" ]; then

		if [ "$single" != "absent" ]; then

			singularity exec kat.img kat filter seq -t "$num_cpus" -o kat_filt"$range"u -T 0.5 --seq "$temp2" --seq2 "$temp3" filt"$range".kmer-in.jf27 &
			singularity exec kat.img kat filter seq -t "$num_cpus" -o kat_filt"$range"m -T 0.5 --seq "$temp1" filt"$range".kmer-in.jf27
			wait

			rm "$temp2" "$temp3" && pigz kat_filt"$range"u.in.R1.fastq && pigz kat_filt"$range"u.in.R2.fastq
			rm "$temp1" && pigz kat_filt"$range"m.in.fastq

		else

			singularity exec kat.img kat filter seq -t "$num_cpus" -o kat_filt"$range" -T 0.5 --seq "$temp2" --seq2 "$temp3" filt"$range".kmer-in.jf27

			rm "$temp2" "$temp3" && pigz kat_filt"$range".in.R1.fastq && pigz kat_filt"$range".in.R2.fastq

		fi

	else	# interleaved; not sure if this will work with kat

		singularity exec kat.img kat filter seq -t "$num_cpus" -o kat_filt"$range"i -T 0.5 --seq "$temp2" filt"$range".kmer-in.jf27

		rm "$temp2" && pigz kat_filt"$range"i.in.fastq

	fi


elif [ "$single" != "absent" ]; then

	singularity exec kat.img kat filter seq -t "$num_cpus" -o kat_filt"$range" -T 0.5 --seq "$temp1" filt"$range".kmer-in.jf27

	rm "$temp1" && pigz kat_filt"$range".in.fastq

else

	usage

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

rm kat.img

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
