#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 28 January 2020
# Modified: Oct 2020
# Description: following Brian Bushnell's recommended steps, error correct PE reads before running an initial Tadpole assembly
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 7-00:00:00
#SBATCH -J corr_assemble


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
	echo "input_type = $input_type, infile_1 = $1"
	echo

else

	echo "Starting at $(date) with the following input parameters:"
	echo "input_type = $input_type, infile_1 = $1, infile_2 = $2"
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


# transform to interleaved if paired
echo
echo "Transforming to interleaved if paired input"
echo

if [ "$input_type" = "paired" ]; then

	reformat.sh in1="$file1" in2="$file2" out=temp.fastq.gz

else

	ln -s "$file1" temp.fastq.gz

fi


# error correction
echo
echo "Correcting sequencing errors with bbmerge (for paired), clumpify and tadpole"
echo

if [ "$input_type" = "paired" ]; then

	bbmerge.sh in=temp.fastq.gz out=ecco.fastq.gz ecco mix vstrict ordered prefilter=2 prealloc=t minlength=50 -Xmx"$max_mem" && rm temp.fastq.gz && ln -s ecco.fastq.gz temp.fastq.gz

echo

fi

 clumpify.sh in=temp.fastq.gz out=eccc.fastq.gz ecc passes=4 -Xmx"$max_mem" && rm temp.fastq.gz && ln -s eccc.fastq.gz temp.fastq.gz

echo

 tadpole.sh in=temp.fastq.gz out=ecct.fastq.gz ecc k=62 ordered prefilter=2 prealloc=t -Xmx"$max_mem" errormult1=64 && rm temp.fastq.gz && ln -s ecct.fastq.gz temp.fastq.gz

# I added the errormult1=64 to the tadpole command since the default is fairly aggressive in error correction and might over-correct


if [ "$input_type" = "paired" ]; then

	# merging overlapping reads
	echo
	echo "Merging overlapping reads"
	echo

	bbmerge-auto.sh in=temp.fastq.gz out=temp_merged.fastq.gz outu=temp_unmerged.fastq.gz strict k=93 extend2=80 rem ordered prefilter=2 prealloc=t -Xmx"$max_mem"

fi


# limit the read length of output from merging and from correction, and rename corrected reads split into two files (for paired)
echo
echo "Limiting read length"
echo

if [ "$input_type" = "paired" ]; then

	bbduk.sh in=temp_merged.fastq.gz out=merged.fastq.gz minlength=50

	bbduk.sh in=temp_unmerged.fastq.gz out1=unmerged_1.fastq.gz out2=unmerged_2.fastq.gz minlength=50

	bbduk.sh in=ecct.fastq.gz out1=corrected_1.fastq.gz out2=corrected_2.fastq.gz minlength=50

else

	bbduk.sh in=ecct.fastq.gz out=corrected.fastq.gz minlength=50

fi


# assemble for multiple values of K
echo
echo "Assembling for different values of K to find the optimum"
echo

if [ "$input_type" = "paired" ]; then

	tadwrapper.sh in=merged.fastq.gz,unmerged_1.fastq.gz,unmerged_2.fastq.gz out=tadwrapper_contigs_%.fasta outfinal=tadwrapper_contigs k=40,81,102,124,217 bisect mincountretain=5 prefilter=2 prealloc=t -Xmx"$max_mem"

else

	tadwrapper.sh in=corrected.fastq.gz out=tadwrapper_contigs_%.fasta outfinal=tadwrapper_contigs k=40,81,102,124,217 bisect mincountretain=5 prefilter=2 prealloc=t -Xmx"$max_mem"

fi


# remove raw files and copy remainder to the storage directory
rm "$file1"

if [ "$input_type" = "paired" ]; then
	rm "$file2"
fi

rm ecc* temp*

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
