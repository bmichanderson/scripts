#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 18/25 Oct 2019
# Modified: Oct 2020
# Description: run Kraken2 classification on an input set of reads (args 2 and 3 for paired) using a pre-built Kraken2 database (arg 1)
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J Kraken2_run


# set default input type
input_type="paired"


# define a function for when the script is called incorrectly
usage()
{
	echo "Usage: script.sbatch kraken2_database illumina_1 [illumina_2 if paired]" 1>&2
	exit 1
}

if [ -z "$1" ]; then
	usage
elif [ -z "$2" ]; then
	usage
elif [ -z "$3" ]; then
	input_type="single"
fi


# report command line arguments
echo "Starting a Kraken2 analysis at $(date)."

kraken2db="$1"

input1="$( basename $2 )"

if [ "$input_type" = "single" ]; then

	echo "k2db = $kraken2db input_type = $input_type infile_1 = $input1"

else

	input2="$( basename $3 )"

	echo "k2db = $kraken2db input_type = $input_type infile_1 = $input1, infile_2 = $input2"

fi


# load modules
echo "Loading modules"

module load bioinfo-tools
module load Kraken2 Krona/2.7


# set the number of threads available for Kraken2
export KRAKEN2_NUM_THREADS="$SLURM_NPROCS"


# make a directory in SNIC temporary storage to hold the database
k2db_dir="$SNIC_TMP"/k2db

if [ ! -d "$k2db_dir" ]; then
	mkdir -p "$k2db_dir"
fi


# copy the input reads and Kraken2 database to the SNIC temporary directory
cp "$2" "$SNIC_TMP" &

if [ "$input_type" != "single" ]; then

	cp "$3" "$SNIC_TMP" &

fi

cp -av "$kraken2db"/* "$k2db_dir"
wait


# change to the SNIC temporary working directory and run the Kraken2 analysis
cd "$SNIC_TMP"

if [ "$input_type" = "single" ]; then

	kraken2 --db "$k2db_dir" --threads "$SLURM_NPROCS" --report kraken2.rpt --fastq-input --gzip-compressed "$input1" > output_kraken2.tsv

else

	kraken2 --db "$k2db_dir" --threads "$SLURM_NPROCS" --report kraken2.rpt --fastq-input --gzip-compressed --paired "$input1" "$input2" > output_kraken2.tsv

fi


# run Krona for producing an html file for displaying the results
ktImportTaxonomy <( cut -f2,3 "output_kraken2.tsv" ) -o output_kraken2_krona.html


# copy the results to storage
out_dir=/proj/snic2020-16-240/runs/"$SLURM_JOB_ID"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

cp kraken2.rpt "$out_dir" &
cp output_kraken2* "$out_dir"
wait


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

echo "Finished running a Kraken2 analysis at $(date)"
