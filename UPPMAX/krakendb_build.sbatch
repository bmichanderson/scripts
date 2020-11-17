#!/bin/bash -l

###################
# Author: B. Anderson
# Date: 9 Aug 2019
# Modified: Oct 2020
# Description: build a Kraken2 database for screening for contaminant sequences; name of database is first command arg
###################

#SBATCH -A snic2018-8-336
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J build_Krakendb
#SBATCH -C mem256GB

echo "Starting to build a Kraken2 database at $(date)."

# load modules
echo "Loading modules"
module load bioinfo-tools
module load Kraken2

# set the number of threads available for Kraken2
export KRAKEN2_NUM_THREADS="$SLURM_NPROCS"

# set a name for the database
if [ -z "$1" ]; then
	exit 1
else
	kraken2db_name="$1"
fi

# change to the SNIC temporary working directory
cd "$SNIC_TMP"

# download the taxonomy needed for building the database
kraken2-build --threads "$SLURM_NPROCS" --download-taxonomy --db "$kraken2db_name"

# download and install various reference libraries
kraken2-build --threads "$SLURM_NPROCS" --download-library archaea --db "$kraken2db_name"
kraken2-build --threads "$SLURM_NPROCS" --download-library bacteria --db "$kraken2db_name"
kraken2-build --threads "$SLURM_NPROCS" --download-library viral --db "$kraken2db_name"
kraken2-build --threads "$SLURM_NPROCS" --download-library fungi --db "$kraken2db_name"
kraken2-build --threads "$SLURM_NPROCS" --download-library protozoa --db "$kraken2db_name"
kraken2-build --threads "$SLURM_NPROCS" --download-library plant --db "$kraken2db_name"

# build the database
kraken2-build --threads "$SLURM_NPROCS" --build --db "$kraken2db_name"

# reduce the size by deleting intermediate files
kraken2-build --threads "$SLURM_NPROCS" --clean --db "$kraken2db_name"

# create a folder on the storage drive and copy the database to it
kraken2db_dir=/proj/snic2020-16-240/kraken2db
if [ ! -d "$kraken2db_dir" ]; then
	mkdir -p "$kraken2db_dir"
fi

cp -av "$SNIC_TMP"/"$kraken2db_name" "$kraken2db_dir"

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

out_dir=/proj/snic2020-16-240/runs/"$SLURM_JOB_ID"
if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

cp /sw/share/slurm/"$system"/uppmax_jobstats/"$node"/"$SLURM_JOB_ID" "$out_dir"

echo "Finished creating Kraken2 database at $(date)"
