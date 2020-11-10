#!/bin/bash

##############
# Author: Ben Anderson
# Date: Nov 2020
# Modified:
# Description: assemble Illumina reads after initial_qc.sh and read depth assessment with a reference
##############

# BLAST filter python script
bscript="/home/banderson/scripts/blast_filter.py"


# Define a help function for an incorrect or empty call
help()
{
(	echo -e '\n''\t'"This script will assemble Illumina paired reads following QC, using a reference and KAT filtering"
	echo -e '\t'"It requires NCBI tools (makeblastdb, blastn) and singularity containers for KAT and Unicycler, as well as BBMap scripts"
	echo -e '\t'"It will create two rounds of Unicycler assembly, the second following a tadpole.sh extension"'\n'
	echo -e '\t'"Usage: assembly_ref.sh merged unmerged_1 unmerged_2 ref klower kupper"'\n'
) 1>&2
	exit 1
}


# Check if there are the correct number of command arguments
if [ $# -ne 6 ]; then
	help
fi


# set working directory
work_dir="$(pwd)"


# Create a log file
echo -e "Assembly log file"'\n' > "$work_dir"/assembly.log


# Set full paths for file arguments
merged=$(readlink -f "$1")
unmerged1=$(readlink -f "$2")
unmerged2=$(readlink -f "$3")
ref=$(readlink -f "$4")


# Report starting time and arguments
echo "Starting assembly at $(date) with the following input:" | tee -a "$work_dir"/assembly.log
echo "Working directory: $work_dir" | tee -a "$work_dir"/assembly.log
echo "merged = $merged" | tee -a "$work_dir"/assembly.log
echo "unmerged_1 = $unmerged1" | tee -a "$work_dir"/assembly.log
echo "unmerged_2 = $unmerged2" | tee -a "$work_dir"/assembly.log
echo "ref = $ref" | tee -a "$work_dir"/assembly.log
echo "klower = $5" | tee -a "$work_dir"/assembly.log
echo "kupper = $6" | tee -a "$work_dir"/assembly.log

start=$(date +%s)


# Step 1: KAT filtering
echo -e '\n'"****"'\n'"Filtering reads with KAT"'\n'"****"'\n' | tee -a "$work_dir"/assembly.log

mkdir kat && cd kat
temp="$(mktemp -u --suffix .fastq)"
mkfifo "$temp" && zcat "$merged" "$unmerged1" "$unmerged2" > "$temp" &
singularity exec -B $(pwd) ~/singularity-containers/kat.img kat filter kmer --low_count "$5" --high_count "$6" -t 12 -o filt"$5"-"$6".kmer "$temp" \
|& tee -a "$work_dir"/assembly.log && rm "$temp"

temp1="$(mktemp -u --suffix .fastq)"
mkfifo "$temp1" && zcat "$merged" > "$temp1" &
singularity exec -B $(pwd) ~/singularity-containers/kat.img kat filter seq -t 12 -o kat_filt"$5"-"$6"m -T 0.5 --seq "$temp1" filt"$5"-"$6".kmer-in.jf27 \
|& tee -a "$work_dir"/assembly.log && rm "$temp1" && pigz kat_filt"$5"-"$6"m.in.fastq

temp2="$(mktemp -u --suffix .fastq)"
temp3="$(mktemp -u --suffix .fastq)"
mkfifo "$temp2" && zcat "$unmerged1" > "$temp2" &
mkfifo "$temp3" && zcat "$unmerged2" > "$temp3" &
singularity exec -B $(pwd) ~/singularity-containers/kat.img kat filter seq -t 12 -o kat_filt"$5"-"$6"u -T 0.5 --seq "$temp2" --seq2 "$temp3" \
filt"$5"-"$6".kmer-in.jf27 |& tee -a "$work_dir"/assembly.log && rm "$temp2" "$temp3" && pigz kat_filt"$5"-"$6"u.in.R1.fastq && pigz kat_filt"$5"-"$6"u.in.R2.fastq


# Step 2: First round of Unicycler
echo -e '\n'"****"'\n'"Running first round of Unicycler"'\n'"****"'\n' | tee -a "$work_dir"/assembly.log
echo "See analysis details in unicycler.log" | tee -a "$work_dir"/assembly.log

mkdir unicycler
singularity exec -B $(pwd) ~/singularity-containers/unicycler.img unicycler -1 kat_filt"$5"-"$6"u.in.R1.fastq.gz -2 kat_filt"$5"-"$6"u.in.R2.fastq.gz \
-s kat_filt"$5"-"$6"m.in.fastq.gz -o unicycler/r01 --verbosity 2 --keep 0 --threads 12 --no_correct --no_rotate --no_pilon \
--min_component_size 250 --min_dead_end_size 250


# Step 3: Filter resulting contigs against the reference then extend with tadpole and map
echo -e '\n'"****"'\n'"Filtering assembled contigs against the reference"'\n'"****"'\n' | tee -a "$work_dir"/assembly.log

cd unicycler/r01
"$bscript" -r "$ref" assembly.fasta |& tee -a "$work_dir"/assembly.log


# Step 4: Extending with tadpole and mapping to extended contigs
echo -e '\n'"****"'\n'"Extending with tadpole and mapping to extended contigs"'\n'"****"'\n' | tee -a "$work_dir"/assembly.log

tadpole.sh in=hits.fasta out=temp_ext.fasta el=100000 er=100000 mode=extend extra="$merged","$unmerged1","$unmerged2" k=31 |& tee -a "$work_dir"/assembly.log

cd ..
echo -e '\n'"Mapping merged reads"'\n' | tee -a "$work_dir"/assembly.log
bbmap.sh ref=r01/temp_ext.fasta in="$merged" outm=mapm.fastq.gz nodisk |& tee -a "$work_dir"/assembly.log && \

echo -e '\n'"Mapping unmerged reads"'\n' | tee -a "$work_dir"/assembly.log
bbmap.sh ref=r01/temp_ext.fasta in1="$unmerged1" in2="$unmerged2" outm=mapu#.fastq.gz nodisk |& tee -a "$work_dir"/assembly.log


# Step 5: Second round of Unicycler
echo -e '\n'"****"'\n'"Running second round of Unicycler"'\n'"****"'\n' | tee -a "$work_dir"/assembly.log
echo "See analysis details in unicycler.log" | tee -a "$work_dir"/assembly.log

singularity exec -B $(pwd) ~/singularity-containers/unicycler.img unicycler -1 mapu1.fastq.gz -2 mapu2.fastq.gz \
-s mapm.fastq.gz -o r02 --verbosity 2 --keep 0 --threads 12 --no_correct --no_rotate --no_pilon \
--min_component_size 250 --min_dead_end_size 250


# Remove temporary files and return to working directory
rm map*
cd ..
rm kat* filt*
cd "$work_dir"


# Report time taken for the assemblies
end=$(date +%s)
duration=$(( $end - $start ))
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e '\n'"Finished at $(date), after running for $duration seconds, or $duration_mins minutes, or $duration_hours hours" | tee -a "$work_dir"/assembly.log
