#!/bin/bash
##########################
# Author: B. Anderson
# Date: 11 June 2020
# Description: using a gene list and database, search for hits in specified fasta files
##########################

parse_script=~/scripts/blast_parse.py


# a function for when the script is called incorrectly or without arguments
usage ()
{
	echo "This script will BLAST search databases of specified regions for specified fasta files."
	echo
	echo "Usage: $(basename $0) options(- ... - ...) file1 file2 ..."
	echo
	echo "Options:"
	echo
	echo "	-f|--file	File (required) with list of regions, one per line"
	echo
	echo "	-d|--dir	Directory (required) containing folders of region databases, each folder named by region"
	echo
	echo "	-t|--type	Type of database to search, nucl (default) or prot"
	echo
	echo "	-r|--rna	Run tRNAScan-SE on the files, yes or no (default)"
	echo
	echo "	-m|--max	Maximum number of target sequences returned in the BLAST searches (default = 15)"
	echo
	echo "	-p|--id		Minimum percentage identity filtered during blast result parsing (default = 0)"
	echo
	echo "	-l|--len	Minimum length of a blast result kept (default = 0)"
	echo

	exit 1
}

if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi

FILES=()
regions_present="FALSE"
dir_present="FALSE"
type_specified="FALSE"
rna_specified="FALSE"

max_targs=15
pid=0
len=0

while [[ $# -gt 0 ]]		# while the number of args is greater than 0
do
key="$1"

case $key in
	-f|--file)
	regions_list="$2"
	regions_present="TRUE"
	shift
	shift
	;;
	-d|--dir)
	directory="$2"
	dir_present="TRUE"
	shift
	shift
	;;
	-t|--type)
	dbtype="$2"
	type_specified="TRUE"
	shift
	shift
	;;
	-r|--rna)
	run_trna="$2"
	rna_specified="TRUE"
	shift
	shift
	;;
	-m|--max)
	max_targs="$2"
	shift
	shift
	;;
	-p|--id)
	pid="$2"
	shift
	shift
	;;
	-l|--len)
	len="$2"
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

if [ $# -eq 0 ]; then		# this would be true if there were no files specified
	usage
fi

if [ "$regions_present" = "FALSE" ]; then
	usage
fi

if [ "$dir_present" = "FALSE" ]; then
	usage
fi

if [ "$type_specified" = "FALSE" ]; then
	dbtype="nucl"
elif [ "$dbtype" != "prot" ] && [ "$dbtype" != "nucl" ]; then
	usage
fi

if [ "$rna_specified" = "FALSE" ]; then
	run_trna="no"
elif [ "$run_trna" != "yes" ] && [ "$run_trna" != "no" ]; then
	usage
fi


# determine available processors for running BLAST searches
threads="$(grep -c ^processor /proc/cpuinfo)"


# assign working directory
work_dir="$(pwd)"


# get full path for the regions input file
regions_path="$(readlink -f $regions_list)"


# run searches for each fasta file
for fasta_file in "$@"; do

	# get full path and basename
	fasta_path="$(readlink -f $fasta_file)"
	fasta_base="$(basename $fasta_path)"


	if [ ! -d "$work_dir/${fasta_base/.f*/}" ]; then
		mkdir -p "$work_dir/${fasta_base/.f*/}"
	fi

	cd "${fasta_base/.f*/}"


	# Run BLAST searches and parse results

	echo "**** Running BLAST searches for $fasta_base"
	echo

	for region in $(cat "$regions_path") ; do

		if [ "$dbtype" = "nucl" ]; then

			echo "#### Blasting for $region"

			blastn -task blastn -query "$fasta_path" -db "$directory"/"$region"/"$region"_nucldb -outfmt "6 std qseq sseq qlen slen stitle" -max_target_seqs "$max_targs" -num_threads "$threads" > "$region"_blast.tab

			"$parse_script" -q "$fasta_path" -p "$pid" -l "$len" "$region"_blast.tab

			if [ -s hits.txt ]; then		# if there is something in the file
				mv hits.txt "$region"_"$dbtype"_hits.txt
			elif [ -f hits.txt ]; then
				rm hits.txt
			fi
			rm "$region"_blast.tab

		elif [ "$dbtype" = "prot" ]; then

			echo "#### Blasting for $region"

			blastx -query "$fasta_path" -db "$directory"/"$region"/"$region"_protdb -outfmt "6 std qseq sseq qlen slen stitle" -max_target_seqs "$max_targs" -num_threads "$threads" > "$region"_blast.tab

			"$parse_script" -q "$fasta_path" -p "$pid" -l "$len" -t x "$region"_blast.tab

			if [ -s hits.txt ]; then		# if there is something in the file
				mv hits.txt "$region"_"$dbtype"_hits.txt
			elif [ -f hits.txt ]; then
				rm hits.txt
			fi
			rm "$region"_blast.tab

		fi

	done


	# Run tRNAScan-SE if requested

	if [ "$run_trna" = "yes" ]; then

		echo "Running tRNAscan-SE on the target fasta file"

		tRNAscan-SE -q -Q -O -o tRNAscan_results.txt "$fasta_path"

	fi

	cd "$work_dir"

done
