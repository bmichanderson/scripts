#!/bin/bash

##########################
# Author: B. Anderson
# Date: 10 Mar 2020
# Modified: Oct 2020
# Description: using a region list and dictionary (optional), build nucleotide or protein BLAST databases from specified genbank files
##########################


# python scripts to call
parse_script=~/scripts/genbank_parse.py
comb_script=~/scripts/extract_sort_comb.py


# a function for when the script is called incorrectly or without arguments
usage ()
{
	echo "This script will compile sequences and build BLAST databases from a set of regions in given genbank files."
	echo
	echo "Usage: $(basename $0) options(- ... - ...) gbfile1 gbfile2 ..."
	echo
	echo "Options:"
	echo
	echo "	-f|--file	File (required) with list of regions, one per line, e.g. atp1"
	echo
	echo "	-d|--dict	Dictionary (optional) of alternative region names, tab-delimited, one per line, e.g. atp1-1	atp1"
	echo
	echo "	-t|--type	Type of database to build, nucl (default) or prot"
	echo
	exit 1
}

if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi

FILES=()
dict_present="FALSE"
type_specified="FALSE"

while [[ $# -gt 0 ]]		# while the number of args is greater than 0
do
key="$1"

case $key in
	-f|--file)
	regions_list="$2"
	shift
	shift
	;;
	-d|--dict)
	dict="$2"
	dict_present="TRUE"
	shift
	shift
	;;
	-t|--type)
	dbtype="$2"
	type_specified="TRUE"
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

if [ "$type_specified" = "FALSE" ]; then
	dbtype="nucl"
elif [ "$dbtype" != "prot" ] && [ "$dbtype" != "nucl" ]; then
	usage
fi


# assign working directory
work_dir="$(pwd)"


# create output top directory
out_dir="$work_dir"/databases_"$(date +%b%Y)"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi


# run a series of steps to extract all regions from each of the genbank files
for gbfile in "$@"; do

	echo "$gbfile"

	"$parse_script" -l CDS "$gbfile" > CDS.tmp

	if [ "$dbtype" = "nucl" ]; then

		"$parse_script" -l rRNA "$gbfile" > rRNA.tmp
		"$parse_script" -l tRNA "$gbfile" > tRNA.tmp

	fi

	cat *.tmp > regions.list

	"$parse_script" -f regions.list -t "$dbtype" "$gbfile"

	rm *.tmp regions.list

done


# run a second script to deal with all the extracted fasta files
echo
echo "Combining extracted regions then deleting intermediate files"

if [ "$dict_present" = "TRUE" ]; then
	"$comb_script" -f "$regions_list" -d "$dict" *extract.fasta
else
	"$comb_script" -f "$regions_list" *extract.fasta
fi

rm *extract.fasta


# for each region in the input list, make a directory and a database in that directory based on the corresponding file
echo
echo "Creating folders and BLAST databases"

for region in $(cat "$regions_list"); do

	if [ ! -d "$out_dir"/"$region" ] && [ -f "$region".fasta ]; then
		mkdir -p "$out_dir"/"$region"
	fi

	if [ -f "$region".fasta ]; then

		mv "$region".fasta "$out_dir"/"$region"/"$region"_"$dbtype".fasta

		cd "$out_dir"/"$region"

		makeblastdb -in "$region"_"$dbtype".fasta -dbtype "$dbtype" -out "$region"_"$dbtype"db -logfile temp.log

		rm temp*

		cd "$work_dir"

	fi

done
