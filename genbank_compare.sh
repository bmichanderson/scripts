#!/bin/bash

##########################
# Author: B. Anderson
# Date: 4 Mar 2020
# Updated: Oct 2020
# Description: output a tab-delimited table summarizing what query regions are present in a set of genbank files
##########################


# python script to call
parse_script=~/scripts/genbank_parse.py


# a function for when the script is called incorrectly or without arguments
usage ()
{
	echo "This script will assess region presence for a set of genbank files."
	echo
	echo "Usage: $(basename $0) options(- ... - ...) file1 file2 ..."
	echo
	echo "Options:"
	echo
	echo "	-f|--file	File (required) with list of regions, one per line"
	echo
	echo "	-d|--dict	Dictionary (optional) of entries and standard equivalencies in region names, tab-delimited, one per line"
	echo
	exit 1
}

if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi


# a function to check whether a string is present in an array (see https://stackoverflow.com/questions/3685970/check-if-a-bash-array-contains-a-value)
string_in_array () {
	local i match="$1"
	shift
	for i; do
		if [ "$i" = "$match" ]; then
			return 0
		fi
	done
	return 1
}


# parse arguments (see https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash)

FILES=()
dict_present="FALSE"

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


# create a header for the summary file
out_head="Gene"

for region in $(cat "$regions_list"); do

	out_head="$out_head	$region"

done

echo "$out_head" > "summary.tab"


# output a line for each genbank file indicating whether each of the target regions is present or not
for gbfile in "$@"; do

	"$parse_script" -l CDS "$gbfile" > CDS.tmp
	"$parse_script" -l rRNA "$gbfile" > rRNA.tmp
	"$parse_script" -l tRNA "$gbfile" > tRNA.tmp

	filename="${gbfile%_*.gb}"
	shortened="${filename/_NC/}"

	output="$shortened"

	for region in $(cat "$regions_list"); do

		if string_in_array "$region" $(cat *.tmp); then		# if the general name for the region is somewhere in the listed regions present
			output="$output	yes"

		elif [ "$dict_present" = "TRUE" ]; then
			hit="FALSE"

			for alt in $(grep "$region" "$dict" | cut -f 1); do
				if string_in_array "$alt" $(cat *.tmp); then
					output="$output	yes"
					hit="TRUE"
					break
				fi
			done

			if [ "$hit" = "FALSE" ]; then
				output="$output	no"
			fi

		else
			output="$output	no"

		fi

	done

	echo "$output" >> "summary.tab"

	rm *.tmp

done
