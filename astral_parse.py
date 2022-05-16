#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: May 2022
# Description: parse the output of ASTRAL to make trees with desired node labels
# Note: Use the `--branch-annotate 2` option when running ASTRAL to get all labels
##########################


import sys
import argparse
from Bio import Phylo


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to parse ASTRAL trees with multiple node labels')


# add arguments to parse
parser.add_argument('-t', type = str, dest = 'tree_file', help = 'The ASTRAL tree file')
parser.add_argument('-f', type = str, dest = 'out_format', help = 'The desired format: \"q\" (quartet supports),' + 
	' \"f\" (frequencies), \"p\" (posterior probs; default)')
parser.add_argument('-o', type = str, dest = 'out_pre', help = 'The output prefix [default \"output\"]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

tree_file = args.tree_file
out_format = args.out_format
out_pre = args.out_pre

if not tree_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_format:
	out_format = 'p'

if out_format not in ['p', 'q', 'f']:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# Load the tree
tree = Phylo.read(tree_file, 'newick')


# process the treefile to relabel internal nodes to whichever desired
# a typical label looks like:
# [q1=0.47430465014957124;q2=0.20129744142138076;q3=0.3243979084290481;f1=136.12543459292695;f2=57.77236568793628;f3=93.1021997191368;
# pp1=0.9999976302178388;pp2=3.9950444203194846E-7;pp3=1.9702777190850334E-6;QC=142;EN=287.0]
# q = quartet support
# f = quartet frequency
# pp = posterior probability (for all three alternatives)
# QC = total number of quartets around the branch
# EN = effective number of genes


for node in tree.get_nonterminals():
	new_name = []
	if node.name:
		elements = str(node.name).strip('[').split(';')
		for element in elements:
			el_parts = element.split('=')
			if el_parts[0][0] == out_format:
				new_name.append(el_parts[1])
		if len(new_name) > 0:
			node.name = '/'.join(new_name)

# write the tree
Phylo.write(tree, open(out_pre + '_' + out_format + '.tre', 'w'), 'newick')
