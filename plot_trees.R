#######################
# Author: B. Anderson
# Date: 16 Mar 2020
# Description: plot an input set of unrooted trees from RAxML, potentially rooting each on a specified outgroup taxon set
#######################

library(ape)


# a help function for when the script is called without arguments or incorrectly

help <- function(help_message) {
	cat("A script to plot a set of unrooted trees, or rooted if outgroups are specified.\n\n")
	cat("Usage: Rscript plot_trees.R options(-... -...) tree_file1 tree_file2...\n\n")
	cat("Options:\n")
	cat("	-l	A file with labels of regions, one per line, for naming the input tree set (so same order as the trees).\n")
	cat("	-o	A file with a list of outgroup taxa to use for rooting, with the preferred taxa first, one per line.\n")
	cat("		These taxa names will be searched against tip labels and should ideally be found uniquely once.\n")
	cat("	-t	A file with lists of taxa to colour (max 10 unless using your own colour pallete).\n")
	cat("		The lists should be entered one per line, with each line receiving a different colour.\n")
	cat("		Names in a list should be delimited by commas (',') e.g. Cau,Cca,Cgr\n")
	cat("	-d	A table to be taken as a dataframe for equating taxa (if not tip labels) and tip labels.\n")
	cat("		The table should be tab delimited, with the taxon from '-t' as the first column, and tip label as the second.\n")
	cat("	-c	A file with a colour pallete, listing colours to use (one per line). Optional (defaults to a 10 colour pallete).\n")
	cat("	-b	Bootstrap support level below which labels are not shown (default: 75).\n\n")
	cat(help_message)
}


# parse the command line

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help("<This is where the error message will be printed>\n"), call.=FALSE)
} else {
	catch_args <- vector("list")
	i <- 1

	labels_present <- FALSE
	labels_file <- ""
	outgroup_present <- FALSE
	outgroup_file <- ""
	taxa_present <- FALSE
	taxa_file <- ""
	df_present <- FALSE
	df_file <- ""
	colours_present <- FALSE
	colours_file <- ""
	bootstrap <- 75

	for (index in 1:length(args)) {
		if (args[index] == "-l") {
			labels_present <- TRUE
			labels_file <- args[index + 1]
		} else if (args[index] == "-o") {
			outgroup_present <- TRUE
			outgroup_file <- args[index + 1]
		} else if (args[index] == "-t")  {
			taxa_present <- TRUE
			taxa_file <- args[index + 1]
		} else if (args[index] == "-d")  {
			df_present <- TRUE
			df_file <- args[index + 1]
		} else if (args[index] == "-c")  {
			colours_present <- TRUE
			colours_file <- args[index + 1]
		} else if (args[index] == "-b")  {
			bootstrap <- as.numeric(args[index + 1])
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

tree_file_list <- setdiff(catch_args, c(labels_file, outgroup_file, taxa_file, df_file, colours_file, bootstrap))	# exclude the options that were added to catch_args


# read in and parse the optional files, if present

if (labels_present) {
	labels <- vector("list")
	i <- 1

	con <- file(labels_file, "r")
	while (TRUE) {
		line <- readLines(con, n = 1)
		if (length(line) == 0) { break }
		labels[i] <- line
		i <- i + 1
	}
	close(con)

	cat("Labels are:", unlist(labels), "\n")
}

if (outgroup_present) {
	outgroups <- vector("list")
	i <- 1

	con <- file(outgroup_file, "r")
	while (TRUE) {
		line <- readLines(con, n = 1)
		if (length(line) == 0) { break }
		outgroups[i] <- line
		i <- i + 1
	}
	close(con)

	cat("Outgroups are:", unlist(outgroups), "\n")
}

if (taxa_present) {
	tax_colour <- vector("list")
	i <- 1

	con <- file(taxa_file, "r")
	while (TRUE) {
		line <- readLines(con, n = 1)
		if (length(line) == 0) { break }
		tax_colour[i] <- strsplit(line, ",")
		i <- i + 1
	}
	close(con)

	cat("Taxa to colour are:", unlist(tax_colour), "\n")
}


if (df_present) {
	df_table <- read.table(df_file, sep="\t", stringsAsFactors=FALSE)
}


if (colours_present) {
	colours <- vector("list")
	i <- 1

	con <- file(colours_file, "r")
	while (TRUE) {
		line <- readLines(con, n = 1)
		if (length(line) == 0) { break }
		colours[i] <- line
		i <- i + 1
	}
	close(con)

	cat("Colours to use are:", unlist(colours), "\n")
}


cat("Bootstrap threshold for node label printing is", bootstrap, "\n")


# establish a colour pallete

if (colours_present) {
	col_pal <- colours
} else {
	col_pal <- c('red', 'dodgerblue', 'gold', 'darkorange', 'burlywood', 'forestgreen', 'turquoise', 'saddlebrown', 'blueviolet', 'salmon')
}


# read in the trees

tree_list <- vector("list")
i <- 1

for (tree_file in tree_file_list) { 
	tree_list[[i]] <- read.tree(tree_file)
	i <- i + 1
}

if (labels_present) {
	i <- 1
	if (length(labels) == length(tree_list)) {
		for (tree in tree_list) {
			tree$tree.names[1] <- labels[[i]]
			tree_list[[i]] <- tree
			i <- i + 1
		}
	} else {
		stop(help("The number of trees and labels do not match!\n"), call.=FALSE)
	}
}


# root the trees if outgroups are present

rooted_tree_list <- vector("list")
i <- 1

unrooted_tree_list <- vector("list")
j <- 1

if (outgroup_present) {
	for (tree in tree_list) {
		tree_taxa <- tree$tip.label
		unrooted <- TRUE

		for (outgroup in outgroups) {
			if (length(grep(outgroup, tree_taxa)) > 0) {
				rooted_tree_list[[i]] <- root(tree, tree_taxa[grep(outgroup, tree_taxa)], resolve.root = TRUE, edgelabel = TRUE)
				i <- i + 1
				unrooted <- FALSE
				break
			}
		}
		
		if (unrooted) {
			unrooted_tree_list[[j]] <- tree
			j <- j + 1
		}
	}

	cat("Trees rooted:", length(rooted_tree_list), ", Trees unrooted:", length(unrooted_tree_list), "\n")

} else {
	
	unrooted_tree_list <- tree_list

	cat("Trees all unrooted:", length(unrooted_tree_list), "\n")
	
}


# Remove node labels (bootstrap support) if less than specified bootstrap percentage (default 75)

i <- 1
for (tree in rooted_tree_list) {
	for (index in 1: length(tree$node.label)) {
		if (tree$node.label[index] == "") {
			tree$node.label[index] <- ""
		} else if (tree$node.label[index] == "Root") {
			tree$node.label[index] <- ""
		} else if (as.numeric(tree$node.label[index]) < bootstrap) {
			tree$node.label[index] <- ""
		}
	}

	rooted_tree_list[[i]] <- tree
	i <- i + 1
}

j <- 1
for (tree in unrooted_tree_list) {
	for (index in 1: length(tree$node.label)) {
		if (tree$node.label[index] == "") {
			tree$node.label[index] <- ""
		} else if (tree$node.label[index] == "Root") {
			tree$node.label[index] <- ""
		} else if (as.numeric(tree$node.label[index]) < bootstrap) {
			tree$node.label[index] <- ""
		}
	}

	unrooted_tree_list[[j]] <- tree
	j <- j + 1
}


# create a pdf and plot

pdf("trees.pdf", family="ArialMT", pointsize=12, paper="A4", width=7.5, height=11)

num_trees <- length(unrooted_tree_list) + length(rooted_tree_list)
cat("Plotting", num_trees, "trees to pdf.\n")

#if (num_trees > 8) {
#	layout(matrix(1:8, 4, 2))
#} else {
#	layout(matrix(1:4, 2, 2))
#}

par(cex=0.75)

if (taxa_present) {
	for (tree in unrooted_tree_list) {
		tree_taxa <- tree$tip.label

		tip_col <- rep("black", length(tree_taxa))
		i <- 1

		for (tax_list in tax_colour) {
			for (taxon in tax_list) {
				if (df_present) {
					taxa = df_table[grep(taxon, df_table$V1), 2]
					for (tax in taxa) {
						tip_col[grep(tax, tree_taxa)] <- col_pal[i]
					}
				} else {
					tip_col[grep(taxon, tree_taxa)] <- col_pal[i]
				}
			}
			i <- i + 1
		}

		plot.phylo(tree, "unrooted", lab4ut = "axial", tip.col = tip_col, main = tree$tree.names, show.node.label = TRUE)
		add.scale.bar()

	}


	if (outgroup_present) {
		for (tree in rooted_tree_list) {
			tree_taxa <- tree$tip.label

			tip_col <- rep("black", length(tree_taxa))
			i <- 1

			for (tax_list in tax_colour) {
				for (taxon in tax_list) {
					if (df_present) {
						taxa = df_table[grep(taxon, df_table$V1), 2]
						for (tax in taxa) {
							tip_col[grep(tax, tree_taxa)] <- col_pal[i]
						}
					} else {
						tip_col[grep(taxon, tree_taxa)] <- col_pal[i]
					}
				}
				i <- i + 1
			}

			plot.phylo(tree, tip.col = tip_col, main = tree$tree.names, show.node.label = TRUE)
			add.scale.bar()	

		}

	}

} else {

	for (tree in unrooted_tree_list) {
		plot.phylo(tree, "unrooted", lab4ut = "axial", main = tree$tree.names, show.node.label = TRUE)
		add.scale.bar()

	}


	if (outgroup_present) {
		for (tree in rooted_tree_list) {
			plot.phylo(tree, main = tree$tree.names, show.node.label = TRUE)
			add.scale.bar()	

		}

	}


}

dev.off()

