#######################
# Author: B.M. Anderson
# Date: 16 Mar 2020
# Modified: Nov 2023 (simplified and made more dependent on input text files), Sep 2024
# Description: plot input Newick trees with ape, potentially rooting them
#######################


# load libraries
suppressMessages(library(ape))
suppressMessages(library(phytools))


# a help function for when the script is called without arguments or incorrectly
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to plot a set of unrooted Newick trees, or rooted if outgroups are specified\n\n")
		cat("Usage: Rscript plot_trees.R <options> tree_file1 tree_file2...\n")
		cat("Options:\n")
		cat("\t-b\tBootstrap support level below which node labels are not displayed [optional]\n")
		cat("\t-c\tColours for tips in a text file, listing taxon and colour (tab separated, one per line) [optional]\n")
		cat("\t\tNote: needs a \"taxon\" entry in the samples file to know which to colour\n")
		cat("\t-l\tLabels to title the trees in a text file (one per line, same order as the files) [optional]\n")
		cat("\t-o\tOutgroup sampleIDs to use for rooting in a text file (one per line) [optional]\n")
		cat("\t-s\tSampleIDs, display names, and taxa in a text file (tab separated, one per line) [optional]\n")
		cat("\t-svg\tFlag for whether to output an SVG file for each tree [default: do not]\n")
		cat("\t-w\tFlag for whether to write the tips in the order plotted, from base of tree up\n\n")
	} else {
	cat(help_message)
	}
}


# parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) { # nolint
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	extra <- 1
	catch <- TRUE
	bootstrap <- 0
	colours_present <- FALSE
	colours_file <- ""
	labels_present <- FALSE
	labels_file <- ""
	outgroup_present <- FALSE
	outgroup_file <- ""
	samples_present <- FALSE
	samples_file <- ""
	svg_out <- FALSE
	write_tips <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-b") {
			bootstrap <- as.numeric(args[index + 1])
			catch <- FALSE
		} else if (args[index] == "-c") {
			colours_present <- TRUE
			colours_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-l") {
			labels_present <- TRUE
			labels_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-o") {
			outgroup_present <- TRUE
			outgroup_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-s")  {
			samples_present <- TRUE
			samples_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-svg")  {
			svg_out <- TRUE
		} else if (args[index] == "-w")  {
			write_tips <- TRUE
		} else {
			if (catch) {
				catch_args[extra] <- args[index]
				extra <- extra + 1
			} else {
				catch <- TRUE
			}
		}
	}
}
if (length(catch_args) < 1) {
	stop(help("Missing tree file(s)!\n"), call. = FALSE)
}

cat("Bootstrap threshold for node label printing is", bootstrap, "\n")


# read in and parse the files
if (colours_present) {
	colour_table <- read.table(colours_file, sep = "\t", header = FALSE)
}

if (labels_present) {
	labels <- read.table(labels_file, sep = "\t", header = FALSE)[, 1]
	cat("Labels for tree titles are:", unlist(labels), "\n")
}

if (outgroup_present) {
	outgroups <- read.table(outgroup_file, sep = "\t", header = FALSE)[, 1]
	cat("Outgroups are:", unlist(outgroups), "\n")
}

if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
	if (ncol(sample_table) == 2) {
		colnames(sample_table) <- c("ID", "label")
	} else if (ncol(sample_table) == 3) {
		colnames(sample_table) <- c("ID", "label", "taxon")
	} else {
		stop(help("Samples table formatted improperly!\n"), call. = FALSE)
	}
}

tree_list <- vector("list")
for (index in seq_len(length(catch_args))) {
	tree <- read.tree(catch_args[[index]])
	tree_list[[index]] <- tree
}


# check that trees have enough titles, if provided, or generate based on file name
if (labels_present) {
	if (length(labels) != length(tree_list)) {
		stop(help("The number of trees and labels do not match!\n"), call. = FALSE)
	}
} else {
	labels <- vector("list")
	for (index in seq_len(length(tree_list))) {
		labels[[index]] <- basename(catch_args[[index]])
	}
}


# root the trees if outgroups are present
if (outgroup_present) {
	for (index in seq_len(length(tree_list))) {
		if (sum(outgroups %in% tree_list[[index]]$tip.label) > 0) {
			these_outgroups <- outgroups[outgroups %in% tree_list[[index]]$tip.label]
			if (is.monophyletic(tree_list[[index]], as.character(these_outgroups))) {
				rootnode <- getMRCA(tree, as.character(these_outgroups))
				position <- 0.5 * tree$edge.length[which(tree$edge[, 2] == rootnode)]
				rooted_tree <- reroot(tree, rootnode, position, edgelabel = TRUE)
				tree_list[[index]] <- rooted_tree
			} else {
				cat("Tree", basename(catch_args[[index]]),
					"does not have monophyletic outgroups, so it is not rooted\n")
			}
		} else {
			cat("Tree", basename(catch_args[[index]]),
				"has no outgroups, so it is not rooted\n")
		}
	}
}


# Remove node labels if less than specified
# also determine the tree with the most taxa for dimensions
max_tips <- 0
for (index in seq_len(length(tree_list))) {
	for (index2 in seq_len(length(tree_list[[index]]$node.label))) {
		if (tree_list[[index]]$node.label[index2] == "") {
			tree_list[[index]]$node.label[index2] <- ""
		} else if (tree_list[[index]]$node.label[index2] == "Root") {
			tree_list[[index]]$node.label[index2] <- ""
		} else if (as.numeric(tree_list[[index]]$node.label[index2]) < bootstrap) {
			tree_list[[index]]$node.label[index2] <- ""
		}
	}
	if (length(tree_list[[index]]$tip.label) > max_tips) {
		max_tips <- length(tree_list[[index]]$tip.label)
	}
}


# Set sample colours, modifying them if both a colours file is specified AND there is a taxon column in the samples file
if (samples_present) {
	sample_table$colour <- rep("black", length(nrow(sample_table)))
	if (all(c(colours_present, ncol(sample_table) == 4))) {
		for (taxon in unique(sample_table$taxon)) {
			if (taxon %in% colour_table$V1) {
				colour <- colour_table$V2[match(taxon, colour_table$V1)]
				sample_table$colour[sample_table$taxon == taxon] <- colour
			}
		}
	}
}


# Substitute labels and colours if a samples file is present
tip_colours <- vector("list")
if (samples_present) {
	for (index in seq_len(length(tree_list))) {
		tips <- tree_list[[index]]$tip.label
		new_tips <- tips
		tip_col <- rep("black", length(tips))
		for (tind in seq_len(length(tips))) {
			if (tips[tind] %in% sample_table$ID) {
				new_tips[tind] <- sample_table$label[match(tips[tind], sample_table$ID)]
				tip_col[tind] <- sample_table$colour[match(tips[tind], sample_table$ID)]
			}
		}
		tree_list[[index]]$tip.label <- new_tips
		tip_colours[[index]] <- tip_col
	}
} else {
	for (index in seq_len(length(tree_list))) {
		tip_colours[[index]] <- rep("black", length(tree_list[[index]]$tip.label))
	}
}


# Plot trees
max_height <- max(c(max_tips / 5, 12))

pdf("trees.pdf", family = "ArialMT", width = (2 * max_height / 3), height = max_height)
cat("Plotting", length(tree_list), "trees to pdf\n")
for (index in seq_len(length(tree_list))) {
	plot.phylo(ladderize(tree_list[[index]], right = FALSE),
		no.margin = FALSE,
		font = 1,
		edge.width = 2,
		label.offset = max(nodeHeights(tree)) / 200,
		tip.col = tip_colours[[index]],
		main = labels[[index]])
	add.scale.bar(x = mean(par("usr")[1:2]),
		y = par("usr")[3] + 1,
		font = 1, lwd = 2)
	drawSupportOnEdges(tree_list[[index]]$node.label,
		adj = c(0.5, -0.5),
		frame = "none")

	if (write_tips) {
		lad_tree <- ladderize(tree_list[[index]], right = FALSE)
		# determine which edges are tips and get the order
		is_tip <- lad_tree$edge[, 2] <= length(lad_tree$tip.label)
		ordered_tips <- lad_tree$edge[is_tip, 2]
		# get the tips in order of plotting
		output_ordered_ids <- lad_tree$tip.label[ordered_tips]
		# write to a file
		connection <- file(paste0("tips_", index, ".txt"))
		writeLines(output_ordered_ids, connection)
		close(connection)
	}
}
invisible(dev.off())

if (svg_out) {
	for (index in seq_len(length(tree_list))) {
		svg(paste0("tree_", index, ".svg"), family = "ArialMT",
			width = (2 * max_height / 3), height = max_height)
		cat("Plotting Tree", basename(catch_args[[index]]), "to svg\n")
		plot.phylo(ladderize(tree_list[[index]], right = FALSE),
			no.margin = FALSE,
			font = 1,
			edge.width = 2,
			label.offset = max(nodeHeights(tree)) / 200,
			tip.col = tip_colours[[index]],
			main = labels[[index]])
		add.scale.bar(x = mean(par("usr")[1:2]),
			y = par("usr")[3] + 1,
			font = 1, lwd = 2)
		drawSupportOnEdges(tree_list[[index]]$node.label,
			adj = c(0.5, -0.5),
			frame = "none")
		invisible(dev.off())
	}
}
