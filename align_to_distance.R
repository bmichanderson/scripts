
##########
# Author: Ben Anderson
# Date: Nov 2022
# Modified: April 2023, May 2023, July 2023
# Description: convert a DNA alignment (fasta) or multiple alignments into a distance matrix (Nexus format)
# Note: the resulting distance matrix will be saved as "dist_out.nex"
##########


# load required library
suppressMessages(library(ape))


# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to convert a DNA alignment (fasta) or multiple alignments into a distance matrix (Nexus)\n")
		cat("Usage: Rscript align_to_distance.R [-m model] [-p pofadinr method] ",
			"[-s samples] alignment.fasta [alignment2.fasta ...]\n")
		cat("Options:\n")
		cat("\t-m\tThe ape DNA distance model [default F84]\n")
		cat("\t-a\tResolve ambiguous bases randomly before calculating distance in ape (\"y\" or \"n\" [default])\n")
		cat("\t-p\tThe pofadinr nucleotide distance method (\"g\" for GENPOFAD, \"m\" for MATCHSTATES)\n")
		cat("\t\tIf set, this will be the method used rather than ape dist.dna\n")
		cat("\t-s\tA file with tab-separated sample IDs and corresponding tip labels for conversion, ",
			"one per line [optional]\n")
	} else {
		cat(help_message)
	}
}


# parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	extra <- 1
	catch <- TRUE
	model <- "F84"
	ambig <- "n"
	samples_present <- FALSE
	pofad_method_set <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-m") {
			model <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-a") {
			ambig <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-p") {
			pofad_method_set <- TRUE
			pofad_method <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
			catch <- FALSE
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
	stop(help("Missing argument for alignment!\n"), call. = FALSE)
}


# read in the input
if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
}

alignments <- vector("list")
for (index in seq_len(length(catch_args))) {
	alignment <- read.FASTA(catch_args[[index]], type = "DNA")
	alignments[[index]] <- alignment
}


# capture taxa present; if samples file present, change the labels
taxa <- c()
if (samples_present) {
	renamed_alignments <- vector("list")
	big_index <- 1
	for (alignment in alignments) {
		renamed_alignment <- alignment
		for (index in seq_len(length(names(alignment)))) {
			name <- names(alignment[index])
			if (name %in% sample_table$V1) {
				names(renamed_alignment)[index] <- sample_table$V2[match(name, sample_table$V1)]
			}
		}
		renamed_alignments[[big_index]] <- renamed_alignment
		big_index <- big_index + 1
		taxa <- unique(c(taxa, names(renamed_alignment)))
	}
	alignments <- renamed_alignments
} else {
	for (alignment in alignments) {
		taxa <- unique(c(taxa, names(alignment)))
	}
}
cat(paste0("\nRead in ", length(alignments), " alignments\n\n"))


# Set the method
if (pofad_method_set) {
	suppressMessages(library(pofadinr))
	if (pofad_method == "m") {
		cat(paste0("Calculating distances using pofadinr dist.snp and the MATCHSTATES method\n"))
		model <- "MATCHSTATES"
	} else if (pofad_method == "g") {
		cat(paste0("Calculating distances using pofadinr dist.snp and the GENPOFAD method\n"))
		model <- "GENPOFAD"
	} else {
		stop(help("Please specify pofadinr method as \"m\" or \"g\"\n"), call. = FALSE)
	}
} else {
	if (ambig == "y") {
		cat("Replacing ambiguous bases in the alignment randomly and relative to frequency in the column\n")
	}
	cat(paste0("Calculating distances using ape dist.dna and the ", model, " model\n"))
}


# for each alignment, calculate a distance matrix using the preferred method
dist_list <- vector("list")
index <- 1
cat("\nProcessing alignment:")
for (alignment in alignments) {
	cat(paste0(" ", index))
	# calculate the distances
	if (pofad_method_set) {
		# convert missing data to "?"
		temp <- as.character(alignment)
		for (subindex in seq_len(length(temp))) {
			temp[[subindex]][temp[[subindex]] == "n"] <- "?"
		}
		alignment <- as.DNAbin(temp)
		distances <- dist.snp(alignment, model = model)
	} else {
		if (ambig == "y") {
			alignment <- solveAmbiguousBases(alignment, method = "columnwise", random = TRUE)
		}
		distances <- dist.dna(alignment, model = model, pairwise.deletion = TRUE)
	}
	dist_list[[index]] <- as.matrix(distances)
	index <- index + 1
}


# if there was more than one alignment, calculate the average distance matrix
if (length(dist_list) > 1) {
	cat(paste0("\n\nCalculating the average distance matrix for ", length(taxa), " taxa...\n"))
	dimensions <- length(taxa)
	avg_matrix <- matrix(0, nrow = dimensions, ncol = dimensions)
	rownames(avg_matrix) <- taxa
	colnames(avg_matrix) <- taxa
	for (index in seq_len(length(taxa))) {
		cat(paste0(" ", index))
		taxon <- taxa[index]
		for (taxon_2 in tail(taxa, -index)) {
			distances <- c()
			found <- FALSE
			for (dist_mat in dist_list) {
				if (taxon %in% rownames(dist_mat) && taxon_2 %in% rownames(dist_mat)) {
					distances <- c(distances, dist_mat[taxon, taxon_2])
					found <- TRUE
				} else {
					distances <- c(distances, NA)
				}
			}
			avg_dist <- mean(as.numeric(distances), na.rm = TRUE)
			if (! found) {
				cat(paste0("\nWarning: did not find a distance comparison for ", taxon,
					" and ", taxon_2, "!\n"))
			}
			avg_matrix[taxon, taxon_2] <- avg_dist
			avg_matrix[taxon_2, taxon] <- avg_dist
		}
	}
	distances <- avg_matrix
} else {
	distances <- as.matrix(dist_list[[1]])
}


# if the rownames have spaces (e.g. "G. alba"), need to quote for Nexus output
space_present <- FALSE
for (rowname in rownames(distances)) {
	if (length(strsplit(rowname, " ")[[1]]) > 1) {
		space_present <- TRUE
		break
	}
}
if (space_present) {
	rownames(distances) <- paste0("'", rownames(distances), "'")
	taxa <- paste0("'", taxa, "'")
}


# output the distance matrix in Nexus format
taxa_block <- paste0("BEGIN TAXA;\n\tDIMENSIONS NTAX=", length(taxa), ";\n\t",
	"TAXLABELS ", paste(taxa, collapse = " "), ";\nEND;\n")
dist_block <- paste0("BEGIN DISTANCES;\n\tFORMAT\n\t\tTRIANGLE=BOTH\n\t\tDIAGONAL\n\t\t",
	"LABELS=LEFT\n\t;\n\tMATRIX\n")
outfile <- file("dist_out.nex", open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(distances, file = outfile, col.names = FALSE, append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)
cat("\n\nDone!\n")
