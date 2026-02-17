
##########
# Author: B.M. Anderson
# Date: Nov 2022
# Modified: April 2023, May 2023, July 2023,
# Apr 2025 (changed calculations and reporting to be more efficient for many alignments),
# Feb 2026 (adjusted so distance averaging occurs by chunk to avoid memory issues with very large datasets >10,000 loci)
# Description: convert a DNA alignment (fasta) or multiple alignments into a distance matrix (Nexus format)
# Note: the resulting distance matrix will be saved as "dist_out.nex"
##########


# load required library
suppressMessages(library(ape))


# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to convert a DNA alignment (fasta) or multiple alignments into a distance matrix (Nexus)\n")
		cat("Usage: Rscript align_to_distance.R [-m model] [-p pofadinr method]",
			"[-s samples] alignment.fasta [alignment2.fasta ...]\n")
		cat("Options:\n")
		cat("\t-m\tThe ape DNA distance model [default \"F84\"]\n")
		cat("\t-a\tResolve ambiguous bases randomly before calculating distance in ape (\"y\" or \"n\" [default])\n")
		cat("\t-p\tThe pofadinr nucleotide distance method (\"g\" for GENPOFAD, \"m\" for MATCHSTATES)\n")
		cat("\t\tIf set, this will be the method used rather than ape dist.dna\n")
		cat("\t-s\tA file with tab-separated sample IDs and corresponding tip labels for conversion,",
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


# capture taxa present
taxa <- c()
for (alignment in alignments) {
	taxa <- unique(c(taxa, names(alignment)))
}

cat(paste0("\nRead in ", length(alignments), " alignments with ", length(taxa), " taxa\n\n"))


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
big_data <- FALSE
batch_index <- 1
batch_list <- vector("list")
if (length(alignments) > 10000) {
	increment <- 5000
	big_data <- TRUE
} else if (length(alignments) > 1000) {
	increment <- 500
} else if (length(alignments) > 100) {
	increment <- 50
} else {
	increment <- 1
}
cat("\nProcessed alignment:")
for (alignment in alignments) {
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

	# adjust the matrix to ensure all taxa are present
	# if there are missing taxa, create columns and rows of NAs for them
	dist_mat <- as.matrix(distances)
	for (taxon in taxa) {
		if (! taxon %in% rownames(dist_mat)) {
			curr_names <- rownames(dist_mat)
			dist_mat <- rbind(cbind(dist_mat, rep(NA, nrow(dist_mat))), rep(NA, ncol(dist_mat) + 1))
			rownames(dist_mat) <- c(curr_names, taxon)
			colnames(dist_mat) <- c(curr_names, taxon)
		}
	}

	# reorder to the order of taxa (ensuring same order across all matrices)
	dist_mat <- dist_mat[taxa, taxa]

	# add to the list
	dist_list[[index]] <- dist_mat

	# report progress and increment index
	if (index %% increment == 0) {
		if (big_data) {
			# if there are many alignments, try to save memory usage and run in batches of 1000
			# see: https://stackoverflow.com/a/3321659
			chunks <- split(dist_list, ceiling(seq_along(dist_list) / 1000))
			rm(dist_list)
			gc()		# garbage cleanup
			for (chunk in chunks) {
				# since the matrices have the same dimensions and order, calculate means across matrices
				# see: https://stackoverflow.com/a/19220503
				mat_array <- array(unlist(chunk), c(length(taxa), length(taxa), length(chunk)))
				chunk_dist <- as.matrix(rowMeans(mat_array, dims = 2, na.rm = TRUE))
				batch_list[[batch_index]] <- chunk_dist
				batch_index <- batch_index + 1
			}
			index <- 0
			dist_list <- vector("list")
			cat(paste0(" ", (batch_index - 1) * 1000))
		} else {
			cat(paste0(" ", index))
		}
	}
	index <- index + 1
}

# report final index if not already reported
# if many alignments, measure the remaining that haven't been added to the batch_list
# note: this will be slightly inaccurate for big data because the last chunk will have a smaller length
if (index %% increment > 1) {
	if (big_data) {
		chunks <- split(dist_list, ceiling(seq_along(dist_list) / 1000))
		rm(dist_list)
		for (chunk in chunks) {
			mat_array <- array(unlist(chunk), c(length(taxa), length(taxa), length(chunk)))
			chunk_dist <- as.matrix(rowMeans(mat_array, dims = 2, na.rm = TRUE))
			batch_list[[batch_index]] <- chunk_dist
			batch_index <- batch_index + 1
		}
		dist_list <- batch_list
		cat(paste0(" ", (batch_index - 1) * 1000 + (index %% 1000) - 1, "\n"))
	} else {
		cat(paste0(" ", index - 1, "\n"))
	}
} else {
	if (big_data) {
		dist_list <- batch_list
	}
}


# if there was more than one alignment, calculate the average distance matrix
if (length(dist_list) > 1) {
	# since the matrices have the same dimensions and order, calculate means across matrices
	# see: https://stackoverflow.com/a/19220503
	cat(paste0("\nTransforming into an array..."))
	mat_array <- array(unlist(dist_list), c(length(taxa), length(taxa), length(dist_list)))
	cat(paste0("\nCalculating the average distance matrix for ", length(taxa), " taxa...\n"))
	distances <- as.matrix(rowMeans(mat_array, dims = 2, na.rm = TRUE))
	rownames(distances) <- taxa
} else {
	distances <- as.matrix(dist_list[[1]])
}


# if there was a samples table provided, rename the taxa to match
if (samples_present) {
	for (index in seq_len(length(rownames(distances)))) {
		name <- rownames(distances)[index]
		taxon <- taxa[index]
		if (name != taxon) {
			cat(paste0("\nWARNING: coding error (shouldn't happen). The taxon ", name, " doesn't match the sort!\n"))
		}
		if (name %in% sample_table$V1) {
			rownames(distances)[index] <- sample_table$V2[match(name, sample_table$V1)]
			taxa[index] <- sample_table$V2[match(taxon, sample_table$V1)]
		} else {
			cat(paste0("\nWARNING: taxon ", name, " is missing from samples table and won't be renamed!\n"))
		}
	}
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
