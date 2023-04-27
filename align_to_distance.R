
##########
# Author: Ben Anderson
# Date: Nov 2022
# Description: convert a DNA alignment (fasta) into a distance matrix (Nexus format)
# Note: the resulting distance matrix will be saved as "dist_out.nex"
##########


# load required library
suppressMessages(library(ape))


# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to convert a DNA alignment (fasta) to distances (Nexus)\n")
		cat("Usage: Rscript align_to_distance.R [-m model] [-p pofadinr method] [-s samples] alignment.fasta\n")
		cat("Options:\n")
		cat("\t-m\tThe ape DNA distance model [default F84]\n")
		cat("\t-p\tThe pofadinr nucleotide distance method (\"g\" for GENPOFAD, \"m\" for MATCHSTATES)\n")
		cat("\t\tIf set, this will be the method used rather than ape dist.dna\n")
		cat("\t-s\tA file with tab-separated sample IDs and desired tip labels, one per line [optional]\n")
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
	samples_present <- FALSE
	pofad_method_set <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-m") {
			model <- args[index + 1]
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
alignment <- read.FASTA(catch_args[[1]], type = "DNA")


# if samples file present, change the labels
if (samples_present) {
	for (index in seq_len(length(names(alignment)))) {
		name <- names(alignment[index])
		if (name %in% sample_table$V1) {
			names(alignment)[index] <- sample_table$V2[match(name, sample_table$V1)]
		}
	}
}


# calculate the distances
if (pofad_method_set) {
	suppressMessages(library(pofadinr))
	# convert missing data to "?"
	temp <- as.character(alignment)
	temp[temp == "n"] <- "?"
	dnabin <- as.DNAbin(temp)
	if (pofad_method == "m") {
		cat(paste0("Calculating distances using pofadinr dist.snp and the MATCHSTATES method...\n"))
		distances <- dist.snp(dnabin, model = "MATCHSTATES")
	} else if (pofad_method == "g") {
		cat(paste0("Calculating distances using pofadinr dist.snp and the GENPOFAD method...\n"))
		distances <- dist.snp(dnabin, model = "GENPOFAD")
	} else {
		stop(help("Please specify pofadinr method as \"m\" or \"g\"\n"), call. = FALSE)
	}
} else {
	cat(paste0("Calculating distances using ape dist.dna and the ", model, " model...\n"))
	distances <- dist.dna(alignment, model = model, pairwise.deletion = TRUE)
}


# output the distance matrix in Nexus format
taxa <- names(alignment)
taxa_block <- paste0("BEGIN TAXA;\n\tDIMENSIONS NTAX=", length(taxa), ";\n\t",
	"TAXLABELS ", paste(taxa, collapse = " "), ";\nEND;\n")
dist_block <- paste0("BEGIN DISTANCES;\n\tFORMAT\n\t\tTRIANGLE=BOTH\n\t\tDIAGONAL\n\t\t",
	"LABELS=LEFT\n\t;\n\tMATRIX\n")
outfile <- file("dist_out.nex", open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(as.matrix(distances), file = outfile, col.names = FALSE,
	append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)
cat("Done!\n")
