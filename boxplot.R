
##########
# Author: B.M. Anderson
# Date: Aug 2026
# Description: plot default R stats boxplot(s) for column(s) of a tab-separated file
#	The input file (first argument) should have a header for the title of each column
#	The output will be named as the first "."-demarcated field of the input file name, with added suffix "_plot#.pdf"
##########


# Define a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to plot boxplot(s) for column(s) of a tab-separated file\n")
		cat("The input file should have a header row with the title(s) of the column(s)\n")
		cat("Usage: Rscript boxplot.R input_file\n")
	} else {
		cat(help_message)
	}
}


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	input_file <- args[1]
}


# Read the input file
prefix <- strsplit(basename(input_file), "\\.")[[1]]
mydata <- read.table(input_file, sep = "\t", header = TRUE)
multiple <- FALSE
if (ncol(mydata) > 1) {
	multiple <- TRUE
}


# for each column, plot a boxplot
if (multiple) {
	for (column in seq_len(ncol(mydata))) {
		pdf(paste0(prefix, "_plot", column, ".pdf"), width = 4, height = 12)
		boxplot(mydata[, column],
			ylim = c(0, max(mydata[, column], na.rm = TRUE) * 1.05),
			main = colnames(mydata)[column])
		invisible(dev.off())
	}
} else {
	column <- 1
	pdf(paste0(prefix, "_plot.pdf"), width = 4, height = 12)
	boxplot(mydata[, column],
		ylim = c(0, max(mydata[, column], na.rm = TRUE) * 1.05),
		main = colnames(mydata)[column])
	invisible(dev.off())
}
