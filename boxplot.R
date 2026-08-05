
##########
# Author: B.M. Anderson
# Date: Aug 2026
# Description: plot default R stats boxplot(s) for column(s) of a tab-separated file
#	The input file (first argument) should have a header for the title of each column
#	The output will be named as the first "."-demarcated field of the input file name, with added suffix "_plot#.pdf"
#	If there are categories in the first column (not data), indicate this with "-c y"
##########


# Define a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to plot boxplot(s) for column(s) of a tab-separated file\n")
		cat("The input file should have a header row with the title(s) of the column(s)\n")
		cat("Option:\n")
		cat("\t-c\tWhether the first column contains categories for splitting, \"y\" or \"n\" [default]\n")
		cat("Usage: Rscript boxplot.R input_file [-c \"y\" or \"n\"]\n")
	} else {
		cat(help_message)
	}
}


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	extra <- 1
	catch <- TRUE
	categories <- "n"
	for (index in seq_len(length(args))) {
		if (args[index] == "-c") {
			categories <- args[index + 1]
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
	stop(help("Missing input file!\n"), call. = FALSE)
} else if (extra > 2) {
	stop(help("Too many arguments\n"), call. = FALSE)
}


# Read the input file
input_file <- catch_args[[1]]
prefix <- strsplit(basename(input_file), "\\.")[[1]]
mydata <- read.table(input_file, sep = "\t", header = TRUE)
multiple <- FALSE


# Determine plot dimensions depending on whether there are multiple categories
if (categories == "y") {
	cats <- TRUE
	number_cats <- length(unique(mydata[, 1]))
	pwidth <- number_cats * 2
	pheight <- 12
	first_column <- 2
	if (ncol(mydata) > 2) {
		multiple <- TRUE
	}
} else if (categories == "n") {
	cats <- FALSE
	pwidth <- 4
	pheight <- 12
	first_column <- 1
	if (ncol(mydata) > 1) {
		multiple <- TRUE
	}
} else {
	stop(help("Categories option specified incorrectly!\n"), call. = FALSE)
}


# plot boxplot(s)
if (cats) {
	mylevels <- unique(mydata[, 1])
	mydata[, 1] <- factor(mydata[, 1], levels = mylevels)
	# define a function to plot the boxplot with categories
	plot_boxplot <- function(dframe, column) {
		boxplot(dframe[, column] ~ dframe[, 1],
			ylim = c(0, max(dframe[, column], na.rm = TRUE) * 1.05),
			main = colnames(dframe)[column],
			ylab = "", xlab = "", xaxt = "n")
		yrange <- max(dframe[, column], na.rm = TRUE) * 1.05
		axis(1, at = seq_len(length(mylevels)), labels = FALSE)
		text(seq_len(length(mylevels)), par("usr")[3] - (0.02 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)
	}
} else {
	# define a function to plot the boxplot without categories
	plot_boxplot <- function(dframe, column) {
		boxplot(dframe[, column],
			ylim = c(0, max(dframe[, column], na.rm = TRUE) * 1.05),
			main = colnames(dframe)[column])
	}
}

if (multiple) {
	iter <- 1
	for (column in seq(first_column, ncol(mydata))) {
		pdf(paste0(prefix, "_plot", iter, ".pdf"), width = pwidth, height = pheight)
		plot_boxplot(mydata, column)
		invisible(dev.off())
		iter <- iter + 1
	}
} else {
	pdf(paste0(prefix, "_plot.pdf"), width = pwidth, height = pheight)
	plot_boxplot(mydata, first_column)
	invisible(dev.off())
}
