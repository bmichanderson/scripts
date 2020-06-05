#######################
# Author: B. Anderson; some ideas from Sanchez-Puerta et al.
# Date: 4, 1 June 2020; 31 May 2020; 18 Jun 2019
# Description: plot blastn hits in BED type format (can be obtained from blastn_parse.py)
#######################


# graphical parameters
height <- 5		# plot height
width <- 15		# plot width
margins <- c(3,8,3,8)	# plot margins (bottom, left, top, right)


# a help function for when the script is called without arguments or incorrectly
help <- function(help_message) {
	cat("A script to plot BLASTN hits in BED format with slen\n\n")
	cat("Usage: Rscript blastn_plot.R [-d dict_file] [-c plot_colour] plot_data.tab\n\n")
	cat("Options:\n")
	cat("	-c	The colour scheme of the plot [default Viridis]\n")
	cat("	-d	A dictionary file with format: accession\tequivalent\n")
	cat(help_message)
}


# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help("<This is where the error message will be printed>\n"), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1

	colour <- "Viridis"
	dict_present <- FALSE
	d_file <- ""

	for (index in 1:length(args)) {
		if (args[index] == "-c") {
			colour <- args[index + 1]
		} else if (args[index] == "-d") {
			dict_present <- TRUE
			d_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}


# load the library for plotting (and dependencies)
library(Sushi)


# read in the plot data
input_list <- setdiff(catch_args, c(colour, d_file))		# exclude the options that were added to catch_args

if (length(input_list) == 1) {
	data <- read.table(input_list[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
	stop(help("More than one plotting file or too many command line arguments"), call. = FALSE)
}


# establish a colour palette
if (colour %in% hcl.pals(type = "sequential")) {
	col_pal <- function(n, ...) { hcl.colors(n, colour, ...) }
} else if (colour %in% colors()) {
	col_pal <- colorRampPalette(c("black", "dodgerblue", colour))
} else {
	cat("Colour not recognized; using green\n")
	col_pal <- colorRampPalette(c("black", "green"))
}


# read in the dictionary (if present) and reset labels
if (dict_present) {
	dict <- read.table(d_file, sep = "\t", stringsAsFactors = FALSE)
	data_working <- data
	for (row in (1: nrow(dict))) {
		data_working$query <- sub(paste0(dict[row, 1], ".*"), dict[row, 2], data_working$query)
		data_working$sbjct <- sub(paste0(dict[row, 1], ".*"), dict[row, 2], data_working$sbjct)
	}
	data <- data_working
}


# start making a pdf
cat("Plotting", length(unique(data$sbjct)), "plots to a single plot.pdf.\n")
pdf("plot.pdf", width = width, height = height)
par(mar = margins)


# determine a range of values for plotting colours consistently across plots
col_range <- range(floor(min(data$score)), 100)


# split data into separate subjects
sbjct_list <- split(data, f = data$sbjct)


# plot for each subject
for (sbjct_df in sbjct_list) {
	chromstart <- 0	

	# set chromosome name and BED data format
	chrom <- as.character(sbjct_df$sbjct[1])
	chromend <- as.numeric(sbjct_df$slen[1])
	beddata <- sbjct_df[c("sbjct", "sbjct_start", "sbjct_end", "query", "score", "sstrand")]
	beddata <- beddata[order(beddata$query, -beddata$score, decreasing = TRUE), ]

	# create a new column of values for rows that have the same query
	index <- 1
	beddata$row_num <- 0
	for (query in unique(beddata$query)) {
		beddata$row_num[beddata$query == query] <- index
		index <- index + 1
	}

	# plot so that each query or query bin has its own row
	if (nrow(beddata) > 1) {
		plotBed(	beddata = beddata, chrom = chrom, chromstart = chromstart, chromend = chromend, 
				colorby = beddata$score, colorbycol = col_pal, colorbyrange = col_range, type = "region",	
				row = "supplied", rowlabels = unique(beddata$query), rownumber = beddata$row_num, rowlabelcex = 0.5, 
				plotbg = "white")
	} else {
		plotBed(	beddata = beddata, chrom = chrom, chromstart = chromstart, chromend = chromend, 
				colorby = beddata$score, colorbycol = col_pal, colorbyrange = col_range, type = "region",	
				row = "auto", rowlabels = unique(beddata$query), rowlabelcex = 0.75,
				plotbg = "white")
	}
	labelgenome(	chrom, chromstart, chromend, scale = "Kb")
	addlegend(	col_range, side = "right", palette = col_pal, title = "identity", title.offset = 0.05, 
			bottominset = 0.05, topinset = 0.05, labelside = "right", width = 0.025, tick.num = 8, xoffset = 0.005)
}

# stop making a pdf
dev.off()

