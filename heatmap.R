
##########
# Author: Ben Anderson
# Date: May 2024
# Description: create a heatmap based on an input tab-delimited text file (with header)
# Note: The file should have display sample labels in the first column, values as subsequent columns
##########


# Define functions

## a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to create a heatmap from an input tab-delimited text file\n")
		cat("Usage: Rscript heatmap.R -o out_pref -f input.tab -c color_pal\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-f\tThe input tab-delimited file (with header without a label for first column)\n")
		cat("\t-c\tColour palette to use [default \"greens\"]\n")
	} else {
		cat(help_message)
	}
}

## a function to plot a heatmap
## for making legend, see https://stackoverflow.com/a/13355440 and https://stackoverflow.com/a/70522655
heatmapper <- function(dmat, palette = "greens", ...) {
	# heatmap
	image(t(dmat),
		col = hcl.colors(n = 100, palette = palette, rev = TRUE),
		axes = FALSE,
		main = expression("Heatmap"),
		useRaster = TRUE)
	axis(2, at = seq(0, 1, length.out = nrow(dmat)),
		labels = rownames(dmat), las = 2, lwd.ticks = 0)

	# legend (still working on this)
	subx <- grconvertX(c(1.1, 1.3), from = "user", to = "ndc")
	suby <- grconvertY(c(0.9, 1), from = "user", to = "ndc")
	op <- par(fig = c(subx, suby),
		mar = c(0, 0, 0, 0),
		new = TRUE)
	legend_colours <- as.raster(hcl.colors(n = 100, palette = palette))
	range_min <- min(dmat, na.rm = TRUE)
	range_max <- max(dmat, na.rm = TRUE)
	if (range_max - range_min < 10) {
		# might want to use a decimal place
		mult <- 10
		step <- round((range_max - range_min) / 5, 1)
		rnum <- 1
	} else {
		# use integers
		mult <- 1
		step <- round((range_max - range_min) / 50, 0) * 10
		rnum <- 0
	}
	digits <- floor(log10(range_max)) + 1
	if (digits < 3) {
		adjust <- 1
	} else {
		adjust <- digits - 2
	}
	legend_seq <- seq(ceiling(range_min * mult) / mult,
		floor(range_max * mult) / mult,
		by = step)
	legend_seq <- floor(legend_seq / (10 ^ adjust)) * (10 ^ adjust)
	if (legend_seq[1] < range_min) {
		legend_seq <- legend_seq[2: length(legend_seq)]
	}
	legend_labels <- format(round(legend_seq, rnum), nsmall = rnum)
	plot(x = c(0, 2), y = c(0, 1), type = "n",
		axes = FALSE, xlab = "", ylab = "", main = "")
	axis(side = 4, at = (legend_seq - range_min) / (range_max - range_min), pos = 1, labels = FALSE,
		col = 0, col.ticks = 1)
	mtext(legend_labels, side = 4, line = -0.5, at = (legend_seq - range_min) / (range_max - range_min), las = 2)
	rasterImage(legend_colours, xleft = 0, ybottom = 0,
		xright = 1, ytop = 1)
	par(op)
}


# Read in and format the data

# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	out_pref <- "output"
	input_present <- FALSE
	col_pal <- "greens"
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
		} else if (args[index] == "-f") {
			input_present <- TRUE
			input_file <- args[index + 1]
		} else if (args[index] == "-c") {
			col_pal <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! input_present) {
	stop(help("Missing argument for input file!\n"), call. = FALSE)
}

# read in the input file
dmat <- as.matrix(read.csv(input_file, sep = "\t", header = TRUE, row.names = 1))

# determine graphics parameters to best fit the data
num_samples <- nrow(dmat)
# for Arial 12 point, each letter is ~ 1/6 inch; we want to have a buffer of about 25%, so * 1.25
my_height <- num_samples * 0.21
if (my_height < 11) {
	my_height <- 11
}
# determine width based on how many features are being displayed
num_elements <- ncol(dmat)
my_width <- num_elements * 0.24
if (my_width < 8) {
	my_width <- 8
}


# start creating a pdf
pdf(paste0(out_pref, "_heatmap.pdf"), width = my_width, height = my_height, family = "ArialMT")
# plot the heatmap
par(mar = c(5, 4, 4, 6) + 0.1)
heatmapper(dmat, palette = col_pal)
# stop creating the pdf
invisible(dev.off())


# similarly, a png
png(paste0(out_pref, "_heatmap.png"), width = my_width, height = my_height, units = "in", res = 600)
par(mar = c(5, 4, 4, 6) + 0.1)
heatmapper(dmat, palette = col_pal)
invisible(dev.off())
