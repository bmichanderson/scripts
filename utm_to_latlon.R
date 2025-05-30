
#####################
# Author: B.M. Anderson
# Date: Sep 2022
# Modified: May 2025
# Description: convert UTM coordinates to longitude and latitude (GDA94)
# Note: call this script with the first argument being the UTM coordinates text file
#	The file should have three or four columns: (sampleID) zone/grid# easting northing
#	The output will be printed to screen and also output to a text file ("out_coords.txt")
#####################


suppressMessages(library(sf))


# set whether southern hemisphere
south <- TRUE
crs_text <- "+proj=utm"
if (south) {
	crs_text <- paste0(crs_text, " +south")
}


# set coord system (GDA94 is 4283; WGS84 is 4326)
coord_sys <- "4283"


# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to convert UTM coordinates to longitude and latitude (GDA94)\n")
		cat("Usage: Rscript utm_to_latlon.R coord_file\n")
		cat("The coord_file should have an optional ID colum and three columns = zone easting northing\n")
	} else {
		cat(help_message)
	}
}


# read in the text file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	coord_table <- read.table(args[1], sep = "\t", header = FALSE)
}


# check the format fo the coord_file and add ID if needed
if (ncol(coord_table) == 3) {		# no ID
	id_column <- seq_len(nrow(coord_table))
	coord_table <- cbind(id_column, coord_table)
	colnames(coord_table) <- c("ID", "zone", "easting", "northing")
} else if (ncol(coord_table) == 4) {		# ID present
	colnames(coord_table) <- c("ID", "zone", "easting", "northing")
} else {
	stop(help("Improperly formatted coordinates file!"), call. = FALSE)
}


# create a dataframe for storing the output
mydf <- data.frame(matrix(nrow = nrow(coord_table), ncol = 3))
colnames(mydf) <- c("ID", "Longitude", "Latitude")


# subset to one table per zone and convert
zones <- unique(coord_table$zone)
index <- 1
for (zone in zones) {
	my_subset <- coord_table[coord_table$zone == zone, ]
	# convert to an sf object
	my_sf <- st_as_sf(my_subset, coords = c("easting", "northing"),
		crs = paste0(crs_text, " +zone=", zone))
	# transform to longitude and latitude
	new_sf <- st_transform(my_sf, crs = paste0("+proj=longlat ", coord_sys))
	# store to the combo dataframe
	rows <- length(new_sf$ID)
	mydf[index: (index + rows - 1), 1] <- new_sf$ID
	mydf[index: (index + rows - 1), 2:3] <- matrix(unlist(new_sf$geometry), ncol = 2, byrow = TRUE)
	index <- index + rows
}


# report to screen and save to file
outdf <- mydf[order(mydf$ID), ]
print.data.frame(outdf, row.names = FALSE)
outfile <- file("out_coords.txt", open = "w")
write.table(outdf, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE)
close(outfile)
