
#####################
# Author: B. Anderson
# Date: Sep 2022
# Description: convert UTM coordinates to longitude and latitude (GDA94)
# Note:	call this script with the first argument being the UTM coordinates text
#		file in the tab-delimited format (no header) 4 columns: ID grid# easting northing
#		The output will be printed to screen and also output to a text file ("out_coords.txt")
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


# read in the text file
args <- commandArgs(trailingOnly = TRUE)
coord_table <- read.table(args[1], sep = "\t", header = FALSE)


# create a dataframe for storing the output
mydf <- data.frame(matrix(nrow = nrow(coord_table), ncol = 3))
colnames(mydf) <- c("ID", "long", "lat")


# subset to one table per zone and convert
zones <- unique(coord_table$V2)
index <- 1
for (zone in zones) {
	my_subset <- coord_table[coord_table$V2 == zone, c(1, 3, 4)]
	# convert to an sf object
	my_sf <- st_as_sf(my_subset, coords = c("V3", "V4"),
		crs = paste0(crs_text, " +zone=", zone))
	# transform to longlat (GDA94 is 4283; WGS84 is 4326)
	new_sf <- st_transform(my_sf, crs = paste0("+proj=longlat ", coord_sys))
	# store to the combo dataframe
	rows <- length(new_sf$V1)
	mydf[index: (index + rows - 1), 1] <- new_sf$V1
	mydf[index: (index + rows - 1), 2:3] <- matrix(unlist(new_sf$geometry), ncol = 2, byrow = TRUE)
	index <- index + rows
}


# report to screen and save to file
outdf <- mydf[order(mydf$ID), ]
print.data.frame(outdf, row.names = FALSE)
outfile <- file("out_coords.txt", open = "w")
write.table(outdf, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE)
close(outfile)
