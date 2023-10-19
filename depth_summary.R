
##########
# Author: Ben Anderson
# Date: Nov 2020
# Modified: Oct 2023
# Description: Calculate summary stats from samtools depth output (arg1)
##########

args <- commandArgs()
depths_file <- args[6]


# read in the depths file
data <- read.table(depths_file, header = FALSE, sep = "\t")


# print a summary of the depths column
summary(data$V3)
print(paste0("Standard deviation: ", sd(data$V3)))
