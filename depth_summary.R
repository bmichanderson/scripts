

########
# Author: Ben Anderson
# Date: Nov 2020
# Description: Calculate summary stats for read depths from samtools (arg1)
########

args <- commandArgs()
depths_file <- args[6]


# read in the depths file
data <- read.table(depths_file, header=FALSE, sep='\t')


# print a summary of the depths column
summary(data$V3)
