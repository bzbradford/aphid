# Script to collect dbf files from arcpy output
# and merge into single data frame


library(foreign)
library(tidyverse)

# specify location of folder structure
oldwd <- getwd()
setwd("c:/arc/aphids/buffers/tables")  

# generate file list
files <- list.files(pattern = "dbf$")
numfiles <- length(files)

for (i in 1:numfiles) {
  if(i == 1) {
    merge <- data.frame() # initialize
    }
  infile <- read.dbf(files[i])
  infile <- filter(infile, infile[2] > 0) # remove nonzero counts (count must be 2nd column)
  names(infile) <- c("Class", "Pixels") # rename columns
  parse <- unlist(strsplit(gsub(".dbf", "", files[i]), "_")) # parse SITE, BUFFER, YEAR from filename
  infile <- infile %>%
    add_column(Site = parse[1],
               Buffer_m = as.numeric(gsub("m", "", parse[2])),
               Year = as.numeric(parse[3]),
               .before = "Class")
  infile$Class <- gsub("_", " ", infile$Class) # swap _ for space in class names
  merge <- bind_rows(merge, infile) # add latest block to growing data frame
}
merge$Class = as.factor(merge$Class) # return to factor from character vector

# clean up
setwd(oldwd)
rm(oldwd, files, numfiles, infile, parse)

# output if desired
write.csv(merge, file = "buffers.csv")  ## write output file
