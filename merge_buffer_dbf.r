# Script to collect dbf files from arcpy output
# and merge into single data frame

library(foreign)
library(tidyverse)


# specify location of folder structure
oldwd <- getwd()
setwd("c:/arcgis/aph/tables")

# generate file list
files <- list.files(pattern = "dbf$")
numfiles <- length(files)

for (i in 1:numfiles) {
  if(i == 1) {merge = data.frame()} # initialize temp file
  infile = read.dbf(files[i])
  names(infile) = c("Class", "Pixels") # rename columns
  parse = unlist(strsplit(gsub(".dbf", "", files[i]), "_")) # parse SITE, BUFFER, YEAR from filename
  if (parse[3] %in% c("2013", "2014", "2015")) {  # fix bullshit off by one error in zonal histogram
    infile = na.omit(mutate(infile, Pixels = lag(Pixels)))
  }
  infile <- infile %>%
    add_column(Site = parse[1],
               Buffer_m = as.numeric(gsub("m", "", parse[2])),
               Year = as.numeric(parse[3]),
               .before = "Class")
  infile$Class <- gsub("_", " ", infile$Class) # swap _ for space in class names
  merge <- bind_rows(merge, infile) # add latest block to growing data frame
}
merge$Class = as.factor(merge$Class) # return to factor from character vector
merge2 = filter(merge, Pixels > 0)

# clean up
setwd(oldwd)
rm(oldwd, files, numfiles, infile, parse)

# output if desired
write.csv(merge2, file = "L:/aphids/cdl_extract.csv")  ## write output file
