library(foreign)
library(tidyverse)

setwd("c:/arc/aphids/buffers/tables")   ## location of folder structure


if(file.exists("merge.dbf")){  ## remove previous merge files
  file.remove("merge.dbf")
}

files <- list.files(pattern = "dbf$")  ## get list of files in folder
numfiles <- length(files)


for (i in 1:numfiles) { ## append each dbf to working file
  if(i == 1) {
    all = data.frame()
    }
  infile = read.dbf(files[i])
  infile = filter(infile, infile[2] > 0)
  names(infile) = c("Class", "Pixels")
  vars = unlist(strsplit(gsub(".dbf", "", files[i]), "_"))
  infile = infile %>%
    add_column(Site = vars[1],
               Buffer_m = as.numeric(gsub("m", "", vars[2])),
               Year = as.numeric(vars[3]),
               .before = "Class")
  all = bind_rows(all, infile)
}

all$Class = as.factor(all$Class)

write.csv(all, file = "~merge.csv")  ## write output file