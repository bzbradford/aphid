# Compute GDDs from PRISM data ---------------------------------------------------

# For handling weather data acquired from http://prism.oregonstate.edu/explorer/bulk.php
# Combines prism csvs into single file

library(tidyverse)


# Read existing saved prism data ----

prism =
  read.csv("data/prism.csv", header = TRUE) %>%
  mutate(Date = as.Date(Date))



# OR Generate prism data from scratch ----

# pick prism output sheets
files = choose.files()

# append all prism files into single data frame
for (i in 1:length(files)){
  if (i == 1) {prism_in = NULL} # clear temp file on first loop
  prism_in <- rbind(prism_in, read.csv(files[i], skip = 10)) # append each file
}

# assign column names
names(prism_in) =
  c("SiteID",
    "Longitude",
    "Latitude",
    "ElevFt",
    "Date",
    "pptIn",
    "tminF",
    "tmaxF")

# define GDD function
fn.gdd =
  function(tmin, tmax, lower = 50, upper = 86) {
    pmax(0, (pmax(tmin, lower) + pmin(tmax, upper)) / 2 - lower)
  }

# compute new columns
prism =
  prism_in %>%
  mutate(Date = as.Date(Date),
         Year = format(Date, "%Y")) %>%
  group_by(SiteID, Year) %>%
  arrange(SiteID, Date) %>%
  mutate(GDD39 = cumsum(fn.gdd(tminF, tmaxF, 39, 86)),
         GDD50 = cumsum(fn.gdd(tminF, tmaxF, 50, 86))) %>%
  ungroup() %>%
  mutate(SiteID = as.factor(SiteID),
         Year = as.integer(Year))

# save prism data
prism %>% write.csv("data/prism.csv")


# had to remove underscores from SiteID field
# prism %>%
#   mutate(SiteID = gsub("_","",SiteID)) %>%
#   arrange(SiteID, Date) ->
#   prism

