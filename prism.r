# Compute GDDs from PRISM data ---------------------------------------------------

# For handling weather data acquired from http://prism.oregonstate.edu/explorer/bulk.php
# Combines prism csvs into single file

library(tidyverse)

files = choose.files() # pick prism output sheets

# append all prism files into single data frame
for (i in 1:length(files)){
  if (i == 1) {prism_in = NULL} # clear temp file on first loop
  prism_in <-
    rbind(prism_in, read.csv(files[i], skip = 10)) # append each file
}

# assign column names
names(prism_in) <- c("SiteID",
                     "Longitude",
                     "Latitude",
                     "ElevFt",
                     "Date",
                     "pptIn",
                     "tminF",
                     "tmaxF")

# define GDD function
fn.gdd <- function(tmin, tmax, lower = 50, upper = 86) {
  pmax(0, (pmax(tmin, lower) + pmin(tmax, upper)) / 2 - lower)
}

# compute new columns
prism_join = prism_in %>%
  mutate(Date = as.Date(Date),
         Year = format(Date, "%Y")) %>%
  group_by(SiteID, Year) %>%
  arrange(SiteID, Date) %>%
  mutate(GDD39 = cumsum(fn.gdd(tminF, tmaxF, 39, 86)),
         GDD50 = cumsum(fn.gdd(tminF, tmaxF, 50, 86))) %>%
  ungroup()

# add new data to prism
prism <- rbind(prism, prism_join)

# fix column types
prism$SiteID <- as.factor(prism$SiteID)
prism$Year <- as.integer(prism$Year)
str(prism)

# save prism data
prism %>% write.csv("prism.csv")