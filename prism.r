
# For handling weather data acquired from http://prism.oregonstate.edu/explorer/bulk.php
# Combines prism csvs into single file

library(tidyverse)


# GDD functions -----------------------------------------------------------

# simple min/max
gdd_simple <- function(tmin, tmax, lower = 50, upper = 86) {
  pmax(0, (pmax(tmin, lower) + pmin(tmax, upper)) / 2 - lower)
}

gdd_sine <- function(min, max, base, upper) {
  if (length(min) > 1 || length(max) > 1) stop("Function must be called with mapply for use with vectors")
  if (min > max) stop("Min cannot be greater than Max")
  if (base > upper) stop("Base cannot be greater than Upper")
  
  # min and max > upper
  if (min >= upper) return(upper - base)
  
  # min and max < lower
  if (max <= base) return(0)
  
  average = (min + max) / 2
  
  # min and max between base and upper
  if (max <= upper && min >= base) return(average - base)
  
  alpha = (max - min) / 2
  
  # min < base, max between base and upper
  if (max <= upper && min < base) {
    base_radians = asin((base - average) / alpha)
    a = average - base
    b = pi / 2 - base_radians
    c = alpha * cos(base_radians)
    return((1 / pi) * (a * b + c))
  }
  
  # max > upper and min between base and upper
  if (max > upper && min >= base) {
    upper_radians = asin((upper - average) / alpha)
    a = average - base
    b = upper_radians + pi / 2
    c = upper - base
    d = pi / 2 - upper_radians
    e = alpha * cos(upper_radians)
    return((1 / pi) * (a * b + c * d - e))
  }
  
  # max > upper and min < base
  if (max > upper && min < base) {
    base_radians = asin((base - average) / alpha)
    upper_radians = asin((upper - average) / alpha)
    a = average - base
    b = upper_radians - base_radians
    c = alpha * (cos(base_radians) - cos(upper_radians))
    d = upper - base
    e = pi / 2 - upper_radians
    return((1 / pi) * ((a * b + c) + (d * e)))
  }
}



# Load climate data -------------------------------------------------------

# read in PRISM climate data
prism_in <- 
  list.files("prism", full.names = T) %>%
  map_dfr(read_csv, skip = 10) %>%
  rename(
    "TminF" = `tmin (degrees F)`,
    "TmaxF" = `tmax (degrees F)`
  ) %>%
  select(c(Name, Longitude, Latitude, Date, TminF, TmaxF))



# Compute GDD accumulations -----------------------------------------------

# sine method
prism_gdds <-
  prism_in %>%
  mutate(
    Year = format(Date, "%Y"),
    Day = as.numeric(format(Date, "%j"))) %>%
  group_by(Name, Year) %>%
  arrange(Date) %>%
  mutate(GDD39 = mapply(gdd_sine, TminF, TmaxF, 39, 86)) %>%
  mutate(GDD39_cum = cumsum(GDD39)) %>%
  arrange(Name, Date)



# Plot GDDs ---------------------------------------------------------------

prism_gdds %>%
  ggplot(aes(x = Day, y = GDD39_cum, color = Year, group = Year)) +
  geom_line() +
  facet_wrap(~ Name) +
  labs(
    x = "Day of year", y = "GDD accumulation, simple (top) vs sine (bottom)",
    title = "Degree-day accumulation")



