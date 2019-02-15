# Load packages ----

library(tidyverse)
library(ggmap)
library(Cairo)

# update ggmap from github
devtools::install_github("dkahle/ggmap")


# set google api key
register_google(key = "AIzaSyBOFuziwPbmD9qFt0EwomTSLyUpo5Guqx0")



# examples ----
get_googlemap("iowa", zoom = 6, maptype = "osm") %>% ggmap() 


us <- c(left = -125, bottom = 25.75, right = -67, top = 49)
get_stamenmap(us, zoom = 5, maptype = "toner-lite") %>% ggmap() 


qmplot(Lon, Lat, data = aphid_sites, maptype = "toner-lite", color = I("red"))
