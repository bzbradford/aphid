# Load packages ----

library(tidyverse)
library(ggmap)
library(Cairo)

# update ggmap from github
devtools::install_github("dkahle/ggmap")


# set google api key
register_google(key = Sys.getenv(map_key))



# examples ----
get_googlemap("iowa", zoom = 6, maptype = "osm") %>% ggmap() 


us <- c(left = -125, bottom = 25.75, right = -67, top = 49)
get_stamenmap(us, zoom = 5, maptype = "toner-lite") %>% ggmap() 


qmplot(Lon, Lat, data = aphid_sites, maptype = "toner-background", color = I("red"))


glycines_jul_oct %>%
  qmplot(Lon, Lat, data = ., geom = "point", maptype = "toner-background", zoom = 7, darken = .7, color = "red") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .3, color = NA) +
  scale_fill_gradient2("Mean\nGlycines", low = "white", mid = "yellow", high = "red", midpoint = 10)


glycines_jul_oct %>%
  qmplot(Lon, Lat, data = ., geom = "point", maptype = "toner-background", zoom = 7, darken = .7, color = "red") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .3, color = NA) +
  scale_fill_gradient2("Mean\nGlycines", low = "white", mid = "yellow", high = "red", midpoint = 10)

get_stamenmap(c(left = -90, bottom = 40, right = -70, top = 50),
              zoom = 5, maptype = "terrain") %>%
  ggmap() +
  geom_point()

mapframe = round(c(
  min(aphid_sites$Lon) - 4,
  min(aphid_sites$Lat) - 2,
  max(aphid_sites$Lon) + 4,
  max(aphid_sites$Lat) + 2
))
mapframe = c(-100, 35, -80, 50)
mymap = get_stamenmap(mapframe, zoom = 6, maptype = "terrain-background")
mymap %>% ggmap() +
  stat_density_2d(
    aes(
      x = Lon,
      y = Lat,
      size = MeanCount,
      fill = ..level..
    ),
    data = filter(glycines_jul_oct, Year == 2005),
    geom = "polygon",
    alpha = .3,
    color = NA
  ) +
  geom_point(
    aes(x = Lon, y = Lat),
    data = aphid_sites,
    color = "black",
    size = 2
  )
