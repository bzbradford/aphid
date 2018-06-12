## elissa code

##Function needed for calculating sine wave heat units
calcHeat <- function(fk1, tsum, diff) {
  twopi = 6.283185308
  pihlf = 1.570796327
  d2 = fk1 - tsum
  theta = atan(d2 / sqrt(diff * diff - d2 * d2))
  if ((d2 < 0) & (theta > 0)) {
    theta = theta - pi
  }
  return((diff * cos(theta) - d2 * (pihlf - theta)) / twopi)
}

Degday_sine <- function(species) {
  if (species == "plant") {
    TH = 85
    TL = 41
  }
  if (species == "SFW") {
    TH = 86
    TL = 50
  }
  if (species == "CFW") {
    TH = 87
    TL = 44
  }
  GDD <-
    c() # empty vector for degree day accumulations for each location
  for (i in 1:46) {
    dd <- c() # empty vector for each degree day total
    filename <-
      paste("./location data/loc", i, ".txt", sep = "") #read in data from location data folder
    dat <- read.table(filename, sep = "\t")
    #print(sum(is.na(dat)))
    dat <-
      na.omit(dat) # looking to find a way to impute missing values in above function
    for (j in 1:nrow(dat)) {
      tmax = dat[j, 3]
      tmin = dat[j, 2]
      heat = 0
      fk1 = 2 * TL
      diff = tmax - tmin
      tsum = tmax + tmin
      
      #return 0 if invalid inputs or max below TL
      if ((tmin > tmax) | (tmax <= tmin) | (tmax <= TL))
        #return(0)
      {
        DDay <- 0
      }
      else if (tmin >= TL)
      {
        DDay <- (tsum - fk1) / 2
      }
      else if (tmin < TL)
      {
        DDay <- calcHeat(fk1, tsum, diff)
      }
      else if (tmax > TH) {
        fk1 = 2 * TH
        zheat = heat
        heat = calcHeat(fk1, tsum, diff)
        #if (cutoff == "intermediate")
        #{
        DDay <- zheat - 2 * heat
        #}
        #if (cutoff == "horizontal")
        #{
        #       DDay <- zheat - heat
        #}
      }
      dd <- append(dd, DDay)
    }
    DDsum <- sum(dd)
    GDD <- append(GDD, DDsum) #vector of all GDD accumulations
  }
  #write file that will contain location data (from id_list) and GDD accumuations
  id_list <- read.csv("StationIDs.csv") #read in file with ID names
  df <-
    data.frame(
      ID = id_list$ID,
      lat = id_list$Lat,
      long = id_list$Long,
      location = id_list$Station,
      deg = GDD
    )
  fname <- paste0("GDDdata", species, ".txt")
  write.table(df, fname, row.names = FALSE, col.names = TRUE)
}


map_gdd_static <- function(filename, date, organism) {
  require(ggplot2) #load libraries
  require(rgdal) #load libraries
  dat <- read.table(filename, header = TRUE) #read in GDD data
  #read in WI shape file
  counties = readOGR(
    dsn = path.expand("shape files/WI_County_Boundaries/WI_County_Bnds.shp"),
    layer = "WI_County_Bnds"
  )
  coordinates(dat) <-  ~ long + lat
  proj4string(dat) <-
    CRS("+proj=longlat +datum=NAD83") #match projections between shape file and point file
  dat <-
    spTransform(dat, CRS(proj4string(counties))) #match projections between shape file and point file
  dat <- data.frame(dat)
  if (organism == "plant") {
    chart_title <-
      paste("Cranberry Growing Degree Days", date, sep = ": ")
    # make plot
    ggplot() +
      geom_polygon(
        data = counties,
        aes(x = long, y = lat, group = group),
        fill = "grey20",
        colour = "grey100",
        alpha = 1
      ) +
      labs(x = "", y = "", title = chart_title) + #labels
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        # get rid of x ticks/text
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        # get rid of y ticks/text
        plot.title = element_text(
          size = 16,
          lineheight = .8,
          face = "bold",
          vjust = 1
        )
      ) + # make title bold and add space
      geom_point(
        aes(x = long, y = lat, color = deg),
        data = dat,
        alpha = 1,
        size = 8
      ) +
      scale_colour_gradientn("GDD", colours = c("green", "red")) + # change color scale
      coord_equal(ratio = 1) # square plot to avoid the distortion
    ggsave(
      file = paste("Cran_", date, ".png", sep = ""),
      width = 5.8,
      height = 4.7
    )
  }
  if (organism == "sparg") {
    chart_title <- paste("Sparganothis Degree Days", date, sep = ": ")
    # make plot
    ggplot() +
      geom_polygon(
        data = counties,
        aes(x = long, y = lat, group = group),
        fill = "grey20",
        colour = "grey100",
        alpha = 1
      ) +
      labs(x = "", y = "", title = chart_title) + #labels
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        # get rid of x ticks/text
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        # get rid of y ticks/text
        plot.title = element_text(
          size = 16,
          lineheight = .8,
          face = "bold",
          vjust = 1
        )
      ) + # make title bold and add space
      geom_point(
        aes(x = long, y = lat, color = deg),
        data = dat,
        alpha = 1,
        size = 8
      ) +
      scale_colour_gradientn("GDD", colours = c("purple", "orange")) + # change color scale
      coord_equal(ratio = 1) # square plot to avoid the distortion
    ggsave(
      file = paste("Sparg_", date, ".png", sep = ""),
      width = 5,
      height = 4.7
    )
  }
  if (organism == "cfw") {
    chart_title <-
      paste("Cranbery Fruitworm Degree Days", date, sep = ": ")
    # make plot
    ggplot() +
      geom_polygon(
        data = counties,
        aes(x = long, y = lat, group = group),
        fill = "grey20",
        colour = "grey100",
        alpha = 1
      ) +
      labs(x = "", y = "", title = chart_title) + #labels
      theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        # get rid of x ticks/text
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        # get rid of y ticks/text
        plot.title = element_text(
          size = 16,
          lineheight = .8,
          face = "bold",
          vjust = 1
        )
      ) + # make title bold and add space
      geom_point(
        aes(x = long, y = lat, color = deg),
        data = dat,
        alpha = 1,
        size = 8
      ) +
      scale_colour_gradientn("GDD", colours = c("blue", "yellow")) + # change color scale
      coord_equal(ratio = 1) # square plot to avoid the distortion
    ggsave(
      file = paste("CFW_", date, ".png", sep = ""),
      width = 5.8,
      height = 4.7
    )
  }
}
