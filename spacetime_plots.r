

####################################################################
### Bring in data file
####################################################################

wd <- getwd()
and.dat <- read.csv(file.choose(), header = TRUE, na = 'n.s.')

####################################################################
### Functions to make time plots
####################################################################

library(plyr)
library(ggplot2)
and.dat2 <-
  subset(and.dat, and.dat$Count > 0) ## Eliminate data when no insects were found

### Time plot
### Loc might be pulled in calc of mn # weeks

time <- function(ll, ind) {
  totyrs <- nlevels(factor(ll$Year[ll[, ind] > 0]))
  yrs <-
    mean(aggregate(ll$Year[ll[, ind] > 0], list(ll$Location[ll[, ind] > 0]), function(x)
      nlevels(factor(x)))$x)
  mnwks <-
    mean(aggregate(ll$Date[ll[, ind] > 0], list(ll$Year[ll[, ind] > 0], ll$Location[ll[, ind] > 0]), function(x)
      nlevels(factor(x)))$x)
  mnabun <-
    mean(aggregate(log1p(ll[, ind]), list(ll$Year, ll$Location), mean)$x)
  data.frame(totyrs, yrs, mnwks, mnabun)
}
Time_out <- ddply(and.dat2, .(Aphid.Species), function(x) {
  time(x, 4)
})
names(Time_out) <-
  c("Pest", "TotYrs", "YearsinLoc", "WeekinYrLoc", "lAbun")

### Plot time

(
  timeP1 <-
    qplot(
      YearsinLoc,
      WeekinYrLoc,
      label = Pest,
      data = Time_out,
      size = lAbun,
      colour = lAbun,
      xlab = "# of Years Detected",
      ylab = "Mean # of Weeks in Year Detected"
    ) +
    scale_size_continuous(range = c(1, 10)) + scale_colour_gradient(
      limits = c(0, 2.4),
      low = "blue",
      high = "orange"
    ) +
    geom_text(
      hjust = 0,
      vjust = 0.5,
      colour = 'black',
      size = 2
    )
)


### Space plot
### Loc might be pulled in calc of mn # weeks

space <- function(ll, ind) {
  Locs <- nlevels(factor(ll$Location[ll[, ind] > 0]))
  mnlocs <-
    mean(aggregate(ll$Location[ll[, ind] > 0], list(ll$Year[ll[, ind] > 0]), function(x)
      nlevels(factor(x)))$x)
  mnlocswk <-
    mean(aggregate(ll$Location[ll[, ind] > 0], list(ll$Year[ll[, ind] > 0], ll$Date[ll[, ind] > 0]), function(x)
      nlevels(factor(x)))$x)
  mnabun <-
    mean(aggregate(log1p(ll[, ind]), list(ll$Year, ll$Location), mean)$x)
  data.frame(Locs, mnlocs, mnlocswk, mnabun)
}

Space_out <-
  ddply(and.dat2, .(Aphid.Species), function(x) {
    space(x, 4)
  })
names(Space_out) <- c("Pest", "Locs", "LocinYr", "LocinWk", "lAbun")

### Plot space

(
  spaceP1 <-
    qplot(
      LocinYr,
      LocinWk,
      label = Pest,
      data = Space_out,
      size = lAbun,
      colour = lAbun,
      xlab = "# of Location Detected",
      ylab = "Mean # of Locations in a Week Detected"
    ) +
    scale_size_continuous(range = c(3, 13)) + scale_colour_gradient(
      limits = c(0.001, 2.4),
      low = "green",
      high = "red"
    ) +
    geom_text(
      hjust = 0,
      vjust = 0.5,
      colour = 'black',
      size = 2
    )
)

#comb<-cbind(Space_out, Time_out)

#(spaceP2<-qplot(LocinYr, WeekinYr, label = Pest, data = comb, size = lAbun, colour = lAbun,
#xlab = "Mean # of Locations in a Year Detected", ylab = "Mean # of Weeks in Year Detected")+
#scale_size_continuous(range = c(3,13))+ scale_colour_gradient(limits=c(0.001, 2.4), low="green", high="red")+
#geom_text(hjust=0, vjust=0.5,colour = 'black', size = 2))


pdf('bubbles.pdf', h = 8.5, w = 11)
timeP1
spaceP1
dev.off()
