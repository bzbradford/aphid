# Load packages -----------------------------------------------------------

library(lme4)
library(mgcv)
library(ggplot2)
library(tidyverse)

# Read in data ------------------------------------------------------------

aphid_in = read.csv(file.choose(), header = TRUE, na = c('','.')) # aphid data in dummied long format
aphid_in$Date = as.Date(aphid_in$Date)
aphid_sp = read.csv(file.choose(), header = TRUE, na = c('','.')) # load species


# Expand dataset with zeros -------------------------------------------------------

fn.expand = function(df) {
  cols = ncol(df)
  df %>%
    spread(SpeciesName,
           Count,
           fill = 0) %>%
    gather((cols - 1):ncol(.),
           key = "SpeciesName",
           value = "Count",
           factor_key = TRUE) %>%
    filter(SpeciesName != "_dummy_") %>%
    droplevels()
}
aphid_long <- fn.expand(aphid_in)

# File export -------------------------------------------------------------

aphid_long %>%
  filter(Count > 0) %>%
#  select(-c(GDD, isFound)) %>%
  arrange(SiteName, Date) %>%
  write.csv(file.choose())

# Aphid data subset options ---------------------------------------------------

# Species list 1 (from Frost code)
aphid <-
  aphid_long %>%
  filter(
    SpeciesName %in% c(
      "Aphis glycines",
      "Rhopalosiphum padi",
      "Rhopalosiphum maidis",
      "Macrosiphum euphorbiae",
      "Myzus persicae",
      "Acyrthosiphon pisum",
      "Capitophorus elaeagni",
      "Schizaphis graminum",
      "Aphis craccivora",
      "Sitobion avenae",
      "Therioaphis trifolii"
    )
  ) %>%
  droplevels()

# aphid list 2: PVY-relevant species
aphid <-
  aphid_long %>%
  filter(
    SpeciesName %in% c(
      "Myzus persicae",
      "Macrosiphum euphorbiae",
      "Aphis glycines",
      "Acyrthosiphon pisum",
      "Rhopalosiphum padi",
      "Capitophorus elaeagni",
      "Schizaphis graminum",
      "Aphis craccivora",
      "Sitobion avenae",
      "Rhopalosiphum maidis",
      "Brachycaudus helichrysi",
      "Brevicoryne brassicae",
      "Hayhurstia atriplicis",
      "Lipaphis pseudobrassicae"
    )
  ) %>%
  droplevels()

# shortlist
aphid <-
  aphid_long %>%
  filter(
    SpeciesName %in% c(
      "Aphis glycines",
      "Rhopalosiphum padi",
      "Rhopalosiphum maidis",
      "Therioaphis trifolii"
    )
  ) %>%
  droplevels()

# Filter dataset by top n species
fn.topsp <- function(df, n) {
  top <-
    df %>%
    group_by(SpeciesName) %>%
    summarise(totcount = sum(Count)) %>%
    arrange(desc(totcount)) %>%
    top_n(n, totcount)
  
  df %>%
    filter(SpeciesName %in% top$SpeciesName) %>%
    droplevels()
}

aphid <- fn.topsp(aphid_long, 9)

# Compute GDDs from PRISM data ---------------------------------------------------

# For handling weather data acquired from http://prism.oregonstate.edu/explorer/bulk.php
# Combine prism csv into single file
files = choose.files()
for (i in 1:length(files)){
  if (i ==1) {prism_in = NULL}
  prism_in <- rbind(prism_in, read.csv(files[i], skip = 10))
}
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

# create columns
prism_join = prism_in %>%
  mutate(Date = as.Date(Date),
         Year = format(Date, "%Y")) %>%
  group_by(SiteID, Year) %>%
  arrange(SiteID, Date) %>%
  mutate(GDD39 = cumsum(fn.gdd(tminF, tmaxF, 39, 86)),
         GDD50 = cumsum(fn.gdd(tminF, tmaxF, 50, 86))) %>%
  ungroup()

# remove existing data for 2017
prism <- prism %>%
  filter(Year != 2017)

# add new data to prism
prism <- rbind(prism, prism_join)

# fix column types
prism$SiteID <- as.factor(prism$SiteID)
prism$Year <- as.integer(prism$Year)
str(prism)


# Join PRISM GDDs with Aphids ----------------------------------------------
aphid$Date <- as.Date(aphid$Date)
aphid <-
  left_join(aphid,
            prism[, c("SiteID", "Date", "GDD39", "GDD50")],
            by = c("SiteID", "Date"))
aphid$SiteID <- as.factor(aphid$SiteID)

aphid$Date <- 
  as.Date(aphid$Date)

cbind(sort(levels(prism$SiteID)),sort(levels(aphid$SiteID)))

# Quick data summaries ----------------------------------------------------

# summarise by total count
aphid_long %>%
  group_by(SpeciesName) %>%
  summarise(TotalCount = sum(Count)) %>%
  arrange(desc(TotalCount))

# unique species per state
aphid_long %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  mutate(SpeciesName = as.character(SpeciesName)) %>%
  summarise(N_Taxa = length(unique(SpeciesName)))

# unique sites per state
aphid_long %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  mutate(SiteID = as.character(SiteID)) %>%
  summarise(N_Sites = length(unique(SiteID)))

### summaries of species diversity ###
aphid_joined <-
  aphid_long %>%
  filter(Count > 0) %>%
  left_join(aphid_sp, by = "SpeciesName") %>%
  mutate(SpeciesName = as.factor(SpeciesName))

f <- function(x) {length(unique(as.character(x)))} # returns number of unique factors

aphid_joined %>%
  group_by(State) %>%
  summarise(Sites = f(SiteID),
            Years = f(Year),
            Genera = f(Genus),
            Species = f(Species),
            TotCt = sum(Count),
            Samples = f(SampleID)
  ) %>%
  add_column(MeanCt = .$TotCt/.$Samples)


# summarise by frequency
aphid_long <-
  aphid_long %>%
  mutate(isFound = case_when(Count == 0 ~ 0, Count > 0 ~ 1))
aphid_long %>%
  group_by(SpeciesName) %>%
  summarise(Freq = sum(isFound)/n()) %>%
  arrange(desc(Freq))

# return number of observations
aphid %>%
  filter(Count > 0) %>%
  group_by(SpeciesName) %>%
  summarise(Obs = n())

# weeks per year
aphid_long %>%
  group_by(Year) %>%
  summarise(N_Weeks = length(unique(Week)),
            N_Captures = sum(Count))


# Quick data vizualization -------------------------------------------

# group by state
p1 <-
  aphid %>%
  group_by(State, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(
    x = reorder(SpeciesName, totcount),
    y = totcount,
    color = SpeciesName
  )) +
  facet_grid(~ State) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  coord_flip() +
  labs(
    aes(y = "Mean count",
        x = "",
        title = "NCR Aphid Suction Trap Network ('05-'17): Captures by state")
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
p1
pdf("NCR_CaptureByState.pdf", h = 8.5, w = 11)
p1
dev.off()

# group by year
p2 <-
  aphid %>%
  group_by(Year, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = as.factor(Year), y = totcount)) +
  facet_wrap(~ reorder(SpeciesName, -totcount)) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  labs(
    aes(y = "Mean count",
        x = "",
        title = "NCR Aphid Suction Trap Network ('05-'17): Captures by year")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")
p2
pdf("NCR_CaptureByYear.pdf", h = 8.5, w = 11)
p2
dev.off()

# group by state and year
p3 <- 
  aphid %>%
  group_by(State, Year, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = Year, y = totcount)) +
  facet_grid(State ~ SpeciesName, labeller = label_wrap_gen(10)) +
  geom_bar(stat = "identity", aes(fill = State)) +
  scale_y_sqrt() +
  labs(
    aes(y = "Mean count",
        x = "Year",
        title = "NCR Aphid Suction Trap Network ('05-'16): Top captures by year and state")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none",
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(angle = 0)
  )
p3
pdf("NCR_CaptureByYearAndState.pdf", h = 8.5, w = 11)
p3
dev.off()

# CM: YWL (multi-year by year + week + loc) ----------------------------------------------

# create variables
aphid$Location <- aphid$SiteID
aphid$YW <- interaction(aphid$Year, aphid$Week)
aphid$YW <- factor(aphid$YW)
aphid$YL <- interaction(aphid$Year, aphid$Location)
aphid$YL <- factor(aphid$YL)
aphid$YWL <- interaction(aphid$Year, aphid$Week, aphid$Location)
aphid$YWL <- factor(aphid$YWL)

# define function
rem.YWL <- function(df) {
  fmod <-
    glmer(
      Count ~ 1 +
        (1 | Year) +
        (1 | Week) +
        (1 | Location) +
        (1 | YW) +
        (1 | YL) +
        (1 | YWL),
      data = df,
      family = "poisson",
      control = glmerControl(optimizer ="Nelder_Mead")
    )
  data.frame(Week = as.numeric(rownames(ranef(fmod)$Week)),
             CMs = ranef(fmod)$Week[, 1])
}

# generate CMs
CM.YWL <-
  aphid %>%
  group_by(Species = SpeciesName) %>%
  do(rem.YWL(.))

# plot with smooth
p <-
  ggplot(CM.YWL, aes(x = Week, y = CMs)) +
  facet_wrap( ~ Species, scales = "free_y") +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x, k = 15),
    fill = "grey50",
    colour = "Black",
    size = 1.5
  ) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point() +
  coord_cartesian(xlim = c(0, 52)) +
  labs(aes(x = "Week",
           y = "CMs (Week)",
           title = "Calendar fits for NCR suction traps (2005-2017)"))
p

# save to pdf
pdf("CMsByWeek_NCR2005-2016.pdf", h = 8.5, w = 11)
p
dev.off()

# predict intercepts
YWLgam <- CM.YWL %>%
  group_by(Species) %>%
  do(
    data.frame(
      Week = 1:52,
      pred = predict.gam(gam(CMs ~ s(Week, k = 20), data = .))
    )
  ) %>%
  filter(sign(lag(pred)) != sign(pred) | sign(lead(pred)) != sign(pred))
YWLgam

YWLgam <- CM.YWL %>%
  group_by(Species) %>%
  do(
    data.frame(
      Week = 1:52,
      pred = predict.gam(gam(CMs ~ s(Week, k = 20), data = .))
    )
  ) %>%
  filter(sign(lag(pred)) != sign(pred) | sign(lead(pred)) != sign(pred))
YWLgam




# CM: WL (single year by week + loc) ---------------------------------------------

### define function ###
rem.WL <- function(df) {
  fmod <-
    glmer(
      Count ~ 1 +
        (1 | Week) +
        (1 | Location),
      data = df,
      family = "poisson",
      control = glmerControl(optimizer ="Nelder_Mead") #fix 'degenerate hessian' warning
    )
  data.frame(Week = as.numeric(rownames(ranef(fmod)$Week)),
             CMs = ranef(fmod)$Week[, 1])
}

# generate CMs
CM.WL <-
  aphid %>%
  group_by(Species = SpeciesName) %>%
  do(rem.WL(.))

# plot results
p <- ggplot(CM.WL, aes(x = Week, y = CMs)) +
  facet_wrap( ~ Species, scales = "free_y") +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = "grey50",
    colour = "Black",
    size = 1.5
  ) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point() +
  labs(aes(x = "Week",
           y = "CMs (Week)",
           title = "Calendar fits for Oregon Aphids (2016)"))
p

pdf("CMs_OR2016.pdf", h = 8.5, w = 11)
p
dev.off()

# prediction function
WLgam <- CM.WL %>%
  group_by(Species) %>%
  do(
    data.frame(
      Week = 1:52,
      pred = predict.gam(gam(CMs ~ s(Week, k = 20), data = .))
    )
  ) %>%
  filter(sign(lag(pred)) != sign(pred) | sign(lead(pred)) != sign(pred))
WLgam

# CMs by GDD ---------------------------------------------------------------

# define location variable
aphid$Location <- aphid$SiteID

# choose GDD model
aphid$GDD <- aphid$GDD39

# define function
rem.GDD <- function(df) {
  fmod <-
    glmer(
      Count ~ 1 +
        (1 | GDD) +
        (1 | Year) +
        (1 | SiteID), # location
      data = df,
      family = "poisson",
      control = glmerControl(optimizer = "Nelder_Mead") #fix 'degenerate hessian' warning
    )
  data.frame(GDD = as.numeric(rownames(ranef(fmod)$GDD)),
             CMs = ranef(fmod)$GDD[, 1])
}

# generate CMs from GDD
CM.GDD <-
  aphid %>%
  group_by(Species = SpeciesName) %>%
  do(rem.GDD(.))

# plot CMs
p <-
  ggplot(CM.GDD, aes(x = GDD, y = CMs)) +
  facet_wrap( ~ Species, scales = "free") +
  scale_x_sqrt() +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = "grey50",
    colour = "Black",
    size = 1.5) +
  geom_abline(intercept = 0, slope = 0) +
  labs(aes(x = "GDDs",
           y = "CMs (GDD39/86)",
           title = "GDD39/86 fits for NCR Suction Traps (2005-2017)"))
p

# export to pdf
pdf("NCR_CMsByGDD50.pdf", h = 8.5, w = 11)
p
dev.off()

# prediction function
GDDgam <- CM.GDD %>%
  group_by(Species) %>%
  do(
    data.frame(
      GDD = seq(10, 7500, by = 10),
      pred = predict.gam(gam(CMs ~ s(GDD, k = 20), data = .), GDD)
    )
  ) %>%
  filter(sign(lag(pred)) != sign(pred) | sign(lead(pred)) != sign(pred))
GDDgam


### prediction of intercepts ----------------------------------------
require(mgcv)

gdds <- seq(100, 3000, by = 10)

# incercepts by loess smoothing
weeklist <- 1:52
YWLpreds <- CM.YWL %>%
  group_by(Species) %>%
  do(data.frame(Week = weeklist, pred = predict(loess(CMs ~ Week, .), weeklist))) %>%
  filter(pred <= 0 & lag(pred) >= 0 | pred >= 0 & lead(pred) <= 0) 
YWLpreds

YWLpreds <- CM.YWL %>%
  group_by(Species) %>%
  do(
    data.frame(
      Week = weeklist,
      pred = predict(loess(CMs ~ Week, .), weeklist)
      )
    ) %>%
  filter(sign(lag(pred)) != sign(pred) | sign(lead(pred)) != sign(pred))
YWLpreds

YWLgam <- CM.YWL %>%
  group_by(Species) %>%
  do(
    data.frame(
      Week = weeklist,
      pred = predict.gam(gam(CMs ~ s(Week, k = 20), data = .))
      )
    ) %>%
  filter(sign(lag(pred)) != sign(pred) | sign(lead(pred)) != sign(pred))
YWLgam

# try gam method
YWLgam <- 
  data.frame(weeklist, predict(gam(CM.YWL, CMs ~ s(Week)), CM.YWL$Week))


smooth <- gam(CMs ~ s(Week), data = CM.ag, drop.unused.levels = FALSE)

predict.gam(smooth, 26, type = "response")
str(smooth)

### Plots -------------------------------------------------------------------

# counts by GDD; facets by species
a <- ggplot(aphid, aes(x = GDD, y = log10(Count + 1), color = as.factor(Year))) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = NA,
    size = .5
  ) +
  geom_abline(intercept = 0, slope = 0) +
  labs(aes(x = "Growing Degree Days (50F/86F)",
           y = "Log10 Count of Aphids",
           title = "Multi-year comparison of aphid phenologies, NCR Suction Trap data 2005-2016",
           legend = "Year"))
a
pdf("multiyear NCR aphid phenology by GDD.pdf", h = 8.5, w = 11)
a
dev.off()

# multi-year comparison of counts by week
b <- ggplot(aphid, aes(x = Week, y = log10(Count + 1), color = as.factor(Year))) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = NA,
    size = .5
  ) +
  geom_abline(intercept = 0, slope = 0) +
  labs(aes(x = "Week",
           y = "Log10 Count of Aphids",
           title = "Multi-year comparison of aphid phenologies, NCR Suction Trap data 2005-2016",
           legend = "Year"))
b
pdf("multiyear NCR aphid phenology by Week.pdf", h = 8.5, w = 11)
b
dev.off()

b <- ggplot(aphid, aes(x = Week, y = log10(Count + 1), color = Year)) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = "grey50",
    size = 1
  ) +
  stat_density() +
  geom_abline(intercept = 0, slope = 0) +
  labs(aes(x = "Week",
           y = "Log10 Count of Aphids",
           title = "Multi-year comparison of aphid phenologies, NCR Suction Trap data 2005-2016",
           legend = "Year"))
b

pdf("multiyear NCR aphid phenology by Week.pdf", h = 8.5, w = 11)
b
dev.off()



b <- ggplot(aphid, aes(x = as.factor(Year), y = log10(Count + 1), group = Year, color = as.factor(Year))) +
  geom_boxplot() +
  facet_grid(SpeciesName ~ .) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x) + 1,
    fill = "grey50",
    size = 1
  )
b

b <- ggplot(subset(aphid, Year == 2016), aes(x = Julian, y = log10(Count + 1))) +
  facet_grid(SpeciesName ~ .) +
  geom_ribbon(aes(ymin = 0, ymax = 2))
b




pdf("SpCountsByYear.pdf", h = 11, w = 8.5)
b
dev.off()



### KF random effects model ------------------------------------------------------------

require(lme4)
require(plyr) # only until update code for dplyr

### generate weeknums if not present ###
#aphid$WkoYr <- ceiling(aphid$Julian / 7)

### specify Location ###
aphid$Location <- aphid$SiteID

### generate interaction terms ###
aphid$YW <- interaction(aphid$Year, aphid$Week); aphid$YW <- factor(aphid$YW)
aphid$YL <- interaction(aphid$Year, aphid$Location); aphid$YL <- factor(aphid$YL)
aphid$YA <- interaction(aphid$Year, aphid$SpeciesName); aphid$YA <- factor(aphid$YA)
aphid$LA <- interaction(aphid$Location, aphid$SpeciesName); aphid$LA <- factor(aphid$LA)
aphid$YWL <- interaction(aphid$Year, aphid$Week, aphid$Location); aphid$YWL <- factor(aphid$YWL)
aphid$YWA <- interaction(aphid$Year, aphid$Week, aphid$SpeciesName); aphid$YWA <- factor(aphid$YWA)
aphid$YWLA <- interaction(aphid$Year, aphid$Week, aphid$Location, aphid$SpeciesName); aphid$YWLA <- factor(aphid$YWLA)

### run the model (takes a long time!) ###
flme1 <-
  glmer(
    Count ~ 1 +
      (1 | Year) +
      (1 | Week) +
      (1 | Location) +
      (1 | SpeciesName) +
      (1 | YW) +
      (1 | YL) +
      (1 | YA) +
      (1 | LA) +
      (1 | YWL) +
      (1 | YWA) +
      (1 | YWLA),
    data = aphid,
    family = "poisson",
    verbose = TRUE
  )

summary(flme1)

### plot model results ###
xx <- as.numeric(rownames(ranef(flme1)$Week))
yy <- ranef(flme1)$Week[, 1]
plot(xx, yy)


### total variance ###
totVAR <-
  sum((AS <- VarCorr(flme1)$Aphid.Species[, 1]),
      (AY <- VarCorr(flme1)$YA[, 1]),
      (AL <- VarCorr(flme1)$LA[, 1]),
      (Y <- VarCorr(flme1)$Year[, 1]),
      (W <-  VarCorr(flme1)$WkoYr[, 1]),
      (YW <- VarCorr(flme1)$YW[, 1]),
      (YL <- VarCorr(flme1)$YL[, 1]),
      (L <- VarCorr(flme1)$Location[, 1]),
      (YWL <- VarCorr(flme1)$YWL[, 1]),
      (YWA <- VarCorr(flme1)$YWA[, 1]),
      (YWLA <- VarCorr(flme1)$YWLA[, 1])
  )

vcs <- c(AS, AY, AL, Y, W, YW, YL, L, YWL, YWA, YWLA)

(vcs / totVAR) * 100


### relative variance ###
# Examine relative importance of sources of
# variability for individual aphid species

remodels2 <- function(df) {
  fmod <-
    glmer(
      Count ~ 1 +
        (1 | Year) +
        (1 | Week) +
        (1 | Location) +
        (1 | YW) +
        (1 | YL) +
        (1 | YWL),
      data = df,
      family = "poisson"
    )
  totVAR <-
    sum(c(
      (Y <- VarCorr(fmod)$Year[, 1]),
      (W <- VarCorr(fmod)$Week[, 1]),
      (L <- VarCorr(fmod)$Location[, 1]),
      (YW <- VarCorr(fmod)$YW[, 1]),
      (YL <- VarCorr(fmod)$YL[, 1]),
      (YWL <- VarCorr(fmod)$YWL[, 1])
    ))
  est <- c(Y, W, L, YW, YL, YWL)
  pest <- (est / totVAR) * 100
  varcs <-
    data.frame(round(est, digits = 3),
               round(pest, digits = 3),
               round(totVAR, digits = 3))
  varcs
}

vacs <- ddply(aphid, .(Aphid.Species), remodels2) # try to convert to dplyr code
names(vacs) <- c("Aphid", "Variance", "PropVar", "totVar")
vacs



# Making plots for individual aphid species spatial components ------------

library(gstat)
library(maps)

remodels4 <- function(df) {
  fmod <-
    glmer(
      Count ~ 1 + (1 |
                     Year) + (1 | Location) + (1 | YW) +  (1 | YL) + (1 | YWL),
      data = df,
      family = "poisson"
    )
  
  lBLUPs <- ranef(fmod)$Location[, 1]
  Location <- rownames(ranef(fmod)$Location)
  pred2 <- data.frame(Location, lBLUPs)
  pred2
  
}
turds2 <- ddply(aphid, .(Aphid.Species), remodels4)
names(turds2) <- c("Aphidspecies", "Location", "CMs")
turds2

Loc <- levels(turds2$Location)
pp2 <- c()
for (i in 1:length(Loc)) {
  pp <- subset(turds2, turds2$Location == Loc[i])
  ppLat <- subset(aphid$Latitude, aphid$Location == Loc[i])
  ppLon <- subset(aphid$Longitude, aphid$Location == Loc[i])
  pp$Lat <- ppLat[1]
  pp$Lon <- ppLon[1]
  pp2 <- rbind(pp2, pp)
}

midwest <-
  data.frame(map(
    'state',
    region = c(
      'wisconsin',
      'minn',
      'iowa',
      'illinois',
      'missouri',
      'indiana',
      'kentucky',
      'michigan',
      'south dakota',
      'kansas'
    ),
    plot = FALSE
  )[c("x", "y")])

(
  mwmap <-
    qplot(
      x,
      y,
      data = midwest,
      geom = "path",
      xlab = "Longitude",
      ylab = "Latitude"
    ) +
    geom_point(aes(
      -Lon, Lat, size = exp(CMs), colour = exp(CMs)
    ), data = pp2) +
    facet_wrap( ~ Aphidspecies)
)
pdf("rlgmap_prototype.pdf")
mwmap
dev.off()


fmod <- gam(Count ~ factor(Year) + s(WkoYr, by = Location),
            data = aphid[aphid$Aphid.Species == "Aphis glycines", ],
            family = "poisson")


# older version??? --------------------------------------------------------

library(gstat)
library(maps)

midwest <-
  data.frame(map(
    'state',
    region = c(
      'wisconsin',
      'minn',
      'iowa',
      'illinois',
      'missouri',
      'indiana',
      'kentucky',
      'michigan',
      'south dakota',
      'kansas'
    ),
    plot = FALSE
  )[c("x", "y")])

yrs <- as.numeric(levels(as.factor(aphid$Year)))
wks <- as.numeric(levels(as.factor(aphid$YroWk)))


for (j in 18:44) {
  pdf(paste(j, ".pdf", sep = ""))
  aphsubs <- subset(aphid, aphid$Year == 2007 & aphid$YroWk == j)
  
  print(
    mwplot <-
      qplot(
        x,
        y,
        data = midwest,
        geom = "path",
        xlab = "Longitude",
        ylab = "Latitude"
      ) +
      geom_point(aes(
        -Longitude, Latitude, size = Count, colour = Count
      ), data = aphsubs) +
      facet_wrap( ~ Aphid.Species)
  )
  
  dev.off()
}



# Examining correlations among aphid species ------------------------------

flme1c <- glmer(
  Count ~ 0 + Aphid.Species + (0 + Aphid.Species | Year),
  data = aphid,
  family = "poisson",
  verbose = TRUE
)

flme2c <- glmer(
  Count ~ 0 + Location + (0 + Location | Aphid.Species),
  data = aphid,
  family = "poisson",
  verbose = TRUE
)
