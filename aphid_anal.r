# Load and prepare data ------------------------------------------------

# required packages
library(lme4)
library(mgcv)
library(ggplot2)
library(tidyverse)

# read count (including dummy counts)
aphidin <-
  read.csv("aphid.csv", header = TRUE, na = c('','.')) %>%
  mutate(Date = as.Date(Date))

# read species names
aphidsp <- read.csv("aphidsp.csv", header = TRUE, na = c('','.'))

# read prism data
prism <-
  read.csv("prism.csv", header = TRUE) %>%
  mutate(Date = as.Date(Date))


### Expand dataset with zeros ###
aph.expand <- function(df) {
  cols <- ncol(df)
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
aphidlong <- aph.expand(aphidin)


### Join PRISM GDDs with aphid data ###
aphidlong <-
  aphidlong %>%
  left_join(prism[, c("SiteID", "Date", "GDD39", "GDD50")],
            by = c("SiteID", "Date")) %>%
  mutate(SiteID = as.factor(SiteID))

### Join aphid names to dataset ###
aphidlong <-
  aphidlong %>%
  left_join(aphidsp, by = "SpeciesName") %>%
  mutate(SpeciesName = as.factor(SpeciesName))


### optional file export ###
aphidlong %>%
  filter(Count > 0) %>%
  arrange(SiteName, Date) %>%
  write.csv("aphidlong.csv", na = "")



# Aphid data subset options -------------------------------------

# Species list 1 (from Frost code)
aphid <-
  aphidlong %>%
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

# Species list 2: PVY-relevant species
aphid <-
  aphidlong %>%
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

# Just four species
aphid <-
  aphidlong %>%
  filter(
    SpeciesName %in% c(
      "Aphis glycines",
      "Rhopalosiphum padi",
      "Rhopalosiphum maidis",
      "Therioaphis trifolii"
    )
  ) %>%
  droplevels()

# Select top n species
aph.topsp <- function(df, n) {
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
aphid <- aph.topsp(aphidlong, 9)



# Quick data summaries of subset ------------------------------------------------

# summarise by total count
aphid %>%
  group_by(SpeciesName) %>%
  summarise(TotalCount = sum(Count), n = n()) %>%
  arrange(desc(TotalCount))

# unique species per state
aphidlong %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  mutate(SpeciesName = as.character(SpeciesName)) %>%
  summarise(N_Taxa = length(unique(SpeciesName)))

# unique sites per state
aphidlong %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  mutate(SiteID = as.character(SiteID)) %>%
  summarise(N_Sites = length(unique(SiteID)))


# Species diversity from full aphid dataset -----------------------------

f <- function(x) {length(unique(as.character(x)))} # returns number of unique factors
aphidlong %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  summarise(Sites = f(SiteID),
            Years = f(Year),
            Samples = f(SampleID),
            Families = f(Family),
            Genera = f(Genus),
            Species = f(Species),
            TotCt = sum(Count)
  ) %>%
  add_column(MeanCt = .$TotCt/.$Samples)

# summarise by frequency
aphidlong %>%
  mutate(isFound = case_when(Count == 0 ~ 0, Count > 0 ~ 1)) %>%
  group_by(SpeciesName) %>%
  summarise(Freq = sum(isFound)/n()) %>%
  arrange(desc(Freq))

# return number of observations
aphid %>%
  filter(Count > 0) %>%
  group_by(SpeciesName) %>%
  summarise(Obs = n())

# weeks per year
aphidlong %>%
  group_by(Year) %>%
  summarise(N_Weeks = length(unique(Week)),
            N_Captures = sum(Count))



# Quick data vizualization -------------------------------------------

### group by state ###
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
  labs(x = "",
       y = "Mean Count",
       title = paste("Aphid Suction Trap Network (",
                     min(aphid$Year),
                     "-",
                     max(aphid$Year),
                     "): Captures by State",
                     sep = '')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
p1

# write pdf
pdf("out/STN_CountsByState.pdf", h = 8.5, w = 11); p1; dev.off()

### Counts by year ###
p2 <-
  aphid %>%
  group_by(Year, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = as.factor(Year), y = totcount)) +
  facet_wrap(~ reorder(SpeciesName, -totcount)) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  labs(x = "",
       y = "Mean Count",
       title = paste("Aphid Suction Trap Network (",
                     min(aphid$Year),
                     "-",
                     max(aphid$Year),
                     "): Captures by Year",
                     sep = '')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none")
p2

# Write pdf
pdf("out/STN_CountsByYear.pdf", h = 8.5, w = 11); p2; dev.off()


### group by state and year ###
p3 <- 
  aphid %>%
  group_by(State, Year, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = Year, y = totcount)) +
  facet_grid(State ~ SpeciesName, labeller = label_wrap_gen(10)) +
  geom_bar(stat = "identity", aes(fill = State)) +
  scale_y_sqrt() +
  labs(x = "",
       y = "Mean Count",
       title = paste("Aphid Suction Trap Network (",
                     min(aphid$Year),
                     "-",
                     max(aphid$Year),
                     "): Captures by Year and State",
                     sep = '')) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none",
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(angle = 0)
  )
p3

# Write pdf
pdf("out/STN_CaptureByYearAndState.pdf", h = 8.5, w = 11); p3; dev.off()



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
remYWL <- function(df) {
  require(lme4)
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
  data.frame(Week = as.numeric(rownames(ranef(fmod)$Week)),
             CMs = ranef(fmod)$Week[, 1])
}

# generate CMs
CMYWL <-
  aphid %>%
  group_by(Species = SpeciesName) %>%
  do(remYWL(.))

# plot with smooth
p <-
  ggplot(CMYWL, aes(x = Week, y = CMs)) +
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

# Write to pdf
pdf("CMsByWeek_NCR2005-2016.pdf", h = 8.5, w = 11); p; dev.off()

# predict intercepts
YWLgam <- CMYWL %>%
  group_by(Species) %>%
  do(
    data.frame(
      Week = 1:52,
      CM = predict.gam(gam(CMs ~ s(Week, k = 20), data = .))
    )
  ) %>%
  filter(sign(lag(CM)) != sign(CM) | sign(lead(CM)) != sign(CM))
YWLgam



# CM: WL (single year by week + loc) ---------------------------------------------

# define function
remWL <- function(df) {
  require(lme4)
  fmod <-
    glmer(
      Count ~ 1 +
        (1 | Week) +
        (1 | Location),
      data = df,
      family = "poisson"
    )
  data.frame(Week = as.numeric(rownames(ranef(fmod)$Week)),
             CMs = ranef(fmod)$Week[, 1])
}

# generate CMs
CMWL <-
  aphid %>%
  group_by(Species = SpeciesName) %>%
  do(rem.WL(.))

# plot results
p <- ggplot(CMWL, aes(x = Week, y = CMs)) +
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
WLgam <- CMWL %>%
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

### define functions ###
remfn <- function(df) {
  require(lme4)
  fmod <-
    glmer(
      Count ~ 1 +
        (1 | GDD39) +
        (1 | Year) +
        (1 | SiteID), # location
      data = df,
      family = "poisson"
    )
  data.frame(GDD = as.numeric(rownames(ranef(fmod)$GDD39)),
             CM = ranef(fmod)$GDD39[, 1])
}
gampredfn <- function(df){
  df %>%
    group_by(SpeciesName) %>%
    do(
      data.frame(
        GDD = 1:7500,
        CM = predict.gam(gam(CM ~ s(GDD), data = .),
                         data.frame(GDD = 1:7500))
      )
    )
}
gamptfn <- function(df){
  df %>%
    mutate(
      type = case_when(
        sign(lag(CM)) < sign(CM) ~ "RI",         # rising intercept
        sign(lag(CM)) > sign(CM) ~ "FI",        # falling intercept
        lag(CM) < CM & CM > lead(CM) ~ "max", # local maxima
        lag(CM) > CM & CM < lead(CM) ~ "min"  # local minima
      )
    ) %>%
    na.omit()
}

# generate CMs from GDD
aphid %>%
  group_by(SpeciesName) %>%
  do(remfn(.)) ->
  CMGDD

# generate predictions from gam fit
CMGDD %>%
  do(gampredfn(.)) ->
  gampreds

# identify critical points
gampreds %>%
  do(gamptfn(.)) ->
  gampts

# plot phenology curves
p <-
  ggplot(CMGDD, aes(x = GDD, y = CM)) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  scale_x_sqrt() +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = "grey50",
    colour = "Black",
    size = 1.5
  ) +
  geom_abline(intercept = 0, slope = 0) +
  labs(x = "GDDs",
       y = "CMs (GDD39/86)",
       title = paste("GDD39/86 fits for Suction Traps",
                     min(aphid$Year), "-", max(aphid$Year))
  )
p

# add and label points of interest from gam prediction
p2 <- p +
  geom_line(data = gampreds, color = "red", size = 1.5) +
  geom_point(data = gampts, aes(x = GDD, y = CM)) +
  geom_label(data = gampts, aes(x = GDD, y = CM, label = GDD))
p2

# write pdf
pdf("out/STN_CMsByGDD.pdf", h = 8.5, w = 11); p2; dev.off()



# GDD cross-validation ----------------------------------------------------

# Cross-validation subsetting #
CV.aphid <- list()
CV.aphid[[1]] <- aphid %>% filter(Year < 2014) # First two-thirds
CV.aphid[[2]] <- aphid %>% filter(Year > 2008) # Recent two-thirds
CV.aphid[[3]] <- aphid %>% filter(Year < 2009 | Year > 2013) # Outside two-thirds
CV.aphid[[4]] <- aphid %>% filter(Year > 2006 & Year < 2016) # middle two-thirds

# generate CMs from each subset
CV.CMGDD <- list()
for (i in 1:4) {
  CV.aphid[[i]] %>%
    group_by(SpeciesName) %>%
    do(remfn(.)) ->
    CV.CMGDD[[i]]
}

# compute gam predictions from each subset
CV.gampreds <- list()
for (i in 1:4) {
  CV.CMGDD[[i]] %>%
    do(gampredfn(.)) ->
    CV.gampreds[[i]]
}

# identify critical points on curve from each subset
CV.gampts <- list()
for (i in 1:4) {
  CV.gampreds[[i]] %>%
    do(gamptfn(.)) ->
    CV.gampts[[i]]
}

# cv plots
pdf("out/STN_CMsByGDD_CV.pdf", h = 8.5, w = 11)
p <- list()
for (i in 1:4){
  ggplot(CV.CMGDD[[i]], aes(x = GDD, y = CM)) +
    facet_wrap(~ SpeciesName) +
    scale_x_sqrt() +
    scale_y_continuous(limits = c(-2, 2)) +
    geom_abline(intercept = 0, slope = 0) +
    labs(
      x = "GDDs",
      y = "CMs (GDD39/86)",
      title = paste("GDD39/86 fits for Suction Traps (CV block ", i, ")", sep = "")
    ) +
    geom_line(data = CV.gampreds[[i]], color = "red", size = 1.5) +
    geom_point(data = CV.gampts[[i]], aes(x = GDD, y = CM)) +
    geom_label(data = CV.gampts[[i]], aes(x = GDD, y = CM, label = GDD)) ->
    p[[i]]
  print(p[[i]])
}
dev.off()



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
