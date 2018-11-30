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
