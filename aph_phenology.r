###
### Conditional mode plots of aphid phenology
###

library(tidyverse)
library(ggplot2)
library(lme4)

# CM: Year + Week + Loc ----------------------------------------------

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



# CM: Week + Loc (single year) ---------------------------------------------

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



# CM by GDD ---------------------------------------------------------------

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

p3 <-
  ggplot(CMGDD, aes(x = GDD, y = CM)) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  scale_x_sqrt() +
  geom_abline(intercept = 0, slope = 0) +
  geom_line(data = gampreds, color = "red", size = 1.5) +
  geom_point(data = gampts, aes(x = GDD, y = CM)) +
  geom_label(data = gampts, aes(x = GDD, y = CM, label = GDD)) +
  labs(x = "GDDs",
       y = "CMs (GDD39/86)",
       title = paste("GDD39/86 fits for Suction Traps (Full dataset)",
                     min(aphid$Year), "-", max(aphid$Year))
  )
p3


# CM GDD cross-validation ----------------------------------------------------

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


### generate crossvalidation plots ###


plt_full <-
  ggplot(CMGDD, aes(x = GDD, y = CM)) +
  facet_wrap(~ SpeciesName) +
  scale_x_sqrt() +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_line(data = gampreds, color = "red", size = 1.5) +
  geom_point(data = gampts, aes(x = GDD, y = CM)) +
  geom_label(data = gampts, aes(x = GDD, y = CM, label = GDD)) +
  labs(
    x = "GDDs",
    y = "CMs (GDD39/86)",
    title = paste(
      "GDD39/86 fits for Suction Traps (Full dataset)",
      min(aphid$Year),
      "-",
      max(aphid$Year)
    )
  )
plt_full

plt_cvs <- list()
for (i in 1:4) {
  plt_cvs[[i]] <-
    ggplot(CV.CMGDD[[i]], aes(x = GDD, y = CM)) +
    facet_wrap( ~ SpeciesName) +
    scale_x_sqrt() +
    scale_y_continuous(limits = c(-2, 2)) +
    geom_abline(intercept = 0, slope = 0) +
    geom_line(data = CV.gampreds[[i]], color = "red", size = 1.5) +
    geom_point(data = CV.gampts[[i]], aes(x = GDD, y = CM)) +
    geom_label(data = CV.gampts[[i]], aes(x = GDD, y = CM, label = GDD)) +
    labs(
      x = "GDDs",
      y = "CMs (GDD39/86)",
      title = paste(
        "GDD39/86 fits for Suction Traps (Crossvalidation block ", i, ")", sep = "")
    )
}

plt_overlay <-
  ggplot(CMGDD, aes(x = GDD, y = CM)) +
  facet_wrap(~ SpeciesName) +
  scale_x_sqrt() +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_line(data = gampreds, color = "red", size = 1.5) +
  geom_line(data = CV.gampreds[[1]], color = "black", size = .75, alpha = 0.5) +
  geom_line(data = CV.gampreds[[2]], color = "black", size = .75, alpha = 0.5) +
  geom_line(data = CV.gampreds[[3]], color = "black", size = .75, alpha = 0.5) +
  geom_line(data = CV.gampreds[[4]], color = "black", size = .75, alpha = 0.5) +
  labs(
    x = "GDDs",
    y = "CMs (GDD39/86)",
    title = paste(
      "GDD39/86 fits for Suction Traps (Full dataset & CV subsets)",
      min(aphid$Year),
      "-",
      max(aphid$Year)
    )
  )
plt_overlay


# generate pdf
pdf("out/STN_CMsByGDD_CV.pdf", h = 8.5, w = 11)
plt_full
plt_cvs
plt_overlay
dev.off()
