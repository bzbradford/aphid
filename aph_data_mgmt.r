library(tidyverse)

# Read in data ------------------------------------------------------------

# read aphid counts (including dummy counts)
aphid_in <-
  read.csv(file.choose(), header = TRUE, na = c('','.')) %>%
  mutate(Date = as.Date(Date))

# read species names
aphid_spp <- read.csv("data/aphidsp.csv", header = TRUE, na = c('','.'))

# read prism data
prism <-
  read.csv("data/prism.csv", header = TRUE) %>%
  mutate(Date = as.Date(Date))



# Data preparation ---------------------------------------------------------

# Expand dataset with zeros, drop dummy variable #
expandfn <- function(df) {
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
aphid_full <- expandfn(aphid_in)

# Join GDD data with aphid data
aphid_full <-
  aphid_full %>%
  left_join(prism[, c("SiteID", "Date", "GDD39", "GDD50")],
            by = c("SiteID", "Date")) %>%
  mutate(SiteID = as.factor(SiteID))

# Join aphid names to dataset
aphid_full <-
  aphid_full %>%
  left_join(aphid_spp, by = "SpeciesName") %>%
  mutate(SpeciesName = as.factor(SpeciesName))

# Export file (optional)
aphid_full %>%
  filter(Count > 0) %>%
  arrange(SiteName, Date) %>%
  write.csv("aphidlong.csv", na = "")



# Aphid data subset options -------------------------------------

# Species list from Frost code
aphid <-
  aphid_full %>%
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


# Just four species
aphid <-
  aphid_full %>%
  filter(
    SpeciesName %in% c(
      "Aphis glycines",
      "Rhopalosiphum padi",
      "Rhopalosiphum maidis",
      "Therioaphis trifolii"
    )
  ) %>%
  droplevels()


# Select by top n species
topSpFn <- function(df, n) {
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
aphid <- topSpFn(aphid_full, 10)

# Top spp in wisconsin
aphid <- aphid_full %>%
  filter(State == "WI") %>%
  topSpFn(10)

# Wisconsin only
aphid <- aphid %>% 
  filter(State == "WI")


# Numerical summaries of subset ------------------------------------------------

# summarise by total count
aphid %>%
  group_by(SpeciesName) %>%
  summarise(TotalCount = sum(Count), n = n()) %>%
  arrange(desc(TotalCount)) %>%
  write.csv("out/totalCount.csv")


# generate annual count summaries and a total count by species, then save
a <- aphid %>%
  select(SpeciesName, Year, Count) %>%
  group_by(SpeciesName, Year) %>%
  summarise(TotalCount = sum(Count)) %>%
  ungroup(Year) %>%
  spread(Year, TotalCount)
b <- aphid %>%
  group_by(SpeciesName) %>%
  summarise(TotalCount = sum(Count))
left_join(a, b) %>%
  write.csv("out/totalCount.csv")


# piped method of the above
totCt <-
  left_join({
    aphid %>%
      select(SpeciesName, Year, Count) %>%
      group_by(SpeciesName, Year) %>%
      summarise(TotalCount = sum(Count)) %>%
      ungroup(Year) %>%
      spread(Year, TotalCount)
  },
  {
    aphid %>%
      group_by(SpeciesName) %>%
      summarise(TotalCount = sum(Count))
  })
totCt %>% write.csv("out/totalCount.csv")

# generate weekly mean counts
aphid %>%
  filter(State == "WI") %>%
  select(SpeciesName, Week, Count) %>%
  group_by(SpeciesName, Week) %>%
  summarise(MeanCount = mean(Count)) %>%
  write.csv("out/wi_wkly.csv")



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



# Summary plots -------------------------------------------
library(ggplot2)

### group by state ###
pltByState <-
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
pltByState
pdf("out/STN_CountsByState.pdf", h = 8.5, w = 11); pltByState; dev.off()


### Counts by year ###
pltByYear <-
  aphid %>% filter(State == "WI") %>%
  group_by(Year, SpeciesName) %>%
  summarise(totcount = sum(Count)) %>%
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
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")
pltByYear
pdf("out/STN_CountsByYear.pdf", h = 8.5, w = 11); pltByYear; dev.off()


### Counts by state and year ###
pltByStateYear <- 
  aphid %>%
  group_by(State, Year, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = Year, y = totcount)) +
  facet_grid(State ~ SpeciesName, labeller = label_wrap_gen(10)) +
  geom_bar(stat = "identity", aes(fill = State)) +
  scale_y_sqrt() +
  labs(x = "",
       y = "Mean Count/Week",
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
pltByStateYear
pdf("out/STN_CaptureByYearAndState.pdf", h = 8.5, w = 11); pltByStateYear; dev.off()

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



