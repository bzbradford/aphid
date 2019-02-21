# Load packages ----

library(tidyverse)
library(Cairo)



# Read in suction data ----

# read aphid counts (including dummy counts)
aphid_in =
  read.csv(file.choose(), header = TRUE, na = c('','.')) %>%
  mutate(Date = as.Date(Date))
aphid_in = filter(aphid_in, SpeciesName != 'Phylloxeridae')


# site info
aphid_sites = read.csv("data/aphid_sites.csv", header = T, na = '')


# read species names
aphid_spp = read.csv(file.choose(), header = TRUE, na = c('','.'))



# Expand dataset ----

# Expand dataset with zeros, drop dummy variable
expandFn = function(df) {
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
aphid_full = expandFn(aphid_in)



# SCRI data ----
SCRI_ID = read.csv(file.choose(), header = TRUE, na = c('', '.')) %>%
  mutate(Date = as.Date(Date))
SCRI_ID_full = 
  expandFn(SCRI_ID) %>%
  left_join(aphid_spp,
            by = "SpeciesName") %>%
  mutate(SpeciesName = as.factor(SpeciesName))
SCRI_ID_full %>% write.csv("out/scri_id_aphids.csv")


# Joining datasets ----

# Join aphid taxonomy
aphid_full =
  aphid_full %>%
  left_join(aphid_spp,
            by = "SpeciesName") %>%
  mutate(SpeciesName = as.factor(SpeciesName))


# Join GDD data
aphid_full =
  aphid_full %>%
  left_join(prism[, c("SiteID", "Date", "GDD39", "GDD50")],
            by = c("SiteID", "Date")) %>%
  mutate(SiteID = as.factor(SiteID))


# Join site locations
aphid_full =
  aphid_full %>%
  left_join(aphid_sites[, c('SiteID', 'Lat', 'Lon')],
            by = 'SiteID') %>%
  mutate(SiteID = as.factor(SiteID))


# add month
aphid_full =
  aphid_full %>%
  mutate(Month = as.numeric(format.Date(Date, format = "%m")))


# Export file (optional)
aphid_full %>%
  arrange(SiteName, Date, SpeciesName) %>%
  write.csv("data/aphid_full.csv", na = '')



# Data subset ----

# Species list from Frost code
aphid =
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
aphid =
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
topSpFn = function(df, n) {
  top =
    df %>%
    group_by(SpeciesName) %>%
    summarise(totcount = sum(Count)) %>%
    arrange(desc(totcount)) %>%
    top_n(n, totcount)
  df %>%
    filter(SpeciesName %in% top$SpeciesName) %>%
    droplevels()
}

# Select top spp by count in entire dataset
aphid = topSpFn(aphid_full, 10)

# Select top spp in wisconsin
aphid = aphid_full %>%
  filter(State == "WI") %>%
  topSpFn(10)

# Filter Wisconsin only
aphid = aphid %>% 
  filter(State == "WI")



# Numerical summaries ----

# summarise by total count
aphid %>%
  group_by(SpeciesName) %>%
  summarise(TotalCount = sum(Count), n = n()) %>%
  arrange(desc(TotalCount)) %>%
  write.csv("out/totalCount.csv")


# generate annual count summaries and a total count by species, then save
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
}) %>%
  write.csv("out/totalCount.csv")


# generate weekly mean counts for wi
aphid %>%
  filter(State == "WI") %>%
  select(SpeciesName, Week, Count) %>%
  group_by(SpeciesName, Week) %>%
  summarise(MeanCount = mean(Count)) %>%
  write.csv("out/wi_wkly.csv")


# unique species per state
aphid_full %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  mutate(SpeciesName = as.character(SpeciesName)) %>%
  summarise(N_Taxa = length(unique(SpeciesName)))


# unique sites per state
aphid_full %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  summarise(N_Sites = length(levels(droplevels(SiteID))))



# Species diversity ----

f = function(x) {length(unique(as.character(x)))} # returns number of unique factors
aphid_summary =
  aphid_full %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  summarise(
    Sites = f(SiteID),
    Years = f(Year),
    Samples = f(SampleID),
    Families = f(Family),
    Genera = f(Genus),
    Species = f(Species),
    TotCt = sum(Count)
  ) %>%
  add_column(MeanCt = .$TotCt / .$Samples)
aphid_summary %>% write.csv("out/STN_diversity_allspp.csv")


# species diversity counts for whole dataset
aphid_full %>%
  filter(Count > 0) %>%
  summarise(
    States = f(State),
    Sites = f(SiteID),
    Years = f(Year),
    Samples = f(SampleID),
    Families = f(Family),
    Genera = f(Genus),
    Species = f(Species),
    TotCt = sum(Count)
  ) %>%
  add_column(MeanCt = .$TotCt / .$Samples)


# summarise by frequency
aphid_full %>%
  mutate(isFound = case_when(Count == 0 ~ 0, Count > 0 ~ 1)) %>%
  group_by(SpeciesName) %>%
  summarise(Freq = sum(isFound)/n()) %>%
  arrange(desc(Freq))


# return number of observations by species name
aphid %>%
  filter(Count > 0) %>%
  group_by(SpeciesName) %>%
  summarise(Obs = n(), Count = sum(Count))


# number of observations by full taxonomy
aphid_full %>%
  filter(Count > 0) %>%
  group_by(Family, Subfamily, Genus, Species, SpeciesName) %>%
  summarise(Obs = n(), Count = sum(Count)) %>%
  write.csv("out/spp_counts.csv")


# weeks per year
aphid_full %>%
  group_by(Year) %>%
  summarise(N_Weeks = length(unique(Week)),
            N_Captures = sum(Count))



# Plot counts by state ----
pltByState =
  aphid %>%
  group_by(State, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(
    x = reorder(SpeciesName, totcount),
    y = totcount
  )) +
  facet_grid( ~ State) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  coord_flip() +
  labs(
    x = "",
    y = "Mean Count",
    title = paste(
      "Aphid Suction Trap Network (",
      min(aphid$Year),
      "-",
      max(aphid$Year),
      "): Captures by State",
      sep = ''
    )
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
pltByState
pdf("out/STN_CountsByState.pdf", h = 8.5, w = 11); pltByState; dev.off()


# Plot counts by year ----
pltByYear =
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


# Plot counts by state and year ----
p = 
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
       title = paste0("Aphid Suction Trap Network (",
                      min(aphid$Year),"-",max(aphid$Year),
                     "): Captures by Year and State")) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(angle = 0, size = 12, face = "bold")
  )
p
CairoPNG("out/STN_CaptureByYearAndState.png", w = 1200, h = 900); p; dev.off()


# Plot counts by GDD; facets by species ----
p = aphid %>%
  ggplot(aes(
    x = GDD,
    y = log10(Count + 1),
    color = as.factor(Year)
  )) +
  facet_wrap(~ SpeciesName, scales = "free") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = NA,
    size = .5
  ) +
  geom_abline(intercept = 0, slope = 0) +
  labs(
    x = "Growing Degree Days (50F/86F)",
    y = "Log10 Count of Aphids",
    title = "Multi-year comparison of aphid phenologies, NCR Suction Trap data 2005-2017",
    legend = "Year"
  )
p
pdf("multiyear NCR aphid phenology by GDD.pdf", h = 8.5, w = 11); p; dev.off()


# Plot multi-year comparison of counts by week ----
p = aphid %>%
  ggplot(aes(x = Week, y = log10(Count + 1), color = as.factor(Year))) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  stat_smooth(method = "gam", formula = y ~ s(x), fill = NA, size = .5) +
  geom_abline(intercept = 0, slope = 0) +
  labs(aes(x = "Week",
           y = "Log10 Count of Aphids",
           title = "Multi-year comparison of aphid phenologies, NCR Suction Trap data 2005-2017",
           legend = "Year"))
p
pdf("multiyear NCR aphid phenology by Week.pdf", h = 8.5, w = 11); p; dev.off()


p = aphid %>%
  ggplot(aes(x = as.factor(Year), y = log10(Count + 1), group = Year, color = as.factor(Year))) +
  geom_boxplot() +
  facet_grid(SpeciesName ~ .) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  geom_smooth( method = "gam", formula = y ~ s(x) + 1, fill = "grey50", size = 1)
p

p = subset(aphid, Year == 2016) %>%
  ggplot(aes(x = Julian, y = log10(Count + 1))) +
  facet_grid(SpeciesName ~ .) +
  geom_ribbon(aes(ymin = 0, ymax = 2))
p



# Russ charts Jan 2019 ----

# Top 5 wis aphids
russ_aphids = c(
  "Aphis glycines",
  "Acyrthosiphon pisum",
  "Rhopalosiphum padi",
  "Rhopalosiphum maidis",
  "Myzus persicae"
  )

#generate dataset
aphid =
  aphid_full %>%
  filter(SpeciesName %in% russ_aphids) %>%
  droplevels()


# all species same graph
p =
  aphid %>%
  group_by(Year, Week, SpeciesName) %>%
  summarize(TotCount = log(sum(Count)+1)) %>%
  ggplot(aes(x = as.Date("2017-01-01")+(Week-1)*7, y = TotCount)) +
  geom_area(aes(fill = SpeciesName)) +
  facet_grid(Year ~ SpeciesName, scales = "free_x") +
  scale_x_date(date_labels = "%b") +
  labs(x = "Week of Year",
       y = "Log abundance") +
  theme(panel.spacing = unit(.1, "lines"),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
p
CairoPNG("out/russ aphids same graph.png", w=1200, h=900); p; dev.off()


# all species same graph
p =
  aphid %>%
  group_by(Year, Week, SpeciesName) %>%
  summarize(TotCount = mean(Count)) %>%
  ggplot(aes(x = as.Date("2017-01-01")+(Week-1)*7, y = TotCount)) +
  geom_area(aes(fill = SpeciesName)) +
  scale_x_date(date_labels = "%b") +
  labs(x = "",
       y = "Mean Count") +
  theme(panel.spacing = unit(.1, "lines"),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.position = "none")
p
p + scale_y_sqrt()

CairoPNG("out/mean_counts_freey_sqrty.png", w=1200, h=900)
p + scale_y_sqrt() + facet_grid(Year ~ SpeciesName, scales = "free_y")
dev.off()

CairoPNG("out/mean_counts_freey_plainy.png", w=1200, h=900)
p + facet_grid(Year ~ SpeciesName, scales = "free_y")
dev.off()

CairoPNG("out/mean_counts_samey_sqrty.png", w=1200, h=900)
p + scale_y_sqrt() + facet_grid(Year ~ SpeciesName)
dev.off()

CairoPNG("out/mean_counts_samey_plainy.png", w=1200, h=900)
p + facet_grid(Year ~ SpeciesName)
dev.off()


# individual graphs for each species
for (s in russ_aphids) {
  print(s)
  p =
  aphid %>%
    filter(SpeciesName == s) %>%
    group_by(Year, Week, SpeciesName) %>%
    summarize(TotCount = log(sum(Count) + 1)) %>%
    ggplot(aes(x = as.Date("2017-01-01") + (Week - 1) * 7, y = TotCount)) +
    geom_area(aes(fill = SpeciesName)) +
    facet_grid(Year ~ SpeciesName, scales = "free_x") +
    scale_x_date(date_labels = "%b") +
    labs(x = "Week of Year",
         y = "Log abundance") +
    theme(
      panel.spacing = unit(.1, "lines"),
      strip.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
  ggsave(p, file=paste0("./out/Abund ", s, ".png"), width = 8, height = 6)
}



# glycines IDW export for Emily ----


# filter
glycines = aphid_full %>%
  filter(SpeciesName == "Aphis glycines") %>%
  mutate(MonthName = as.factor(format.Date(Date, "%b"))) %>%
  droplevels()


# export
glycines %>% write.csv('out/glycines.csv')


# mean glycines counts by site and year for Jul & Aug
glycines_jul_aug =
  glycines %>%
  filter(Month %in% c(7, 8)) %>%
  group_by(Year, SiteID) %>%
  summarize(MeanCount = mean(Count))
glycines_jul_aug %>% write.csv("out/glycines-jul-aug.csv")


# mean glycines counts by site and year for Sep & Oct
glycines_sep_oct =
  glycines %>%
  filter(Month %in% c(9, 10)) %>%
  group_by(Year, SiteID) %>%
  summarize(MeanCount = mean(Count))
glycines_sep_oct %>% write.csv("out/glycines-sep-oct.csv")


# mean glycines counts by site and year for Jul thru Oct
glycines_jul_oct =
  glycines %>%
  filter(Month %in% c(7, 8, 9, 10)) %>%
  group_by(Year, SiteID) %>%
  summarize(MeanCount = mean(Count), Lon = min(Lon), Lat = min(Lat)) %>%
  ungroup()
glycines_jul_oct %>% write.csv("out/glycines-jul-oct.csv")
