
# libraries ---------------------------------------------------------------

library(tidyverse)




# prism handler -----------------------------------------------------------

# load prism csvs
prism_in <- list.files("./prism", pattern = glob2rx("PRISM_*.csv"), full.names = T) %>%
  map_df(~read_csv(., skip = 10))


# define GDD function
gdd_fn <- function(tmin, tmax, lower, upper) {
  pmax(0, (pmax(tmin, lower) + pmin(tmax, upper)) / 2 - lower)
}


# fix names and add gdd
prism <- prism_in %>%
  rename(
    SiteID = Name,
    ElevFt = `Elevation (ft)`,
    tminF = `tmin (degrees F)`,
    tmaxF = `tmax (degrees F)`
  ) %>%
  mutate(
    State = substr(SiteID, 1, 2),
    Year = as.numeric(format(Date, "%Y")),
    Julian = lubridate::yday(Date)) %>%
  group_by(SiteID, Year) %>%
  arrange(SiteID, Date) %>%
  mutate(
    GDD39 = cumsum(gdd_fn(tminF, tmaxF, 39, 86)),
    GDD50 = cumsum(gdd_fn(tminF, tmaxF, 50, 86))) %>%
  ungroup()


# GDD accumulation by state and year
prism %>%
  ggplot(aes(x = Date, y = GDD39, color = SiteID)) +
  geom_line() +
  facet_grid(State ~ Year, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = "none")



# load data files ---------------------------------------------------------

# read aphid counts (including dummy counts)
aph_in <- read_csv("data/stn_data_20191016.csv") %>%
  filter(!SpeciesName %in% c('Phylloxeridae', 'Adelgidae')) %>%
  mutate(
    Month = as.numeric(format.Date(Date, format = "%m")),
    Week = lubridate::week(Date)) %>%
  mutate_if(is.character, as.factor)

# site info
aph_sites <- read_csv("data/stn_sites.csv")

# species names
aph_spp <- read_csv("data/aphid_species.csv", na = c('', '.'))



# Expand and join dataset --------------------------------------------------

expandFn = function(df) {
  df %>%
    pivot_wider(
      names_from = "SpeciesName",
      values_from = "Count",
      values_fill = list(Count = 0)
    ) %>%
    select(-dummy) %>%
    pivot_longer(
      cols = 11:ncol(.),
      names_to = "SpeciesName",
      values_to = "Count") %>%
    mutate_if(is.character, as.factor)
}

## Expand dataset with zeros, drop dummy variable ##
aph_exp <- 
  aph_in %>%
  pivot_wider(
    names_from = "SpeciesName",
    values_from = "Count",
    values_fill = list(Count = 0)) %>%
  select(-"_dummy_") %>%
  pivot_longer(
    cols = 11:ncol(.),
    names_to = "SpeciesName",
    values_to = "Count") %>%
  mutate(SpeciesName = as.factor(SpeciesName))


## create state name and centroid dataset ##
states <- 
  tibble(
    State = state.abb,
    StateName = state.name,
    StateLon = state.center$x,
    StateLat = state.center$y
  ) %>%
  arrange(desc(StateLat)) %>%
  mutate(
    State = fct_inorder(State),
    StateName = fct_inorder(StateName)
  )

# get a vector of state abbreviations sorted by northiness
StateOrder <- states %>%
  arrange(desc(StateLon)) %>%
  .$State


# join month num, state names, degree days, site info, and taxonomic info
aph_full <- 
  aph_exp %>%
  mutate_if(is.factor, as.character) %>%
  left_join(prism[, c('SiteID', 'Date', 'GDD39', 'GDD50')],
            by = c('SiteID', 'Date')) %>%
  left_join(aphid_sites[, c('SiteID', 'Lat', 'Lon')],
            by = 'SiteID') %>%
  left_join(aphid_spp,
            by = 'SpeciesName') %>%
  left_join(states) %>%
  mutate_if(is.character, as.factor)


# Export file (optional)
aph_full %>%
  write_csv("data/aphid stn data export.csv.gz")

# # export Minnesota data
# aphid_full %>%
#   filter(State == "MN") %>%
#   write_csv(paste0("mn_aphids_", min(.$Year), "-", max(.$Year), ".csv"))



# Subset aphid data -------------------------------------------------------

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
  aph_full %>%
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

# Select top n spp by count in entire dataset
aphid <- topSpFn(aph_full, 10)

# Select top n spp in wisconsin
aphid <- 
  aph_full %>%
  filter(State == "WI") %>%
  topSpFn(10)

# Filter Wisconsin only
aphid <- filter(aphid, State == "WI")




# Numerical summaries -----------------------------------------------------

# summarise by total count
aphid %>%
  group_by(SpeciesName) %>%
  summarise(TotalCount = sum(Count), n = n()) %>%
  arrange(desc(TotalCount)) %>%
  write_csv("out/totalCount.csv")


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
  write_csv("out/totalCount.csv")


# generate weekly mean counts for wi
aphid %>%
  filter(State == "WI") %>%
  select(SpeciesName, Week, Count) %>%
  group_by(SpeciesName, Week) %>%
  summarise(MeanCount = mean(Count)) %>%
  write_csv("out/wi_wkly.csv")


# unique species per state
aph_full %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  mutate(SpeciesName = as.character(SpeciesName)) %>%
  summarise(N_Taxa = length(unique(SpeciesName)))


# unique sites per state
aph_full %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  summarise(N_Sites = length(unique(SiteID)))


# summary by state
aph_summary_states <-
  aph_full %>%
  filter(Count > 0) %>%
  group_by(State) %>%
  summarise(
    Sites = n_distinct(SiteID),
    Years = n_distinct(Year),
    Samples = n_distinct(SampleID),
    Families = n_distinct(Family),
    Genera = n_distinct(Genus),
    Species = n_distinct(Species),
    TotCt = sum(Count)
  ) %>%
  mutate(MeanCt = TotCt / Samples)

# summary for whole dataset
aph_summary_all <-
  aph_full %>%
  filter(Count > 0) %>%
  summarise(
    State = n_distinct(State),
    Sites = n_distinct(SiteID),
    Years = n_distinct(Year),
    Samples = n_distinct(SampleID),
    Families = n_distinct(Family),
    Genera = n_distinct(Genus),
    Species = n_distinct(Species),
    TotCt = sum(Count)
  ) %>%
  mutate(MeanCt = TotCt / Samples)

# bind and save summary
rbind(aph_summary_states, aph_summary_all) %>%
  write_csv("out/STN summary - diversity all spp.csv")


# summarise by frequency
aph_summary_freq <- 
  aph_full %>%
  mutate(isFound = case_when(Count == 0 ~ 0, Count > 0 ~ 1)) %>%
  group_by(SpeciesName) %>%
  summarise(Freq = sum(isFound)/n()) %>%
  arrange(desc(Freq))
write_csv(aph_summary_freq,
          "out/STN summary - detection frequency by species.csv")


# return number of observations by species name
aph_summary_spcounts <- 
  aph_full %>%
  filter(Count > 0) %>%
  group_by(Family, Subfamily, Genus, Species, SpeciesName) %>%
  summarise(Obs = n(),
            Count = sum(Count)) %>%
  mutate(Mean = Count / Obs)
aph_summary_spcounts %>%
  write_csv("out/STN summary - spp counts.csv")


# weeks per year
aph_full %>%
  group_by(Year) %>%
  summarise(N_Weeks = length(unique(Week)),
            N_Captures = sum(Count))





# Summary plots -----------------------------------------------------------

## Plot of counts by state ##
aphid %>%
  group_by(State, SpeciesName) %>%
  summarise(totcount = mean(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = reorder(SpeciesName, totcount), y = totcount)) +
  facet_grid( ~ State) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  coord_flip() +
  labs(
    x = "",
    y = "Mean Count",
    title = paste0("Aphid Suction Trap Network (", min(aphid$Year), "-", max(aphid$Year), "): Captures by State")) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5),
    legend.position = "none")

ggsave("out/Plot of counts by state.svg", h = 6, w = 10)



## Plot counts by year ##
aphid %>%
  filter(State == "WI") %>%
  group_by(Year, SpeciesName) %>%
  summarise(totcount = sum(Count)) %>%
  arrange(SpeciesName, desc(totcount)) %>%
  ggplot(aes(x = as.factor(Year), y = totcount)) +
  facet_wrap(~ reorder(SpeciesName, -totcount)) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  labs(
    x = "", y = "Mean Count",
    title = paste("Aphid Suction Trap Network (", min(aphid$Year), "-", max(aphid$Year), "): Captures by Year", sep = '')) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(size = 12, face = "italic"),
    legend.position = "none")

ggsave("out/STN plot - counts by year.svg", h = 6, w = 10)



## Plot counts by state and year ##
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
    strip.text.x = element_text(size = 10, face = "italic"),
    strip.text.y = element_text(angle = 0, size = 12, face = "bold")
  )
plt
ggsave("out/STN plot - counts by state and year.svg", plt, h = 10, w = 16)



## Plot of aphid abundance by GDD and year ##
plt <- 
  aphid %>%
  mutate(Year = as.character(Year),
         Count = log(Count + 1)) %>%
  ggplot(aes(x = GDD50, y = Count)) +
  facet_wrap(~ SpeciesName, scales = "free") +
  scale_x_sqrt() +
  stat_smooth(
    aes(color = Year),
    method = "gam",
    formula = y ~ s(x),
    fill = NA,
    size = .5
  ) +
  stat_smooth(
    color = "red",
    method = "gam",
    formula = y ~ s(x),
    size = 1.25
  ) +
  geom_abline(intercept = 0, slope = 0) +
  labs(
    x = "Growing Degree Days (50F/86F)",
    y = "Log-normalized aphid count",
    title = "Aphid abundance by year and GDD (Suction Trap data, 2005-2018)"
  ) +
  theme(strip.text = element_text(angle = 0, size = 12, face = "italic"))
plt
ggsave("out/STN plot - aphid abundance by year and GDD.svg", h = 10, w = 16)





## Plot multi-year comparison of counts by week ##
plt <- 
  aphid %>%
  mutate(Year = as.character(Year),
         Count = log(Count + 1)) %>%
  ggplot(aes(x = Week, y = Count)) +
  facet_wrap( ~ SpeciesName, scales = "free") +
  stat_smooth(aes(color = Year),
              method = "gam",
              formula = y ~ s(x),
              fill = NA,
              size = .5) +
  stat_smooth(color = "red",
              method = "gam",
              formula = y ~ s(x),
              size = 1.25) +
  geom_abline(intercept = 0, slope = 0) +
  labs(aes(x = "Week",
           y = "Log-normalized aphid abundance",
           title = "Aphid abundance by week and year (Suction Trap data, 2005-2018)",
           legend = "Year")) +
  theme(strip.text = element_text(angle = 0, size = 12, face = "italic"))
plt
ggsave("out/STN plot - aphid abundance by year and week.svg", h = 10, w = 16)


# boxplots by year and species, not very useful
p = aphid %>%
  ggplot(aes(x = as.factor(Year), y = log10(Count + 1), group = Year, color = as.factor(Year))) +
  geom_boxplot() +
  facet_grid(SpeciesName ~ .) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  geom_smooth( method = "gam", formula = y ~ s(x) + 1, fill = "grey50", size = 1)
p




# Russ charts Jan 2019 ----------------------------------------------------

# Top 5 wis aphids
russ_aphids = c(
  "Aphis glycines",
  "Acyrthosiphon pisum",
  "Rhopalosiphum padi",
  "Rhopalosiphum maidis",
  "Myzus persicae"
  )

#generate dataset
aphid <- 
  aph_full %>%
  filter(SpeciesName %in% russ_aphids) %>%
  droplevels()


## Violin plots by species and year ##
plt <-
  aphid %>%
  group_by(Year, Week, SpeciesName, SiteID) %>%
  summarise(Count = sum(Count)) %>%
  summarise(MeanCount = log(mean(Count) + 1)) %>%
  ggplot(aes(x = as.Date("2020-01-01") + (Week - 1) * 7,
             fill = SpeciesName)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = .25, color = "grey") +
  geom_area(aes(y = MeanCount)) +
  geom_area(aes(y = -MeanCount)) +
  facet_grid(Year ~ SpeciesName, scales = "free_x") +
  scale_x_date(date_labels = "%b") +
  labs(title = "Aphid abundance over time by year (Suction trap data 2005-2018)",
       x = "",
       y = "Log abundance") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = .5),
    panel.spacing = unit(.25, "lines"),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text.x = element_text(size = 10, face = "italic"),
    strip.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
plt
ggsave("out/Russ aphids 2019 violin plots by year.png", plt, h = 10, w = 16)


## violin plots by species and state ##
plt <-
  aphid %>%
  mutate(State = factor(State, levels = StateOrder)) %>%
  group_by(State, Week, SpeciesName, SiteID) %>%
  summarise(Count = sum(Count)) %>% # get site totals
  summarise(MeanCount = log(mean(Count) + 1)) %>% # average across sites
  ggplot(aes(x = as.Date("2017-01-01") + (Week - 1) * 7,
             fill = SpeciesName)) +
  geom_hline(yintercept = 0) +
  geom_area(aes(y = MeanCount)) +
  geom_area(aes(y = -MeanCount)) +
  facet_grid(State ~ SpeciesName, scales = "free_x") +
  scale_x_date(date_labels = "%b") +
  labs(title = "Aphid abundance at trap sites (Suction trap data, 2005-2018)",
       x = "",
       y = "Log abundance") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = .5),
    panel.spacing = unit(.25, "lines"),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text.x = element_text(size = 10, face = "italic"),
    strip.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
plt
ggsave("out/Russ aphids 2019 violin plots by state.png", plt, h = 10, w = 16)



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




# Glycines IDW export for Emily -------------------------------------------

# filter
glycines = aph_full %>%
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




# Aphis glycines counts by state ------------------------------------------

aph_full %>%
  filter(SpeciesName == "Aphis glycines") %>%
  mutate(State = factor(State, levels = StateOrder),
         Year = as.character(Year)) %>%
  group_by(State, Year) %>%
  summarise(Count = mean(Count) + 1) %>%
  arrange(desc(Count)) %>%
  ggplot(aes(x = Year, y = Count)) +
  geom_hline(yintercept = 1) +
  facet_grid(State ~ .) +
  geom_bar(stat = "identity", aes(fill = State)) +
  scale_y_log10() +
  labs(title = "Aphis glycinces abundance by year and state",
       x = "",
       y = "Mean count per week") +
  theme_gray() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(angle = 0, size = 12, face = "bold")
  )


## bar plots by year for specific species ##
sp = "Aphis glycines"
sp = "Rhopalosiphum maidis"
sp = "Rhopalosiphum padi"
sp = "Sitobion avenae"
sp = "Schizaphis graminum"
plt <- 
  aph_full %>%
  filter(SpeciesName == sp) %>% # edit species name here
  mutate(Year = as.Date(paste0(Year, "-01-01"))) %>%
  group_by(StateName, Year) %>%
  summarise(mean = mean(Count) + 1,
            sd = sd(Count)) %>%
  arrange(desc(StateName)) %>%
  ggplot(aes(x = Year, y = mean, fill = StateName)) +
  geom_hline(yintercept = 1) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean,
                    ymax = mean + sd)) +
  scale_y_log10() +
  facet_wrap( ~ StateName, nrow = 2) +
  labs(title = bquote(bold(bolditalic(.(sp)) ~ "captures, 2005-2018")),
       x = "",
       y = bquote(bold("Mean captures per trap per week"))) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plt
ggsave("out/padi captures.png", p, w = 12, h = 5, u = "in")



# Total insect abundance over time ----------------------------------------

aph_full %>%
  group_by(Year, SiteID, SampleID) %>%
  summarise(totwkly = sum(Count)) %>%
  summarise(meanbysite = mean(totwkly)) %>%
  ggplot(aes(x = Year,
             y = meanbysite)) +
  geom_smooth() +
  geom_point() +
  scale_y_sqrt(expand = c(0, 0)) +
  scale_x_continuous(breaks = 2005:2018, minor_breaks = NULL) +
  labs(y = "Mean weekly capture per site",
       title = "Aphid abundance per site in the Aphid Suction Trap Network") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))
