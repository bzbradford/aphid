

# read pvy infectivity data
pvy_inf <-
  read.csv("data/pvy-inf.csv", header = TRUE)
names(pvy_inf) <- c("SpeciesName", "InfLow", "InfUp", "InfAvg")



# Filters -----------------------------------------------------------------

# all potential pvy vectors
aphid <-
  aphid_full %>%
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


# PVY species with known infectivity

# keep aphid species with PVY risk values
aphid_pvy <- aphid_full %>%
  left_join(pvy_inf) %>%
  mutate(SpeciesName = as.factor(SpeciesName)) %>%
  drop_na(InfAvg) %>%
  mutate(PVYCount = Count*asin(sqrt(InfAvg)),
         PVYRisk = log(Count + 1)*asin(sqrt(InfAvg))) %>%
  droplevels()

# keep wisconsin only
aphid_pvy <- filter(aphid_pvy, State == "WI")

# collapse to aggregate counts
aphid_pvy_aggr <-
  aphid_pvy %>%
  group_by(SampleID) %>%
  summarise(
    SiteID = SiteID[1],
    SpeciesName = "Aggregate",
    Year = Year[1],
    GDD39 = GDD39[1],
    Count = sum(PVYRisk)
  )

# collapse to aggregate counts, intergers only
aphid_pvy_aggr <-
  aphid_pvy %>%
  group_by(SampleID) %>%
  summarise(
    SiteID = SiteID[1],
    SpeciesName = "Aggregate",
    Year = Year[1],
    GDD39 = GDD39[1],
    Count = as.integer(sum(PVYCount))
  )



# Summaries ---------------------------------------------------------------

# generate weekly mean risk scores
aphid_pvy %>%
  select(SpeciesName, Week, Count, PVYRisk) %>%
  group_by(SpeciesName, Week) %>%
  summarise(MeanCount = mean(Count), MeanRisk = mean(PVYRisk)) %>%
  write.csv("out/wi_pvy_wkly.csv")


# Phenology code ----------------------------------------------------------

# CMs from PVY Risk
cm_gdd_pvycount =
  aphid_pvy_aggr %>%
  group_by(SpeciesName) %>%
  do(remfn(.))

# generate predictions from gam fit
gampreds = gampredfn(cm_gdd_pvycount)

# identify critical points
gampts = gamptfn(gampreds)




# Subset into groups ------------------------------------------------------

# subsets for phenology plots
aphid1 <- aphid_full %>%
  filter(SpeciesName %in% c(
    "Aphis glycines",
    "Rhopalosiphum padi",
    "Acyrthosiphon pisum",
    "Myzus persicae"
  )) %>% droplevels()

aphid2 <- aphid_full %>%
  filter(SpeciesName %in% c(
    "Rhopalosiphum maidis",
    "Capitophorus elaeagni",
    "Macrosiphum euphorbiae",
    "Aphis craccivora"
  )) %>% droplevels()

aphid3 <- aphid_full %>%
  filter(SpeciesName %in% c(
    "Sitobion avenae",
    "Brachycaudus helichrysi"
  )) %>% droplevels()
