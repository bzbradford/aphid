

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
      "Aphis craccivora",
      "Sitobion avenae",
      "Rhopalosiphum maidis",
      "Brachycaudus helichrysi",
      "Brevicoryne brassicae"
    )
  ) %>% droplevels()


# PVY-relevant species in wisconsin only
aphidpvy <-
  aph_full %>%
  filter(
    SpeciesName %in% c(
      "Myzus persicae",
      "Macrosiphum euphorbiae",
      "Aphis glycines",
      "Acyrthosiphon pisum",
      "Rhopalosiphum padi",
      "Capitophorus elaeagni",
      "Aphis craccivora",
      "Sitobion avenae",
      "Rhopalosiphum maidis",
      "Brachycaudus helichrysi",
      "Brevicoryne brassicae"
    )
  ) %>%
  filter(State == "WI") %>%
  left_join(pvy_inf) %>%
  mutate(SpeciesName = as.factor(SpeciesName)) %>%
  droplevels() %>%
  mutate(PVYRisk = log(Count + 1)*asin(sqrt(InfAvg)))


# generate weekly mean risk scores
aphidpvy %>%
  select(SpeciesName, Week, Count, PVYRisk) %>%
  group_by(SpeciesName, Week) %>%
  summarise(MeanCount = mean(Count), MeanRisk = mean(PVYRisk)) %>%
  write.csv("out/wi_pvy_wkly.csv")


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
