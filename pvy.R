### PVY-relevant aphid species and infection risk ###

# Load packages ----
library(tidyverse)
library(lme4)
library(mgcv)

# Read pvy infectivity data ----

pvy_inf = read.csv("data/pvy-inf.csv", header = TRUE)
names(pvy_inf) = c("SpeciesName", "InfLo", "InfHi", "InfAvg")


# Subsets of the aphid dataset ----

# all potential pvy vectors
aphid_pvy =
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


# keep aphid species with PVY risk values
aphid_pvy = 
  aphid_full %>%
  left_join(pvy_inf) %>%
  mutate(SpeciesName = as.factor(SpeciesName)) %>%
  drop_na(InfAvg) %>%
  mutate(PVYCount = Count*asin(sqrt(InfAvg)),
         PVYRisk = log(Count + 1)*asin(sqrt(InfAvg))) %>%
  droplevels()


# collapse to aggregate counts
aphid_pvy_aggr_risk =
  aphid_pvy %>%
  group_by(SampleID) %>%
  summarise(
    State = State[1],
    SiteID = SiteID[1],
    Year = Year[1],
    GDD39 = GDD39[1],
    Risk = as.integer(sum(PVYRisk)*10)
  )

aphid_pvy_aggr_risk_wi =
  aphid_pvy_aggr_risk %>%
  filter(State == "WI")

# collapse to aggregate counts, integers only
aphid_pvy_aggr_count =
  aphid_pvy %>%
  group_by(SampleID) %>%
  summarise(
    State = State[1],
    SiteID = SiteID[1],
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



# generate annual count summaries and a total count by species, then save

aphid_pvy %>%
  filter(State == "WI") %>%
  group_by(SpeciesName, Year) %>%
  summarise(TotalAnnual = sum(Count)) %>%
  summarise(MeanAnnual = mean(TotalAnnual)) %>%
  write.csv("out/pvy wi total count.csv")

test = aphid_pvy %>%
  filter(State == "WI") %>%
  droplevels()
levels(test$SiteName)

# CMs ---------------------------------------------------------------------

# CMs from PVY Risk
pvy_cm =
  aphid_pvy %>%
  {
    fmod =
      ranef(glmer(
        Risk ~ 1 + (1 | GDD39) + (1 | Year) + (1 | SiteID),
        data = .,
        family = "poisson"
      ))
    data.frame(GDD = as.numeric(rownames(fmod$GDD39)),
               CM = fmod$GDD39[, 1])
  }

# CMs for particular species
pvy_cm =
  aphid_pvy %>%
  filter(State == "WI",
         SpeciesName == "Rhopalosiphum padi") %>%
  {
    fmod =
      ranef(glmer(
        Count ~ 1 + (1 | GDD39) + (1 | Year) + (1 | SiteID),
        data = .,
        family = "poisson"
      ))
    data.frame(GDD = as.numeric(rownames(fmod$GDD39)),
               CM = fmod$GDD39[, 1])
  }

pvy_gampreds =
  cbind(data.frame(GDD = 1:7500),
        CM = predict.gam(gam(CM ~ s(GDD), data = pvy_cm),
                    data.frame(GDD = 1:7500)))

pvy_gampts =
  pvy_gampreds %>%
  mutate(
    type = case_when(
      sign(lag(CM)) < sign(CM) ~ "RI",
      sign(lag(CM)) > sign(CM) ~ "FI",
      lag(CM) < CM & CM > lead(CM) ~ "max",
      lag(CM) > CM & CM < lead(CM) ~ "min"
    )
  ) %>%
  na.omit()


# plots -------------------------------------------------------------------

# risk index
p =
  pvy_cm %>%
  ggplot(aes(x = GDD, y = CM)) +
  scale_x_sqrt(limits = c(500, 5000)) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_abline(intercept = 0, slope = 0) +
  labs(x = "Degree Days (39/86)",
       y = "Conditional mode (spline fit deviance)",
       title = paste("Aggregate PVY risk index for Wisconsin",
                     min(aphid$Year), "-", max(aphid$Year))
  ) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  geom_line(data = pvy_gampreds, color = "red", size = 1.5) +
  geom_point(data = pvy_gampts, aes(x = GDD, y = CM)) +
  geom_label(data = pvy_gampts, aes(x = GDD, y = CM, label = GDD))
p

ggsave("out/pvy risk cm wi.png",
       plot = p,
       type = "cairo-png",
       h = 4,
       w = 6,
       dpi = 300)


# species plots
p =
  pvy_cm %>%
  ggplot(aes(x = GDD, y = CM)) +
  scale_x_sqrt() +
  geom_abline(intercept = 0, slope = 0) +
  labs(x = "Degree Days (39/86)",
       y = "Conditional mode (spline fit deviance)",
       title = paste("Rhopalosiphum padi derived flight phenology")
  ) +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  geom_line(data = pvy_gampreds, color = "red", size = 1.5) +
  geom_point(data = pvy_gampts, aes(x = GDD, y = CM)) +
  geom_label(data = pvy_gampts, aes(x = GDD, y = CM, label = GDD))
p

ggsave("out/wi r padi model fit.png",
       plot = p,
       type = "cairo-png",
       h = 4,
       w = 6,
       dpi = 300)

# soybean aphid counts
aphid_pvy %>%
  filter(State == "WI", SpeciesName == "Aphis glycines") %>%
  mutate(Week = as.factor(Week)) %>%
  ggplot(aes(x = Week, y = Count)) +
  geom_boxplot() +
  scale_y_sqrt()
p




# test by state ----
pvy_cm =
  aphid_pvy_aggr %>%
  group_by(State) %>%
  do({
    fmod =
      ranef(glmer(
        Risk ~ 1 + (1 | GDD39) + (1 | Year),
        data = .,
        family = "poisson"
      ))
    data.frame(GDD = as.numeric(rownames(fmod$GDD39)),
               CM = fmod$GDD39[, 1])
  })
p =
  pvy_cm %>%
  ggplot(aes(x = GDD, y = CM)) +
  geom_hline(yintercept = 0) +
  stat_smooth(
    method = "gam",
    formula = y ~ s(x),
    fill = "grey50",
    colour = "Black",
    size = 1.5
  ) +
  facet_wrap( ~ State, scales = "free") +
  labs(title = "PVY Risk estimates by degree day (39/86)") +
  theme(strip.text = element_text(face = "bold"))
p
CairoPNG("out/pvy_risk_cm.png", h = 800, w = 1200); p; dev.off()

# predicted values from a smooth fit of CMs
pvy_gampreds =
  pvy_cm %>%
  group_by(State) %>%
  do(
  cbind(data.frame(GDD = 0:as.integer(max(.$GDD))),
        CM = predict.gam(gam(CM ~ s(GDD), data = .),
                         data.frame(GDD = 0:as.integer(max(.$GDD)))))
  )

# inflection points and intercepts
pvy_gampts =
  pvy_gampreds %>%
  group_by(State) %>%
  mutate(
    type = case_when(
      sign(lag(CM)) < sign(CM) ~ "RI",
      sign(lag(CM)) > sign(CM) ~ "FI",
      lag(CM) < CM & CM > lead(CM) ~ "max",
      lag(CM) > CM & CM < lead(CM) ~ "min"
    )
  ) %>%
  na.omit()


# graph
p =
  pvy_gampreds %>%
  group_by(State) %>%
  ggplot(aes(x = GDD, y = CM)) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = 0) +
  geom_line(color = "red", size = 1) +
  facet_wrap(~ State, scales = "free_x") +
  geom_point(data = pvy_gampts, aes(x = GDD, y = CM)) +
  geom_vline(data = pvy_gampts, aes(xintercept = GDD)) +
  labs(title = "PVY Risk estimates by degree day (39/86)") +
  theme(strip.text = element_text(face = "bold"))
p
CairoPNG("out/pvy risk by state.png", h = 800, w = 1200); p; dev.off()


# # Subset into groups ------------------------------------------------------
# 
# # subsets for phenology plots
# aphid1 <- aphid_full %>%
#   filter(SpeciesName %in% c(
#     "Aphis glycines",
#     "Rhopalosiphum padi",
#     "Acyrthosiphon pisum",
#     "Myzus persicae"
#   )) %>% droplevels()
# 
# aphid2 <- aphid_full %>%
#   filter(SpeciesName %in% c(
#     "Rhopalosiphum maidis",
#     "Capitophorus elaeagni",
#     "Macrosiphum euphorbiae",
#     "Aphis craccivora"
#   )) %>% droplevels()
# 
# aphid3 <- aphid_full %>%
#   filter(SpeciesName %in% c(
#     "Sitobion avenae",
#     "Brachycaudus helichrysi"
#   )) %>% droplevels()
