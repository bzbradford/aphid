
# Libraries ---------------------------------------------------------------

library(tidyverse)


# PVY-relevant aphid species and infection risk ---------------------------

# Read pvy infectivity data
pvy_inf <- read_csv("data/pvy-inf.csv")


# filter by all potential pvy vectors
aph_pvy <- 
  aph_full %>%
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


# Or filter by aphid species with PVY risk values
aph_pvy <- 
  aph_full %>%
  filter(SpeciesName %in% pvy_inf$SpeciesName) %>%
  mutate_if(is.factor, as.character) %>%
  left_join(pvy_inf) %>%
  mutate(
    PVYCount = Count * asin(sqrt(InfRtAvg)),
    PVYRisk = log(Count + 1) * asin(sqrt(InfRtAvg))) %>%
  mutate_if(is.character, as.factor)




# Summaries ---------------------------------------------------------------

# generate weekly mean risk scores
pvy_risk_wkly_wi <- 
  aph_pvy %>%
  filter(State == "WI") %>%
  select(SpeciesName, Week, Count, PVYRisk) %>%
  group_by(SpeciesName, Week) %>%
  summarise(MeanCount = mean(Count), MeanRisk = mean(PVYRisk))
pvy_risk_wkly_wi %>%
  write_csv("out/pvy risk wkly wisc.csv")



# generate annual count summaries and a total count by species, then save

pvy_risk_yrly_wi <- 
  aph_pvy %>%
  filter(State == "WI") %>%
  group_by(SpeciesName, Year) %>%
  summarise(TotalRiskYrly = sum(Count)) %>%
  summarise(MeanRiskYrly = mean(TotalRiskYrly))
pvy_risk_yrly_wi %>%
  write.csv("out/pvy wi total count.csv")

test = aphid_pvy %>%
  filter(State == "WI") %>%
  droplevels()
levels(test$SiteName)




# PVY phenology, all spp -----------------------------------------------------

# collapse to aggregate counts
pvy_aggr_risk <- 
  aph_pvy %>%
  filter(State == "MI", Year >= 2010) %>% # MICHIGAN
  group_by(SampleID) %>%
  summarise(
    State = State[1],
    SiteID = SiteID[1],
    Year = Year[1],
    GDD39 = GDD39[1],
    Risk = as.integer(sum(PVYRisk)*10)
  )

# Generate conditional modes from aggregate risk scores
pvy_cm <-
  pvy_aggr_risk %>%
  {
    fmod <- lme4::ranef(lme4::glmer(
      Risk ~ 1 +
        (1 | Year) +
        # (1 | State) +
        (1 | SiteID) +
        (1 | GDD39),
      data = .,
      family = "poisson"))
    tibble(GDD = as.numeric(rownames(fmod$GDD39)),
           CM = fmod$GDD39[, 1])
  }

pvy_gampreds <-
  tibble(
    GDD = 1:7500,
    CM = mgcv::predict.gam(
      mgcv::gam(CM ~ s(GDD), data = pvy_cm),
      data.frame(GDD = 1:7500)))

pvy_gampts <- 
  pvy_gampreds %>%
  mutate(
    type = case_when(
      sign(lag(CM)) < sign(CM) ~ "RI",
      sign(lag(CM)) > sign(CM) ~ "FI",
      lag(CM) < CM & CM > lead(CM) ~ "max",
      lag(CM) > CM & CM < lead(CM) ~ "min")) %>%
  na.omit()


# risk index
pvy_cm %>%
  ggplot(aes(x = GDD, y = CM)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_line(data = pvy_gampreds, color = "red", size = 1.5) +
  geom_point(data = pvy_gampts, aes(x = GDD, y = CM)) +
  geom_label(data = pvy_gampts, aes(x = GDD, y = CM, label = GDD)) +
  scale_x_sqrt(limits = c(500, 5000)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Degree Days (39/86)",
    y = "Conditional mode (spline fit deviance)",
    title = paste("Aggregate PVY risk index derived from Suction Trap data",
      min(aphid$Year), "-", max(aphid$Year))
  ) +
  theme(strip.text = element_text(size = 14, face = "bold"))

ggsave("out/pvy risk curve, all spp.png", h = 4, w = 6)





# plots -------------------------------------------------------------------




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



# Michigan ----------------------------------------------------------------


aph_sites %>%
  filter(State == "MI") %>%
  leaflet::leaflet() %>%
  leaflet::addTiles() %>%
  leaflet::addMarkers


# collapse to aggregate counts
pvy_aggr_risk <- 
  aph_pvy %>%
  filter(State == "MI", Year >= 2010) %>% # MICHIGAN
  group_by(SampleID) %>%
  summarise(
    State = State[1],
    SiteID = SiteID[1],
    Year = Year[1],
    GDD39 = GDD39[1],
    Risk = as.integer(sum(PVYRisk)*10)
  )

pvy_aggr_risk %>%
  ggplot(aes(x = GDD39, y = Risk)) +
  geom_point(aes(color = factor(Year)), alpha = .75) + 
  geom_smooth(method = "loess", span = 1/4, se = F) +
  scale_color_viridis_d() +
  labs(
    title = "PVY risk score by GDD for Michigan",
    color = "Year")

ggsave("out/mi/raw PVY risk by GDD.png", type = "cairo")


# GDD peaks @ 1258 and 3882

mi_prism <- prism %>%
  filter(State == "MI") %>%
  mutate(DOY = lubridate::yday(Date)) %>%
  group_by(Year, DOY) %>%
  summarise(GDD39 = mean(GDD39)) %>%
  mutate(Date = as.Date("2021-01-01") + DOY - 1)

mi_risk_ranges <- mi_prism %>%
  group_by(Year) %>%
  summarise(
    {
      approx(GDD39, DOY, xout = c(1258, 3882), ties = "ordered", rule = 2) %>%
        setNames(c("GDD39", "DOY")) %>%
        as_tibble() %>%
        drop_na()
    },
    .groups = "drop"
  ) %>%
  mutate(Date = as.Date("2021-01-01") + DOY - 1) %>%
  group_by(GDD39) %>%
  summarise(
    RiskStart = min(Date),
    RiskPeak = mean(Date),
    RiskEnd = max(Date)
  )

mi_risk_ranges %>%
  ggplot() +
  geom_hline(aes(yintercept = GDD39), color = "#c5050c", size = 0.5, linetype = "dashed") +
  geom_text(aes(x = as.Date("2021-1-1"), y = GDD39, label = GDD39), vjust = -0.25, color = "#c5050c") +
  geom_rect(aes(xmin = RiskStart, xmax = RiskEnd, ymin = -Inf, ymax = GDD39), alpha = 0.5) +
  geom_segment(aes(x = RiskPeak, y = -Inf, xend = RiskPeak, yend = GDD39), alpha = 0.5) +
  geom_text(aes(x = RiskStart, y = -Inf, label = format(RiskStart, format = "%b %d")), angle = 90, hjust = -.1) +
  geom_text(aes(x = RiskPeak, y = -Inf, label = format(RiskPeak, format = "%b %d")), angle = 90, hjust = -.1) +
  geom_text(aes(x = RiskEnd, y = -Inf, label = format(RiskEnd, format = "%b %d")), angle = 90, hjust = -.1) +
  geom_line(data = mi_prism, aes(x = Date, y = GDD39, color = factor(Year))) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_color_viridis_d() +
  labs(
    title = "Degree day (base 39F) accumulation in Michigan, 2005-2021",
    x = NULL,
    y = "Cumulative degree days (base 39F)",
    color = "Year") +
  theme_classic() +
  theme(panel.grid.major = element_line(linetype = "dotted"))

ggsave("out/mi/mi gdd plot.png", type = "cairo")




# Generate conditional modes from aggregate risk scores
pvy_cm <- {
    fmod <- lme4::ranef(lme4::glmer(
      Risk ~ 1 +
        (1 | Year) +
        # (1 | State) +
        (1 | SiteID) +
        (1 | GDD39),
      data = pvy_aggr_risk,
      family = "poisson"))
    tibble(GDD = as.numeric(rownames(fmod$GDD39)),
      CM = fmod$GDD39[, 1])
  }

pvy_gampreds <-
  tibble(
    GDD = 1:7500,
    CM = mgcv::predict.gam(
      mgcv::gam(CM ~ s(GDD), data = pvy_cm),
      data.frame(GDD = 1:7500)))

pvy_gampts <- 
  pvy_gampreds %>%
  mutate(
    type = case_when(
      sign(lag(CM)) < sign(CM) ~ "RI",
      sign(lag(CM)) > sign(CM) ~ "FI",
      lag(CM) < CM & CM > lead(CM) ~ "max",
      lag(CM) > CM & CM < lead(CM) ~ "min")) %>%
  na.omit()


max(mi_aphids$GDD39)


# risk index
pvy_cm %>%
  ggplot(aes(x = GDD, y = CM)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_line(data = pvy_gampreds, color = "red", size = 1.5) +
  geom_point(data = pvy_gampts, aes(x = GDD, y = CM)) +
  geom_label(data = pvy_gampts, aes(x = GDD, y = CM, label = GDD)) +
  scale_x_sqrt(limits = c(500, 6000), breaks = scales::pretty_breaks()) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Degree day accumulation (base 39F, upper 86F)",
    y = "Conditional mode (~relative abundance)",
    title = paste("Michigan aggregate PVY risk index derived from Suction Trap data 2010-2021")
  ) +
  theme(strip.text = element_text(size = 14, face = "bold"))

ggsave("out/mi/pvy risk curve.png", h = 4, w = 6)
