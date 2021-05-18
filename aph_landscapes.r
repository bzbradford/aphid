
library(sf)
library(raster)
library(landscapemetrics)
library(tidyverse)

# required from aph_data_mgmt
aph_full
aph_sites


# create sf object
aph_sites_sf <- aph_sites %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)
aph_sites_sf %>% ggplot() + geom_sf()

# transform to aea
aph_sites_sf_albers <- aph_sites_sf %>%
  st_transform(crs = crs(nlcd))
aph_sites_sf_albers %>% ggplot() + geom_sf()



# NCDL extract and plot ----

# NLCD raster
nlcd <- raster("E:/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img")

plot(nlcd)
plot(aph_sites_sf_albers$geometry, add = T)

nlcd_classes <- read_csv("data/nlcd_classes.csv")


# get landscape
aph_nlcd_lsm <- 
  sample_lsm(
    landscape = nlcd,
    y = aph_sites_sf_albers,
    plot_id = aph_sites_sf_albers$SiteID,
    shape = "circle",
    size = 10000,
    what = "lsm_c_pland",
    progress = T
  )

aph_nlcd <- 
  aph_nlcd_lsm %>%
  left_join(nlcd_classes) %>%
  select(plot_id, class, class_name, value) %>%
  mutate(class_label = paste(class, class_name)) %>%
  rename(area = value)


# nlcd data
aph_nlcd %>%
  ggplot(aes(x = plot_id, y = area, fill = class_label)) +
  geom_col(position = position_stack(reverse = T), width = 1, color = "white") +
  coord_flip() +
  scale_y_continuous(expand = expansion()) +
  labs(
    x = "Suction trap site",
    y = "Percent of landscape",
    title = "NLCD-derived landscape class composition within 10km of suction trap sites",
    fill = "NLCD class") +
  theme_classic()

ggsave("land/nlcd class comparison.png", type = "cairo")



# aphid dataset -----------------------------------------------------------

# samples per site
aph_full %>%
  group_by(SiteID, Lat, Lon) %>%
  summarise(n = n())

aph_by_week <- aph_full %>%
  filter(Count > 0) %>%
  arrange(Lat) %>%
  mutate(SiteID = fct_inorder(SiteID)) %>%
  group_by(SiteID, Lat, Lon, Year, Week) %>%
  summarise(
    abundance = sum(Count),
    log_abund = log1p(abundance),
    richness = n_distinct(Species),
    .groups = "drop"
  )


plt.abund <- 
  aph_by_week %>%
  group_by(SiteID, Week) %>%
  summarise(abundance = log1p(mean(abundance))) %>%
  ggplot(aes(x = SiteID, y = Week, fill = abundance)) +
  geom_tile() +
  coord_flip() +
  scale_fill_viridis_c() +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  labs(
    x = "<-S      Site Name      N->",
    y = "Week of year",
    title = "Aphid abundance by week of year"
  )

  
plt.richness <- 
  aph_by_week %>%
  group_by(SiteID, Week) %>%
  summarise(richness = mean(species_richness)) %>%
  ggplot(aes(x = SiteID, y = Week, fill = richness)) +
  geom_tile() +
  coord_flip() +
  scale_fill_viridis_c() +
  scale_y_continuous(expand = expansion()) +
  theme_classic() +
  labs(
    x = "<-S      Site Name      N->",
    y = "Week of year",
    title = "Species richness"
  )

plt <- gridExtra::grid.arrange(plt.abund, plt.richness, ncol = 2)
ggsave("land/abundance and richness.png", plt, type = "cairo")




# NLCD models ----

nlcd_classes

aph_nlcd_wide <- 
  aph_nlcd %>%
  select(plot_id, class_name, area) %>%
  pivot_wider(
    names_from = class_name,
    values_from = area,
    values_fill = 0
  )

df <- aph_by_week %>%
  filter(between(Week, 20, 45)) %>%
  left_join(aph_nlcd, by = c("SiteID" = "plot_id")) %>%
  mutate(
    Year = factor(Year),
    Week = factor(Week),
    Class = factor(class_name),
    area = area / 100
  )



## abundance model ----

fit <- lm(log_abund ~ Year + Week + Class, weights = area, data = df)
# plot(fit)
hist(df$log_abund)
# try a shapiro-wilkes test


# plot
enframe(coef(fit)) %>%
  mutate(var = case_when(
    name == "(Intercept)" ~ "(Intercept)",
    grepl("Year", name) ~ "Year",
    grepl("Week", name) ~ "Week",
    grepl("Class", name) ~ "Class"
  )) %>%
  mutate(name = str_remove(name, var)) %>%
  ggplot(aes(x = reorder(name, desc(name)), y = value, fill = value)) +
  geom_col(width = 1, color = "white") +
  coord_flip() +
  labs(
    y = "Coefficient",
    x = "Variable",
    title = "Effect on aphid abundance of land cover within 10km of suction trap sites",
    subtitle = paste("Model:", fit$call[2])
  ) +
  facet_wrap(~ var, scales = "free")

ggsave("land/model coefs - abundance.png", type = "cairo")



## richness model ----

hist(df$richness)

fit <- lm(richness ~ Year + Week + Class, weights = area, data = df)
plot(fit)

# deal with the richness data shape, modelled by poisson
fit <- glm(richness ~ Year + Week + Class, weights = area, data = df, family = "poisson")
plot(fit)

# plot
enframe(coef(fit)) %>%
  mutate(var = case_when(
    name == "(Intercept)" ~ "(Intercept)",
    grepl("Year", name) ~ "Year",
    grepl("Week", name) ~ "Week",
    grepl("Class", name) ~ "Class"
  )) %>%
  mutate(name = str_remove(name, var)) %>%
  ggplot(aes(x = reorder(name, desc(name)), y = value, fill = value)) +
  geom_col(width = 1, color = "white") +
  coord_flip() +
  labs(
    y = "Coefficient",
    x = "Variable",
    title = "Effect on aphid species richness of land cover within 10km of suction trap sites",
    subtitle = paste("Model:", fit$call[2])
  ) +
  facet_wrap(~ var, scales = "free")

ggsave("land/model coefs - richness.png", type = "cairo")




## stats ----

library(lme4)

fit <- 
  aph_by_week %>%
  filter(between(Week, 20, 45)) %>%
  mutate(Year = factor(Year), Week = factor(Week)) %>%
  lmer(log_abund ~ Week + (1 | SiteID) + (1 | Year), data = ., REML = F)
summary(fit)


fortify.merMod(fit) %>%
  mutate(across(c(Year, Week), as.numeric)) %>%
  ggplot(aes(Week, log_abund, color = SiteID)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(aes(y = .fitted), fun = mean, geom = "line")



# CDL models ----

cdl_classes <- read_csv("data/cdl_classes.csv") %>%
  mutate(class_label = paste(class, class_name)) %>%
  mutate_if(is.character, fct_inorder)

# load extract, I had already done this apparently
aph_cdl <- read_csv("data/cdl_extract.csv") %>%
  mutate(
    buffer_sqm = pi * buffer^2,
    area_sqm = area_pxls * 30^2,
    area = area_sqm / buffer_sqm * 100
  ) %>%
  left_join(cdl_classes)

top_cdl <- aph_cdl %>%
  group_by(class_label) %>%
  summarise(max_area = max(area), mean_area = mean(area)) %>%
  filter(mean_area > 1)

# plot cdl class breakdown
aph_cdl %>%
  filter(buffer == 10000) %>%
  group_by(site, class_label) %>%
  summarise(area = mean(area)) %>%
  ggplot(aes(x = site, y = area, fill = class_label)) +
  geom_col(position = position_stack(reverse = T), width = 1, color = "white", size = .5) +
  coord_flip() +
  scale_y_continuous(expand = expansion()) +
  labs(
    x = "Suction trap site",
    y = "Percent of landscape",
    title = "CDL landscape class composition within 10km of suction trap sites",
    fill = "CDL class") +
  theme_classic() +
  guides(fill = guide_legend(ncol = 2, keywidth = unit(.25, "cm"), keyheight = unit(.25, "cm")))


# grid color
aph_cdl %>%
  filter(buffer == 10000) %>%
  group_by(site, class_label) %>%
  summarise(area = mean(area)) %>%
  complete(site, class_label, fill = list(area = 0)) %>%
  group_by(class_label) %>%
  mutate(mean_area = mean(area)) %>%
  filter(mean_area > .1) %>%
  ggplot(aes(x = site, y = class_label, fill = area)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    x = "Suction trap site",
    y = "Landscape class",
    title = "CDL landscape class composition within 10km of suction trap sites",
    fill = "Percent of buffer") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))



## Coef plots ----

# join specific CDL data
aph_cdl_joined <- aph_by_week %>%
  filter(between(Week, 20, 45)) %>%
  left_join({
    aph_cdl %>%
      filter(buffer == 10000) %>%
      group_by(class) %>%
      mutate(mean_area = mean(area)) %>%
      filter(mean_area > .1) %>%
      rename(SiteID = site, Year = year)
    }) %>%
  mutate(
    Year = factor(Year),
    Week = factor(Week),
    Class = factor(class_name)
  ) %>%
  filter(!is.na(area))


# what classes remain
aph_cdl_joined %>%
  group_by(class_name) %>%
  summarise(
    mean = mean(area),
    max = max(area)
  ) %>%
  arrange(max)


# model
hist(aph_cdl_joined$log_abund)
fit <- lm(log_abund ~ Year + Week + Class, weights = area, data = aph_cdl_joined)


# plot
enframe(coef(fit)) %>%
  mutate(var = case_when(
    name == "(Intercept)" ~ "(Intercept)",
    grepl("Year", name) ~ "Year",
    grepl("Week", name) ~ "Week",
    grepl("Class", name) ~ "Class"
  )) %>%
  mutate(name = str_remove(name, var)) %>%
  ggplot(aes(x = reorder(name, desc(name)), y = value, fill = value)) +
  geom_col(width = 1, color = "white") +
  coord_flip() +
  labs(
    y = "Coefficient",
    x = "Variable",
    title = "Effect on aphid abundance of land cover within 10km of suction trap sites",
    subtitle = paste("Model:", fit$call[2])
  ) +
  facet_wrap(~ var, scales = "free")

ggsave("land/model coefs - richness.png", type = "cairo")
