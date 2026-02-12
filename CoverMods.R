#Step 1. Load packages
library(tidyverse)
library(car)
library(ggeffects)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(mgcv)
library(emmeans)
library(splines)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Clean_Data.csv", sep = ";")

df$Fire_Int_Groups <- factor(
  df$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

#Figure 1

DataPlot <- ggplot(df, aes(x = Years_since_fire)) +
  geom_bar(aes(fill = Fire_Int_Groups), 
           position = position_dodge2(preserve = "single")) +
  facet_wrap(~Continent, 
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  theme_bw() +
  labs(
    x = "Years since fire",
    y = "Datapoints",
    fill = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(0,22), n.breaks= 12) +
  scale_y_continuous(limits = c(0, 10), n.breaks = 6) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.direction = "vertical",
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", colour = "NA")
  )

df_studies <- df %>%
  group_by(Years_since_fire, Fire_Int_Groups, Continent) %>%
  summarise(n_unique = n_distinct(StudyID), .groups = "drop")

StudyPlot <- ggplot(df_studies, aes(x = Years_since_fire,
                                    y = n_unique)) +
  geom_col(aes(fill = Fire_Int_Groups), 
           position = position_dodge2(preserve = "single")) +
  facet_wrap(~Continent, 
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  theme_bw() +
  labs(
    x = "Years since fire",
    y = "Studies",
    fill = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(0,22), n.breaks= 12) +
  scale_y_continuous(limits = c(0, 10), n.breaks = 6) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.direction = "vertical",
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", colour = "NA")
  )

StudyPlot

#Size investigated
df_sizes <- df %>%
  mutate(
    studysize = Plot_size * Sample_size
  )

df_sizes <- df_sizes %>%
  group_by(Years_since_fire, Fire_Int_Groups, Continent) %>%
  summarise(SumSize = sum(studysize), .groups = "drop")

df_sizes$Fire_Int_Groups <- factor(
  df_sizes$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

SizePlot <- ggplot(df_sizes, aes(x = Years_since_fire,
                                               y = log10(SumSize))) +
  geom_col(aes(fill = Fire_Int_Groups), 
           position = position_dodge2(preserve = "single")) +
  facet_wrap(~Continent, 
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  theme_bw() +
  labs(
    x = "Years since fire",
    y = expression(log10(Summed~area~studied~(m^2))),
    fill = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(0,22), n.breaks= 12) +
#  scale_y_continuous(limits = c(0, 10), n.breaks = 6) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.direction = "vertical",
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", colour = "NA")
  )

SizePlot

df_sizes2$Fire_Int_Groups <- factor(
  df_sizes2$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

IndSizePlot <- ggplot(df_sizes2, aes(x = Years_since_fire,
                                 y = log10(studysize))) +
  geom_jitter(aes(color = Fire_Int_Groups), size = 2.5, 
             height = 0.15, width = 0) +
  facet_wrap(~Continent, 
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  theme_bw() +
  labs(
    x = "Years since fire",
    y = expression(log10(Area~studied~(m^2))),
    fill = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(0,22), n.breaks= 12) +
  #  scale_y_continuous(limits = c(0, 10), n.breaks = 6) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.direction = "vertical",
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", colour = "NA")
  )

IndSizePlot

Figure1Plots <- (StudyPlot | DataPlot) / (SizePlot | IndSizePlot)

ggsave(Figure1Plots, filename = "Fig1.png",
       dpi = 300, height = 10.52, width = 17.3)

#####################################

metadata <- df %>%
  select(-contains("postfire"))

metadata <- metadata[,1:20]

species_raw <- df %>%
  select(StudyID, RowID, contains("postfire")) %>%
  select(
    StudyID, RowID,
    matches("_postfire$|_postfire_cover$")
  )


species_names <- species_raw %>%
  select(RowID, StudyID, matches("_postfire$")) %>%
  pivot_longer(
    cols = matches("_postfire$"),
    names_to = "base",
    names_pattern = "(.*)_postfire$",
    values_to = "species"
  )

species_cover <- species_raw %>%
  select(RowID, StudyID, matches("_postfire_cover$")) %>%
  pivot_longer(
    cols = matches("_postfire_cover$"),
    names_to = "base",
    names_pattern = "(.*)_postfire_cover$",
    values_to = "cover"
  )

species_long <- left_join(
  species_names,
  species_cover,
  by = c("RowID", "StudyID", "base")
)

species_long <- species_long %>%
  mutate(cover = as.numeric(unlist(cover)))

#Drop zeros
species_long <- subset(species_long, !species == 0)

#Add metadata
species_long_meta <- species_long %>%
  left_join(metadata, by = c("RowID", "StudyID"))


#Step 3. Load climate data

Perc_Data <- read.csv("site_climatological_percipitation_WorldClim_v2.csv",
                      sep = ";")

Temp_Data <- read.csv("site_climatological_temperature_WorldClim_v2.csv",
                      sep = ";")

species_long_Temp <- species_long_meta %>%
  left_join(Temp_Data, by = "Title")

species_long_Perc <- species_long_Temp %>%
  left_join(Perc_Data, by = "Title")

species_long_Perc <- species_long_Perc %>% 
  mutate(across(c("Fire_Int_Groups", "Continent", "StudyID"), as.factor))

#Create a plantgroup column 
df_long <- species_long_Perc %>%
  mutate(
    PlantGroup = str_extract(base, "(?<=_)[^_]+(?=_)")
  )

#Scale and center percipitation and temperature
df_long$Temp_sc <- as.numeric(scale(df_long$Avg_Temp, center = TRUE, scale = TRUE))
df_long$Per_sc <- as.numeric(scale(df_long$AvgPer, center = TRUE, scale = TRUE))


#Create a total size of study column
df_long <- df_long %>%
  mutate(
    studysize = Plot_size * Sample_size
  )

#
#Divide by 100 to get a 0-1 range but
#avoid 0 and 1 one in the dataset (beta must be >0 and <1)

df_long$coverstd <- (df_long$cover + 0.01) / 101

df_long$base <- as.factor(df_long$base)

#Use only first 10 years after fire
df_long_sub <- df_long %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

df_long_sub$species <- as.factor(df_long_sub$species)

####What species have data?

#

###################################################
#Dominant species per plant group plots. Herbs

herbs_df <- df_long_sub %>%
  filter(PlantGroup == "herb") %>%
  droplevels()

#Filter out all species that appear less than five times, in less than five timesteps,
#and less than two studies
herbs_filtered <- herbs_df %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 5, 
n_distinct(Years_since_fire) >= 5,
  n_distinct(StudyID) >=2)%>%
  ungroup() %>%
  droplevels()

herbs_EU <- herbs_filtered %>%
  filter(Continent == "Eurasia")

species <- herbs_EU %>%
  filter(Fire_Int_Groups == "Low")

ggplot(species, aes(x = Years_since_fire)) +
         geom_histogram(aes(fill = species))

mod <- glmmTMB(
  coverstd ~ Years_since_fire * species +
    Temp_sc+
    (1 | StudyID.x),
  dispformula = ~ Years_since_fire,
  weights = log(studysize),
  family = beta_family(),
  data = species
)

sim_res <- DHARMa::simulateResiduals(fittedModel = mod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, species$Years_since_fire)


mod <- gam(
  coverstd ~ 
    s(Years_since_fire, k = 5) +
    s(Temp_sc, k = 5) +
    s(StudyID, bs = "re"),   # random effect
  family = betar(link = "logit"),
  weights = log10(studysize),
  data = species,
  method = "REML"
)

gam.check(mod)
summary(mod)
sim <- simulateResiduals(mod)
plot(sim)
testUniformity(sim)
testDispersion(sim)

prediction_plot <- function(herbs_filtered, continent_level, fire_level) {
  
  # Filter data
  df_sub <- herbs_filtered %>%
    filter(Continent == continent_level,
           Fire_Int_Groups == fire_level) %>%
    mutate(species = as.factor(species)) %>%
    droplevels()
  
  # Number of species
  n_sp <- nlevels(df_sub$species)
  
  # Fit model depending on number of species
  if (n_sp > 1) {
    mod <- gam(
      coverstd ~ 
        species +
        s(Years_since_fire, species, bs = "fs", k = 5) +
        s(Temp_sc, k =5) +
        s(StudyID, bs = "re"),
      family = betar(link = "logit"),
      weights = log10(studysize),
      data = df_sub,
      method = "REML"
    )
  } else {
    # Single-species model
    mod <- gam(
      coverstd ~ 
        s(Years_since_fire, k = 5) +   # no species interaction
        s(Temp_sc, k = 5) +
        s(StudyID, bs = "re"),
      family = betar(link = "logit"),
      weights = log10(studysize),
      data = df_sub,
      method = "REML"
    )
  }
  
#  pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
 #   df_sp <- df_sub[df_sub$species == sp, ]
#    data.frame(
#      Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
#      species = sp,
#      Per_sc = mean(df_sp$Per_sc, na.rm = TRUE),
#      Temp_sc = mean(df_sp$Temp_sc, na.rm = TRUE),
#      StudyID.x = NA
#      studysize <- median(df_sub$studysize, na.rm = TRUE)
#    )
#  }))
 
  # Prediction grid
  pred_grid <- expand.grid(
    Years_since_fire = seq(
      min(df_sub$Years_since_fire, na.rm = TRUE),
      max(df_sub$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    species = levels(df_sub$species),
    Temp_sc = mean(df_sub$Temp_sc, na.rm = TRUE),
    StudyID = levels(df_sub$StudyID)[1],
    studysize = median(df_sub$studysize)
  )
   
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    exclude = "s(StudyID)"  # ignore random effect
  )
  
  pred_grid$fit <- mod$family$linkinv(pred$fit)
  pred_grid$upper <- mod$family$linkinv(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- mod$family$linkinv(pred$fit - 1.96 * pred$se.fit)

  
  # Plot
  p <- ggplot(pred_grid,
              aes(x = Years_since_fire,
                  y = fit,
                  color = species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = species),
                alpha = 0.2, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.4, aes(color = species)) +
    labs(
      x = "Years since fire",
      y = "Predicted cover",
      color = "Species",
      title = paste("Herbs", continent_level, fire_level, "intensity")
    ) +
    theme_bw() +
    scale_x_continuous(limits= c(1, 10), n.breaks = 6) +
    scale_color_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    scale_fill_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    theme(
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.direction = "vertical",
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text = element_text(size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5)
    )
  
  return(p)
}

herbsEUhigh_plot <- prediction_plot(herbs_filtered, "Eurasia", "High")
herbsEUmed_plot  <- prediction_plot(herbs_filtered, "Eurasia", "Medium")
herbsEUlow_plot  <- prediction_plot(herbs_filtered, "Eurasia", "Low")
herbsNAhigh_plot  <- prediction_plot(herbs_filtered, "North_America", "High")
herbsNAmed_plot  <- prediction_plot(herbs_filtered, "North_America", "Medium")
herbsNAlow_plot  <- prediction_plot(herbs_filtered, "North_America", "Low")

EUherbplots <- herbsEUhigh_plot / herbsEUmed_plot / herbsEUlow_plot
EUherbplots

NAherbplots <- herbsNAhigh_plot / herbsNAmed_plot / herbsNAlow_plot
NAherbplots
####
#Dwarfshrubs

dwarf_df <- df_long_sub %>%
  filter(PlantGroup == "dwarfshrub")

#Filter out all species that appear less than three times
dwarf_filtered <- dwarf_df %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 5, 
         n_distinct(Years_since_fire) >= 5,
         n_distinct(StudyID) >=2)%>%
  ungroup() %>%
  droplevels()

prediction_plot <- function(dwarf_filtered, continent_level, fire_level) {
  
  # Filter data
  df_sub <- dwarf_filtered %>%
    filter(Continent == continent_level,
           Fire_Int_Groups == fire_level) %>%
    mutate(species = as.factor(species)) %>%
    droplevels()
  
  # Fit model
  mod <- gam(
    coverstd ~ 
      species +
      s(Years_since_fire, species, bs = "fs", k = 5) +
      Temp_sc+
      s(StudyID, bs = "re"),   # random effect
    family = betar(link = "logit"),
    data = df_sub,
    method = "REML"
  )
  
  #  pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
  #   df_sp <- df_sub[df_sub$species == sp, ]
  #    data.frame(
  #      Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
  #      species = sp,
  #      Per_sc = mean(df_sp$Per_sc, na.rm = TRUE),
  #      Temp_sc = mean(df_sp$Temp_sc, na.rm = TRUE),
  #      StudyID.x = NA
  #      studysize <- median(df_sub$studysize, na.rm = TRUE)
  #    )
  #  }))
  
  # Prediction grid
  pred_grid <- expand.grid(
    Years_since_fire = seq(
      min(df_sub$Years_since_fire, na.rm = TRUE),
      max(df_sub$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    species = levels(df_sub$species),
    Temp_sc = mean(df_sub$Temp_sc, na.rm = TRUE),
    StudyID = levels(df_sub$StudyID)[1]
  )
  
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    exclude = "s(StudyID)"  # ignore random effect
  )
  
  pred_grid$fit <- mod$family$linkinv(pred$fit)
  pred_grid$upper <- mod$family$linkinv(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- mod$family$linkinv(pred$fit - 1.96 * pred$se.fit)
  
  
  # Plot
  p <- ggplot(pred_grid,
              aes(x = Years_since_fire,
                  y = fit,
                  color = species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = species),
                alpha = 0.2, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.2, aes(color = species)) +
    labs(
      x = "Years since fire",
      y = "Predicted cover",
      color = "Species",
      title = paste("Dwarfshrub", continent_level, fire_level, "intensity")
    ) +
    theme_bw() +
    scale_x_continuous(limits= c(1, 10), n.breaks = 6) +
    scale_color_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    scale_fill_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    theme(
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.direction = "vertical",
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text = element_text(size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5)
    )
  
  return(p)
}

dwarfEUhigh_plot <- prediction_plot(dwarf_filtered, "Eurasia", "High")
dwarfEUmed_plot  <- prediction_plot(dwarf_filtered, "Eurasia", "Medium")
dwarfEUlow_plot  <- prediction_plot(dwarf_filtered, "Eurasia", "Low")
dwarfNAhigh_plot  <- prediction_plot(dwarf_filtered, "North_America", "High")
dwarfNAmed_plot  <- prediction_plot(dwarf_filtered, "North_America", "Medium")
dwarfNAlow_plot  <- prediction_plot(dwarf_filtered, "North_America", "Low")

EUdwarfplots <- dwarfEUhigh_plot / dwarfEUmed_plot / dwarfEUlow_plot
EUdwarfplots

NAdwarfplots <- dwarfNAhigh_plot / dwarfNAmed_plot / dwarfNAlow_plot
NAdwarfplots

####
#Grasses

grass_df <- df_long_sub %>%
  filter(PlantGroup == "graminoid")

#Filter out all species that appear less than three times
grass_filtered <- grass_df %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 5, 
         n_distinct(Years_since_fire) >= 5,
         n_distinct(StudyID) >=2)%>%
  ungroup() %>%
  droplevels()

prediction_plot <- function(grass_filtered, continent_level, fire_level) {
  
  # Filter data
  df_sub <- grass_filtered %>%
    filter(Continent == continent_level,
           Fire_Int_Groups == fire_level) %>%
    mutate(species = as.factor(species)) %>%
    droplevels()
  
  # Number of species
  n_sp <- nlevels(df_sub$species)
  
  # Fit model depending on number of species
  if (n_sp > 1) {
    mod <- gam(
      coverstd ~ 
        species +
        s(Years_since_fire, species, bs = "fs", k = 5) +
        Temp_sc +
        s(StudyID, bs = "re"),
      family = betar(link = "logit"),
      data = df_sub,
      method = "REML"
    )
  } else {
    # Single-species model
    mod <- gam(
      coverstd ~ 
        s(Years_since_fire, k = 5) +   # no species interaction
        Temp_sc +
        s(StudyID, bs = "re"),
      family = betar(link = "logit"),
      data = df_sub,
      method = "REML"
    )
  }
  
  #  pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
  #   df_sp <- df_sub[df_sub$species == sp, ]
  #    data.frame(
  #      Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
  #      species = sp,
  #      Per_sc = mean(df_sp$Per_sc, na.rm = TRUE),
  #      Temp_sc = mean(df_sp$Temp_sc, na.rm = TRUE),
  #      StudyID.x = NA
  #      studysize <- median(df_sub$studysize, na.rm = TRUE)
  #    )
  #  }))
  
  # Prediction grid
  pred_grid <- expand.grid(
    Years_since_fire = seq(
      min(df_sub$Years_since_fire, na.rm = TRUE),
      max(df_sub$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    species = levels(df_sub$species),
    Temp_sc = mean(df_sub$Temp_sc, na.rm = TRUE),
    StudyID = levels(df_sub$StudyID)[1]
  )
  
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    exclude = "s(StudyID)"  # ignore random effect
  )
  
  pred_grid$fit <- mod$family$linkinv(pred$fit)
  pred_grid$upper <- mod$family$linkinv(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- mod$family$linkinv(pred$fit - 1.96 * pred$se.fit)
  
  
  # Plot
  p <- ggplot(pred_grid,
              aes(x = Years_since_fire,
                  y = fit,
                  color = species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = species),
                alpha = 0.2, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.2, aes(color = species)) +
    labs(
      x = "Years since fire",
      y = "Predicted cover",
      color = "Species",
      title = paste("Graminoids", continent_level, fire_level, "intensity")
    ) +
    theme_bw() +
    scale_x_continuous(limits= c(1, 10), n.breaks = 6) +
    scale_color_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    scale_fill_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    theme(
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.direction = "vertical",
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text = element_text(size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5)
    )
  
  return(p)
}

grassEUhigh_plot <- prediction_plot(grass_filtered, "Eurasia", "High")
grassEUmed_plot  <- prediction_plot(grass_filtered, "Eurasia", "Medium")
grassEUlow_plot  <- prediction_plot(grass_filtered, "Eurasia", "Low")
grassNAhigh_plot  <- prediction_plot(grass_filtered, "North_America", "High")
grassNAmed_plot  <- prediction_plot(grass_filtered, "North_America", "Medium")
grassNAlow_plot  <- prediction_plot(grass_filtered, "North_America", "Low")

EUgrassplots <- grassEUhigh_plot / grassEUmed_plot / grassEUlow_plot
EUgrassplots

NAgrassplots <- grassNAhigh_plot / grassNAmed_plot / grassNAlow_plot
NAgrassplots

####
#Bryophytes
moss_df <- df_long_sub %>%
  filter(PlantGroup == "bryophyte")

#Filter out all species that appear less than three times
moss_filtered <- moss_df %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 5, 
         n_distinct(Years_since_fire) >= 5,
         n_distinct(StudyID) >=2)%>%
  ungroup() %>%
  droplevels()

prediction_plot <- function(moss_filtered, continent_level, fire_level) {
  
  # Filter data
  df_sub <- moss_filtered %>%
    filter(Continent == continent_level,
           Fire_Int_Groups == fire_level) %>%
    mutate(species = as.factor(species)) %>%
    droplevels()
  
  # Number of species
  n_sp <- nlevels(df_sub$species)
  
  # Fit model depending on number of species
  if (n_sp > 1) {
    mod <- gam(
      coverstd ~ 
        species +
        s(Years_since_fire, species, bs = "fs", k = 5) +
        Temp_sc +
        s(StudyID, bs = "re"),
      family = betar(link = "logit"),
      data = df_sub,
      method = "REML"
    )
  } else {
    # Single-species model
    mod <- gam(
      coverstd ~ 
        s(Years_since_fire, k = 5) +   # no species interaction
        Temp_sc +
        s(StudyID, bs = "re"),
      family = betar(link = "logit"),
      data = df_sub,
      method = "REML"
    )
  }
  
  # Prediction grid
  pred_grid <- expand.grid(
    Years_since_fire = seq(
      min(df_sub$Years_since_fire, na.rm = TRUE),
      max(df_sub$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    species = levels(df_sub$species),
    Temp_sc = mean(df_sub$Temp_sc, na.rm = TRUE),
    StudyID = levels(df_sub$StudyID)[1]
  )
  
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    exclude = "s(StudyID)"  # ignore random effect
  )
  
  pred_grid$fit <- mod$family$linkinv(pred$fit)
  pred_grid$upper <- mod$family$linkinv(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- mod$family$linkinv(pred$fit - 1.96 * pred$se.fit)
  
  
  # Plot
  p <- ggplot(pred_grid,
              aes(x = Years_since_fire,
                  y = fit,
                  color = species)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = species),
                alpha = 0.2, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.4, aes(color = species)) +
    labs(
      x = "Years since fire",
      y = "Predicted cover",
      color = "Species",
      title = paste("Bryophytes", continent_level, fire_level, "intensity")
    ) +
    theme_bw() +
    scale_x_continuous(limits= c(1, 10), n.breaks = 6) +
    scale_color_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    scale_fill_manual(values = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    ))+
    theme(
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.direction = "vertical",
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text = element_text(size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5)
    )
  
  return(p)
}

mossEUhigh_plot <- prediction_plot(moss_filtered, "Eurasia", "High")
mossEUmed_plot  <- prediction_plot(moss_filtered, "Eurasia", "Medium")
mossEUlow_plot  <- prediction_plot(moss_filtered, "Eurasia", "Low")
mossNAhigh_plot  <- prediction_plot(moss_filtered, "North_America", "High")
mossNAmed_plot  <- prediction_plot(moss_filtered, "North_America", "Medium")
mossNAlow_plot  <- prediction_plot(moss_filtered, "North_America", "Low")

EUmossplots <- mossEUhigh_plot / mossEUmed_plot / mossEUlow_plot
EUmossplots

NAmossplots <- mossNAhigh_plot / mossNAmed_plot
NAmossplots


ggsave(EUherbplots, filename = "EUHerbs.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(NAherbplots, filename = "NAHerbs.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(NAdwarfplots, filename = "NAdwarf.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(EUdwarfplots, filename = "EUDwarf.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(NAgrassplots, filename = "NAGrass.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(EUgrassplots, filename = "EUGrasss.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(NAmossplots, filename = "NAMoss.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(EUmossplots, filename = "EUMoss.png",
       dpi = 300, height = 10.56, width = 13)






####################################################
#Most dominant species mods
#Remake to include all three species by summing the cover of all

#Normalize studysize for better model fit
df_long$weight_sc <- df_long$studysize / mean(summed_cover$studysize)


summed_cover <- df_long %>%
  group_by(across(-c(base, species, cover))) %>%
  summarize(
    TotalCover = sum(cover),
    .groups = "drop"
  )


#What if we include all three levels,

#Try beta regression by turning summed cover to 0-1

summed_cover$CoverStd <- summed_cover$TotalCover/300

df_long$base <- as.factor(df_long$base)

#Herbs
herbs_dom <- df_long %>%
  filter(PlantGroup == "herb")

# Scale Cover to proportion (0–1) and avoid exact 0/1
herbs_dom <- herbs_dom %>%
  mutate(
    Cover_prop = cover / 100.01       # adjust 0/100 for Beta regression
  )

range(herbs_dom$Cover_prop)

#Max 10 years 

herbs_dom <- herbs_dom %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

herbmod <- glmmTMB(
  Cover_prop ~ 
    ns(Years_since_fire, df = 4) * Fire_Int_Groups +
    ns(Temp_sc, df = 3) +
    Per_sc +
    (1 | base) +   # random intercept for dominance rank
    (1 | StudyID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = herbs_dom,
  weights = weight_sc
)

AIC_vals <- AIC(herbmod, herbmod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

res <- resid(herbmod, type = "pearson")
fitted_vals <- fitted(herbmod)
qqnorm(res)
qqline(res, col = "red")

plot(fitted(herbmod), resid(herbmod, type="pearson"))
abline(h=0, col="red", lty=2)

pearson_resid <- resid(herbmod, type="pearson")
df_resid <- df.residual(herbmod)
sum(pearson_resid^2) / df_resid

ranef(herbmod)

summary(herbmod)
Anova(herbmod, type = 'III')


#Plot predictions!
herb_pred <- emmeans(
  herbmod,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(herbs_dom$Years_since_fire, na.rm = TRUE),
      max(herbs_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(herbs_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(herbs_dom$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

###NEW####
pred_grid <- as.data.frame(herb_pred) %>%
  mutate(
    fit_cover   = emmean * 100,   # back-transform to cover %
    lower_cover = lower.CL * 100,
    upper_cover = upper.CL * 100
  )
#####

#########NEW###

# Sum predicted cover for top 3 species per combination
pred_summed <- pred_grid %>%
  group_by(Years_since_fire, Fire_Int_Groups) %>%
  summarise(
    fit = sum(fit_cover),
    lower = sum(lower_cover),
    upper = sum(upper_cover),
    .groups = "drop"
  )

############

pred_grid_herb$Fire_Int_Groups <- factor(
  pred_grid_herb$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predherbplot <- ggplot(pred_summed,
       aes(x = Years_since_fire,
           y = fit,
           color = Fire_Int_Groups,
           group = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted summed cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Herbs")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 


predherbplot

ggsave(plot = predherbplot, filename = "Pred_herb_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#DWARFSHRUBS

dwarf_dom <- summed_cover %>%
  filter(PlantGroup == "dwarfshrub")

ggplot(dwarf_dom, aes(x = Years_since_fire)) +
  geom_histogram(aes(fill = Fire_Int_Groups)) +
  facet_wrap(~Continent)

dwarfmod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups * Continent +
    poly(Temp_sc, 2) +
    poly(Per_sc, 3),
  data = dwarf_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varPower(form = ~ weight_sc),
  method      = "REML"
)

AIC_vals <- AIC(dwarfmod, dwarfmod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(dwarfmod, weights = varFixed(~ 1 / weight_sc), method = "ML")
mod_power <- update(dwarfmod, weights = varPower(form = ~ weight_sc), method = "ML")
mod_exp   <- update(dwarfmod, weights = varExp(form = ~ weight_sc), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(dwarfmod, type = "normalized"))
qqline(resid(dwarfmod, type = "normalized"))
plot(dwarfmod, resid(., type = "normalized") ~ fitted(.))
summary(dwarfmod)
Anova(dwarfmod, type = 'III')

#Plot predictions!
dwarf_pred <- emmeans(
  dwarfmod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(dwarf_dom$Years_since_fire, na.rm = TRUE),
      max(dwarf_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(dwarf_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(dwarf_dom$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_dwarf <- as.data.frame(dwarf_pred)

pred_grid_dwarf <- pred_grid_dwarf %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_dwarf$Fire_Int_Groups <- factor(
  pred_grid_dwarf$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

preddwarfplot <- ggplot(pred_grid_dwarf,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted summed cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Dwarfshrubs")+
  scale_y_continuous(limits = c(0,100)) +
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 


preddwarfplot

ggsave(plot = preddwarfplot, filename = "Pred_dwarf_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#Graminoids

grass_dom <- summed_cover %>%
  filter(PlantGroup == "graminoid")

ggplot(grass_dom, aes(x = Years_since_fire)) +
  geom_histogram(aes(fill = Fire_Int_Groups)) +
  facet_wrap(~Continent)

grassmod <- gls(
  log(TotalCover) ~
    Years_since_fire * Fire_Int_Groups * Continent +
    Temp_sc +
    poly(Per_sc, 3),
  data = grass_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varPower(form = ~ weight_sc),
  method      = "REML"
)

AIC_vals <- AIC(grassmod, grassmod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(grassmod, weights = varFixed(~ 1 / weight_sc), method = "ML")
mod_power <- update(grassmod, weights = varPower(form = ~ weight_sc), method = "ML")
mod_exp   <- update(grassmod, weights = varExp(form = ~ weight_sc), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(grassmod, type = "normalized"))
qqline(resid(grassmod, type = "normalized"))
plot(grassmod, resid(., type = "normalized") ~ fitted(.))
summary(grassmod)
Anova(grassmod, type = 'III')

#Plot predictions!
grass_pred <- emmeans(
  grassmod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(grass_dom$Years_since_fire, na.rm = TRUE),
      max(grass_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(grass_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(grass_dom$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_grass <- as.data.frame(grass_pred)

pred_grid_grass <- pred_grid_grass %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_grass$Fire_Int_Groups <- factor(
  pred_grid_grass$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predgramplot <- ggplot(pred_grid_grass,
                        aes(x = Years_since_fire,
                            y = fit,
                            color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted summed cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Graminoids")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

predgramplot

ggsave(plot = predgramplot, filename = "Pred_gram_plot.png", dpi =300,
       height = 4.2, width = 6.5)

####
#Trees

tree_dom <- summed_cover %>%
  filter(PlantGroup == "tree")

tree_dom <- tree_dom %>%
  filter(!is.na(TotalCover))

ggplot(tree_dom, aes(x = Years_since_fire)) +
  geom_histogram(aes(fill = Fire_Int_Groups)) +
  facet_wrap(~Continent)

treemod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 2) +
    Fire_Int_Groups * Continent +
    poly(Temp_sc, 3),
  data = tree_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varExp(form = ~ weight_sc),
  method      = "REML"
)


AIC_vals <- AIC(treemod, treemod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(treemod, weights = varFixed(~ 1 / weight_sc), method = "ML")
mod_power <- update(treemod, weights = varPower(form = ~ weight_sc), method = "ML")
mod_exp   <- update(treemod, weights = varExp(form = ~ weight_sc), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(treemod, type = "normalized"))
qqline(resid(treemod, type = "normalized"))
plot(treemod, resid(., type = "normalized") ~ fitted(.))
summary(treemod)
Anova(treemod, type = 'III')

#Plot predictions!
tree_pred <- emmeans(
  treemod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(tree_dom$Years_since_fire, na.rm = TRUE),
      max(tree_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
#    Per_sc = mean(tree_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(tree_dom$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_tree <- as.data.frame(tree_pred)

pred_grid_tree <- pred_grid_tree %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_tree$Fire_Int_Groups <- factor(
  pred_grid_tree$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predtreeplot <- ggplot(pred_grid_tree,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted summed cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Trees")+
  scale_y_continuous(limits = c(0,200)) +
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA"))

predtreeplot

ggsave(plot = predtreeplot, filename = "Pred_tree_plot.png", dpi =300,
       height = 4.2, width = 6.5)


#####
#Shrubs

shrubs_dom <- df_long %>%
  filter(base == "Dominant_shrub_1")

shrubs_dom <- shrubs_dom %>%
  filter(!is.na(coverstd))

#Only data from one study in Europe. Use only NA for shrub model
shrubs_dom <- shrubs_dom %>%
  filter(! Continent %in% "Eurasia")

shrubmod <- glmmTMB(
  coverstd ~
    Fire_Int_Groups * Years_since_fire,
  family = beta_family(),
  # dispformula = ~ Fire_Int_Groups,
  data = shrubs_dom
)

sim_res <- simulateResiduals(fittedModel = shrubmod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, shrubs_dom$Years_since_fire)
plotResiduals(sim_res, shrubs_dom$Fire_Int_Groups)


summary(shrubmod)
Anova(shrubmod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(shrubs_dom$Years_since_fire, na.rm = TRUE),
    max(shrubs_dom$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(shrubs_dom$Fire_Int_Groups)
) 

#type = response ok?
pred <- predict(
  shrubmod,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE,
  re.form = NA,
  allow.new.levels = TRUE
)

eta <- pred$fit
se  <- pred$se.fit

lower_eta <- eta - 1.96 * se
upper_eta <- eta + 1.96 * se

# inverse link = plogis for beta regression
pred_grid$fit   <- plogis(eta)
pred_grid$lower <- plogis(lower_eta)
pred_grid$upper <- plogis(upper_eta)

pred_grid$Fire_Int_Groups <- factor(
  pred_grid$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predshrubplot<- ggplot(pred_grid,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted shrub cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Shrubs")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

predshrubplot

ggsave(plot = predshrubplot, filename = "Pred_shrub_plot.png", dpi =300,
       height = 4.2, width = 6.5)

#####
#Mosses
bryo_dom <- df_long %>%
  filter(base == "Dominant_bryophyte_1")

bryo_dom <- bryo_dom %>%
  filter(!is.na(coverstd))

bryomod <- glmmTMB(
  coverstd ~
    Fire_Int_Groups * Years_since_fire +
    Fire_Int_Groups * Continent +
    Years_since_fire * Continent,
  family = beta_family(),
  dispformula = ~ Years_since_fire,
  data = bryo_dom
)

sim_res <- simulateResiduals(fittedModel = bryomod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, bryo_dom$Years_since_fire)
plotResiduals(sim_res, bryo_dom$Fire_Int_Groups)
plotResiduals(sim_res, bryo_dom$Continent)

summary(bryomod)
Anova(bryomod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(bryo_dom$Years_since_fire, na.rm = TRUE),
    max(bryo_dom$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(bryo_dom$Fire_Int_Groups),
  Continent = levels(bryo_dom$Continent)
) 

#type = response ok?
pred <- predict(
  bryomod,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE,
  re.form = NA,
  allow.new.levels = TRUE
)

eta <- pred$fit
se  <- pred$se.fit

lower_eta <- eta - 1.96 * se
upper_eta <- eta + 1.96 * se

# inverse link = plogis for beta regression
pred_grid$fit   <- plogis(eta)
pred_grid$lower <- plogis(lower_eta)
pred_grid$upper <- plogis(upper_eta)

pred_grid$Fire_Int_Groups <- factor(
  pred_grid$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predmossplot<- ggplot(pred_grid,
                      aes(x = Years_since_fire,
                          y = fit,
                          color = Fire_Int_Groups)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted moss cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Bryophytes")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

predmossplot

ggsave(plot = predmossplot, filename = "Pred_moss_plot.png", dpi =300,
       height = 4.2, width = 6.5)

#Combinedplot

combinedplot <- (predtreeplot|predshrubplot)/(predherbplot|preddwarfplot)/(predgramplot|predmossplot)
combinedplot

ggsave(plot=combinedplot, filename = "coverplots_combined.png", dpi =300,
       height = 10.52, width = 13)


#Where models come to die

#Best gls models on summed cover 0-100 for all three dominance levels
herbmod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups * Continent +
    poly(Temp_sc, 3) +
    poly(Per_sc, 2),
  data = herbs_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varExp(form = ~ weight_sc),
  method      = "REML"
)