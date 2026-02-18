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
    y = expression(Summed~area~studied~(m^2)),
    fill = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(0,22), n.breaks= 12) +
  scale_y_continuous(
    labels = scales::label_math(10^.x)
  )+
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
    y = expression(Area~studied~(m^2)),
    fill = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(0,21), n.breaks= 12) +
  scale_y_continuous(
    labels = scales::label_math(10^.x)
  )+
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



Figure1Plots <- StudyPlot / DataPlot / SizePlot / IndSizePlot
Figure1Plots
ggsave(Figure1Plots, filename = "Fig1.png",
       dpi = 300, height = 12.62, width = 13)

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

#Use only first 10 years after fire
df_long <- df_long %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

#
#Divide by 100 to get a 0-1 range but
#avoid 0 and 1 one in the dataset (beta must be >0 and <1)

df_long$coverstd <- (df_long$cover + 0.01) / 101

df_long$base <- as.factor(df_long$base)


df_long$species <- as.factor(df_long$species)

####What species have data?

#

###################################################
#Dominant species per plant group plots. Herbs

herbs_df <- df_long %>%
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
  weights = weights_sc,
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
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire * species +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
      )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
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
    weight_sc = median(df_sub$weight_sc)
  )
   
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    re.form = NA  # ignore random effect
  )
  
  pred_grid$fit   <- plogis(pred$fit)
  pred_grid$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  
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
  
  # Number of species
  n_sp <- nlevels(df_sub$species)
  
  # Fit model depending on number of species
  if (n_sp > 1) {
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire * species +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
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
    weight_sc = median(df_sub$weight_sc)
  )
  
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    re.form = NA  # ignore random effect
  )
  
  pred_grid$fit   <- plogis(pred$fit)
  pred_grid$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  
  
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
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire * species +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
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
    weight_sc = median(df_sub$weight_sc)
  )
  
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    re.form = NA  # ignore random effect
  )
  
  pred_grid$fit   <- plogis(pred$fit)
  pred_grid$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  
  
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
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire * species +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        Temp_sc +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = weight_sc,
      data = df_sub
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
    weight_sc = median(df_sub$weight_sc)
  )
  
  # Predictions
  pred <- predict(
    mod,
    newdata = pred_grid,
    type = "link",
    se.fit = TRUE,
    re.form = NA  # ignore random effect
  )
  
  pred_grid$fit   <- plogis(pred$fit)
  pred_grid$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  pred_grid$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  
  
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
df_long$area_sc <- df_long$studysize / mean(df_long$studysize)

df_long$base <- as.factor(df_long$base)

#Herbs

herbs_dom1 <- df_long %>%
  filter(base == "Dominant_herb_1")
herbs_dom2 <- df_long %>%
  filter(base == "Dominant_herb_2")
herbs_dom3 <- df_long %>%
  filter(base == "Dominant_herb_3")

#herbs_dom <- herbs_dom %>%
#  mutate(across(where(is.factor), droplevels)) %>%   # remove unused factor levels
#  mutate(species = as.character(species)) %>%
#  group_by(RowID) %>%
#  complete(
#    base = base,
#    fill = list(cover = 0, coverstd = 0, species = "NULL")
#  ) %>%
#  ungroup()

ggplot(herbs_dom1, aes(x = Years_since_fire, y = coverstd)) +
  stat_smooth(method = 'lm', aes(color = Fire_Int_Groups)) +
  geom_jitter(aes(color = Fire_Int_Groups)) + 
  facet_wrap(~ Continent)

herbmod_1 <- glmmTMB(
  coverstd ~ 
    ns(Years_since_fire, df = 2) * Continent +
    Fire_Int_Groups +
    ns(Temp_sc, df = 3) +
    (1 | StudyID/RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = herbs_dom1
)

herbmod2_1 <- glmmTMB(
  coverstd ~ 
    ns(Years_since_fire, df = 2) * Continent +
    Fire_Int_Groups +
    ns(Temp_sc, df = 3) +
    (1 | StudyID/RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = herbs_dom1
)


AIC_vals <- AIC(herbmod_1, herbmod2_1)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(herbmod)
plot(simres)

herbs_dom1$pred <- predict(herbmod_1, type = "response")
plot(herbs_dom1$pred, herbs_dom1$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(herbmod_1)
Anova(herbmod_1, type = 'III')


#Plot predictions!
herb_pred_1 <- emmeans(
  herbmod_1,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(herbs_dom1$Years_since_fire, na.rm = TRUE),
      max(herbs_dom1$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(herbs_dom$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_herb_1 <- as.data.frame(summary(herb_pred_1, infer = TRUE))

pred_grid_herb_1 <- pred_grid_herb_1 %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
  )

pred_grid_herb_1$Fire_Int_Groups <- factor(
  pred_grid_herb_1$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predherbplot <- ggplot(pred_grid_herb_1,
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
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c("Eurasia" = "Eurasia",
                             "North_America" = "North America"))) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Herbs")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_blank(),
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
dwarf_dom <- df_long %>%
  filter(PlantGroup == "dwarfshrub")

ggplot(dwarf_dom, aes(x = Years_since_fire, y = coverstd)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) + 
  facet_grid(Continent ~base)

dwarf_mod <- glmmTMB(
  coverstd ~ 
    Fire_Int_Groups * base +
    Years_since_fire +
    Continent +
    ns(Per_sc, df = 3) +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = dwarf_dom,
  dispformula = ~ Fire_Int_Groups + base + Years_since_fire,
  weights = weight_sc
)

dwarf_mod2 <- glmmTMB(
  coverstd ~ 
    Fire_Int_Groups *  base +
    Years_since_fire +
    Continent +
    ns(Per_sc, df = 3) +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = dwarf_dom,
    dispformula = ~ Fire_Int_Groups + base + Years_since_fire,
#  weights = weight_sc
)

AIC_vals <- AIC(dwarf_mod, dwarf_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(dwarf_mod)
plot(simres)

dwarf_dom$pred_marginal <- predict(
  dwarf_mod,
  type = "response",
  re.form = NA   # removes random effects
)

dwarf_dom$pred <- predict(dwarf_mod, type = "response")
plot(dwarf_dom$pred_marginal, dwarf_dom$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(dwarf_mod)
Anova(dwarf_mod, type = 'III')

predict(dwarf_mod, type="response", re.form=NULL)

predict(dwarf_mod, type="response", re.form=NA)

#Plot predictions!
dwarf_pred <- emmeans(
  dwarf_mod,
  ~ Years_since_fire | base * Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(dwarf_dom$Years_since_fire, na.rm = TRUE),
      max(dwarf_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(dwarf_dom$Per_sc, na.rm = TRUE)
  )
)

pred_grid_dwarf <- as.data.frame(summary(dwarf_pred, infer = TRUE))

pred_grid_dwarf <- pred_grid_dwarf %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
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
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  facet_grid(Continent ~ base,
              labeller = labeller(
                base = c("Dominant_dwarfshrub_1" = "Most dominant sp.",
                         "Dominant_dwarfshrub_2" = "Second most dominant sp.",
                         "Dominant_dwarfshrub_3" = "Third most dominant sp."))) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Dwarfshrubs")+
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_blank(),
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

grass_dom <-df_long %>%
  filter(PlantGroup == "graminoid")

#Very few 2nd & 3rd dominant. Combine to one

grass_dom$Dominance <- fct_collapse(
  grass_dom$base,
  Dominant_graminoid_1 = "Dominant_graminoid_1",
  Dominant_graminoid_2 = c("Dominant_graminoid_2", "Dominant_graminoid_3")
)

ggplot(grass_dom, aes(x = Years_since_fire, y = coverstd)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) + 
  facet_grid(Continent~Dominance)

grass_mod <- glmmTMB(
  coverstd ~ 
    Years_since_fire +
    Fire_Int_Groups * Dominance +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  dispformula = ~ Fire_Int_Groups + Dominance + Years_since_fire,
  data = grass_dom,
  weights = weight_sc
)

grass_mod2 <- glmmTMB(
  coverstd ~ 
    Years_since_fire +
    Fire_Int_Groups * Dominance +
    Temp_sc +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = grass_dom,
  dispformula = ~ Fire_Int_Groups + Dominance + Years_since_fire,
  weights = weight_sc
)


AIC_vals <- AIC(grass_mod, grass_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(grass_mod)
plot(simres)

grass_dom$pred <- predict(grass_mod, type = "response")
plot(grass_dom$pred, grass_dom$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(grass_mod)
Anova(grass_mod, type = 'III')


#Plot predictions!
grass_pred <- emmeans(
  grass_mod,
  ~ Years_since_fire | Dominance* Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(grass_dom$Years_since_fire, na.rm = TRUE),
      max(grass_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
   # Per_sc = mean(grass_dom$Per_sc, na.rm = TRUE)
        Temp_sc = mean(grass_dom$Temp_sc, na.rm = TRUE)
  )
)

pred_grid_grass <- as.data.frame(summary(grass_pred, infer = TRUE))

pred_grid_grass <- pred_grid_grass %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
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
    scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
    scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Years since fire",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  facet_wrap(~ Dominance,
             labeller = labeller(
               Dominance = c("Dominant_graminoid_1" = "Most dominant sp.",
                        "Dominant_graminoid_2" = "Second most dominant sp.")
)) +
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

tree_dom <- df_long %>%
  filter(PlantGroup == "tree")

tree_dom <- tree_dom %>%
  filter(!is.na(cover))

ggplot(tree_dom, aes(x = Years_since_fire, y = coverstd)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) + 
  facet_grid(Continent ~ base)



#Europe has shit data for eerything but most dominant in low and medium int.
#Can we subset the dataset somehow to account for this?
#Let's start with just NA

Tree_NA <- tree_dom %>%
  filter(Continent == "North_America")

tree_modNA <- glmmTMB(
  coverstd ~ 
    Years_since_fire * Fire_Int_Groups +
    Years_since_fire * base +
    Fire_Int_Groups * base +
    Per_sc +
    Temp_sc +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  dispformula = ~ Years_since_fire + base,
  data = Tree_NA,
  weights = weight_sc
)

AIC_vals <- AIC(tree_mod, tree_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(tree_mod)
plot(simres)

Tree_NA$pred <- predict(tree_modNA, type = "response")
plot(Tree_NA$pred, Tree_NA$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(tree_mod)
Anova(tree_modNA, type = 'III')


#Plot predictions!
tree_pred_NA <- emmeans(
  tree_modNA,
  ~ Years_since_fire | base * Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(Tree_NA$Years_since_fire, na.rm = TRUE),
      max(Tree_NA$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
     Per_sc = mean(Tree_NA$Per_sc, na.rm = TRUE),
    Temp_sc = mean(Tree_NA$Temp_sc, na.rm = TRUE)
  )
)

pred_grid_tree_NA <- as.data.frame(summary(tree_pred_NA, infer = TRUE))

pred_grid_tree_NA <- pred_grid_tree_NA %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
  )

pred_grid_tree_NA$Fire_Int_Groups <- factor(
  pred_grid_tree_NA$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predtreeplot_NA <- ggplot(pred_grid_tree_NA,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  facet_wrap(~ base,
             labeller = labeller(
               base = c("Dominant_tree_1" = "Most dominant sp.",
                             "Dominant_tree_2" = "Second most dominant sp.",
                        "Dominant_tree_3" = "Third most dominant sp.",
                        "Dominant_tree_4" = "Fourth most dominant sp."))) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Trees North America")+
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA"))

predtreeplot_NA


Tree_EU <- tree_dom %>%
  filter(Continent == "Eurasia")

#Use only two levels of dominance
Tree_EU$Dominance <- fct_collapse(
  Tree_EU$base,
  Dominant_tree_1 = "Dominant_tree_1",
  Dominant_tree_2 = c("Dominant_tree_2", "Dominant_tree_3",
                      "Dominant_tree_4")
)

#Combine low and medium intensity
Tree_EU$Fire_Int_2 <- fct_collapse(
  Tree_EU$Fire_Int_Groups,
  High = "High",
  `Medium-Low` = c("Medium", "Low")
)

ggplot(Tree_EU, aes(x = Years_since_fire, y = coverstd)) +
  geom_smooth(aes(color = Fire_Int_2)) +
  geom_point(aes(color = Fire_Int_2)) + 
  facet_wrap(~ Dominance)

tree_mod <- glmmTMB(
  coverstd ~ 
    Years_since_fire +
    Fire_Int_2 * Dominance +
    ns(Per_sc, df = 3) +
    ns(Temp_sc,df = 3) +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  dispformula = ~ Dominance + Per_sc + Temp_sc,
  data = Tree_EU,
  weights = weight_sc
)

tree_mod2 <- glmmTMB(
  coverstd ~ 
    Years_since_fire +
    Fire_Int_2 * Dominance +
    ns(Per_sc, df = 3) +
    ns(Temp_sc, df = 3) +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  dispformula = ~ Dominance + Per_sc + Temp_sc,
  data = Tree_EU,
  weights = weight_sc
)

AIC_vals <- AIC(tree_mod, tree_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(tree_mod)
plot(simres)

Tree_EU$pred <- predict(tree_mod, type = "response")
plot(Tree_EU$pred, Tree_EU$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(tree_mod)
Anova(tree_mod, type = 'III')


#Plot predictions!
tree_pred <- emmeans(
  tree_mod,
  ~ Years_since_fire | Dominance * Fire_Int_2,
  at = list(
    Years_since_fire = seq(
      min(Tree_EU$Years_since_fire, na.rm = TRUE),
      max(Tree_EU$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(Tree_EU$Per_sc, na.rm = TRUE),
    Temp_sc = mean(Tree_EU$Temp_sc, na.rm = TRUE)
  )
)

pred_grid_tree <- as.data.frame(summary(tree_pred, infer = TRUE))

pred_grid_tree <- pred_grid_tree %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
  )

pred_grid_tree$Fire_Int_2 <- factor(
  pred_grid_tree$Fire_Int_2,
  levels = c("High", "Medium-Low")
)

predtreeplot_EU <- ggplot(pred_grid_tree,
                          aes(x = Years_since_fire,
                              y = fit,
                              color = Fire_Int_2)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_2),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "goldenrod")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  facet_wrap(~ Dominance,
             labeller = labeller(
               Dominance = c("Dominant_tree_1" = "Most dominant sp.",
                        "Dominant_tree_2" = "Second most dominant sp."))) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Trees Europe")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA"))

predtreeplot_EU


ggsave(plot = predtreeplot_EU, filename = "Pred_tree_plot_EU.png", dpi =300,
       height = 4.2, width = 6.5)


#####
#Shrubs
#Model is shit, can't be used.

shrubs_dom <- df_long %>%
  filter(PlantGroup == "shrub")

shrubs_dom <- shrubs_dom %>%
  filter(!is.na(coverstd))

#Only data from one study in Europe. Use only NA for shrub model
shrubs_dom <- shrubs_dom %>%
  filter(! Continent %in% "Eurasia")

ggplot(shrubs_dom, aes(x = Years_since_fire, y = coverstd)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  facet_wrap(~ Dominance)

#Remove Medium intensity, we don't have enough data
shrubs_dom <- shrubs_dom %>%
  filter(!Fire_Int_Groups %in% "Medium")

#Very few 2nd & 3rd dominant. Combine to one

shrubs_dom$Dominance <- fct_collapse(
  shrubs_dom$base,
  Dominant_shrub_1 = "Dominant_shrub_1",
  Dominant_shrub_2 = c("Dominant_shrub_2", "Dominant_shrub_3")
)

#Actually just use most dominant
shrubs_dom <- shrubs_dom %>%
  filter(Dominance == "Dominant_shrub_1")

shrub_mod <- glmmTMB(
  coverstd ~ 
  Years_since_fire +
   Temp_sc +
    Per_sc +
    (1 | StudyID),      # random intercept for study
  family = beta_family(link = "logit"),
#  dispformula = ~ Fire_Int_Groups,
  data = shrubs_dom,
  weights = weight_sc
)

shrub_mod2 <- glmmTMB(
  coverstd ~ 
    Years_since_fire +
    Temp_sc +
    Per_sc +
    (1 | StudyID),      # random intercept for study
  family = beta_family(link = "logit"),
  #  dispformula = ~ Fire_Int_Groups,
  data = shrubs_dom,
  weights = weight_sc
)

AIC_vals <- AIC(shrub_mod, shrub_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(shrub_mod)
plot(simres)

shrubs_dom$pred <- predict(shrub_mod, type = "response")
plot(shrubs_dom$pred, shrubs_dom$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(shrub_mod)
Anova(shrub_mod, type = 'III')


#Plot predictions!
shrub_pred <- emmeans(
  shrub_mod,
  ~ Years_since_fire,
  at = list(
    Years_since_fire = seq(
      min(shrubs_dom$Years_since_fire, na.rm = TRUE),
      max(shrubs_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
        Per_sc = mean(shrubs_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(shrubs_dom$Temp_sc, na.rm = TRUE)
  )
)

pred_grid_shrub <- as.data.frame(summary(shrub_pred, infer = TRUE))

pred_grid_shrub <- pred_grid_shrub %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
  )

predshrubplot <- ggplot(pred_grid_shrub,
                       aes(x = Years_since_fire,
                           y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Time since fire (years)",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Most dominant shrub in North America")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 


predshrubplot

ggsave(plot = predshrubplot, filename = "Pred_shrub_plot.png", dpi =300,
       height = 4.2, width = 6.5)

#####
#Mosses
bryo_dom <- df_long %>%
  filter(PlantGroup == "bryophyte")

bryo_dom <- bryo_dom %>%
  filter(!is.na(coverstd))

ggplot(bryo_dom, aes(x = Years_since_fire, y = coverstd)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  facet_grid(Continent ~ Dominance)

#Very few 2nd & 3rd dominant. Combine to one

bryo_dom$Dominance <- fct_collapse(
  bryo_dom$base,
  Dominant_bryophyte_1 = "Dominant_bryophyte_1",
  Dominant_bryophyte_2 = c("Dominant_bryophyte_2", "Dominant_bryophyte_3")
)

bryo_mod <- glmmTMB(
  coverstd ~ 
    ns(Years_since_fire, df= 2) * Fire_Int_Groups * Dominance * Continent +
    Temp_sc + 
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = bryo_dom,
  dispformula = ~ Continent,
  weights = weight_sc
)

bryo_mod2 <- glmmTMB(
  coverstd ~ 
    Years_since_fire * Fire_Int_Groups * Dominance +
    (1 | StudyID / RowID),      # random intercept for study
  family = beta_family(link = "logit"),
  data = bryo_dom,
  dispformula = ~ Dominance + Fire_Int_Groups,
  weights = weight_sc
)

AIC_vals <- AIC(bryo_mod, bryo_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

simres <- simulateResiduals(bryo_mod)
plot(simres)

bryo_dom$pred <- predict(bryo_mod, type = "response")
plot(bryo_dom$pred, bryo_dom$coverstd,
     xlab = "Predicted coverstd",
     ylab = "Observed coverstd")
abline(0,1, col="red")

summary(bryo_mod)
Anova(bryo_mod, type = 'III')


#Plot predictions!
bryo_pred <- emmeans(
  bryo_mod,
  ~ Years_since_fire | Dominance* Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(bryo_dom$Years_since_fire, na.rm = TRUE),
      max(bryo_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    # Per_sc = mean(bryo_dom$Per_sc, na.rm = TRUE)
    Temp_sc = mean(bryo_dom$Temp_sc, na.rm = TRUE)
  )
)

pred_grid_bryo <- as.data.frame(summary(bryo_pred, infer = TRUE))

pred_grid_bryo <- pred_grid_bryo %>%
  mutate(
    fit   = plogis(emmean),
    lower = plogis(asymp.LCL),
    upper = plogis(asymp.UCL)
  )

pred_grid_bryo$Fire_Int_Groups <- factor(
  pred_grid_bryo$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predmossplot <- ggplot(pred_grid_bryo,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Years since fire",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  facet_grid(Continent ~ Dominance,
             labeller = labeller(
               Dominance = c("Dominant_bryophyte_1" = "Most dominant sp.",
                             "Dominant_bryophyte_2" = "Second most dominant sp.",
               "Dominant_bryophyte_3" = "Third most dominant sp."),
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"
               ))) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Bryophytes")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA"))

predmossplot

ggsave(plot = predmossplot, filename = "Pred_moss_plot.png", dpi =300,
       height = 4.2, width = 6.5)

#Combinedplot

combinedplot_tree <- predtreeplot_NA/(predtreeplot_EU|predshrubplot) +
  plot_layout(heights = c(2, 2, 1, 1))
combinedplot_tree

combinedplot_ground <-  predherbplot/preddwarfplot/predgramplot/predmossplot +
  plot_layout(heights = c(2, 1, 1, 2))
combinedplot_ground

ggsave(plot=combinedplot_ground, filename = "coverplots_ground.png", dpi =300,
       height = 14, width = 13)
ggsave(plot=combinedplot_tree, filename = "coverplots_tree.png", dpi =300,
       height = 10.52, width = 16.25)


ggplot(subset(df_long, PlantGroup == "herb"), aes(x = Years_since_fire, y = coverstd)) +
  geom_point(aes(color = StudyID)) +
  facet_grid(Continent ~Fire_Int_Groups)
       