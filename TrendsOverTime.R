#Step 1. Load packages
library(tidyverse)
library(car)
library(ggeffects)
library(DHARMa)
library(glmmTMB)
library(patchwork)
library(lme4)
library(performance)
library(nlme)
library(emmeans)

rm(list=ls())

#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Clean_Data.csv", sep = ";")

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
  left_join(metadata, by = "RowID")


#Step 3. Load climate and trait data
TRY_Traits <- read.csv("TRY_Clean.csv", sep = ";")

Perc_Data <- read.csv("site_climatological_percipitation_WorldClim_v2.csv",
                      sep = ";")

Temp_Data <- read.csv("site_climatological_temperature_WorldClim_v2.csv",
                      sep = ";")

species_long_Temp <- species_long_meta %>%
  left_join(Temp_Data, by = "Title")

species_long_Perc <- species_long_Temp %>%
  left_join(Perc_Data, by = "Title")

traits_wide <- TRY_Traits %>%
  select(-c("n_obs", "trait_SE", "trait_SD")) %>%
  pivot_wider(
    names_from = TraitName,
    values_from = "trait_mean"
  ) 

traits_wide <- traits_wide[,1:9]

df_long <- species_long_Perc %>%
  left_join(traits_wide, by = "species")

#Drop 0 in species column for now

df_long <- subset(df_long, !species == 0)

df_long <- df_long %>% 
  mutate(across(c("Fire_Int_Groups", "Continent", "StudyID.x"), as.factor))

#Create a plantgroup column 
df_long <- df_long %>%
  mutate(
    PlantGroup = str_extract(base, "(?<=_)[^_]+(?=_)")
  )

#Filter out bryophytes
df_traits <- df_long %>%
  filter(! PlantGroup %in% "bryophyte")

# Names of species in each dataset
species_try <- unique(TRY_Traits$species)
species_df <- unique(df_long$species)

# Species in TRY_clean that are not in df_long
missing_species <- setdiff(species_try, species_df)

missing_species

#Scale and center percipitation and temperature
df_traits$Temp_sc <- as.numeric(scale(df_traits$Avg_Temp, center = TRUE, scale = TRUE))
df_traits$Per_sc <- as.numeric(scale(df_traits$AvgPer, center = TRUE, scale = TRUE))

#Add average genera level trait values to species_sp
#1. Define traits
traits <- c(
  "Seed_dry_mass",
  "Seed_longevity",
  "Plant_height_vegetative",
  "Leaf_nitrogen",
  "Leaf_area_PI"
)


#2. Define genus/species info
df_traits <- df_traits %>%
  mutate(
    genus = str_extract(species, "^[^_]+"),
    species2 = str_extract(species, "(?<=_).+"),
    is_genus_level = species2 %in% c("sp.", "sp")
  )

#3. Calculate genus level mean trait values for all genera
#that has data for more than 1 species
genus_traits <- df_traits %>%
  filter(!is_genus_level) %>%
  group_by(genus) %>%
  summarise(
    across(
      all_of(traits),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_genus_mean"
    ),
    n_species = n_distinct(species2),
    .groups = "drop"
  ) %>%
  filter(n_species >= 2)


#4. Add to dataset
df_filled <- df_traits %>%
  left_join(genus_traits, by = "genus") %>%
  mutate(
    across(
      all_of(traits),
      ~ if_else(
        is_genus_level & is.na(.x),
        get(paste0(cur_column(), "_genus_mean")),
        .x
      )
    )
  ) %>%
  select(-ends_with("_genus_mean"))

df_filled$species <- as.factor(df_filled$species)

df_filled <- df_filled %>%
  mutate(
    studysize = Plot_size * Sample_size
  )

#Group into ground/treeshrub layer
Ground_df <-df_filled %>%
  filter(PlantGroup %in% c("dwarfshrub" ,"herb", "graminoid"))

#Group into ground layer
Tree_df <-df_filled %>%
  filter(PlantGroup %in% c("tree", "shrub"))

#Calculate community weighted means
ground.cwm <-   # New dataframe where we can inspect the result
  Ground_df %>%   # First step in the next string of statements
  group_by(studysize, StudyID.x, RowID, Temp_sc, Per_sc, Fire_Int_Groups,
           Years_since_fire, Continent) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    Height_cwm = weighted.mean(Plant_height_vegetative, cover, na.rm = TRUE),   # Actual calculation of CWMs
    Mass_cwm = weighted.mean(Seed_dry_mass, cover, na.rm = TRUE),
    Nitrogen_cwm = weighted.mean(Leaf_nitrogen, cover, na.rm = TRUE),
    Area_cwm = weighted.mean(Leaf_area_PI, cover, na.rm = TRUE),
    Long_cwm = weighted.mean(Seed_longevity, cover, na.rm = TRUE),
    Disp_cwm = weighted.mean(Dispersal_unit_dry_mass, cover, na.rm = TRUE)
  )

#Calculate community weighted means
tree.cwm <-   # New dataframe where we can inspect the result
  Tree_df %>%   # First step in the next string of statements
  group_by(StudyID.x, RowID, Temp_sc, Per_sc, Fire_Int_Groups,
           Years_since_fire, Continent, studysize) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    Height_cwm = weighted.mean(Plant_height_vegetative, cover, na.rm = TRUE),   # Actual calculation of CWMs
    Mass_cwm = weighted.mean(Seed_dry_mass, cover, na.rm = TRUE),
    Nitrogen_cwm = weighted.mean(Leaf_nitrogen, cover, na.rm = TRUE),
    Area_cwm = weighted.mean(Leaf_area_PI, cover, na.rm = TRUE),
    Long_cwm = weighted.mean(Seed_longevity, cover, na.rm = TRUE),
    Disp_cwm = weighted.mean(Dispersal_unit_dry_mass, cover, na.rm = TRUE)
  )

#use only first 

#Ground layer
ggplot(ground.cwm, aes(x = Years_since_fire, y=Mass_cwm, by = Fire_Int_Groups))+
  geom_point(aes(color = Fire_Int_Groups)) + 
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

#Use only first 13 years after fire
ground_sub_cwm <- ground.cwm %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 13)

#Fins best model by comparing model AICs using method = "ML"
options(contrasts = c("contr.sum", "contr.poly"))

Seed_mass_mod <- gls(
  log(Mass_cwm) ~
    Years_since_fire *
    Fire_Int_Groups,
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varExp(form = ~ studysize),
  method      = "REML"
)

Seed_mass_mod2 <- gls(
  log(Mass_cwm) ~
    Years_since_fire *
    Fire_Int_Groups,
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varExp(form = ~ studysize),
  method      = "ML"
)

AIC_vals <- AIC(Seed_mass_mod, Seed_mass_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals
plot(Seed_mass_mod, resid(., type = "normalized") ~ fitted(.))
qqnorm(resid(Seed_mass_mod, type = "normalized"))
qqline(resid(Seed_mass_mod, type = "normalized"))
plot(
  log(ground_sub_cwm$studysize),
  abs(resid(Seed_mass_mod, type = "normalized"))
)
intervals(Seed_mass_mod)$corStruct #Shows we should keep studyId.x in the model


mod_fixed <- update(Seed_mass_mod, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Seed_mass_mod, weights = varPower(form = ~ studysize), method = "ML")
mod_exp   <- update(Seed_mass_mod, weights = varExp(form =~ studysize), method = "ML")

AIC(mod_fixed, mod_power, mod_exp) #Exp model is best

summary(Seed_mass_mod)
Anova(Seed_mass_mod, type = 'III')


#Plot predictions!
emm_sm <- emmeans(
  Seed_mass_mod,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(ground.cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground.cwm$Per_sc, na.rm = TRUE)
  ),
  weights = "proportional"
)


pred_grid_sm <- as.data.frame(emm_sm)

pred_grid_sm <- pred_grid_sm %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )



pred_grid_sm$Fire_Int_Groups <- factor(
  pred_grid_sm$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

seedmassplot<- ggplot(pred_grid_sm,
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
    y = "mg",
    color = "Fire intensity"
  ) +
  theme_bw() +
 # scale_y_continuous(limits = c(0,50)) +
  scale_x_continuous(limits = c(1,13), 
                     breaks = c(1, 3, 5, 7, 9, 11, 13)) +
  ggtitle("Seed dry mass")+
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
        plot.title = element_text(size = 18, hjust = 0.5)) 


seedmassplot

ggsave(plot = seedmassplot, filename = "Seedmass_ground.png", dpi =300,
       height = 4.2, width = 6.5)

##PLANT HEIGHT MODEL

ggplot(ground.cwm, aes(x = Years_since_fire, y = log(Height_cwm), 
                       by = Fire_Int_Groups))+
  geom_point(aes(color= Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

Plant_height_mod <- gls(
  log(Height_cwm) ~
    Years_since_fire * Fire_Int_Groups +
    Temp_sc,
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ 1 / log(studysize)),
  method      = "REML"
)


AIC(Plant_height_mod)

AIC_vals <- AIC(Plant_height_mod, Plant_height_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Plant_height_mod, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Plant_height_mod, weights = varFixed(~ log(studysize)), method = "ML")
mod_exp   <- update(Plant_height_mod, weights = varExp(form =~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Plant_height_mod, type = "normalized"))
qqline(resid(Plant_height_mod, type = "normalized"))
plot(Plant_height_mod, resid(., type = "normalized") ~ fitted(.))

summary(Plant_height_mod)
Anova(Plant_height_mod, type = 'III')


#Plot predictions!
emm_ph <- emmeans(
  Plant_height_mod,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 30
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  ), weights = "proportional"
)


pred_grid_ph <- as.data.frame(emm_ph)

pred_grid_ph <- pred_grid_ph %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_ph$Fire_Int_Groups <- factor(
  pred_grid_ph$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

plantheightplot<- ggplot(pred_grid_ph,
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
    y = "m",
    color = "Fire intensity"
  ) +
  theme_bw() +
#  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(1,13), 
                     breaks = c(1, 3, 5, 7, 9, 11, 13))+
  ggtitle("Plant height")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

plantheightplot

ggsave(plot = plantheightplot, filename = "Plantheight_ground.png", dpi =300,
       height = 4.2, width = 6.5)

#Leaf nitrogen
Leaf_nitrogen_mod <- gls(
  log(Nitrogen_cwm) ~
    poly(Years_since_fire, 2) +
    Fire_Int_Groups +
    poly(Temp_sc, 3),
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ studysize),
  method      = "REML"
)


AIC_vals <- AIC(Leaf_nitrogen_mod, Leaf_nitrogen_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Leaf_nitrogen_mod, weights = varPower(form = ~ studysize), method = "ML")
mod_power <- update(Leaf_nitrogen_mod, weights = varPower(form = ~ log(studysize)), method = "ML")
mod_exp   <- update(Leaf_nitrogen_mod, weights = varExp(form = ~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Leaf_nitrogen_mod, type = "normalized"))
qqline(resid(Leaf_nitrogen_mod, type = "normalized"))
plot(Leaf_nitrogen_mod, resid(., type = "normalized") ~ fitted(.))

summary(Leaf_nitrogen_mod)
Anova(Leaf_nitrogen_mod, type = 'III')
#Plot predictions!
emm_ln <- emmeans(
  Leaf_nitrogen_mod,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 30
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)


pred_grid_ln <- as.data.frame(emm_ln)

pred_grid_ln <- pred_grid_ln %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_ln$Fire_Int_Groups <- factor(
  pred_grid_ln$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafnitrogenplot<- ggplot(pred_grid_ln,
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
    y = "mg/g",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Leaf nitrogen")+
  scale_x_continuous(limits = c(1,13), 
                     breaks = c(1, 3, 5, 7, 9, 11, 13)) +
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

leafnitrogenplot

ggsave(plot = leafnitrogenplot, filename = "LEafnitrogen_ground.png", dpi =300,
       height = 4.2, width = 6.5)

#Leaf_area

Leaf_area_mod <- gls(
  log(Area_cwm) ~
    Years_since_fire * Fire_Int_Groups *
    Continent +
    poly(Temp_sc, 3) +
    Per_sc,
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ studysize),
  method      = "REML"
)


AIC_vals <- AIC(Leaf_area_mod, Leaf_area_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Leaf_area_mod, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Leaf_area_mod, weights = varPower(form = ~ studysize), method = "ML")
mod_exp   <- update(Leaf_area_mod, weights = varPower(form = ~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Leaf_area_mod, type = "normalized"))
qqline(resid(Leaf_area_mod, type = "normalized"))
plot(Leaf_area_mod, resid(., type = "normalized") ~ fitted(.))
summary(Leaf_area_mod)
Anova(Leaf_area_mod, type = 'III')


#Plot predictions!
emm_la <- emmeans(
  Leaf_area_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 30
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)


pred_grid_la <- as.data.frame(emm_la)

pred_grid_la <- pred_grid_la %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_la$Fire_Int_Groups <- factor(
  pred_grid_la$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafareaplot<- ggplot(pred_grid_la,
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
    y = expression(mm^2/mg),
    color = "Fire intensity"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(1,13), 
                     breaks = c(1, 3, 5, 7, 9, 11, 13))+
  ggtitle("SLA") +
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

leafareaplot

ggsave(plot = leafareaplot, filename = "Leafarea_ground.png", dpi =300,
       height = 4.2, width = 6.5)

#Longevity

longevity_df <- ground_sub_cwm %>%
  filter(!(is.na(Long_cwm)))

Longevity_mod <- gls(
  log(Long_cwm) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups *
    Continent +
    poly(Temp_sc, 3) +
    Per_sc,
  data = longevity_df,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ studysize),
  method      = "REML"
)

Longevity_mod2 <- gls(
  log(Long_cwm) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups +
    Continent +
    poly(Temp_sc, 3) +
    Per_sc,
  data = longevity_df,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ studysize),
  method      = "REML"
)

AIC_vals <- AIC(Longevity_mod, Longevity_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Longevity_mod, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Longevity_mod, weights = varPower(form = ~ studysize), method = "ML")
mod_exp   <- update(Longevity_mod, weights = varPower(form = ~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Longevity_mod, type = "normalized"))
qqline(resid(Longevity_mod, type = "normalized"))
plot(Longevity_mod, resid(., type = "normalized") ~ fitted(.))
summary(Longevity_mod)
Anova(Longevity_mod, type = 'III')


#Plot predictions!
emm_la <- emmeans(
  Leaf_area_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 30
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)


pred_grid_la <- as.data.frame(emm_la)

pred_grid_la <- pred_grid_la %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_la$Fire_Int_Groups <- factor(
  pred_grid_la$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafareaplot<- ggplot(pred_grid_la,
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
    y = expression(mm^2/mg),
    color = "Fire intensity"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(1,13), 
                     breaks = c(1, 3, 5, 7, 9, 11, 13))+
  ggtitle("SLA") +
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

leafareaplot

ggsave(plot = leafareaplot, filename = "Leafarea_ground.png", dpi =300,
       height = 4.2, width = 6.5)


###########
combinedtraitplot_ground <- (plantheightplot|leafnitrogenplot)/(leafareaplot|seedmassplot)

combinedtraitplot_ground

ggsave(plot = combinedtraitplot_ground, filename = "traitplots_ground.png",
       dpi = 300, height = 7.39, width = 13)





##################
#Tree layer
Seed_mass_mod_T <- lm(log(Mass_cwm) ~
                      Years_since_fire*Continent +
                      Years_since_fire * Fire_Int_Groups +
                      Continent * Fire_Int_Groups +
                      I(Temp_sc^2)+
                      Temp_sc +
                      Per_sc,
                    data = tree.cwm)

plot(Seed_mass_mod_T)
shapiro.test(resid(Seed_mass_mod_T)) #QQ

summary(Seed_mass_mod_T)
Anova(Seed_mass_mod_T, type = 'III')

#Plot predictions!
pred_grid_sm_T <- expand.grid(
  Years_since_fire = seq(
    min(tree.cwm$Years_since_fire, na.rm = TRUE),
    max(tree.cwm$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(tree.cwm$Fire_Int_Groups),
  Continent     = levels(tree.cwm$Continent)
) %>%
  mutate(
    Temp_sc   = mean(tree.cwm$Temp_sc, na.rm = TRUE),
    Per_sc = mean(tree.cwm$Per_sc, na.rm = TRUE)
  )


pred_sm_T <- predict(
  Seed_mass_mod_T,
  newdata = pred_grid_sm_T,
  type = "response",      
  se.fit = TRUE
)

eta <- pred_sm_T$fit
se  <- pred_sm_T$se.fit

lower_eta <- eta - 1.96 * se
upper_eta <- eta + 1.96 * se

# back-transform from log scale
pred_grid_sm_T$fit   <- exp(eta)
pred_grid_sm_T$lower <- exp(lower_eta)
pred_grid_sm_T$upper <- exp(upper_eta)

pred_grid_sm_T$Fire_Int_Groups <- factor(
  pred_grid_sm_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

seedmassplot_T <- ggplot(pred_grid_sm_T,
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
    y = "mg",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Seed dry mass")+
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
        plot.title = element_text(size = 18, hjust = 0.5)) 


seedmassplot_T

ggsave(plot = seedmassplot_T, filename = "Seedmass_tree.png", dpi =300,
       height = 4.2, width = 6.5)

Leaf_nitrogen_mod_T <- lm(log(Nitrogen_cwm) ~
                        Years_since_fire*Continent +
                        Years_since_fire * Fire_Int_Groups +
                        Continent * Fire_Int_Groups +
                        I(Temp_sc^2)+
                        Temp_sc +
                        Per_sc,
                      data = tree.cwm)

plot(Leaf_nitrogen_mod_T)
shapiro.test(resid(Leaf_nitrogen_mod_T)) #QQ

summary(Leaf_nitrogen_mod_T)
Anova(Leaf_nitrogen_mod_T, type = 'III')

#Plot predictions!
pred_grid_ln_T <- expand.grid(
  Years_since_fire = seq(
    min(tree.cwm$Years_since_fire, na.rm = TRUE),
    max(tree.cwm$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(tree.cwm$Fire_Int_Groups),
  Continent     = levels(tree.cwm$Continent)
) %>%
  mutate(
    Temp_sc   = mean(tree.cwm$Temp_sc, na.rm = TRUE),
    Per_sc = mean(tree.cwm$Per_sc, na.rm = TRUE)
  )


pred_ln_T <- predict(
  Leaf_nitrogen_mod_T,
  newdata = pred_grid_ln_T,
  type = "response",      
  se.fit = TRUE
)

eta <- pred_ln_T$fit
se  <- pred_ln_T$se.fit

lower_eta <- eta - 1.96 * se
upper_eta <- eta + 1.96 * se

# back-transform from log scale
pred_grid_ln_T$fit   <- exp(eta)
pred_grid_ln_T$lower <- exp(lower_eta)
pred_grid_ln_T$upper <- exp(upper_eta)

pred_grid_ln_T$Fire_Int_Groups <- factor(
  pred_grid_ln_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafnitrogenplot_T <- ggplot(pred_grid_ln_T,
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
    y = "mg/g",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Leaf nitrogen")+
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
        plot.title = element_text(size = 18, hjust = 0.5)) 


leafnitrogenplot_T

ggsave(plot = leafnitrogenplot_T, filename = "Leafnitrogen_tree.png", dpi =300,
       height = 4.2, width = 6.5)

Leaf_area_mod_T <- lm(log(Area_cwm) ~
                            Years_since_fire*Continent +
                            Years_since_fire * Fire_Int_Groups +
                            Continent * Fire_Int_Groups +
                            I(Temp_sc^2)+
                            Temp_sc +
                            Per_sc,
                          data = tree.cwm)

plot(Leaf_area_mod_T)
shapiro.test(resid(Leaf_area_mod_T)) #QQ

summary(Leaf_area_mod_T)
Anova(Leaf_area_mod_T, type = 'III')

#Plot predictions!
pred_grid_la_T <- expand.grid(
  Years_since_fire = seq(
    min(tree.cwm$Years_since_fire, na.rm = TRUE),
    max(tree.cwm$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(tree.cwm$Fire_Int_Groups),
  Continent     = levels(tree.cwm$Continent)
) %>%
  mutate(
    Temp_sc   = mean(tree.cwm$Temp_sc, na.rm = TRUE),
    Per_sc = mean(tree.cwm$Per_sc, na.rm = TRUE)
  )


pred_la_T <- predict(
  Leaf_area_mod_T,
  newdata = pred_grid_ln_T,
  type = "response",      
  se.fit = TRUE
)

eta <- pred_la_T$fit
se  <- pred_la_T$se.fit

lower_eta <- eta - 1.96 * se
upper_eta <- eta + 1.96 * se

# back-transform from log scale
pred_grid_la_T$fit   <- exp(eta)
pred_grid_la_T$lower <- exp(lower_eta)
pred_grid_la_T$upper <- exp(upper_eta)

pred_grid_la_T$Fire_Int_Groups <- factor(
  pred_grid_la_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafareaplot_T <- ggplot(pred_grid_la_T,
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
    y = expression(mm^2/mg),
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("SLA")+
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
        plot.title = element_text(size = 18, hjust = 0.5)) 


leafareaplot_T

ggsave(plot = leafareaplot_T, filename = "Leafarea_tree.png", dpi =300,
       height = 4.2, width = 6.5)

combinedtraitplot_tree <- leafareaplot_T / leafnitrogenplot_T / seedmassplot_T

combinedtraitplot_tree

ggsave(plot = combinedtraitplot_tree, filename = "traitplot_tree.png",
       dpi = 300, height = 10.52, width =13)

