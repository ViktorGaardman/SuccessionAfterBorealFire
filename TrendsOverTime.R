#Step 1. Load packages
library(tidyverse)
library(car)
library(ggeffects)
library(patchwork)
library(nlme)
library(emmeans)
library(splines)
library(MASS)
library(glmmTMB)
library(DHARMa)
library(ordinal)
library(nnet)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Clean_Data.csv", sep = ";")

table(subset(df, Continent == "Eurasia")$Fire_Int_Groups, subset(df, Continent == "Eurasia")$Years_since_fire)

table(subset(df, Continent == "North_America")$Fire_Int_Groups, subset(df, Continent == "North_America")$Years_since_fire)

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

#Create a plantgroup column 
species_long_Perc <- species_long_Perc %>%
  mutate(
    PlantGroup = str_extract(base, "(?<=_)[^_]+(?=_)")
  )

traits_bryo <- species_long_Perc %>%
  filter(PlantGroup %in% "bryophyte")

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

#Filter out bryophytes
df_traits <- df_long %>%
  filter(! PlantGroup %in% "bryophyte")

#Check missing data for traits

trait_cols <- c("Dispersal_unit_dry_mass", "Seed_longevity")

species_summary <- df_traits %>%
  group_by(PlantGroup, species) %>%
  summarise(
    freq = n(),
    across(all_of(trait_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(PlantGroup, desc(freq))

plantgroup_dfs <- species_summary %>%
  select(-freq) %>%   # optional: remove frequency column
  group_split(PlantGroup, keep = FALSE)

species_sp_freq <- df_traits %>%
  filter(str_detect(species, "_sp\\.$")) %>%
  count(species, sort = TRUE)

library(dplyr)
library(stringr)

species_sp_freq <- df_traits %>%
  # keep only genus-level species
  filter(str_detect(species, "_sp\\.$")) %>%
  
  # extract genus
  mutate(genus = str_remove(species, "_sp\\.$")) %>%
  
  # count frequency
  count(species, genus, sort = TRUE) %>%
  
  # for each genus, find matching full species
  rowwise() %>%
  mutate(
    species_names = paste(
      unique(
        df_traits$species[
          str_detect(df_traits$species, paste0("^", genus, "_")) &
            !str_detect(df_traits$species, "_sp\\.$")
        ]
      ),
      collapse = "; "
    )
  ) %>%
  ungroup() %>%
  select(-genus)


writexl::write_xlsx(species_sp_freq, "Species_sp_freq.xlsx")

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

#Ground layer
ggplot(ground_sub_cwm, aes(x = Years_since_fire, y=Mass_cwm, by = Fire_Int_Groups))+
  geom_point(aes(color = Fire_Int_Groups)) + 
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

#Use only first 10 years after fire
ground_sub_cwm <- ground.cwm %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

#Fins best model by comparing model AICs using method = "ML"

Seed_mass_mod <- gls(
  log(Mass_cwm) ~
    Years_since_fire *
    Fire_Int_Groups,
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varExp(form = ~ studysize),
  method      = "REML"
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
  scale_x_continuous(limits = c(1,10), 
                     breaks = c(2, 4, 6, 8, 10)) +
  ggtitle("Seed dry mass")+
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


seedmassplot

ggsave(plot = seedmassplot, filename = "Seedmass_ground.png", dpi =300,
       height = 4.2, width = 6.5)

##PLANT HEIGHT MODEL

ggplot(ground_sub_cwm, aes(x = Years_since_fire, y = log(Height_cwm), 
                       by = Fire_Int_Groups))+
  geom_point(aes(color= Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

Plant_height_mod <- gls(
  log(Height_cwm) ~
    Years_since_fire * Fire_Int_Groups *
    Continent +
    poly(Temp_sc,3),
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ 1 / log(studysize)),
  method      = "REML"
)

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
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE)
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
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
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
  scale_x_continuous(limits = c(1,10), 
                     n.breaks = 6)+
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

plantheightplot

ggsave(plot = plantheightplot, filename = "Plantheight_ground.png", dpi =300,
       height = 4.2, width = 6.5)

#Leaf nitrogen
Leaf_nitrogen_mod <- gls(
  log(Nitrogen_cwm) ~
    Fire_Int_Groups * Continent +
    Years_since_fire * Continent +
    poly(Temp_sc, 3) +
    Per_sc,
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
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
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
  scale_x_continuous(limits = c(1,10), 
                     n.breaks = 6) +
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

leafnitrogenplot

ggsave(plot = leafnitrogenplot, filename = "LEafnitrogen_ground.png", dpi =300,
       height = 4.2, width = 6.5)

#Leaf_area


Leaf_area_mod <- gls(
  log(Area_cwm) ~
    Years_since_fire +
    Fire_Int_Groups * Continent +
    poly(Temp_sc, 3),
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
      length.out = 40
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE)
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
  scale_x_continuous(limits = c(1,10), 
                     n.breaks = 6)+
  ggtitle("SLA") +
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

leafareaplot

ggsave(plot = leafareaplot, filename = "Leafarea_ground.png", dpi =300,
       height = 4.2, width = 6.5)

#Longevity

longevity_df <- ground_sub_cwm %>%
  filter(!(is.na(Long_cwm)))

Longevity_mod <- gls(
  log(Long_cwm) ~
    Years_since_fire +
    Continent + 
    Temp_sc,
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
mod_exp   <- update(Longevity_mod, weights = varExp(form = ~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Longevity_mod, type = "normalized"))
qqline(resid(Longevity_mod, type = "normalized"))
plot(Longevity_mod, resid(., type = "normalized") ~ fitted(.))
summary(Longevity_mod)
Anova(Longevity_mod, type = 'III')


#Plot predictions!
emm_lo <- emmeans(
  Longevity_mod,
  ~ Years_since_fire | Continent,
  at = list(
    Years_since_fire = seq(
      min(longevity_df$Years_since_fire, na.rm = TRUE),
      max(longevity_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(longevity_df$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)


pred_grid_lo <- as.data.frame(emm_lo)

pred_grid_lo <- pred_grid_lo %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_lo$Fire_Int_Groups <- factor(
  pred_grid_lo$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

longevityplot<- ggplot(pred_grid_lo,
                      aes(x = Years_since_fire,
                          y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA, show.legend = FALSE
                ) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             scale = "free_y",
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Years",
    color = "Longevity"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), 
                     n.breaks = 6)+
  ggtitle("Seed longevity") +
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

longevityplot

ggsave(plot = longevityplot, filename = "Longevity_ground.png", dpi =300,
       height = 4.2, width = 6.5)


#Dispersal unit mass

disp_df <- ground_sub_cwm %>%
  filter(!(is.na(Disp_cwm)))

Dispersal_mod <- gls(
  log(Disp_cwm) ~
    Years_since_fire +
    Fire_Int_Groups * Continent +
    Temp_sc,
  data = disp_df,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ log(studysize)),
  method      = "REML"
)

AIC_vals <- AIC(Dispersal_mod, Dispersal_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Dispersal_mod, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Dispersal_mod, weights = varPower(form = ~ log(studysize)), method = "ML")
mod_exp   <- update(Dispersal_mod, weights = varExp(form = ~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Dispersal_mod, type = "normalized"))
qqline(resid(Dispersal_mod, type = "normalized"))
plot(Dispersal_mod, resid(., type = "normalized") ~ fitted(.))
summary(Dispersal_mod)
Anova(Dispersal_mod, type = 'III')


#Plot predictions!
emm_di <- emmeans(
  Dispersal_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(disp_df$Years_since_fire, na.rm = TRUE),
      max(disp_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(disp_df$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)


pred_grid_di <- as.data.frame(emm_di)

pred_grid_di <- pred_grid_di %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_di$Fire_Int_Groups <- factor(
  pred_grid_di$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

dispplot<- ggplot(pred_grid_di,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             scale = "free_y",
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "mg",
    color = "Longevity"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), 
                     n.breaks = 6)+
  ggtitle("Dispersal unit dry mass") +
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

dispplot

ggsave(plot = dispplot, filename = "Disp_ground.png", dpi =300,
       height = 4.2, width = 6.5)


###########
combinedtraitplot_ground <- (plantheightplot|leafnitrogenplot)/
  (leafareaplot|seedmassplot)/
  (longevityplot|dispplot)

combinedtraitplot_ground

ggsave(plot = combinedtraitplot_ground, filename = "traitplots_ground.TIFF",
       dpi = 450, height = 10.56, width = 13)


####################




##################
#Tree layer

ggplot(tree_sub_cwm, aes(x = Years_since_fire, y=log(Mass_cwm), by = Fire_Int_Groups))+
  geom_point(aes(color = Fire_Int_Groups)) + 
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

#Use only first 10 years after fire
tree_sub_cwm <- tree.cwm %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)


tree_sub_cwm <- tree_sub_cwm %>%
  filter(!(is.na(Long_cwm)))

Seed_mod_T <- gls(
  log(Mass_cwm) ~
    Years_since_fire * Fire_Int_Groups * Continent +
    Per_sc,
  data = tree_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ log(studysize)),
  method      = "REML"
)

AIC_vals <- AIC(Seed_mod_T, Seed_mod_T2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Seed_mod_T, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Seed_mod_T, weights = varPower(form = ~ log(studysize)), method = "ML")
mod_exp   <- update(Seed_mod_T, weights = varPower(form = ~ studysize), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Seed_mod_T, type = "normalized"))
qqline(resid(Seed_mod_T, type = "normalized"))
plot(Seed_mod_T, resid(., type = "normalized") ~ fitted(.))
summary(Seed_mod_T)
Anova(Seed_mod_T, type = 'III')


#Plot predictions!
emm_mass_T <- emmeans(
  Seed_mod_T,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(tree_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(tree_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(tree_sub_cwm$Per_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_mass_T <- as.data.frame(emm_mass_T)

pred_grid_mass_T <- pred_grid_mass_T %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_mass_T$Fire_Int_Groups <- factor(
  pred_grid_mass_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

seedmassplot_T <- ggplot(pred_grid_mass_T,
                      aes(x = Years_since_fire,
                          y = fit,
                          color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             scales = "free_y",
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

seedmassplot_T

ggsave(plot = seedmassplot_T, filename = "Seedmass_tree.png", dpi =300,
       height = 4.2, width = 6.5)


#Leaf nitrogen
Leaf_nitrogen_T <- gls(
  log(Nitrogen_cwm) ~
    poly(Years_since_fire, 2) * Fire_Int_Groups * Continent +
    Temp_sc,
  data = tree_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varPower(form = ~ studysize),
  method      = "REML"
)

AIC_vals <- AIC(Leaf_nitrogen_T, Leaf_nitrogen_T2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

mod_fixed <- update(Leaf_nitrogen_T, weights = varFixed(~ 1 / log(studysize)), method = "ML")
mod_power <- update(Leaf_nitrogen_T, weights = varPower(form = ~ studysize), method = "ML")
mod_exp   <- update(Leaf_nitrogen_T, weights = varPower(form = ~ log(studysize)), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(Leaf_nitrogen_T, type = "normalized"))
qqline(resid(Leaf_nitrogen_T, type = "normalized"))
plot(Leaf_nitrogen_T, resid(., type = "normalized") ~ fitted(.))
summary(Leaf_nitrogen_T)
Anova(Leaf_nitrogen_T, type = 'III')


#Plot predictions!
emm_nitro_T <- emmeans(
  Leaf_nitrogen_T,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(tree_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(tree_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(tree_sub_cwm$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_nitro_T <- as.data.frame(emm_nitro_T)

pred_grid_nitro_T <- pred_grid_nitro_T %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_nitro_T$Fire_Int_Groups <- factor(
  pred_grid_nitro_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafnitrogenplot_T <- ggplot(pred_grid_nitro_T,
                         aes(x = Years_since_fire,
                             y = fit,
                             color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent,
             scales = "free_y",
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 


leafnitrogenplot_T

combinedtraitplot_tree <- leafnitrogenplot_T / seedmassplot_T

combinedtraitplot_tree

ggsave(plot = combinedtraitplot_tree, filename = "traitplot_tree.TIFF",
       dpi = 450, height = 10.52, width =13)



#####BRYOPHYTE TRAITS
traits_bryo$Temp_sc <- as.numeric(scale(traits_bryo$Avg_Temp, center = TRUE, scale = TRUE))
traits_bryo$Per_sc <- as.numeric(scale(traits_bryo$AvgPer, center = TRUE, scale = TRUE))

Bryo_trait_data <- read.csv("bryoATT_Cleaned.csv", sep = ";")
Bryo_trait_data <- Bryo_trait_data[,1:4]

traits_bryo_df <- traits_bryo %>%
  left_join(Bryo_trait_data, by = "species")

traits_bryo_df <- traits_bryo_df %>%
  filter(!is.na(Len))

traits_bryo_df <- traits_bryo_df %>%
  mutate(
    studysize = Plot_size * Sample_size
  )


#Normalize studysize for better model fit
traits_bryo_df$weight_sc <- traits_bryo_df$studysize / mean(traits_bryo_df$studysize)

traits_bryo_df <- traits_bryo_df %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

#Calculate community weighted means
bryo.cwm <-   # New dataframe where we can inspect the result
  traits_bryo_df %>%   # First step in the next string of statements
  group_by(StudyID.x, RowID, Temp_sc, Per_sc, Fire_Int_Groups,
           Years_since_fire, Continent, weight_sc) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    Length_cwm = weighted.mean(Len, cover, na.rm = TRUE),   # Actual calculation of CWMs
  )

ggplot(bryo.cwm, aes( x = log(Length_cwm))) + 
  geom_histogram()

#Trait models 

#Length

ggplot(bryo.cwm, aes( x = Years_since_fire, y = Length_cwm)) +
  geom_point(aes(color = Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups)) +
  facet_wrap(~Continent)

Length_mod <- gls(
  log(Length_cwm) ~
    ns(Years_since_fire, df = 2) * Fire_Int_Groups +
    Continent + 
    Temp_sc,
  data = bryo.cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ I(1 / weight_sc)),
  method      = "ML"
)

Length_mod <- gls(
  log(Length_cwm) ~
    Fire_Int_Groups * Continent +
    ns(Years_since_fire, df = 3) * Continent +
    ns(Per_sc, df = 3),
  data = bryo.cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ I(1 / weight_sc)),
  method      = "REML"
)


AIC_vals <- AIC(Length_mod, Length_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

qqnorm(resid(Length_mod, type = "normalized"))
qqline(resid(Length_mod, type = "normalized"))
plot(Length_mod, resid(., type = "normalized") ~ fitted(.))
summary(Length_mod)
Anova(Length_mod, type = 'III')


#Plot predictions!
emm_length <- emmeans(
  Length_mod,
  ~ Years_since_fire | Continent * Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(bryo.cwm$Years_since_fire, na.rm = TRUE),
      max(bryo.cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(bryo.cwm$Temp_sc, na.rm = TRUE)),
  weights = 'proportional'
)

pred_grid_length <- as.data.frame(emm_length)

pred_grid_length <- pred_grid_length %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_length$Fire_Int_Groups <- factor(
  pred_grid_length$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

mosslength_plot <- ggplot(pred_grid_length,
                             aes(x = Years_since_fire,
                                 y = fit,
                                 color = Fire_Int_Groups)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  labs(
    x = "Years since fire",
    y = "mm"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0,330))+
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  ggtitle("Leafy shoot length")+
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


mosslength_plot

#Sporophyte

ggplot(traits_bryo_df, aes(x = Sporophyte_frequency, fill = Primary_lifeform)) +
  geom_bar(position = "dodge")

traits_bryo_df$Sporophyte_frequency <- ordered(
  traits_bryo_df$Sporophyte_frequency
)

ord_mod <- clmm(
  Sporophyte_frequency ~
    Years_since_fire * Continent * Fire_Int_Groups +
    Temp_sc +
    (1 | StudyID.x),
  data = traits_bryo_df,
  Hess = TRUE
)


AIC_vals <- AIC(ord_mod, ord_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

nominal_test(ord_mod)

summary(ord_mod)



#Plot predictions!
emm_ord <- emmeans(
  ord_mod,
  ~ Years_since_fire | Continent * Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(traits_bryo_df$Years_since_fire, na.rm = TRUE),
      max(traits_bryo_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Continent = unique(traits_bryo_df$Continent),
    Fire_Int_Groups = c("High", "Medium", "Low"),
    Temp_sc = mean(traits_bryo_df$Temp_sc, na.rm = TRUE),
    StudyID.x = NA)  # set to NA or leave random effect at zero
)

pred_grid_ord <- as.data.frame(emm_ord)


# 2. Extract thresholds from clmm model
alpha <- ord_mod$alpha
alpha_names <- names(alpha)

# 3. Convert emmean (linear predictor) to cumulative probabilities
cumprob <- sapply(pred_grid_ord$emmean, function(eta) plogis(alpha - eta))
cumprob <- t(cumprob)  # rows = observations, cols = thresholds
colnames(cumprob) <- alpha_names

# 4. Convert cumulative probabilities to category probabilities
# P1 = cumprob1, P2 = cumprob2 - cumprob1, P3 = 1 - cumprob2
prob_mat <- cbind(
  cumprob[,1],
  cumprob[,2] - cumprob[,1],
  1 - cumprob[,2]
)
colnames(prob_mat) <- c("Rare", "Occasional", "Frequent")


# 5. Combine with prediction grid
pred_long <- cbind(pred_grid_ord, prob_mat) %>%
  pivot_longer(
    cols = c("Rare", "Occasional", "Frequent"),
    names_to = "Sporophyte_frequency",
    values_to = "Probability"
  )

# 6. Plot
mossord_plot <- ggplot(pred_long,
       aes(x = Years_since_fire,
           y = Probability,
           color = Sporophyte_frequency,
           fill = Sporophyte_frequency)) +
  geom_ribbon(aes(ymin = 0, ymax = Probability), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_grid(Continent ~ Fire_Int_Groups,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"
               )
             )) +
  scale_color_manual(values = c("black", "#E69F00", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#009E73")) +
  labs(x = "Years since fire", y = "Predicted probability") +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  ggtitle("Sporophyte frequency")+
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

mossord_plot



##Primary lifeform

ggplot(traits_bryo_df, aes(x = Primary_lifeform, y = Years_since_fire)) +
  geom_point() +
  facet_grid(Fire_Int_Groups~Continent)

nom_mod <- multinom(
  Primary_lifeform ~
    Years_since_fire +
    Fire_Int_Groups +
    Per_sc,
  data = traits_bryo_df,
  Hess = TRUE,
  trace = FALSE
)

nom_mod2 <- multinom(
  Primary_lifeform ~
    Years_since_fire +
    Fire_Int_Groups +
    Per_sc,
  data = traits_bryo_df,
  Hess = TRUE, 
  trace = FALSE
)

nom_brms <- brm(
  formula = Primary_lifeform ~ 
    Years_since_fire + 
    Fire_Int_Groups + 
    Per_sc + 
    (1 | StudyID.x),
  data = traits_bryo_df,
  family = categorical(link = "logit"),  # nominal multinomial
  cores = 4,      # adjust to your CPU
  iter = 4000,    # increase if model is slow to converge
  warmup = 1000,
  chains = 4,
  seed = 123
)

nom_brms2 <- brm(
  formula = Primary_lifeform ~ 
    Years_since_fire * 
    Fire_Int_Groups + 
    Per_sc + 
    (1 | StudyID.x),
  moment_match = TRUE,
  data = traits_bryo_df,
  family = categorical(link = "logit"),  # nominal multinomial
  cores = 4,      # adjust to your CPU
  iter = 4000,    # increase if model is slow to converge
  warmup = 1000,
  chains = 4,
  seed = 123
)

loo1 <- loo(nom_brms)
loos2 <- loo(nom_brms2)

loo_compare(loo1, loos2)
pp_check(nom_brms, type = "bars")
pp_check(nom_brms2, type = "bars")

summary(nom_brms)



#Plot predictions!
pred_grid_nom <- expand.grid(
  Years_since_fire = seq(
    min(traits_bryo_df$Years_since_fire, na.rm = TRUE),
    max(traits_bryo_df$Years_since_fire, na.rm = TRUE),
    length.out = 40
  ),
  Continent = unique(traits_bryo_df$Continent),
  Fire_Int_Groups = c("High", "Medium", "Low"),
  Per_sc = mean(traits_bryo_df$Per_sc, na.rm = TRUE),  # optional covariates
  StudyID.x = NA  # population-level predictions; remove random effect
)

# Population-level predicted probabilities
pred_probs_nom <- fitted(
  nom_brms,
  newdata = pred_grid_nom,
  re_formula = NA  # ignores random effects
)

pred_probs_df <- as.data.frame(pred_probs_nom)
colnames(pred_probs_df)

pred_long_life <- cbind(pred_grid_nom, pred_probs_nom) %>%
  pivot_longer(
    cols = c("Estimate.P(Y = Mat)", "Estimate.P(Y = Tuft)",
             "Estimate.P(Y = Weft)"),
    names_to = "Primary_lifeform",
    values_to = "Probability"
  )


library(ggplot2)

ggplot(pred_long_life,
       aes(x = Years_since_fire,
           y = Probability,
           color = Primary_lifeform,
           fill = Primary_lifeform)) +
  geom_ribbon(aes(ymin = 0, ymax = Probability), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent + Fire_Int_Groups) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) +
  labs(x = "Years since fire", y = "Predicted probability") +
  theme_bw()


# 6. Plot
mossnom_plot <- ggplot(pred_long,
                       aes(x = Years_since_fire,
                           y = Probability,
                           color = Sporophyte_frequency,
                           fill = Sporophyte_frequency)) +
  geom_ribbon(aes(ymin = 0, ymax = Probability), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_grid(Continent ~ Fire_Int_Groups,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"
               )
             )) +
  scale_color_manual(values = c("black", "#E69F00", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#009E73")) +
  labs(x = "Years since fire", y = "Predicted probability") +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  ggtitle("Sporophyte frequency")+
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

mossnom_plot