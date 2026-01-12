#Step 1. Load packages
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(ggeffects)
library(performance)


#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Raw_Data.csv", sep = ";")

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


#Model what affects community 
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

#Is percipitation/temp and SWI correlated?

ggplot(species_long_Perc, aes(x = SWI, y=log(AvgPer))) +
  geom_point()+
  stat_smooth(method='lm')

ggplot(species_long_Perc, aes(x = SWI, y=Avg_Temp)) +
  geom_point()+
  stat_smooth(method='lm')

#Yes both are correlated. Make a correlation plot to check for sure?

traits_wide <- traits_wide[,1:9]

df_long <- species_long_Perc %>%
  left_join(traits_wide, by = "species")

#Drop 0 in species column for now

df_long <- subset(df_long, !species == 0)

df_long <- df_long %>% 
  mutate(across(c("Fire_Int_Groups", "Continent", "StudyID.x"), as.factor))

#Standardize plot size for looking at dominant traits
df_long_s <- df_long %>%
  group_by(RowID) %>%
  mutate(rel_cover = cover / sum(cover, na.rm = TRUE)) %>%
  ungroup()

#Multiple trait value by cover. Question: What is the relative dominance
#of each trait across time, continent, and fire intensity. Is it affected
#by mean temp, perc, or soil water content?
CWM_df <- df_long_s %>%
  group_by(RowID) %>%
  summarise(
    CWM_Leaf_area_PE = sum(rel_cover * Leaf_area_PE, na.rm = TRUE),
    CWM_Leaf_area_PI = sum(rel_cover * Leaf_area_PI, na.rm = TRUE),
    CWM_Leaf_nitrogen = sum(rel_cover * Leaf_nitrogen, na.rm = TRUE),
    CWM_Plant_height_vegetative = sum(rel_cover * Plant_height_vegetative, na.rm = TRUE),
    CWM_Plant_height_generative = sum(rel_cover * Plant_height_generative, na.rm = TRUE),
    CWM_Seed_dry_mass = sum(rel_cover * Seed_dry_mass, na.rm = TRUE),
    .groups = "drop"
  )

#Look at how traits are affected by successional pathways
#(no dominance, just baisc traits)

# Names of species in each dataset
species_try <- unique(TRY_Traits$species)
species_df <- unique(df_long$species)

# Species in TRY_clean that are not in df_long
missing_species <- setdiff(species_try, species_df)

missing_species

Seed_mass_mod <- lmer(Seed_dry_mass ~
                        Years_since_fire +
                        Fire_Int_Groups +
                        Continent +
                        SWI +
                        Avg_Temp +
                        AvgPer +
                        (1|StudyID.x),
                      data = df_long) 

car::Anova(Seed_mass_mod, type ='II')

check_normality(Seed_mass_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(Seed_mass_mod)  # Checks if residual variance is consistent
check_collinearity(Seed_mass_mod)
check_outliers(Seed_mass_mod)
sim_res <- simulateResiduals(fittedModel = Seed_mass_mod, n = 500)
plot(sim_res)

ggplot(df_long, aes(x=Years_since_fire, y = log(Seed_dry_mass))) +
  geom_point(aes(color = Fire_Int_Groups))+
  stat_smooth(method='lm', aes(color = Fire_Int_Groups))

Plant_height_mod <- lmer(Plant_height_vegetative ~
                           Years_since_fire +
                           Fire_Int_Groups +
                           Continent +
                           SWI +
                           Avg_Temp +
                           AvgPer +
                           (1|StudyID.x),
                         data = df_long) 

car::Anova(Plant_height_mod, type ='II')

check_normality(Plant_height_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(Plant_height_mod)  # Checks if residual variance is consistent
check_collinearity(Plant_height_mod)
check_outliers(Plant_height_mod)
sim_res <- simulateResiduals(fittedModel = Plant_height_mod, n = 500)
plot(sim_res)

ggplot(df_long, aes(x=Years_since_fire, y = log(Seed_dry_mass))) +
  geom_point(aes(color = Fire_Int_Groups))+
  stat_smooth(method='lm', aes(color = Fire_Int_Groups))