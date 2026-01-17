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

#Question: What is the relative dominance
#of each trait across time, continent, and fire intensity. Is it affected
#by mean temp, perc, or soil water content?
#Look at how traits are affected by successional pathways
#(no dominance, just baisc traits)

# Names of species in each dataset
species_try <- unique(TRY_Traits$species)
species_df <- unique(df_long$species)

# Species in TRY_clean that are not in df_long
missing_species <- setdiff(species_try, species_df)

missing_species

Seed_mass_mod <- lmer(log(Seed_dry_mass) ~
                        Years_since_fire +
                        Avg_Temp+
                        AvgPer+
                        Fire_Int_Groups +
                        Latitude*Continent +
                        (1|RowID/StudyID.x),
                      data = df_long) 

car::Anova(Seed_mass_mod, type ='III')

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
  stat_smooth(method='lm', aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

##################
#Community cover tests using a beta-regression
#1. Should we change the ysf category?
#2. Are any effects non-linear?

library(glmmTMB)

#New variable for YSF with more data per step
df_long <- df_long %>%
  mutate(
    YSF_interval = case_when(
      Years_since_fire <= 1 ~ "1",
      Years_since_fire <= 2 ~ "2",
      Years_since_fire <= 3 ~ "3",
      Years_since_fire <= 4 ~ "4",
      Years_since_fire <= 5 ~ "5",
      Years_since_fire <= 7 ~ "6",
      Years_since_fire <= 9 ~ "7",
      Years_since_fire <= 22 ~ "8"
    )
  )

df_long$YSF_interval <- as.factor(df_long$YSF_interval)
  

#Create a plantgroup column 
df_long <- df_long %>%
  mutate(
    PlantGroup = str_extract(base, "(?<=_)[^_]+(?=_)")
  )

df_long <- df_long %>%
  mutate(
    studysize = Plot_size * Sample_size
  )

#Divide by 100 to get a 0-1 range but
#avoid 0 and 1 one in the dataset (beta must be >0 and <1)

df_long$coverstd <- (df_long$cover + 0.01) / 101

#Scale and center Latitude, percipitation, and temperature
df_long$Temp_sc <- scale(df_long$Avg_Temp, center = TRUE, scale = TRUE)
df_long$Per_sc <- scale(df_long$AvgPer, center = TRUE, scale = TRUE)
df_long$Latitude_sc <- scale(df_long$Latitude, center = TRUE, scale = TRUE)



#Variables to add as non-linear (atleast try)
#Temperature (herbs, dwarfshrub, graminoid, shrub, tree)
#Percipitation (herbs, dwarfshrub, graminoid, shrub)
#YSF_interval (herbs, trees)
#Latitude (trees)
#Continent (no shrubs in Eurasia basically, so don't divide by continent)

#Current issues:
#Model convergence failure when we add non-linear terms
#High uncertainties at later years for NA.
#Should we cut the data at e.g., 15 years?


herbs <- df_long %>%
  filter(PlantGroup == "herb")



herbmod <- glmmTMB(
  coverstd ~
    Years_since_fire*Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire +
    Avg_Temp
    AvgPer +
    Latitude +
    (1 | RowID/StudyID.x) +
    (1 | species),
  family = beta_family(),
  weights = studysize,
  data = herbs
)


summary(herbmod)
Anova(herbmod, type = 'III')


#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(herbs$Years_since_fire, na.rm = TRUE),
    max(herbs$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(herbs$Fire_Int_Groups),
  Continent     = levels(herbs$Continent)
) %>%
  mutate(
    Avg_Temp   = mean(herbs$Avg_Temp, na.rm = TRUE),
    AvgPer = mean(herbs$AvgPer, na.rm = TRUE),
    Latitude      = mean(herbs$Latitude, na.rm = TRUE),
    studysize     = mean(herbs$studysize, na.rm = TRUE)
  )

#type = response ok?
pred <- predict(
  herbmodlin,
  newdata = pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA,
  allow.new.levels = TRUE
)

pred_grid <- pred_grid %>%
  mutate(
    fit   = pred$fit,
    se    = pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

pred_grid$Fire_Int_Groups <- factor(
  pred_grid$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predherbplot<- ggplot(pred_grid,
       aes(x = Years_since_fire,
           y = fit,
           color = Fire_Int_Groups)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted herb cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.direction='vertical',
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) 

ggsave(plot = predherbplot, filename = "Pred_herb_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#DWARFSHRUBS

dwarfs <- df_long %>%
  filter(PlantGroup == "dwarfshrub")

dwarfmod <- glmmTMB(
  coverstd ~
    Years_since_fire*Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire +
    Avg_Temp +
    AvgPer +
    Latitude +
    (1 | RowID/StudyID.x) +
    (1 | species),
  family = beta_family(),
  weights = studysize,
  data = dwarfs
)

Anova(herbmod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(dwarfs$Years_since_fire, na.rm = TRUE),
    max(dwarfs$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(dwarfs$Fire_Int_Groups),
  Continent     = levels(dwarfs$Continent)
) %>%
  mutate(
    Avg_Temp   = mean(dwarfs$Avg_Temp, na.rm = TRUE),
    AvgPer = mean(dwarfs$AvgPer, na.rm = TRUE),
    Latitude      = mean(dwarfs$Latitude, na.rm = TRUE),
    studysize     = mean(dwarfs$studysize, na.rm = TRUE)
  )

#type = response ok?
pred <- predict(
  dwarfmod,
  newdata = pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA,
  allow.new.levels = TRUE
)

pred_grid <- pred_grid %>%
  mutate(
    fit   = pred$fit,
    se    = pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

pred_grid$Fire_Int_Groups <- factor(
  pred_grid$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

preddwarfplot<- ggplot(pred_grid,
                      aes(x = Years_since_fire,
                          y = fit,
                          color = Fire_Int_Groups)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted dwarfshrub cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.direction='vertical',
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) 

ggsave(plot = preddwarfplot, filename = "Pred_dwarf_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#Graminoids

graminoid <- df_long %>%
  filter(PlantGroup == "graminoid")

grammod <- glmmTMB(
  coverstd ~
    Years_since_fire*Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire +
    Avg_Temp +
    AvgPer +
    Latitude +
    (1 | RowID/StudyID.x) +
    (1 | species),
  family = beta_family(),
  weights = studysize,
  data = graminoid
)

Anova(grammod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(graminoid$Years_since_fire, na.rm = TRUE),
    max(graminoid$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(graminoid$Fire_Int_Groups),
  Continent     = levels(graminoid$Continent)
) %>%
  mutate(
    Avg_Temp   = mean(graminoid$Avg_Temp, na.rm = TRUE),
    AvgPer = mean(graminoid$AvgPer, na.rm = TRUE),
    Latitude      = mean(graminoid$Latitude, na.rm = TRUE),
    studysize     = mean(graminoid$studysize, na.rm = TRUE)
  )

#type = response ok?
pred <- predict(
  grammod,
  newdata = pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA,
  allow.new.levels = TRUE
)

pred_grid <- pred_grid %>%
  mutate(
    fit   = pred$fit,
    se    = pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

pred_grid$Fire_Int_Groups <- factor(
  pred_grid$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predgramplot<- ggplot(pred_grid,
                       aes(x = Years_since_fire,
                           y = fit,
                           color = Fire_Int_Groups)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Continent) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Fire_Int_Groups),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted graminoid cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.direction='vertical',
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) 

ggsave(plot = predgramplot, filename = "Pred_gram_plot.png", dpi =300,
       height = 4.2, width = 6.5)