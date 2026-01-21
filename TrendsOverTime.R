#Step 1. Load packages
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(ggeffects)
library(performance)
library(DHARMa)
library(glmmTMB)
library(patchwork)

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

#Create a plantgroup column 
df_long <- df_long %>%
  mutate(
    PlantGroup = str_extract(base, "(?<=_)[^_]+(?=_)")
  )

#Filter out bryophytes
df_traits <- df_long %>%
  filter(! PlantGroup %in% "bryophyte")

#Look at the data
ggplot(data = df_traits, aes(x = log(Seed_dry_mass))) +
  geom_histogram() #Log

ggplot(data = df_traits, aes(x = log(Plant_height_vegetative))) +
  geom_histogram() #Log

ggplot(data = df_traits, aes(x = log(Leaf_nitrogen))) +
  geom_histogram() #log

ggplot(data = df_traits, aes(x = log(Leaf_area_PI))) +
  geom_histogram() #log

ggplot(data = df_traits, aes(x = log(Seed_longevity+0.01))) +
  geom_histogram() 

#Decide on leaf area PE vs PI and vegetative vs generative plant height

#PI has more data than PE
#vegetative has more data than generative

#Dispersal unit dry mass has 474 observations (1326 missing). Skip
#this trait
#Seedlongevity has 1119 observations (681 missing). We can try it
#but also much data lacking.

#Question: What is the relative dominance
#of each trait across time, continent, and fire intensity. Is it affected
#by mean temp, perc, or soil water content?
#Look at how traits are affected by successional pathways
#(no dominance, just baisc traits)
#Remember to add average genera-level trait values for species_sp.

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

Seed_mass_mod <- lm(log(Seed_dry_mass) ~
                        Years_since_fire*Fire_Int_Groups +
                        Continent +
                        Per_sc +
                        Temp_sc,
                      weights = log(cover+1),
                      data = df_filled) 


car::Anova(Seed_mass_mod, type ='III')

check_normality(Seed_mass_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(Seed_mass_mod)  # Checks if residual variance is consistent
check_collinearity(Seed_mass_mod)
check_outliers(Seed_mass_mod)
sim_res <- simulateResiduals(fittedModel = Seed_mass_mod, n = 500)
plot(sim_res)

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(df_filled$Years_since_fire, na.rm = TRUE),
    max(df_filled$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(df_filled$Fire_Int_Groups),
  Continent     = levels(df_filled$Continent)
) %>%
  mutate(
    Temp_sc   = mean(df_filled$Temp_sc, na.rm = TRUE),
    Per_sc = mean(df_filled$Per_sc, na.rm = TRUE),
    cover     = mean(df_filled$cover, na.rm = TRUE)
  )


pred <- predict(
  Seed_mass_mod,
  newdata = pred_grid,
  type = "response",      
  se.fit = TRUE,
  re.form = NA
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

seedmassplot<- ggplot(pred_grid,
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
    y = "Predicted value",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Seed dry mass")+
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


seedmassplot

ggsave(plot = seedmassplot, filename = "Seedmass_plot.png", dpi =300,
       height = 4.2, width = 6.5)

temppred <- ggpredict(Seed_mass_mod, terms = "Temp_sc")

#Do we need to predict for temp and precipitation? I think not.
#One example here kept in case we want it later.
plottemp <- ggplot(temppred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 0.8, color = "firebrick") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.3) +
  theme_minimal() +
  ylab("Predicted value") +
  xlab("Temperature (scaled)") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

plottemp


##PLANT HEIGHT MODEL

Plant_height_mod <- lm(log(Plant_height_vegetative) ~
                             Years_since_fire +
                         Fire_Int_Groups +
                             Continent +
                         Temp_sc +
                         Per_sc,
                           weights = log(cover+1),
                           data = df_filled) 

car::Anova(Plant_height_mod, type ='II')

check_normality(Plant_height_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(Plant_height_mod)  # Checks if residual variance is consistent
check_collinearity(Plant_height_mod)
check_outliers(Plant_height_mod)
sim_res <- simulateResiduals(fittedModel = Plant_height_mod, n = 500)
plot(sim_res)

#Leaf nitrogen

leaf_nig_mod <- lm(log(Leaf_nitrogen) ~
        Years_since_fire *
        Fire_Int_Groups +
          Years_since_fire * Temp_sc +
          Temp_sc +
          Per_sc,
      weights = log(cover+1),
      data = df_filled) 

car::Anova(leaf_nig_mod, type ='II')

check_normality(leaf_nig_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(leaf_nig_mod)  # Checks if residual variance is consistent
check_collinearity(leaf_nig_mod)
check_outliers(leaf_nig_mod)
sim_res <- simulateResiduals(fittedModel = leaf_nig_mod, n = 500)
plot(sim_res)


#Leaf_area

leaf_area_mod <- lm(log(Leaf_area_PI) ~
                     Years_since_fire +
                     Fire_Int_Groups +
                     Continent +
                     Temp_sc +
                     Per_sc,
                   weights = log(cover+1),
                   data = df_filled) 

car::Anova(leaf_nig_mod, type ='II')

check_normality(leaf_nig_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(leaf_nig_mod)  # Checks if residual variance is consistent
check_collinearity(leaf_nig_mod)
check_outliers(leaf_nig_mod)
sim_res <- simulateResiduals(fittedModel = leaf_area_mod, n = 500)
plot(sim_res)

#Seed longevity
seed_long_mod <- lm(Seed_longevity ~
                      Years_since_fire +
                      Fire_Int_Groups +
                      Continent +
                      Temp_sc +
                      Per_sc,
                    data = df_filled) 

car::Anova(leaf_nig_mod, type ='II')

check_normality(leaf_nig_mod)  # Tests if residuals are normally distributed
check_heteroscedasticity(leaf_nig_mod)  # Checks if residual variance is consistent
check_collinearity(leaf_nig_mod)
check_outliers(leaf_nig_mod)
sim_res <- simulateResiduals(fittedModel = leaf_area_mod, n = 500)
plot(sim_res)


##################

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

herbs_dom <- df_long %>%
  filter(base == "Dominant_herb_1")

#Normalize studysize for better model fit
herbs_dom$studysize_norm <- sqrt(herbs_dom$studysize) / mean(sqrt(herbs_dom$studysize), na.rm = TRUE)


herbmod <- glmmTMB(
  coverstd ~
    Years_since_fire * Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire,
  dispformula = ~ Fire_Int_Groups,
  family = beta_family(),
  data = herbs_dom
)



sim_res <- simulateResiduals(fittedModel = herbmod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, herbs_dom$Years_since_fire)
plotResiduals(sim_res, herbs_dom$Fire_Int_Groups)
plotResiduals(res0, herbs_dom$Avg_Temp)
plotResiduals(res0, herbs_dom$AvgPer)
plotResiduals(sim_res, herbs_dom$Continent)



summary(herbmod)
Anova(herbmod, type = 'III')


#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(herbs_dom$Years_since_fire, na.rm = TRUE),
    max(herbs_dom$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(herbs_dom$Fire_Int_Groups),
  Continent     = levels(herbs_dom$Continent)
) 
#%>%
#  mutate(
#    Avg_Temp   = mean(herbs_dom$Avg_Temp, na.rm = TRUE),
#    AvgPer = mean(herbs_dom$AvgPer, na.rm = TRUE),
#    studysize     = mean(herbs_dom$studysize, na.rm = TRUE)
#  )


pred <- predict(
  herbmod,
  newdata = pred_grid,
  type = "link",      
  se.fit = TRUE,
  re.form = NA
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

predherbplot<- ggplot(pred_grid,
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
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Herbs")+
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


predherbplot

ggsave(plot = predherbplot, filename = "Pred_herb_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#DWARFSHRUBS

dwarfs_dom <- df_long %>%
  filter(base == "Dominant_dwarfshrub_1")

dwarfmod <- glmmTMB(
  coverstd ~
    Years_since_fire*Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire,
  family = beta_family(),
  dispformula = ~ Continent,
  data = dwarfs_dom
)

sim_res <- simulateResiduals(fittedModel = dwarfmod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, dwarfs_dom$Years_since_fire)
plotResiduals(sim_res, dwarfs_dom$Fire_Int_Groups)
plotResiduals(sim_res, dwarfs_dom$Continent)



summary(dwarfmod)
Anova(dwarfmod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(dwarfs_dom$Years_since_fire, na.rm = TRUE),
    max(dwarfs_dom$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(dwarfs_dom$Fire_Int_Groups),
  Continent     = levels(dwarfs_dom$Continent)
) 

#type = response ok?
pred <- predict(
  dwarfmod,
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

preddwarfplot<- ggplot(pred_grid,
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
    y = "Predicted dwarfshrub cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Dwarfshrubs")+
  theme(legend.position="right",
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

preddwarfplot

ggsave(plot = preddwarfplot, filename = "Pred_dwarf_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#Graminoids

gram_dom <- df_long %>%
  filter(base == "Dominant_graminoid_1")

grammod <- glmmTMB(
  coverstd ~
    Years_since_fire*Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire,
  family = beta_family(),
#  dispformula = ~Years_since_fire,
  data = gram_dom
)

sim_res <- simulateResiduals(fittedModel = grammod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, gram_dom$Years_since_fire)
plotResiduals(sim_res, gram_dom$Fire_Int_Groups)
plotResiduals(sim_res, gram_dom$Continent)



summary(grammod)
Anova(grammod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(gram_dom$Years_since_fire, na.rm = TRUE),
    max(gram_dom$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(gram_dom$Fire_Int_Groups),
  Continent     = levels(gram_dom$Continent)
) 

#type = response ok?
pred <- predict(
  grammod,
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

predgramplot<- ggplot(pred_grid,
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
    y = "Predicted graminoid cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Graminoids")+
  theme(legend.position="none",
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.direction='vertical',
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

predgramplot

ggsave(plot = predgramplot, filename = "Pred_gram_plot.png", dpi =300,
       height = 4.2, width = 6.5)

####
#Trees

tree_dom <- df_long %>%
  filter(base == "Dominant_tree_1")

tree_dom <- tree_dom %>%
  filter(!is.na(coverstd))

treemod <- glmmTMB(
  coverstd ~
    Years_since_fire*Continent +
    Fire_Int_Groups * Continent +
    Fire_Int_Groups * Years_since_fire,
  family = beta_family(),
    dispformula = ~ Fire_Int_Groups,
  data = tree_dom
)

sim_res <- simulateResiduals(fittedModel = treemod, n = 500)
plot(sim_res)
testDispersion(sim_res)
testQuantiles(sim_res)
plotResiduals(sim_res, tree_dom$Years_since_fire)
plotResiduals(sim_res, tree_dom$Fire_Int_Groups)
plotResiduals(sim_res, tree_dom$Continent)



summary(treemod)
Anova(treemod, type = 'III')

#Plot predictions!
pred_grid <- expand.grid(
  Years_since_fire = seq(
    min(tree_dom$Years_since_fire, na.rm = TRUE),
    max(tree_dom$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(tree_dom$Fire_Int_Groups),
  Continent     = levels(tree_dom$Continent)
) 

#type = response ok?
pred <- predict(
  treemod,
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

predtreeplot<- ggplot(pred_grid,
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
    y = "Predicted tree cover",
    color = "Fire intensity"
  ) +
  theme_bw() +
  ggtitle("Trees")+
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
   dispformula = ~ Continent,
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
