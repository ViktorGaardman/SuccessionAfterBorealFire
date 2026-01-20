#Step 1. Load packages
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(ggeffects)
library(performance)
library(DHARMa)
library(glmmTMB)

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
#1. Find a new way to describe three species
#2. Are any effects non-linear?

#Create a plantgroup column 
df_long <- df_long %>%
  mutate(
    PlantGroup = str_extract(base, "(?<=_)[^_]+(?=_)")
  )

#Create a total size of study column
df_long <- df_long %>%
  mutate(
    studysize = Plot_size * Sample_size
  )

#
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

herbs$base <- as.factor(herbs$base)

herbs_dom <- herbs %>%
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


predherbplot

ggsave(plot = predherbplot, filename = "Pred_herb_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#DWARFSHRUBS

dwarfs <- df_long %>%
  filter(PlantGroup == "dwarfshrub")

dwarfs$base <- as.factor(dwarfs$base)

dwarfs_dom <- dwarfs %>%
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

preddwarfplot

ggsave(plot = preddwarfplot, filename = "Pred_dwarf_plot.png", dpi =300,
       height = 4.2, width = 6.5)


####
#Graminoids

graminoid <- df_long %>%
  filter(PlantGroup == "graminoid")

graminoid$base <- as.factor(graminoid$base)

gram_dom <- graminoid %>%
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

predgramplot

ggsave(plot = predgramplot, filename = "Pred_gram_plot.png", dpi =300,
       height = 4.2, width = 6.5)

####
#Trees

trees <- df_long %>%
  filter(PlantGroup == "tree")

trees$base <- as.factor(trees$base)

tree_dom <- trees %>%
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

predtreeplot

ggsave(plot = predtreeplot, filename = "Pred_tree_plot.png", dpi =300,
       height = 4.2, width = 6.5)
