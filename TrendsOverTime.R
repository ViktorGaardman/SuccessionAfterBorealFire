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

ggplot(data = df_traits, aes(x = sqrt(Plant_height_vegetative))) +
  geom_histogram() #Log

ggplot(data = df_traits, aes(x = sqrt(Leaf_nitrogen))) +
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

#Use only first 10 years after fire
ground_sub_cwm <- ground.cwm %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

Seed_mass_mod <- gls(
  log(Mass_cwm) ~
    poly(Years_since_fire, 3) * Continent +
    poly(Years_since_fire, 3) * Fire_Int_Groups +
    poly(Temp_sc, 4)+
    Per_sc,
  data = ground_sub_cwm,
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ 1 / log(studysize)),
  method      = "REML"
)

plot(Seed_mass_mod, resid(., type = "normalized") ~ fitted(.))
summary(Seed_mass_mod)
Anova(Seed_mass_mod, type = 'II')


#Plot predictions!
emm_sm <- emmeans(
  Seed_mass_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  )
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
                          y = emmean,
                          color = Fire_Int_Groups)) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Fire_Int_Groups),
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
    y = "log(mg)",
    color = "Fire intensity"
  ) +
  theme_bw() +
 # scale_y_continuous(limits = c(0,50)) +
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
    poly(Years_since_fire, 3) * Continent +
    poly(Years_since_fire, 3) * Fire_Int_Groups +
    poly(Temp_sc, 3)+
    Per_sc,
  
  data = ground_sub_cwm,
  
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ 1 / log(studysize)),
  method      = "REML"
)

plot(Plant_height_mod, resid(., type = "normalized") ~ fitted(.))
summary(Plant_height_mod)
Anova(Plant_height_mod, type = 'II')


#Plot predictions!
emm_ph <- emmeans(
  Plant_height_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 30
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  )
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
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "m",
    color = "Fire intensity"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0,1))+
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

ggplot(ground.cwm, aes(x = Years_since_fire, y = log(Nitrogen_cwm), 
                       by = Fire_Int_Groups))+
  geom_point(aes(color= Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

ggplot(ground_sub_cwm, aes(x = Temp_sc, y = log(Mass_cwm)))+
  geom_point()+
  geom_smooth()
  #facet_wrap(~Continent)

Leaf_nitrogen_mod <- gls(
  log(Nitrogen_cwm) ~
    poly(Years_since_fire, 3) * Continent +
    poly(Years_since_fire, 3) * Fire_Int_Groups +
    poly(Temp_sc, 3) +
    Per_sc,
  
  data = ground_sub_cwm,
  
  correlation = corCompSymm(form = ~ 1 | StudyID.x),
  weights     = varFixed(~ 1 / log(studysize)),
  method      = "REML"
)

plot(Leaf_nitrogen_mod, resid(., type = "normalized") ~ fitted(.))
summary(Leaf_nitrogen_mod)
Anova(Leaf_nitrogen_mod, type = 'II')


#Plot predictions!
emm_ln <- emmeans(
  Leaf_nitrogen_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      max(ground_sub_cwm$Years_since_fire, na.rm = TRUE),
      length.out = 30
    ),
    Temp_sc = mean(ground_sub_cwm$Temp_sc, na.rm = TRUE),
    Per_sc  = mean(ground_sub_cwm$Per_sc, na.rm = TRUE)
  )
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

Leaf_area_mod <- lm(log(Area_cwm) ~
                          Years_since_fire*Continent +
                          Years_since_fire * Fire_Int_Groups +
                          Continent * Fire_Int_Groups +
                          I(Temp_sc^2)+
                          Temp_sc +
                          Per_sc,
                        data = ground.cwm)

plot(Leaf_area_mod)
shapiro.test(resid(Leaf_area_mod)) #QQ

summary(Leaf_area_mod)
Anova(Leaf_area_mod, type = 'III')

#Plot predictions!
pred_grid_la <- expand.grid(
  Years_since_fire = seq(
    min(ground.cwm$Years_since_fire, na.rm = TRUE),
    max(ground.cwm$Years_since_fire, na.rm = TRUE),
    length.out = 88
  ),
  Fire_Int_Groups = levels(ground.cwm$Fire_Int_Groups),
  Continent     = levels(ground.cwm$Continent)
) %>%
  mutate(
    Temp_sc   = mean(ground.cwm$Temp_sc, na.rm = TRUE),
    Per_sc = mean(ground.cwm$Per_sc, na.rm = TRUE)
  )


pred_la <- predict(
  Leaf_area_mod,
  newdata = pred_grid_la,
  type = "response",      
  se.fit = TRUE
)

eta <- pred_la$fit
se  <- pred_la$se.fit

lower_eta <- eta - 1.96 * se
upper_eta <- eta + 1.96 * se

# back-transform from log scale
pred_grid_la$fit   <- exp(eta)
pred_grid_la$lower <- exp(lower_eta)
pred_grid_la$upper <- exp(upper_eta)


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
  dispformula = ~Years_since_fire,
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
    dispformula = ~ Years_since_fire + Fire_Int_Groups,
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

