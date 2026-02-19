#Step 1. Load packages
library(tidyverse)
library(car)
library(ggeffects)
library(patchwork)
#library(nlme)
library(lme4)
library(emmeans)
library(splines)
library(glmmTMB)
library(DHARMa)
library(brms)
library(ordinal)
#library(nnet)
library(brms)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

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
  left_join(metadata, by = c("StudyID", "RowID"))


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
  mutate(across(c("Fire_Int_Groups", "Continent", "StudyID"), as.factor))

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
df_traits$SWI_sc <- as.numeric(scale(df_traits$SWI, center = TRUE, scale = TRUE))
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

#Drop, Corydalis. sp, Rubus. sp as these genera
#are too morphologically variable to use mean values from

df_filled <- df_filled %>%
  filter(! species %in% c("Corydalis.sp", "Rubus.sp"))

df_filled <- df_filled %>%
  mutate(
    studysize = Plot_size * Sample_size ,
    area_sc = studysize / mean(studysize)
  )

#Use only first 10 years after fire
df_filled <- df_filled %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)


#Group into ground/treeshrub layer
Ground_df <-df_filled %>%
  filter(PlantGroup %in% c("dwarfshrub" ,"herb", "graminoid"))

#Group into ground layer
Tree_df <-df_filled %>%
  filter(PlantGroup %in% c("tree", "shrub"))

#Calculate community weighted means
ground.cwm <-   # New dataframe where we can inspect the result
  Ground_df %>%   # First step in the next string of statements
  group_by(area_sc, StudyID, RowID, Temp_sc, Per_sc, SWI_sc, Fire_Int_Groups,
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
  group_by(StudyID, RowID, Temp_sc, Per_sc, SWI_sc, Fire_Int_Groups,
           Years_since_fire, Continent, area_sc) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    Height_cwm = weighted.mean(Plant_height_vegetative, cover, na.rm = TRUE),   # Actual calculation of CWMs
    Mass_cwm = weighted.mean(Seed_dry_mass, cover, na.rm = TRUE),
    Nitrogen_cwm = weighted.mean(Leaf_nitrogen, cover, na.rm = TRUE),
    Area_cwm = weighted.mean(Leaf_area_PI, cover, na.rm = TRUE),
    Long_cwm = weighted.mean(Seed_longevity, cover, na.rm = TRUE),
    Disp_cwm = weighted.mean(Dispersal_unit_dry_mass, cover, na.rm = TRUE)
  )


#Fins best model by comparing model AICs using method = "ML"

ggplot(ground.cwm, aes(x = Years_since_fire, y = Mass_cwm)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  stat_smooth(method = 'lm', aes(color = Fire_Int_Groups)) +
  facet_wrap( ~Continent)

#SEEDMASS

seed_mass_mod_lm <- glmmTMB(log(Mass_cwm) ~
                              Years_since_fire *
                              Fire_Int_Groups  *Continent +
                            (1 | StudyID),
                          data = ground.cwm,
                          weights = area_sc,
                          dispformula = ~ Fire_Int_Groups + Continent+ Years_since_fire
)

seed_mass_mod_lm2 <- glmmTMB(log(Mass_cwm) ~
                              Years_since_fire *
                              Fire_Int_Groups  *Continent +
                               SWI_sc +
                              (1 | StudyID),
                            data = ground.cwm,
                            weights = area_sc,
                            dispformula = ~ Fire_Int_Groups + Continent + Years_since_fire
)

AIC_vals <- AIC(seed_mass_mod_lm, seed_mass_mod_lm2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(seed_mass_mod_lm)
Anova(seed_mass_mod_lm, type = 'III')


#Plot predictions!
emm_sm <- emmeans(
  seed_mass_mod_lm,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground.cwm$Years_since_fire, na.rm = TRUE),
      max(ground.cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
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

seedmassplot <- ggplot(pred_grid_sm,
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
  scale_x_continuous(limits = c(1,10), 
                     breaks = c(2, 4, 6, 8, 10)) +
  ggtitle("Seed dry mass")+
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
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

ggplot(ground.cwm, aes(x = Years_since_fire, y = Height_cwm, 
                       by = Fire_Int_Groups))+
  geom_jitter(aes(color= Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

Height_mod_lm <- glmmTMB(  log(Height_cwm) ~
                             ns(Years_since_fire, df = 2) * Fire_Int_Groups *
                             Continent + 
                         Temp_sc +
                          (1 | StudyID),
                        data = ground.cwm,
                        weights = area_sc
)

Height_mod_lm2 <- glmmTMB(  log(Height_cwm) ~
                             ns(Years_since_fire, df = 2) * Fire_Int_Groups *
                             Continent + 
                             Temp_sc +
                              SWI_sc +
                             (1 | StudyID),
                           data = ground.cwm,
                           weights = area_sc,
)

AIC_vals <- AIC(Height_mod_lm, Height_mod_lm2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(Height_mod_lm)
Anova(Height_mod_lm, type = 'III')


#Plot predictions!
emm_ph <- emmeans(
  Height_mod_lm,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground.cwm$Years_since_fire, na.rm = TRUE),
      max(ground.cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
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

plantheightplot <- ggplot(pred_grid_ph,
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

ggplot(ground.cwm, aes(x = Years_since_fire, y = Nitrogen_cwm, 
                       by = Fire_Int_Groups))+
  geom_jitter(aes(color= Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

Nitro_mod_lm <- glmmTMB(log(Nitrogen_cwm) ~
                       Years_since_fire * Fire_Int_Groups * Continent +
                       (1 | StudyID),
                     weights = area_sc,
                     data = ground.cwm
)

Nitro_mod_lm2 <- glmmTMB(log(Nitrogen_cwm) ~
                        Years_since_fire * Fire_Int_Groups * Continent +
                          SWI_sc +
                        (1 | StudyID),
                      weights = area_sc,
                      data = ground.cwm
)

AIC_vals <- AIC(Nitro_mod_lm, Nitro_mod_lm2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(Nitro_mod_lm)
Anova(Nitro_mod_lm, type = 'III')

#Plot predictions!
emm_ln <- emmeans(
  Nitro_mod_lm,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground.cwm$Years_since_fire, na.rm = TRUE),
      max(ground.cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
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
    x = "Years since fire",
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
        axis.title.x = element_text(size=18),
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

ggplot(ground.cwm, aes(x = Years_since_fire, y = Area_cwm, 
                       by = Fire_Int_Groups))+
  geom_point(aes(color= Fire_Int_Groups))+
  geom_smooth(aes(color = Fire_Int_Groups))+
  facet_wrap(~Continent)

area_mod_lm <- glmmTMB(log(Area_cwm) ~ 
                      Fire_Int_Groups * Continent *
                      Years_since_fire +
                      (1 | StudyID),
                    weights = area_sc,
                    data = ground.cwm,
                    dispformula = ~ Continent
)

area_mod_lm2 <- glmmTMB(log(Area_cwm) ~ 
                         Fire_Int_Groups * Continent *
                         Years_since_fire +
                          SWI_sc +
                         (1 | StudyID),
                       weights = area_sc,
                       data = ground.cwm,
                       dispformula = ~ Continent
)


AIC_vals <- AIC(area_mod_lm, area_mod_lm2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(area_mod_lm)
Anova(area_mod_lm, type = 'III')

#Plot predictions!
emm_la <- emmeans(
  area_mod_lm,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(ground.cwm$Years_since_fire, na.rm = TRUE),
      max(ground.cwm$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
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
    x = "Years since fire",
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
        axis.title.x = element_text(size=18),
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

###########
combinedtraitplot_ground <- (plantheightplot|seedmassplot)/
  (leafareaplot|leafnitrogenplot)

combinedtraitplot_ground

ggsave(plot = combinedtraitplot_ground, filename = "traitplots_ground.png",
       dpi = 450, height = 8.45, width = 13)


####################




##################
#Tree layer

ggplot(tree.cwm, aes(x = Years_since_fire, y=log(Mass_cwm), by = Fire_Int_Groups))+
  geom_point(aes(color = Fire_Int_Groups)) + 
  geom_smooth(aes(color = Fire_Int_Groups))+
  #  scale_y_continuous(limits= c(0,)) +
  facet_wrap(~Continent)

#Use only first 10 years after fire
tree.cwm <- tree.cwm %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

tree.cwm <- tree.cwm %>%
  filter(!is.na(Mass_cwm))

table(tree.cwm$Continent, tree.cwm$Fire_Int_Groups, tree.cwm$Years_since_fire)

treeNA <- tree.cwm %>%
  filter(Continent == "North_America")

Seed_mod_NA <- glmmTMB(
  log(Mass_cwm) ~
    Years_since_fire *
    Fire_Int_Groups +
    Temp_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ Fire_Int_Groups + Years_since_fire,
  data = treeNA
)

Seed_mod_NA2 <- glmmTMB(
  log(Mass_cwm) ~
    Years_since_fire *
    Fire_Int_Groups +
    Temp_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ Fire_Int_Groups + Years_since_fire,
  data = treeNA
)

AIC_vals <- AIC(Seed_mod_NA, Seed_mod_NA2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(Seed_mod_NA)
Anova(Seed_mod_NA, type = 'III')

#Plot predictions!
emm_mass_NA <- emmeans(
  Seed_mod_NA,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(treeNA$Years_since_fire, na.rm = TRUE),
      max(treeNA$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(treeNA$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_mass_NA <- as.data.frame(emm_mass_NA)

pred_grid_mass_NA <- pred_grid_mass_NA %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_mass_NA$Fire_Int_Groups <- factor(
  pred_grid_mass_NA$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

seedmassplot_NA <- ggplot(pred_grid_mass_NA,
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
  #  scale_y_continuous(limits = c(0,50), n.breaks = 6) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Seed dry mass North America")+
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

seedmassplot_NA

treeEU <- tree.cwm %>%
  filter(Continent == "Eurasia")

treeEU$Intensity <- fct_collapse(
  treeEU$Fire_Int_Groups,
  High = "High",
  `Medium/Low` = c("Medium", "Low")
)

Seed_mod_EU <- glmmTMB(
  log(Mass_cwm) ~
    Years_since_fire *
    Intensity +
    (1 | StudyID),
  weights = area_sc,
  data = treeEU
)

Seed_mod_EU2 <- glmmTMB(
  log(Mass_cwm) ~
    Years_since_fire *
    Intensity +
    SWI_sc +
    (1 | StudyID),
  weights = area_sc,
  data = treeEU
)

AIC_vals <- AIC(Seed_mod_EU, Seed_mod_EU2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(Seed_mod_EU)
Anova(Seed_mod_EU, type = 'III')

#Plot predictions!
emm_mass_EU <- emmeans(
  Seed_mod_EU,
  ~ Years_since_fire | Intensity,
  at = list(
    Years_since_fire = seq(
      min(treeEU$Years_since_fire, na.rm = TRUE),
      max(treeEU$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
  ),
  weights = 'proportional'
)

pred_grid_mass_EU <- as.data.frame(emm_mass_EU)

pred_grid_mass_EU <- pred_grid_mass_EU %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL),
    upper_plot = pmin(upper, 60)   # cap at 50 for plotting
  )

pred_grid_mass_EU$Intensity <- factor(
  pred_grid_mass_EU$Intensity,
  levels = c("High", "Medium/Low")
)

seedmassplot_EU <- ggplot(pred_grid_mass_EU,
                          aes(x = Years_since_fire,
                              y = fit,
                              color = Intensity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Intensity),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "grey60")) + 
  scale_fill_manual(values = c("firebrick", "grey60")) + 
  labs(
    x = "Time since fire (years)",
    y = "mg",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Seed dry mass Europe")+
  theme(legend.position="left",
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

seedmassplot_EU


#Leaf nitrogen

nitro_mod_EU <- glmmTMB(
  log(Nitrogen_cwm) ~
    ns(Years_since_fire, df = 2) *
    Intensity +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ ns(Years_since_fire, df = 2),
  data = treeEU
)

nitro_mod_EU2 <- glmmTMB(
  log(Nitrogen_cwm) ~
    ns(Years_since_fire, df = 2) *
    Intensity +
    SWI_sc+
    (1 | StudyID),
  dispformula = ~ ns(Years_since_fire, df = 2),
  weights = area_sc,
  data = treeEU
)

AIC_vals <- AIC(nitro_mod_EU, nitro_mod_EU2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(nitro_mod_EU)
Anova(nitro_mod_EU, type = 'III')

#Plot predictions!
emm_nitro_EU <- emmeans(
  nitro_mod_EU,
  ~ Years_since_fire | Intensity,
  at = list(
    Years_since_fire = seq(
      min(treeEU$Years_since_fire, na.rm = TRUE),
      max(treeEU$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
    #    Temp_sc = mean(treeEU$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_nitro_EU <- as.data.frame(emm_nitro_EU)

pred_grid_nitro_EU <- pred_grid_nitro_EU %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_nitro_EU$Intensity <- factor(
  pred_grid_nitro_EU$Intensity,
  levels = c("High", "Medium/Low")
)

nitrogenplot_EU <- ggplot(pred_grid_nitro_EU,
                          aes(x = Years_since_fire,
                              y = fit,
                              color = Intensity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Intensity),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "grey60")) + 
  scale_fill_manual(values = c("firebrick", "grey60")) + 
  labs(
    x = "Time since fire (years)",
    y = "mg",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  scale_y_continuous(limits = c(10,25), n.breaks = 5)+
  theme_bw() +
  ggtitle("Leaf nitrogen Europe")+
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

nitrogenplot_EU

nitro_mod_NA <- glmmTMB(
  log(Nitrogen_cwm) ~
    Years_since_fire *
    Fire_Int_Groups +
    Temp_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ Fire_Int_Groups,
  data = treeNA
)

nitro_mod_NA2 <- glmmTMB(
  log(Nitrogen_cwm) ~
    Years_since_fire *
    Fire_Int_Groups +
    Temp_sc +
    SWI_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ Fire_Int_Groups,
  data = treeNA
)

AIC_vals <- AIC(nitro_mod_NA, nitro_mod_NA2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(nitro_mod_NA)
Anova(nitro_mod_NA, type = 'III')

#Plot predictions!
emm_nitro_NA <- emmeans(
  nitro_mod_NA,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(treeNA$Years_since_fire, na.rm = TRUE),
      max(treeNA$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(treeNA$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_nitro_NA <- as.data.frame(emm_nitro_NA)

pred_grid_nitro_NA <- pred_grid_nitro_NA %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_nitro_NA$Fire_Int_Groups <- factor(
  pred_grid_nitro_NA$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

leafnitrogenplot_NA <- ggplot(pred_grid_nitro_NA,
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  scale_y_continuous(limits = c(10,25), n.breaks = 5)+
  ggtitle("Leaf nitrogen North America")+
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


leafnitrogenplot_NA


#SLA

treeSLA <- treeEU %>%
  filter(!is.na(Area_cwm))

ggplot(treeSLA, aes(x = Years_since_fire, y = Area_cwm)) +
  geom_jitter(aes(color = Intensity))+
  geom_smooth(aes(color = Intensity))

SLA_mod_EU <- glmmTMB(
  log(Area_cwm) ~
    ns(Years_since_fire, df = 2) *
    Intensity +
    Temp_sc+
    SWI_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ ns(Years_since_fire, df = 2) + Intensity,
  data = treeSLA
)

SLA_mod_EU2 <- glmmTMB(
  log(Area_cwm) ~
    ns(Years_since_fire, df = 2) *
    Intensity +
    Temp_sc +
    SWI_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ns(Years_since_fire, df = 2) + Intensity,
  data = treeSLA
)

AIC_vals <- AIC(SLA_mod_EU, SLA_mod_EU2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(SLA_mod_EU)
Anova(SLA_mod_EU, type = 'III')

#Plot predictions!
emm_SLA_EU <- emmeans(
  SLA_mod_EU,
  ~ Years_since_fire | Intensity,
  at = list(
    Years_since_fire = seq(
      min(treeEU$Years_since_fire, na.rm = TRUE),
      max(treeEU$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(treeEU$Temp_sc, na.rm = TRUE),
    SWI_sc = mean(treeEU$SWI_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_SLA_EU <- as.data.frame(emm_SLA_EU)

pred_grid_SLA_EU <- pred_grid_SLA_EU %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_SLA_EU$Intensity <- factor(
  pred_grid_SLA_EU$Intensity,
  levels = c("High", "Medium/Low")
)

SLAplot_EU <- ggplot(pred_grid_SLA_EU,
                     aes(x = Years_since_fire,
                         y = fit,
                         color = Intensity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Intensity),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("firebrick", "grey60")) + 
  scale_fill_manual(values = c("firebrick", "grey60")) + 
  labs(
    x = "Years since fire",
    y = expression(mm^2/mg),
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  scale_y_continuous(limits = c(0, 17.5), n.breaks = 5) +
  theme_bw() +
  ggtitle("SLA Europe")+
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

SLAplot_EU

SLA_mod_NA <- glmmTMB(
  log(Area_cwm) ~
    Years_since_fire *
    Fire_Int_Groups +
    Temp_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~Fire_Int_Groups + Years_since_fire,
  data = treeNA
)

SLA_mod_NA2 <- glmmTMB(
  log(Area_cwm) ~
    Years_since_fire *
    Fire_Int_Groups +
    Temp_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~Fire_Int_Groups + Years_since_fire,
  data = treeNA
)

AIC_vals <- AIC(SLA_mod_NA, SLA_mod_NA2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(SLA_mod_NA)
Anova(SLA_mod_NA, type = 'III')

#Plot predictions!
emm_SLA_NA <- emmeans(
  SLA_mod_NA,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(treeNA$Years_since_fire, na.rm = TRUE),
      max(treeNA$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(treeNA$Temp_sc, na.rm = TRUE)
  ),
  weights = 'proportional'
)

pred_grid_SLA_NA <- as.data.frame(emm_SLA_NA)

pred_grid_SLA_NA <- pred_grid_SLA_NA %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_SLA_NA$Fire_Int_Groups <- factor(
  pred_grid_SLA_NA$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

SLAplot_NA <- ggplot(pred_grid_SLA_NA,
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
    y = expression(mm^2/mg),
    color = "Fire intensity"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  scale_y_continuous(limits = c(0, 17.5), n.breaks = 5) +
  ggtitle("SLA North America")+
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


SLAplot_NA


combinedtraitplot_tree <- (seedmassplot_EU|seedmassplot_NA)/
  (nitrogenplot_EU|leafnitrogenplot_NA)/
  (SLAplot_EU|SLAplot_NA)

combinedtraitplot_tree

ggsave(plot = combinedtraitplot_tree, filename = "traitplot_tree.TIFF",
       dpi = 450, height = 10.52, width =13)








#####BRYOPHYTE TRAITS
traits_bryo$Temp_sc <- as.numeric(scale(traits_bryo$Avg_Temp, center = TRUE, scale = TRUE))
traits_bryo$Per_sc <- as.numeric(scale(traits_bryo$AvgPer, center = TRUE, scale = TRUE))
traits_bryo$SWI_sc <- as.numeric(scale(traits_bryo$SWI, center = TRUE, scale = TRUE))

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
traits_bryo_df$area_sc <- traits_bryo_df$studysize / mean(traits_bryo_df$studysize)

traits_bryo_df <- traits_bryo_df %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

#Calculate community weighted means
bryo.cwm <-   # New dataframe where we can inspect the result
  traits_bryo_df %>%   # First step in the next string of statements
  group_by(StudyID, RowID, Temp_sc, Per_sc, SWI_sc, Fire_Int_Groups,
           Years_since_fire, Continent, area_sc) %>%   # Groups the summary file by Plot number
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

Length_mod <- glmmTMB(
  log(Length_cwm) ~
    Fire_Int_Groups * Continent *
    Years_since_fire +
    Temp_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ Fire_Int_Groups,
  data = bryo.cwm
)

Length_mod2 <- glmmTMB(
  log(Length_cwm) ~
    Fire_Int_Groups * Continent *
    Years_since_fire +
    Temp_sc +
    SWI_sc +
    (1 | StudyID),
  weights = area_sc,
  dispformula = ~ Fire_Int_Groups,
  data = bryo.cwm
)

AIC_vals <- AIC(Length_mod, Length_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

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
    y = "mm",
    color = "Fire Intensity"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  ggtitle("Leafy shoot length")+
  theme(legend.position="right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
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

#To increase data, we combine rare and occassional sporophytes

traits_bryo_df$Sporophyte_frequency2 <- fct_collapse(
  traits_bryo_df$Sporophyte_frequency,
  Abundant = "Abundant",
  `Occasional/Rare` = c("Occasional", "Rare")
)

spor_mod <- glmmTMB(Sporophyte_frequency2 ~ 
                      Years_since_fire * 
                      Fire_Int_Groups * Continent +
                      Temp_sc +
                      (1|StudyID),
                    family = binomial,
                    data = traits_bryo_df)

spor_mod2 <- glmmTMB(Sporophyte_frequency2 ~ 
                       Years_since_fire *
                       Fire_Int_Groups * Continent +
                       Temp_sc +
                       (1|StudyID),
                     family = binomial,
                     data = traits_bryo_df)

AIC_vals <- AIC(spor_mod, spor_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

#Plot predictions!
emm_ord <- emmeans(
  spor_mod,
  type = "response",
  ~ Years_since_fire | Continent * Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(traits_bryo_df$Years_since_fire, na.rm = TRUE),
      max(traits_bryo_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Continent = unique(traits_bryo_df$Continent),
    Fire_Int_Groups = c("High", "Medium", "Low"),
    Temp_sc = mean(traits_bryo_df$Temp_sc, na.rm = TRUE)
  ))

pred_grid_spor <- as.data.frame(emm_ord)

pred_grid_spor <- pred_grid_spor %>%
  mutate(
    Probability   = prob,
    lower = exp(asymp.LCL),
    upper = exp(asymp.UCL)
  )

pred_grid_spor <- pred_grid_spor %>%
  mutate(
    Probability_Rare = Probability,
    Probability_Abundant = 1 - Probability
  )

pred_grid_long <- pred_grid_spor %>%
  pivot_longer(
    cols = c(Probability_Abundant, Probability_Rare),
    names_to = "Sporophyte_frequency",
    values_to = "Probability_new"
  ) 


pred_grid_long$Fire_Int_Groups <- factor(
  pred_grid_long$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

pred_grid_long <- pred_grid_long %>%
  mutate(
    Sporophyte_frequency = case_when(
      Sporophyte_frequency == "Probability_Abundant" ~ "Abundant",
      Sporophyte_frequency == "Probability_Rare" ~ "Occassional/Rare"
    )
  )

mossspor_plot <- ggplot(pred_grid_long,
                        aes(x = Years_since_fire,
                            y = Probability_new,
                            color = Sporophyte_frequency,
                            fill = Sporophyte_frequency)) +
  geom_ribbon(aes(ymin = 0, ymax = Probability_new), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_grid(Continent ~ Fire_Int_Groups,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"
               )
             )) +
  scale_color_manual(values = c("black", "#E69F00")) +
  scale_fill_manual(values = c("black", "#E69F00")) +
  labs(x = "Years since fire", y = "Predicted probability",
       color = "Sporophyte frequency",
       fill = "Sporophyte frequency") +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  ggtitle("Fire intensity")+
  theme(legend.position="right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        legend.direction='vertical',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA"))

mossspor_plot

##Primary lifeform

ggplot(traits_bryo_df, aes(x = Primary_lifeform, y = Years_since_fire)) +
  geom_point() +
  facet_grid(Fire_Int_Groups~Continent)

nom_brms <- brm(
  formula = Primary_lifeform ~ 
    Years_since_fire + 
    Fire_Int_Groups + 
    Per_sc + 
    (1 | StudyID),
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
  StudyID = NA  # population-level predictions; remove random effect
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

pred_long_life$Primary_lifeform <- factor(
  pred_long_life$Primary_lifeform,
  levels = c("Estimate.P(Y = Mat)", "Estimate.P(Y = Tuft)", "Estimate.P(Y = Weft)"),
  labels = c("Mat", "Tuft", "Weft")
)

primlife_plot <- ggplot(pred_long_life,
                        aes(x = Years_since_fire,
                            y = Probability,
                            color = Primary_lifeform,
                            fill = Primary_lifeform)) +
  geom_ribbon(aes(ymin = 0, ymax = Probability), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap( ~ Fire_Int_Groups) +
  scale_color_manual(values = c("black", "#E69F00", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#009E73")) +
  labs(x = "Years since fire", y = "Predicted probability",
       color = "Primary lifeform",
       fill = "Primary lifeform") +
  theme_bw() +
  scale_x_continuous(limits = c(1,10), n.breaks = 6)+
  ggtitle("Fire intensity")+
  theme(legend.position="right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        legend.direction='vertical',
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

primlife_plot

mossplots <-  mosslength_plot / mossspor_plot / primlife_plot +
  plot_layout(heights = c(1, 2, 1))

mossplots

ggsave(mossplots, filename = "mossplots.TIFF", height = 12.67,
       width = 13, dpi = 450)
