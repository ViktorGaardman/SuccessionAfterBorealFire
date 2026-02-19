#Step 1. Load packages
library(tidyverse)
library(car)
#library(ggeffects)
library(glmmTMB)
library(DHARMa)
library(patchwork)
#library(mgcv)
library(emmeans)
library(splines)
library(lme4)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Clean_Data.csv", sep = ";")

df$Fire_Int_Groups <- factor(
  df$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

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

df_long <- df_long %>%
  mutate(
    area_sc = studysize / mean(studysize)
  )

#Use only first 10 years after fire
df_long <- df_long %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 10)

#Sum cover per plantgroup
summed_df <- df_long %>%
  group_by(StudyID, Years_since_fire, Fire_Int_Groups,
          Continent, PlantGroup, Per_sc, Temp_sc, area_sc, studysize) %>%
  summarize(
    Total_cov = sum(cover),
      .groups = "drop"
  )

#Fit linear mixed models

#Herbs

herb_df <- summed_df %>%
  filter(PlantGroup == "herb")

ggplot(herb_df, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  facet_wrap(~Continent)

herb_mod <- glmmTMB(log(Total_cov) ~
                      Years_since_fire * Continent +
                      Years_since_fire * Fire_Int_Groups +
                      Continent * Fire_Int_Groups +
                      (1 | StudyID),
                    weights = area_sc,
                 data = herb_df,
                 dispformula = ~ Continent)  

herb_mod2 <- glmmTMB(log(Total_cov) ~
                       ns(Years_since_fire, df = 2) * Continent +
                       ns(Years_since_fire, df = 2) * Fire_Int_Groups +
                       Continent * Fire_Int_Groups +
                       (1 | StudyID),
                     weights = area_sc,
                     data = herb_df,
                     dispformula = ~ Continent)  

AIC_vals <- AIC(herb_mod, herb_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(herb_mod)
Anova(herb_mod, type = 'III')

#Plot predictions!
herb_pred <- emmeans(
  herb_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(herb_df$Years_since_fire, na.rm = TRUE),
      max(herb_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
  ),
  weights = 'proportional'
)

pred_grid_herb <- as.data.frame(herb_pred)

pred_grid_herb <- pred_grid_herb %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

pred_grid_herb$Fire_Int_Groups <- factor(
  pred_grid_herb$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predherbplot <- ggplot(pred_grid_herb,
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
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


#Dwarfshrubs

dwarf_df <- summed_df %>%
  filter(PlantGroup == "dwarfshrub")

ggplot(dwarf_df, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  facet_wrap(~Continent)

dwarf_mod <- glmmTMB(log(Total_cov) ~
                       Years_since_fire * Continent *Fire_Int_Groups+
                       Temp_sc +
                      (1 | StudyID),
                    weights = area_sc,
                    dispformula = ~ Fire_Int_Groups,
                    data = dwarf_df)  

dwarf_mod2 <- glmmTMB(log(Total_cov) ~
                       Years_since_fire * Continent *Fire_Int_Groups+
                       Temp_sc +
                       (1 | StudyID),
                     weights = area_sc,
                     dispformula = ~ Fire_Int_Groups,
                     data = dwarf_df)     

AIC_vals <- AIC(dwarf_mod, dwarf_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(dwarf_mod)
Anova(dwarf_mod, type = 'III')

#Plot predictions!
dwarf_pred <- emmeans(
  dwarf_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(herb_df$Years_since_fire, na.rm = TRUE),
      max(herb_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Temp_sc = mean(herb_df$Temp_sc, na.rm = TRUE)
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
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Years since fire",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

preddwarfplot



#Graminoid

grass_df <- summed_df %>%
  filter(PlantGroup == "graminoid")

ggplot(grass_df, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  facet_wrap(~Continent)

grass_mod <- glmmTMB(log(Total_cov) ~
                       Years_since_fire * Continent *Fire_Int_Groups +
                       (1 | StudyID),
                     weights = area_sc,
                    dispformula = ~ Continent,
                     data = grass_df)  

grass_mod2 <- glmmTMB(log(Total_cov) ~
                        Years_since_fire * Continent *Fire_Int_Groups +
                        (1 | StudyID),
                      weights = area_sc,
                      dispformula = ~ Continent + Years_since_fire,
                      data = grass_df)     

AIC_vals <- AIC(grass_mod, grass_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(grass_mod)
Anova(grass_mod, type = 'III')

#Plot predictions!
grass_pred <- emmeans(
  grass_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(herb_df$Years_since_fire, na.rm = TRUE),
      max(herb_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
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

predgrassplot <- ggplot(pred_grid_grass,
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  ggtitle("Graminoids")+
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

predgrassplot


#Shrubs

#Only shrubs from North America and group Med and Low intensity
shrub_df <- summed_df %>%
  filter(PlantGroup == "shrub") %>%
  filter(Continent == "North_America")

shrub_df$Intensity <- fct_collapse(
  shrub_df$Fire_Int_Groups,
  High = "High",
  `Medium/Low` = c("Medium", "Low")
)


ggplot(shrub_df, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Intensity)) +
  geom_jitter(aes(color = Intensity)) +
  scale_y_continuous(limits = c(0,10))

shrub_mod <- glmmTMB(log(Total_cov) ~
                       Years_since_fire +
                       (1 | StudyID),
                     weights = area_sc,
                     dispformula = ~ Years_since_fire,
                     data = shrub_df)  

shrub_mod2 <- glmmTMB(log(Total_cov) ~
                       Years_since_fire +
                       (1 | StudyID),
                     weights = area_sc,
                     dispformula = ~ Years_since_fire,
                     data = shrub_df)    

AIC_vals <- AIC(shrub_mod, shrub_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(shrub_mod)
Anova(shrub_mod, type = 'III')

#Plot predictions!
shrub_pred <- emmeans(
  shrub_mod,
  ~ Years_since_fire,
  at = list(
    Years_since_fire = seq(
      min(herb_df$Years_since_fire, na.rm = TRUE),
      max(herb_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
  ),
  weights = 'proportional'
)

pred_grid_shrub <- as.data.frame(shrub_pred)

pred_grid_shrub <- pred_grid_shrub %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

predshrubplot <- ggplot(pred_grid_shrub,
                        aes(x = Years_since_fire,
                            y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Years since fire",
    y = "Predicted cover"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Shrubs North America")+
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

predshrubplot


#Trees

#Split data by continent due to large difference in data-amounts
tree_NA <- summed_df %>%
  filter(PlantGroup == "tree") %>%
  filter(Continent == "North_America")

ggplot(tree_NA, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_jitter(aes(color = Fire_Int_Groups))

tree_mod_NA <- glmmTMB(log(Total_cov) ~
                         ns(Years_since_fire, df = 2) +
                         Fire_Int_Groups +
                       (1 | StudyID),
                     weights = area_sc,
                     dispformula = ~ ns(Years_since_fire, df = 2) + Fire_Int_Groups,
                     data = tree_NA)  

tree_mod_NA2 <- glmmTMB(log(Total_cov) ~
                          ns(Years_since_fire, df = 2) +
                          Fire_Int_Groups +
                         (1 | StudyID),
                       weights = area_sc,
                       dispformula = ~ ns(Years_since_fire, df = 2),
                       data = tree_NA)   

AIC_vals <- AIC(tree_mod_NA, tree_mod_NA2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(tree_mod_NA)
Anova(tree_mod_NA, type = 'III')

#Plot predictions!
tree_pred_NA <- emmeans(
  tree_mod_NA,
  ~ Years_since_fire | Fire_Int_Groups,
  at = list(
    Years_since_fire = seq(
      min(tree_NA$Years_since_fire, na.rm = TRUE),
      max(tree_NA$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
  ),
  weights = 'proportional'
)

tree_pred_NA <- as.data.frame(tree_pred_NA)

tree_pred_NA <- tree_pred_NA %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

tree_pred_NA$Fire_Int_Groups <- factor(
  tree_pred_NA$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predtreeplot_NA <- ggplot(tree_pred_NA,
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Trees North America")+
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

predtreeplot_NA

ggplot(tree_EU, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_jitter(aes(color = Fire_Int_Groups))

tree_EU <- summed_df %>%
  filter(PlantGroup == "tree") %>%
  filter(Continent == "Eurasia")

tree_EU <- tree_EU %>%
  filter(Years_since_fire >= 1, Years_since_fire <= 5)

#tree_EU$Intensity <- fct_collapse(
#  tree_EU$Fire_Int_Groups,
#  High = "High",
#  `Medium/Low` = c("Medium", "Low")
#)


tree_mod_EU <- glmmTMB(log(Total_cov) ~
                          Years_since_fire +
                         Intensity +
                          (1 | StudyID),
                        weights = area_sc,
                        data = tree_EU)   

tree_mod_EU2 <- glmmTMB(log(Total_cov) ~
                          Years_since_fire + 
                          (1 | StudyID),
                        weights = area_sc,
                        data = tree_EU)   

AIC_vals <- AIC(tree_mod_EU, tree_mod_EU2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(tree_mod_EU)
Anova(tree_mod_EU, type = 'III')

#Plot predictions!
tree_pred_EU <- emmeans(
  tree_mod_EU,
  ~ Years_since_fire | Intensity,
  at = list(
    Years_since_fire = seq(
      min(tree_EU$Years_since_fire, na.rm = TRUE),
      max(tree_EU$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
  ),
  weights = 'proportional'
)

tree_pred_EU <- as.data.frame(tree_pred_EU)

tree_pred_EU <- tree_pred_EU %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )

tree_pred_EU$Intensity <- factor(
  tree_pred_EU$Intensity,
  levels = c("High", "Medium/Low")
)

predtreeplot_EU <- ggplot(tree_pred_EU,
                          aes(x = Years_since_fire,
                              y = fit,
                              color = Intensity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Intensity),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +

  labs(
    x = "Years since fire",
    y = "Predicted cover",
    color = "Fire intensity"
  ) +
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  ggtitle("Trees Europe")+
  theme(legend.position="right",
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

predtreeplot_EU



#Mosses

moss_df <- summed_df %>%
  filter(PlantGroup == "bryophyte")

ggplot(moss_df, aes(x = Years_since_fire, y = Total_cov)) +
  geom_smooth(aes(color = Fire_Int_Groups)) +
  geom_jitter(aes(color = Fire_Int_Groups)) +
  facet_wrap(~Continent)

moss_mod <- glmmTMB(log(Total_cov) ~
                       Years_since_fire * Continent *Fire_Int_Groups +
                       (1 | StudyID),
                     weights = area_sc,
                     dispformula = ~ Continent + Fire_Int_Groups + Years_since_fire,
                     data = moss_df)  

moss_mod2 <- glmmTMB(log(Total_cov) ~
                       Years_since_fire * Continent *Fire_Int_Groups +
                        (1 | StudyID),
                      weights = area_sc,
                      dispformula = ~ Continent + Fire_Int_Groups + Years_since_fire,
                      data = moss_df)     

AIC_vals <- AIC(moss_mod, moss_mod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals

summary(moss_mod)
Anova(moss_mod, type = 'III')

#Plot predictions!
moss_pred <- emmeans(
  moss_mod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(herb_df$Years_since_fire, na.rm = TRUE),
      max(herb_df$Years_since_fire, na.rm = TRUE),
      length.out = 40
    )
  ),
  weights = 'proportional'
)

pred_grid_moss <- as.data.frame(moss_pred)

pred_grid_moss <- pred_grid_moss %>%
  mutate(
    fit   = exp(emmean),
    lower = exp(lower.CL),
    upper = exp(upper.CL)
  )


pred_grid_moss$Fire_Int_Groups <- factor(
  pred_grid_moss$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

predmossplot <- ggplot(pred_grid_moss,
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
  scale_x_continuous(limits = c(1,10), n.breaks = 6) +
  theme_bw() +
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  ggtitle("Mosses")+
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
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = "NA")) 

predmossplot


totalcoverplot <- (predtreeplot_NA|predshrubplot)/
  (predherbplot|preddwarfplot)/(predgrassplot|predmossplot)

totalcoverplot

ggsave(totalcoverplot, filename = "TotalCover.TIFF", 
       dpi = 450, height = 10.52, width =13)




####SPECIESWISE MODELS
#Standardize cover

df_long$coverstd <- (df_long$cover + 0.01) / 101

df_long$base <- as.factor(df_long$base)

df_long$species <- as.factor(df_long$species)


herbs_df_sp <- df_long %>%
  filter(PlantGroup == "herb") %>%
  droplevels()

#Filter out all species that appear less than five times, in less than five timesteps,
#and less than two studies
herbs_filtered <- herbs_df_sp %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 4, 
         n_distinct(Years_since_fire) >= 4,
         n_distinct(StudyID) >=2) %>%
  ungroup() %>%
  droplevels()

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
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  }
  
    pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
     df_sp <- df_sub[df_sub$species == sp, ]
      data.frame(
        Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
        species = sp,
        StudyID = NA,
        area_sc = median(df_sub$area_sc, na.rm = TRUE)
      )
    }))
  
  # Prediction grid
#  pred_grid <- expand.grid(
#    Years_since_fire = seq(
#      min(df_sub$Years_since_fire, na.rm = TRUE),
#      max(df_sub$Years_since_fire, na.rm = TRUE),
#      length.out = 40
#    ),
#    species = levels(df_sub$species),
#    Temp_sc = mean(df_sub$Temp_sc, na.rm = TRUE),
#    StudyID = levels(df_sub$StudyID)[1],
#    weight_sc = median(df_sub$weight_sc)
#  )
  
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

dwarf_df_sp <- df_long %>%
  filter(PlantGroup == "dwarfshrub")

#Filter out all species that appear less than three times
dwarf_filtered <- dwarf_df_sp %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 4, 
         n_distinct(Years_since_fire) >= 4,
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
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  }
  
    pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
     df_sp <- df_sub[df_sub$species == sp, ]
      data.frame(
        Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
        species = sp,
        StudyID = NA,
        area_sc = median(df_sub$area_sc, na.rm = TRUE)
      )
    }))

  
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
      title = paste("Dwarfshrubs", continent_level, fire_level, "intensity")
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

grass_df_sp <- df_long %>%
  filter(PlantGroup == "graminoid")

#Filter out all species that appear less than three times
grass_filtered <- grass_df_sp %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 4, 
         n_distinct(Years_since_fire) >= 4,
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
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  }
  
    pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
     df_sp <- df_sub[df_sub$species == sp, ]
      data.frame(
        Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
        species = sp,
        StudyID = NA,
        area_sc = median(df_sub$area_sc, na.rm = TRUE)
      )
    }))

  
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
moss_df_sp <- df_long %>%
  filter(PlantGroup == "bryophyte")

#Filter out all species that appear less than four times
moss_filtered <- moss_df_sp %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 4, 
         n_distinct(Years_since_fire) >= 4,
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
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  }
  
    pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
     df_sp <- df_sub[df_sub$species == sp, ]
      data.frame(
        Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
        species = sp,
        StudyID = NA,
        area_sc = median(df_sub$area_sc, na.rm = TRUE)
      )
    }))
  
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


####
#Trees
tree_df_sp <- df_long %>%
  filter(PlantGroup == "tree")

#Filter out all species that appear less than four times
tree_filtered <- tree_df_sp %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 4, 
         n_distinct(Years_since_fire) >= 4,
         n_distinct(StudyID) >=2)%>%
  ungroup() %>%
  droplevels()

prediction_plot <- function(tree_filtered, continent_level, fire_level) {
  
  # Filter data
  df_sub <- tree_filtered %>%
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
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  }
  
  pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
    df_sp <- df_sub[df_sub$species == sp, ]
    data.frame(
      Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
      species = sp,
      StudyID = NA,
      area_sc = median(df_sub$area_sc, na.rm = TRUE)
    )
  }))
  
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
      title = paste("Trees", continent_level, fire_level, "intensity")
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

treeEUhigh_plot <- prediction_plot(tree_filtered, "Eurasia", "High")
treeEUmed_plot  <- prediction_plot(tree_filtered, "Eurasia", "Medium")
treeEUlow_plot  <- prediction_plot(tree_filtered, "Eurasia", "Low")
treeNAhigh_plot  <- prediction_plot(tree_filtered, "North_America", "High")
treeNAmed_plot  <- prediction_plot(tree_filtered, "North_America", "Medium")
treeNAlow_plot  <- prediction_plot(tree_filtered, "North_America", "Low")

EUtreeplots <- treeEUhigh_plot / treeEUmed_plot / treeEUlow_plot
EUtreeplots

NAtreeplots <- treeNAhigh_plot / treeNAmed_plot / treeNAlow_plot
NAtreeplots

####
#Shrubs
shrub_df_sp <- df_long %>%
  filter(PlantGroup == "shrub")

#Filter out all species that appear less than four times
shrub_filtered <- shrub_df_sp %>%
  group_by(species, Continent, Fire_Int_Groups) %>%
  filter(n() >= 4, 
         n_distinct(Years_since_fire) >= 4,
         n_distinct(StudyID) >=2)%>%
  ungroup() %>%
  droplevels()

prediction_plot <- function(shrub_filtered, continent_level, fire_level) {
  
  # Filter data
  df_sub <- shrub_filtered %>%
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
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  } else {
    # Single-species model
    mod <- glmmTMB(
      coverstd ~ 
        Years_since_fire +
        (1 | StudyID),
      family = beta_family(link = "logit"),
      weights = area_sc,
      data = df_sub
    )
  }
  
  pred_grid <- do.call(rbind, lapply(levels(df_sub$species), function(sp) {
    df_sp <- df_sub[df_sub$species == sp, ]
    data.frame(
      Years_since_fire = seq(min(df_sp$Years_since_fire), max(df_sp$Years_since_fire), length.out = 40),
      species = sp,
      StudyID = NA,
      area_sc = median(df_sub$area_sc, na.rm = TRUE)
    )
  }))
  
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
      title = paste("Shrubs", continent_level, fire_level, "intensity")
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

shrubNAhigh_plot  <- prediction_plot(shrub_filtered, "North_America", "High")
shrubNAmed_plot  <- prediction_plot(shrub_filtered, "North_America", "Medium")
shrubNAlow_plot  <- prediction_plot(shrub_filtered, "North_America", "Low")

NAshrubplots <- shrubNAhigh_plot / shrubNAmed_plot / shrubNAlow_plot
NAshrubplots


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
ggsave(EUtreeplots, filename = "EUTree.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(NAtreeplots, filename = "NATree.png",
       dpi = 300, height = 10.56, width = 13)
ggsave(NAshrubplots, filename = "NAShrub.png",
       dpi = 300, height = 10.56, width = 13)
