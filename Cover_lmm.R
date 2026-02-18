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
                      Years_since_fire * Continent * Fire_Int_Groups +
                      (1 | StudyID),
                    weights = area_sc,
                 data = herb_df,
                 dispformula = ~ Continent)  

herb_mod2 <- glmmTMB(log(Total_cov) ~
                       Years_since_fire * Continent +
                       Years_since_fire * Fire_Int_Groups +
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
  facet_wrap(~Continent) +
  ggtitle("Herbs")+
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
  facet_wrap(~Continent) +
  ggtitle("Dwarfshrubs")+
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

preddwarfplot
















#Where models come to die

#Best gls models on summed cover 0-100 for sum of all three dominance levels
herbmod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups * Continent +
    poly(Temp_sc, 3) +
    poly(Per_sc, 2),
  data = herbs_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varExp(form = ~ weight_sc),
  method      = "REML"
)

dwarfmod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups * Continent +
    poly(Temp_sc, 2) +
    poly(Per_sc, 3),
  data = dwarf_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varPower(form = ~ weight_sc),
  method      = "REML"
)



#Old prediction code (for gls())


mod_fixed <- update(grassmod, weights = varFixed(~ 1 / weight_sc), method = "ML")
mod_power <- update(grassmod, weights = varPower(form = ~ weight_sc), method = "ML")
mod_exp   <- update(grassmod, weights = varExp(form = ~ weight_sc), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(grassmod, type = "normalized"))
qqline(resid(grassmod, type = "normalized"))
plot(grassmod, resid(., type = "normalized") ~ fitted(.))
summary(grassmod)
Anova(grassmod, type = 'III')

#Plot predictions!
grass_pred <- emmeans(
  grassmod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(grass_dom$Years_since_fire, na.rm = TRUE),
      max(grass_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(grass_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(grass_dom$Temp_sc, na.rm = TRUE)
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

predgramplot <- ggplot(pred_grid_grass,
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
    y = "Predicted summed cover",
    color = "Fire intensity"
  ) +
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

summed_cover <- df_long %>%
  group_by(across(-c(base, species, cover, coverstd))) %>%
  summarize(
    TotalCover = sum(cover),
    .groups = "drop"
  )

herbs_dom <- summed_cover %>%
  filter(PlantGroup == "herb")


herbmod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 3) * Fire_Int_Groups * Continent +
    poly(Temp_sc, 3) +
    poly(Per_sc, 2),
  data = herbs_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varExp(form = ~ weight_sc),
  method      = "ML"
)

AIC_vals <- AIC(herbmod, herbmod2)
AIC_vals$delta <- AIC_vals$AIC - min(AIC_vals$AIC)
AIC_vals


mod_fixed <- update(herbmod, weights = varFixed(~ 1 / weight_sc), method = "ML")
mod_power <- update(herbmod, weights = varPower(form = ~ weight_sc), method = "ML")
mod_exp   <- update(herbmod, weights = varExp(form = ~ weight_sc), method = "ML")
AIC(mod_fixed, mod_power, mod_exp)

qqnorm(resid(herbmod, type = "normalized"))
qqline(resid(herbmod, type = "normalized"))
plot(herbmod, resid(., type = "normalized") ~ fitted(.))
summary(herbmod)
Anova(herbmod, type = 'III')

#Plot predictions!
herb_pred <- emmeans(
  herbmod,
  ~ Years_since_fire | Fire_Int_Groups * Continent,
  at = list(
    Years_since_fire = seq(
      min(herbs_dom$Years_since_fire, na.rm = TRUE),
      max(herbs_dom$Years_since_fire, na.rm = TRUE),
      length.out = 40
    ),
    Per_sc = mean(herbs_dom$Per_sc, na.rm = TRUE),
    Temp_sc = mean(herbs_dom$Temp_sc, na.rm = TRUE)
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
  facet_wrap(~ Continent,
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  scale_fill_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(
    x = "Time since fire (years)",
    y = "Predicted summed cover",
    color = "Fire intensity"
  ) +
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

predherbplot

treemod <- gls(
  log(TotalCover) ~
    poly(Years_since_fire, 2) +
    Fire_Int_Groups * Continent +
    poly(Temp_sc, 3),
  data = tree_dom,
  correlation = corCompSymm(form = ~ 1 | StudyID),
  weights     = varExp(form = ~ weight_sc),
  method      = "REML"
)
