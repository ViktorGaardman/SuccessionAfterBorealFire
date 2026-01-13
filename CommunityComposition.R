#Step 1. Load packages
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(ggeffects)
library(performance)
library(grid)
library(ggcorrplot)
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

#Add all metadata together
TRY_Traits <- read.csv("TRY_Clean.csv", sep = ";")

Perc_Data <- read.csv("site_climatological_percipitation_WorldClim_v2.csv",
                      sep = ";")

Temp_Data <- read.csv("site_climatological_temperature_WorldClim_v2.csv",
                      sep = ";")

metadata_Temp <- metadata %>%
  left_join(Temp_Data, by = "Title")

metadata_clean <- metadata_Temp %>%
  left_join(Perc_Data, by = "Title")

metadata_clean <- metadata_clean %>%
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
  
table(metadata_clean$YSF_interval, useNA = "ifany")


metadata_clean$YSF_interval <- factor(
  metadata_clean$YSF_interval,
  levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
  ordered = TRUE
)

metadata_clean <- metadata_clean %>%
  mutate(across(c("Continent", "Fire_Int_Groups"), as.factor))

#Check collinearity of continuous variables

#Plant_Div due to high correlations
pond_corr <- metadata_clean %>% 
  select(SWI, Avg_Temp, AvgPer, Latitude) %>%
  cor(use = "pairwise.complete.obs")  # Avoids NA issues

# Plot the correlation matrix
ggcorrplot(pond_corr, lab = TRUE, type = "lower", hc.order = TRUE)

#Split dataset into trees/shrubs, groundlayer plants, and mosses

Trees_long <- species_long %>%
  filter(base %in% c("Dominant_tree_1", "Dominant_tree_2",
                     "Dominant_tree_3", "Dominant_tree_4",
                     "Dominant_shrub_1", "Dominant_shrub_2",
                     "Dominant_shrub_3"))

#Simply matrix
Trees_long <- Trees_long %>%
  mutate(species = if_else(species == "Rosa_sp.",
                           "Rosa_acicularis",
                           species))

Trees_long <- Trees_long %>%
  mutate(species = if_else(species == "Sorbus_americana",
                           "Sorbus_sp.",
                           species))

Trees_long <- Trees_long %>%
  mutate(species = if_else(species == "Amelanchier_sp.",
                           "Amelanchier_alnifolia",
                           species))

Trees_long <- Trees_long %>%
  mutate(species = if_else(species == "Populus_sp.",
                           "Populus_tremuloides",
                           species))

  
Ground_long <-species_long %>%
  filter(base %in% c("Dominant_dwarfshrub_1", "Dominant_dwarfshrub_2",
                     "Dominant_dwarfshrub_3", "Dominant_herb_1",
                     "Dominant_herb_2", "Dominant_herb_3",
                     "Dominant_graminoid_1", "Dominant_graminoid_2",
                     "Dominant_graminoid_3"))

Mosses_long <- species_long %>%
  filter(base %in% c("Dominant_bryohpyte_1", "Dominant_bryophyte_2",
                     "Dominant_bryophyte_3"))


#TREES
Tree_matrix <- Trees_long %>%
  pivot_wider(
    id_cols = c(RowID, StudyID),
    names_from = species,      # base = original Dominant_X_Y column
    values_from = cover,
    values_fill = 0
  )

#Create NMDS
Tree_cols <- setdiff(names(Tree_matrix), c("RowID", "StudyID"))

Tree_filter <- Tree_matrix[
  rowSums(Tree_matrix[, Tree_cols] != 0, na.rm = TRUE) > 0,
]

metadata_tree <- Tree_filter %>%
  left_join(metadata_clean, "RowID")

Tree_info <- Tree_filter %>%
  select(c(StudyID, RowID))

Tree_clean <- Tree_filter %>%
  select(- c(StudyID, RowID))

NMDS_tree <- metaMDS(Tree_clean, distance = "bray", k = 2, trymax = 1000)

NMDS_tree <- metaMDS(Tree_clean, distance = "bray", k = 2, trymax = 1000,
                    previous.best = NMDS_tree)

#extract the site scores
datascores_T = as.data.frame(scores(NMDS_tree)$sites)  

#add metadata
datascores_T$Fire_Int_Groups = metadata_tree$Fire_Int_Groups
datascores_T$YSF_interval = metadata_tree$YSF_interval
datascores_T$SWI = metadata_tree$SWI
datascores_T$Avg_Temp = metadata_tree$Avg_Temp
datascores_T$SWI = metadata_tree$SWI
datascores_T$AvgPer = metadata_tree$AvgPer
datascores_T$Continent = metadata_tree$Continent

datascores_T$Fire_Int_Groups <- factor(
  datascores_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

ysf_paths <- datascores_T %>%
  filter(Continent == "Eurasia") %>%
  group_by(Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Fire_Int_Groups, YSF_interval)


#plot by continent, YSF, and Fire intensity
Tree_plot <- ggplot(subset(datascores_T, Continent == "Eurasia"), aes(x = NMDS1, y = NMDS2, 
                                                                      color = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Fire_Int_Groups)) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("Eurasia trees") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  annotate("text", x = max(datascores_T$NMDS1), 
           y = min(datascores_T$NMDS2), 
           label = paste("Stress =", round(NMDS_tree$stress, 3)), 
           hjust = 1.1, vjust = 0.1, size = 5) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(color = "Fire intensity") +
  geom_path(
    data = ysf_paths,
    aes(x = NMDS1, y = NMDS2, group = Fire_Int_Groups, color = Fire_Int_Groups),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    linewidth = 1,
    show.legend = FALSE
  )

Tree_plot

ysf_path_NA <- datascores_T %>%
  filter(Continent == "North_America") %>%
  group_by(Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Fire_Int_Groups, YSF_interval)


#plot by continent, YSF, and Fire intensity
Tree_plot_NA <- ggplot(subset(datascores_T, Continent == "North_America"), aes(x = NMDS1, y = NMDS2, 
                                                                               color = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Fire_Int_Groups)) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("North America trees")+
  theme(legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction='vertical',
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(color = "Fire intensity") +
  geom_path(
    data = ysf_path_NA,
    aes(x = NMDS1, y = NMDS2, group = Fire_Int_Groups, color = Fire_Int_Groups),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    linewidth = 1,
    show.legend = FALSE
  )

Tree_plot_NA

##Permanova with StudyID as a random factor
permutations <- with(metadata_tree, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_tree <- vegdist(Tree_clean, method = "bray")

###Fit permanova model
Permanova_tree <- adonis2(dist_tree ~ Continent*Years_since_fire*Fire_Int_Groups +
                            Avg_Temp +
                            SWI +
                            Latitude +
                            AvgPer,
                          data=metadata_tree,
                         permutations=permutations, method="bray")

Permanova_tree

#Check assumption of homogeneity of multivariate dispersion OK
Assumption <- anova(betadisper(dist_tree, metadata_tree$Fire_Int_Groups))
Assumption2 <- anova(betadisper(dist_tree, metadata_tree$Continent))


###HERB LAYER

Herb_matrix <- Ground_long %>%
  pivot_wider(
    id_cols = c(RowID, StudyID),
    names_from = species,      # base = original Dominant_X_Y column
    values_from = cover,
    values_fill = 0
  )

#Create NMDS
Herb_cols <- setdiff(names(Herb_matrix), c("RowID", "StudyID"))

Herb_filter <- Herb_matrix[
  rowSums(Herb_matrix[, Herb_cols] != 0, na.rm = TRUE) > 0,
]

metadata_herb <- Herb_filter %>%
  left_join(metadata_clean, "RowID")

Herb_info <- Herb_filter %>%
  select(c(StudyID, RowID))

Herb_clean <- Herb_filter %>%
  select(- c(StudyID, RowID))

NMDS_herb <- metaMDS(Herb_clean, distance = "bray", k = 2, trymax = 1000)

NMDS_herb <- metaMDS(Herb_clean, distance = "bray", k = 2, trymax = 1000,
                     previous.best = NMDS_herb)

#extract the site scores
datascores_H = as.data.frame(scores(NMDS_herb)$sites)  

#add metadata
datascores_H$Fire_Int_Groups = metadata_herb$Fire_Int_Groups
datascores_H$YSF_interval = metadata_herb$YSF_interval
datascores_H$SWI = metadata_herb$SWI
datascores_H$Avg_Temp = metadata_herb$Avg_Temp
datascores_H$SWI = metadata_herb$SWI
datascores_H$AvgPer = metadata_herb$AvgPer
datascores_H$Continent = metadata_herb$Continent

datascores_H$Fire_Int_Groups <- factor(
  datascores_H$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

ysf_paths <- datascores_H %>%
  filter(Continent == "Eurasia") %>%
  group_by(Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Fire_Int_Groups, YSF_interval)


#plot by continent, YSF, and Fire intensity
Herb_plot <- ggplot(subset(datascores_H, Continent == "Eurasia"), aes(x = NMDS1, y = NMDS2, 
                                      color = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Fire_Int_Groups)) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("Eurasia herbs") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(color = "Fire intensity") +
  geom_path(
    data = ysf_paths,
    aes(x = NMDS1, y = NMDS2, group = Fire_Int_Groups, color = Fire_Int_Groups),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    linewidth = 1,
    show.legend = FALSE
  )

Herb_plot

ysf_path_NA <- datascores_H %>%
  filter(Continent == "North_America") %>%
  group_by(Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Fire_Int_Groups, YSF_interval)


#plot by continent, YSF, and Fire intensity
Herb_plot_NA <- ggplot(subset(datascores_H, Continent == "North_America"), aes(x = NMDS1, y = NMDS2, 
                                                                      color = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Fire_Int_Groups)) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("North America herbs")+
  theme(legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction='vertical',
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  annotate("text", x = max(datascores_H$NMDS1), 
           y = min(datascores_H$NMDS2), 
           label = paste("Stress =", round(NMDS_herb$stress, 3)), 
           hjust = 1.1, vjust = 0.1, size = 5) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(color = "Fire intensity") +
  geom_path(
    data = ysf_path_NA,
    aes(x = NMDS1, y = NMDS2, group = Fire_Int_Groups, color = Fire_Int_Groups),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    linewidth = 1,
    show.legend = FALSE
  )

Herb_plot_NA



##Permanova with StudyID as a random factor
permutations <- with(metadata_herb, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_herb <- vegdist(Herb_clean, method = "bray")

###Fit permanova model
Permanova_herb <- adonis2(dist_herb ~ Continent*Years_since_fire*Fire_Int_Groups +
                            Avg_Temp*Fire_Int_Groups +
                            AvgPer*Fire_Int_Groups +
                            Latitude*Fire_Int_Groups +
                            Avg_Temp*Years_since_fire +
                            AvgPer* Years_since_fire +
                            Latitude * Years_since_fire +
                            Avg_Temp * Continent +
                            AvgPer * Continent +
                            Latitude * Continent,
                          data=metadata_herb,
                          permutations=permutations, method="bray")

Permanova_herb #IF THREE WAY INTERACTION IS SIGNIFICANT, IS THERE
# A POINT IN USING THEM IN OTHER INTERACTION TERMS?

#Homogeneity OK
Assumption <- anova(betadisper(dist_herb, metadata_herb$Fire_Int_Groups))
Assumption2 <- anova(betadisper(dist_herb, metadata_herb$Continent))

########################
#MOSSES
Moss_matrix <- Mosses_long %>%
  pivot_wider(
    id_cols = c(RowID, StudyID),
    names_from = species,      # base = original Dominant_X_Y column
    values_from = cover,
    values_fill = 0
  )

#Create NMDS
Moss_cols <- setdiff(names(Moss_matrix), c("RowID", "StudyID"))

Moss_filter <- Moss_matrix[
  rowSums(Moss_matrix[, Moss_cols] != 0, na.rm = TRUE) > 0,
]

metadata_moss <- Moss_filter %>%
  left_join(metadata_clean, "RowID")

Moss_info <- Moss_filter %>%
  select(c(StudyID, RowID))

moss_clean <- Moss_filter %>%
  select(- c(StudyID, RowID))

NMDS_moss <- metaMDS(moss_clean, distance = "bray", k = 2, trymax = 1000)

NMDS_moss <- metaMDS(moss_clean, distance = "bray", k = 2, trymax = 1000,
                     previous.best = NMDS_moss)

##Permanova with StudyID as a random factor
permutations <- with(metadata_moss, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_moss <- vegdist(moss_clean, method = "bray")

###Fit permanova model
Permanova_moss <- adonis2(dist_moss ~ Continent*Years_since_fire*Fire_Int_Groups +
                            Avg_Temp*Fire_Int_Groups +
                            AvgPer*Fire_Int_Groups +
                            Latitude*Fire_Int_Groups +
                            Avg_Temp*Years_since_fire +
                            AvgPer* Years_since_fire +
                            Latitude * Years_since_fire +
                          Avg_Temp * Continent +
                            AvgPer * Continent +
                            Latitude * Continent,
                          data=metadata_moss,
                          permutations=permutations, method="bray")

Permanova_moss

#Homogeneity OK
Assumption <- anova(betadisper(dist_moss, metadata_moss$Fire_Int_Groups))
Assumption2 <- anova(betadisper(dist_moss, metadata_moss$Continent))

#extract the site scores
datascores_M = as.data.frame(scores(NMDS_moss)$sites)  

#add metadata
datascores_M$Fire_Int_Groups = metadata_moss$Fire_Int_Groups
datascores_M$YSF_interval = metadata_moss$YSF_interval
datascores_M$SWI = metadata_moss$SWI
datascores_M$Avg_Temp = metadata_moss$Avg_Temp
datascores_M$SWI = metadata_moss$SWI
datascores_M$AvgPer = metadata_moss$AvgPer
datascores_M$Continent = metadata_moss$Continent

datascores_M$Fire_Int_Groups <- factor(
  datascores_M$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

#Arrows for continuous variables
env_vars <- metadata_moss %>%
  select(Avg_Temp, AvgPer, Latitude)

envfit_res <- envfit(NMDS_moss ~ ., data = env_vars, permutations = 999)

# To add arrows to ggplot, first extract coordinates
arrows_moss <- as.data.frame(scores(envfit_res, display = "vectors"))
arrows_moss$Variable <- rownames(arrows_moss)

ysf_paths <- datascores_M %>%
  filter(Continent == "Eurasia") %>%
  group_by(Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Fire_Int_Groups, YSF_interval)


#plot by continent, YSF, and Fire intensity
Moss_plot <- ggplot(subset(datascores_M, Continent == "Eurasia"), aes(x = NMDS1, y = NMDS2, 
                                                                      color = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Fire_Int_Groups)) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("Eurasia mosses") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  annotate("text", x = max(datascores_M$NMDS1), 
           y = min(datascores_M$NMDS2), 
           label = paste("Stress =", round(NMDS_moss$stress, 3)), 
           hjust = 1.1, vjust = 0.1, size = 5) +
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(color = "Fire intensity") +
  geom_path(
    data = ysf_paths,
    aes(x = NMDS1, y = NMDS2, group = Fire_Int_Groups, color = Fire_Int_Groups),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    linewidth = 1,
    show.legend = FALSE
  )

Moss_plot

ysf_path_NA <- datascores_M %>%
  filter(Continent == "North_America") %>%
  group_by(Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Fire_Int_Groups, YSF_interval)


#plot by continent, YSF, and Fire intensity
moss_plot_NA <- ggplot(subset(datascores_M, Continent == "North_America"), aes(x = NMDS1, y = NMDS2, 
                                                                               color = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Fire_Int_Groups)) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("North America mosses")+
  theme(legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction='vertical',
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  scale_color_manual(values = c("firebrick", "goldenrod", "cornflowerblue")) + 
  labs(color = "Fire intensity") +
  geom_path(
    data = ysf_path_NA,
    aes(x = NMDS1, y = NMDS2, group = Fire_Int_Groups, color = Fire_Int_Groups),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    linewidth = 1,
    show.legend = FALSE
  )

moss_plot_NA

ggplot(datascores_M, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = Fire_Int_Groups)) +
  geom_segment(
    data = arrows_moss,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = arrows_moss,
    aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = Variable),
    inherit.aes = FALSE
  )

##NMDS plots combined

TreePlots <- Tree_plot+Tree_plot_NA
HerbPlots <- Herb_plot+Herb_plot_NA
MossPlots <- Moss_plot+moss_plot_NA

TreePlots
HerbPlots
MossPlots

ggsave(plot = TreePlots, filename = "TreeNMDS.png", dpi = 300,
       height = 5.26, width = 13)
ggsave(plot = HerbPlots, filename = "HerbsNMDS.png", dpi = 300,
       height = 5.26, width = 13)
ggsave(plot = MossPlots, filename = "MossesNMDS.png", dpi = 300,
       height = 5.26, width = 13)
