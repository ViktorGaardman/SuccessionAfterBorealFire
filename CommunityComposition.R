#Step 1. Load packages
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(ggeffects)
library(performance)
library(grid)
library(ggcorrplot)

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
  
Ground_long <-species_long %>%
  filter(base %in% c("Dominant_dwarfshrub_1", "Dominant_dwarfshrub_2",
                     "Dominant_dwarfshrub_3", "Dominant_herb_1",
                     "Dominant_herb_2", "Dominant_herb_3",
                     "Dominant_graminoid_1", "Dominant_graminoid_2",
                     "Dominant_graminoid_3"))

Mosses_long <- species_long %>%
  filter(base %in% c("Dominant_bryohpyte_1", "Dominant_bryophyte_2",
                     "Dominant_bryophyte_3"))

#This version makes a normal matrix. Still the only option for NMDS?
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

#plot by time period and date
Tree_plot <- ggplot(datascores_T, aes(x = NMDS1, y = NMDS2, 
               color = Continent, shape = Fire_Int_Groups)) +
  geom_point(size = 2, aes(color = Continent)) +
  coord_fixed() + 
  theme_bw() +
  stat_ellipse(aes(color = Continent), level=0.9) +
  theme(legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction='vertical',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  annotate("text", x = max(datascores_T$NMDS1), 
           y = min(datascores_T$NMDS2), 
           label = paste("Stress =", round(NMDS_tree$stress, 3)), 
           hjust = 0.4, vjust = 1.5, size = 6) +
  scale_color_brewer(palette = "Dark2") + 
  labs(color = "Continent")

Tree_plot

##Permanova with StudyID as a random factor
permutations <- with(metadata_tree, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_tree <- vegdist(Tree_clean, method = "bray")

###Fit permanova model
Permanova_tree <- adonis2(dist_tree ~ Continent*Years_since_fire*Fire_Int_Groups +
                            Avg_Temp +
                            AvgPer,
                          data=metadata_tree,
                         permutations=permutations, method="bray")

Permanova_tree

#Check assumption of homogeneity of multivariate dispersion
Assumption <- anova(betadisper(dist_tree, metadata_tree$Years_since_fire))


###HERB LAYER

#This version makes a normal matrix. Still the only option for NMDS?
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
                     previous.best = NMDS_tree)

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
  theme(legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction='vertical',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  annotate("text", x = max(datascores_H$NMDS1), 
           y = min(datascores_H$NMDS2), 
           label = paste("Stress =", round(NMDS_herb$stress, 3)), 
           hjust = 1.1, vjust = 0.1, size = 5) +
  scale_color_brewer(palette = "Dark2") + 
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
  theme(legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction='vertical',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + 
  annotate("text", x = max(datascores_H$NMDS1), 
           y = min(datascores_H$NMDS2), 
           label = paste("Stress =", round(NMDS_herb$stress, 3)), 
           hjust = 1.1, vjust = 0.1, size = 5) +
  scale_color_brewer(palette = "Dark2") + 
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
                            Avg_Temp +
                            AvgPer +
                            Latitude,
                          data=metadata_herb,
                          permutations=permutations, method="bray")

Permanova_herb

Assumption <- anova(betadisper(dist_herb, metadata_herb$Fire_Int_Groups))
Assumption2 <- anova(betadisper(dist_herb, metadata_herb$Continent))
