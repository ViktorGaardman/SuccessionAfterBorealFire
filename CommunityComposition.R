#Step 1. Load packages
library(tidyverse)
library(vegan)
library(ggcorrplot)
library(patchwork)
library(writexl)

rm(list=ls())

#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Clean_Data.csv", sep = ";")

metadata <- df %>%
  select(-contains("postfire"))

metadata <- metadata[,1:20]

species_raw <- df %>%
  select(StudyID, RowID, Continent, contains("postfire")) %>%
  select(
    StudyID, RowID, Continent,
    matches("_postfire$|_postfire_cover$")
  )


species_names <- species_raw %>%
  select(RowID, StudyID, Continent, matches("_postfire$")) %>%
  pivot_longer(
    cols = matches("_postfire$"),
    names_to = "base",
    names_pattern = "(.*)_postfire$",
    values_to = "species"
  )

species_cover <- species_raw %>%
  select(RowID, StudyID, Continent, matches("_postfire_cover$")) %>%
  pivot_longer(
    cols = matches("_postfire_cover$"),
    names_to = "base",
    names_pattern = "(.*)_postfire_cover$",
    values_to = "cover"
  )

species_long <- left_join(
  species_names,
  species_cover,
  by = c("Continent", "RowID", "StudyID", "base")
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

#Standardize temp and percipitation
metadata_clean$Temp_sc <- as.numeric(scale(metadata_clean$Avg_Temp, center = TRUE, scale = TRUE))
metadata_clean$Per_sc <- as.numeric(scale(metadata_clean$AvgPer, center = TRUE, scale = TRUE))

#Check collinearity of continuous variables
pond_corr <- metadata_clean %>% 
  select(SWI, Temp_sc, Per_sc, Latitude) %>%
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

##Permanova with StudyID as a random factor
permutations <- with(metadata_tree, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_tree <- vegdist(Tree_clean, method = "bray")

###Fit permanova model
Permanova_tree <- adonis2(dist_tree ~
                            Fire_Int_Groups*Years_since_fire +
                            Years_since_fire * Continent +
                            Temp_sc +
                            Per_sc,
                          data=metadata_tree,
                          permutations=permutations, method="bray")

Permanova_tree

anova(betadisper(dist_tree, metadata_tree$Fire_Int_Groups))
anova(betadisper(dist_tree, metadata_tree$Continent))

Permanova_tree_out <- as.data.frame(Permanova_tree)
Permanova_tree_out$Term <- rownames(Permanova_tree_out)
rownames(Permanova_tree_out) <- NULL

write_xlsx(
  Permanova_tree_out,
  path = "PERMANOVA_results_trees&shrubs.xlsx"
)

#NMDS
NMDS_tree <- metaMDS(Tree_clean, distance = "bray", k = 2, trymax = 5000)

NMDS_tree <- metaMDS(Tree_clean, distance = "bray", k = 2, trymax = 5000,
                    previous.best = NMDS_tree)

#extract the site scores
datascores_T = as.data.frame(scores(NMDS_tree)$sites)  

#add metadata
datascores_T$Fire_Int_Groups = metadata_tree$Fire_Int_Groups
datascores_T$YSF_interval = metadata_tree$YSF_interval
datascores_T$Temp_sc = metadata_tree$Temp_sc
datascores_T$Temp_sc = metadata_tree$Temp_sc
datascores_T$Continent = metadata_tree$Continent

datascores_T$Fire_Int_Groups <- factor(
  datascores_T$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

#Add arrows across time
ysf_paths_T <- datascores_T %>%
  group_by(Continent, Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Continent, Fire_Int_Groups, YSF_interval)

#Add a shape to arrow starts
ysf_starts_T <- ysf_paths_T %>%
  group_by(Continent, Fire_Int_Groups) %>%
  slice_min(YSF_interval, n = 1)

#Add shape to arrowhead
ysf_ends_T <- ysf_paths_T %>%
  group_by(Continent, Fire_Int_Groups) %>%
  slice_max(YSF_interval, n = 1)

Tree_plot <- ggplot(
  datascores_T,
  aes(x = NMDS1, y = NMDS2, color = Fire_Int_Groups)
) +
  geom_point(size = 2) +
  geom_path(
    data = ysf_paths_T,
    aes(group = interaction(Continent, Fire_Int_Groups),
        color = Fire_Int_Groups),
    linewidth = 1
  )+
  geom_point(
    data = ysf_ends_T,
    aes(x = NMDS1, y = NMDS2, fill = Fire_Int_Groups),
    shape = 24,
    color = "black",
    size = 3,
    stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(
    data = ysf_starts_T,
    aes(x = NMDS1, y = NMDS2, fill = Fire_Int_Groups),
    shape = 23,
    color = "black",
    size = 2.5,
    stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  facet_wrap(~ Continent, 
             scales = "free_x",
             labeller = labeller(
               Continent = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  theme_bw() +
  scale_color_manual(
    values = c("firebrick", "goldenrod", "cornflowerblue")
  ) +
  scale_fill_manual(
    values = c("firebrick", "goldenrod", "cornflowerblue")
  ) +
  labs(color = "Fire intensity") +
  ggtitle("Trees & shrubs")+
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = max(datascores_T$NMDS1),
    y = min(datascores_T$NMDS2),
    label = paste("Stress =", round(NMDS_tree$stress, 3)),
    hjust = 1.1, vjust = 0.1, size = 5
  )

Tree_plot

ggsave(plot = Tree_plot, filename = "NMDS_trees.png",
       dpi = 300, width = 13, height = 5.26)









###HERB LAYER

Herb_matrix <- Ground_long %>%
  pivot_wider(
    id_cols = c(RowID, StudyID, Continent),
    names_from = species,      # base = original Dominant_X_Y column
    values_from = cover,
    values_fill = 0
  )

#Create NMDS
Herb_cols <- setdiff(names(Herb_matrix), c("RowID", "StudyID", "Continent"))

Herb_filter <- Herb_matrix[
  rowSums(Herb_matrix[, Herb_cols] != 0, na.rm = TRUE) > 0,
]

metadata_herb <- Herb_filter %>%
  left_join(metadata_clean, "RowID")

Herb_info <- Herb_filter %>%
  select(c(StudyID, RowID, Continent))

Herb_clean <- Herb_filter %>%
  select(- c(StudyID, RowID, Continent))

##Permanova with StudyID as a random factor
permutations <- with(metadata_herb, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_herb <- vegdist(Herb_clean, method = "bray")

###Fit permanova model
Permanova_herb <- adonis2(dist_herb ~
                            Fire_Int_Groups*Years_since_fire +
                            Years_since_fire * Continent.x +
                            Temp_sc +
                            Per_sc,
                          data=metadata_herb,
                          permutations=permutations, method="bray")

Permanova_herb

anova(betadisper(dist_herb, metadata_herb$Fire_Int_Groups))
anova(betadisper(dist_herb, metadata_herb$Continent.x))

Permanova_herb_out <- as.data.frame(Permanova_herb)
Permanova_herb_out$Term <- rownames(Permanova_herb_out)
rownames(Permanova_herb_out) <- NULL

write_xlsx(
  Permanova_herb_out,
  path = "PERMANOVA_results_groundlayer.xlsx"
)

NMDS_herb <- metaMDS(Herb_clean, distance = "bray", k = 2, trymax = 1000)

NMDS_herb <- metaMDS(Herb_clean, distance = "bray", k = 2, trymax = 1000,
                     previous.best = NMDS_herb)

#extract the site scores
datascores_H = as.data.frame(scores(NMDS_herb)$sites)  

#add metadata
datascores_H$Fire_Int_Groups = metadata_herb$Fire_Int_Groups
datascores_H$YSF_interval = metadata_herb$YSF_interval
datascores_H$Continent.x = metadata_herb$Continent.x
datascores_H$Temp_sc = metadata_herb$Temp_sc
datascores_H$Per_sc = metadata_herb$Per_sc


datascores_H$Fire_Int_Groups <- factor(
  datascores_H$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

#Add arrows across time
ysf_paths_H <- datascores_H %>%
  group_by(Continent.x, Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Continent.x, Fire_Int_Groups, YSF_interval)

#Add a shape to arrow starts
ysf_starts_H <- ysf_paths_H %>%
  group_by(Continent.x, Fire_Int_Groups) %>%
  slice_min(YSF_interval, n = 1)

#Add shape to arrowhead
ysf_ends_H <- ysf_paths_H %>%
  group_by(Continent.x, Fire_Int_Groups) %>%
  slice_max(YSF_interval, n = 1)

Herb_plot <- ggplot(
  datascores_H,
  aes(x = NMDS1, y = NMDS2, color = Fire_Int_Groups)
) +
  geom_point(size = 2) +
  geom_path(
    data = ysf_paths_H,
    aes(group = interaction(Continent.x, Fire_Int_Groups),
        color = Fire_Int_Groups),
    linewidth = 1
  )+
  geom_point(
    data = ysf_ends_H,
    aes(x = NMDS1, y = NMDS2, fill = Fire_Int_Groups),
    shape = 24,
    color = "black",
    size = 3,
    stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(
    data = ysf_starts_H,
    aes(x = NMDS1, y = NMDS2, fill = Fire_Int_Groups),
    shape = 23,
    color = "black",
    size = 2.5,
    stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  facet_wrap(~ Continent.x, 
             scales = "free_x",
             labeller = labeller(
               Continent.x = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  theme_bw() +
  scale_color_manual(
    values = c("firebrick", "goldenrod", "cornflowerblue")
  ) +
  scale_fill_manual(
    values = c("firebrick", "goldenrod", "cornflowerblue")
  ) +
  labs(color = "Fire intensity") +
  ggtitle("Ground layer")+
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.position = "none"
  ) +
  annotate(
    "text",
    x = max(datascores_H$NMDS1),
    y = min(datascores_H$NMDS2),
    label = paste("Stress =", round(NMDS_herb$stress, 3)),
    hjust = 1.1, vjust = 0.1, size = 5
  )

Herb_plot

ggsave(plot = Herb_plot, filename = "NMDS_herbs.png",
       dpi = 300, width = 13, height = 5.26)

#combined plot
combinedNMDS <- Tree_plot/Herb_plot
combinedNMDS

ggsave(plot = combinedNMDS, filename = "NMDS_plots.png",
       dpi = 300, height = 7.89, width = 8.66)

###MOSSES
##INSUFFICENT DATA TO MAKE ANY ANALYSES

Moss_matrix <- Mosses_long %>%
  pivot_wider(
    id_cols = c(RowID, StudyID, Continent),
    names_from = species,      # base = original Dominant_X_Y column
    values_from = cover,
    values_fill = 0
  )

#Create NMDS
Moss_cols <- setdiff(names(Moss_matrix), c("RowID", "StudyID", "Continent"))

Moss_filter <- Moss_matrix[
  rowSums(Moss_matrix[, Moss_cols] != 0, na.rm = TRUE) > 0,
]

metadata_moss <- Moss_filter %>%
  left_join(metadata_clean, "RowID")

Moss_info <- Moss_filter %>%
  select(c(StudyID, RowID, Continent))

Moss_clean <- Moss_filter %>%
  select(- c(StudyID, RowID, Continent))

#To aid model fit and account for species_sp, we group all mosses to genus level

Moss_genus <- Moss_clean %>%
  tibble::rownames_to_column("row") %>%
  pivot_longer(-row, names_to = "species", values_to = "value") %>%
  mutate(genus = sub("_.*", "", species)) %>%
  group_by(row, genus) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = value) %>%
  column_to_rownames("row")

##Permanova with StudyID as a random factor
permutations <- with(metadata_moss, how(nperm=999, blocks = StudyID.x))

#Calculate distance matrix
dist_moss <- vegdist(Moss_genus, method = "bray")

###Fit permanova model
Permanova_moss <- adonis2(dist_moss ~
                            Fire_Int_Groups*Years_since_fire +
                            Years_since_fire * Continent.x +
                            Temp_sc +
                            Per_sc,
                          data=metadata_moss,
                          permutations=permutations, method="bray")

Permanova_moss

anova(betadisper(dist_moss, metadata_moss$Fire_Int_Groups))
anova(betadisper(dist_moss, metadata_moss$Continent.x))

Permanova_moss_out <- as.data.frame(Permanova_moss)
Permanova_moss_out$Term <- rownames(Permanova_moss_out)
rownames(Permanova_moss_out) <- NULL

write_xlsx(
  Permanova_moss_out,
  path = "PERMANOVA_results_bryophytes.xlsx"
)

NMDS_moss <- metaMDS(Moss_clean, distance = "bray", k = 2, trymax = 1000)

NMDS_moss <- metaMDS(Moss_clean, distance = "bray", k = 2, trymax = 1000,
                     previous.best = NMDS_herb)

NMDS_moss_genus <- metaMDS(Moss_clean, distance = "bray", k = 2, trymax = 1000)

NMDS_moss_genus <- metaMDS(Moss_clean, distance = "bray", k = 2, trymax = 5000,
                     previous.best = NMDS_herb)

#extract the site scores
datascores_H = as.data.frame(scores(NMDS_herb)$sites)  

#add metadata
datascores_H$Fire_Int_Groups = metadata_herb$Fire_Int_Groups
datascores_H$YSF_interval = metadata_herb$YSF_interval
datascores_H$Continent.x = metadata_herb$Continent.x
datascores_H$Temp_sc = metadata_herb$Temp_sc
datascores_H$Per_sc = metadata_herb$Per_sc


datascores_H$Fire_Int_Groups <- factor(
  datascores_H$Fire_Int_Groups,
  levels = c("High", "Medium", "Low")
)

#Add arrows across time
ysf_paths_H <- datascores_H %>%
  group_by(Continent.x, Fire_Int_Groups, YSF_interval) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Continent.x, Fire_Int_Groups, YSF_interval)

#Add a shape to arrow starts
ysf_starts_H <- ysf_paths_H %>%
  group_by(Continent.x, Fire_Int_Groups) %>%
  slice_min(YSF_interval, n = 1)

#Add shape to arrowhead
ysf_ends_H <- ysf_paths_H %>%
  group_by(Continent.x, Fire_Int_Groups) %>%
  slice_max(YSF_interval, n = 1)

Herb_plot <- ggplot(
  datascores_H,
  aes(x = NMDS1, y = NMDS2, color = Fire_Int_Groups)
) +
  geom_point(size = 2) +
  geom_path(
    data = ysf_paths_H,
    aes(group = interaction(Continent.x, Fire_Int_Groups),
        color = Fire_Int_Groups),
    linewidth = 1
  )+
  geom_point(
    data = ysf_ends_H,
    aes(x = NMDS1, y = NMDS2, fill = Fire_Int_Groups),
    shape = 24,
    color = "black",
    size = 3,
    stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(
    data = ysf_starts_H,
    aes(x = NMDS1, y = NMDS2, fill = Fire_Int_Groups),
    shape = 23,
    color = "black",
    size = 2.5,
    stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  facet_wrap(~ Continent.x, 
             scales = "free_x",
             labeller = labeller(
               Continent.x = c(
                 "Eurasia" = "Eurasia",
                 "North_America" = "North America"))) +
  theme_bw() +
  scale_color_manual(
    values = c("firebrick", "goldenrod", "cornflowerblue")
  ) +
  scale_fill_manual(
    values = c("firebrick", "goldenrod", "cornflowerblue")
  ) +
  labs(color = "Fire intensity") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = max(datascores_T$NMDS1),
    y = min(datascores_T$NMDS2),
    label = paste("Stress =", round(NMDS_tree$stress, 3)),
    hjust = 1.1, vjust = 0.1, size = 5
  )

Herb_plot

ggsave(plot = Herb_plot, filename = "NMDS_herbs.png",
       dpi = 300, width = 13, height = 5.26)



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


ggsave(plot = Herb_plot, filename = "EUHerb_NMDS.png", dpi = 300,
       height = 5.26, width = 13)
ggsave(plot = Herb_plot_NA, filename = "NAHerb_NMDS.png", dpi = 300,
       height = 5.26, width = 13)


#Dominating species per group after fire
#Make a plot showing the two/three (?)species with the highest mean cover
#Per plant group

#Use species_long => split by continent =>
#Add metadata => mean values per species and time_since_fire*fire intensity combination
#plot and show as a list in paper
#Run on genus level to make it more readable

species_long_genus <- species_long %>%
  mutate(genus = sub("_.*", "", species)) %>%
  group_by(RowID, genus)

alldata <- species_long_genus %>%
  left_join(metadata_clean, 'RowID')

summed <- alldata %>%
  group_by(Continent.x, Fire_Int_Groups, YSF_interval, genus) %>%
  summarise(
    meancov = mean(cover),
    .groups = "drop"
  ) 

summed_top5 <- summed %>%
  group_by(Continent.x, Fire_Int_Groups, YSF_interval) %>%
  slice_max(meancov, n = 5, with_ties = FALSE) %>%
  ungroup()

levels(summed_top5$YSF_interval) <- c(
  "1",
  "2",
  "3",
  "4",
  "6-7",
  "8-9",
  "10-22"
)


summed_top5$YSF_interval <- as.numeric(summed_top5$YSF_interval)

summed_top5_low <- summed_top5 %>%
  filter(Fire_Int_Groups == "Low")

dominanceplot_Low_EU <- ggplot(subset(summed_top5_low, Continent.x == "Eurasia"), aes(x = YSF_interval, y = meancov, 
                                         color = genus)) +
  geom_point(size = 3) +
  geom_line(aes(color = genus), size = 1) +
  theme_bw() +
  scale_color_discrete() +
  labs(color = "Genus") +
    ggtitle("Europe", subtitle = "Low intensity fire") +
  xlab("Years since fire") +
  ylab("Mean cover") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 

dominanceplot_Low_EU

dominanceplot_Low_NA <- ggplot(subset(summed_top5_low, Continent.x == "North_America"), aes(x = YSF_interval, y = meancov, 
                                                                                      color = genus)) +
  geom_point(size = 3) +
  geom_line(aes(color = genus), size = 1) +
  theme_bw() +
  scale_color_discrete() +
  labs(color = "Genus") +
  ggtitle("North America", subtitle = "Low intensity fire") +
  xlab("Years since fire") +
  ylab("Mean cover") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 

dominanceplot_Low_NA

summed_top5_medium <- summed_top5 %>%
  filter(Fire_Int_Groups == "Medium")

dominanceplot_medium_EU <- ggplot(subset(summed_top5_medium, Continent.x == "Eurasia"), aes(x = YSF_interval, y = meancov, 
                                                                                      color = genus)) +
  geom_point(size = 3) +
  geom_line(aes(color = genus), size = 1) +
  theme_bw() +
  scale_color_discrete() +
  labs(color = "Genus") +
  ggtitle("Europe", subtitle = "Medium intensity fire") +
  xlab("Years since fire") +
  ylab("Mean cover") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 

dominanceplot_medium_EU

dominanceplot_medium_NA <- ggplot(subset(summed_top5_medium, Continent.x == "North_America"), aes(x = YSF_interval, y = meancov, 
                                                                                            color = genus)) +
  geom_point(size = 3) +
  geom_line(aes(color = genus), size = 1) +
  theme_bw() +
  scale_color_discrete() +
  labs(color = "Genus") +
  ggtitle("North America", subtitle = "Medium intensity fire") +
  xlab("Years since fire") +
  ylab("Mean cover") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 

dominanceplot_medium_NA


summed_top5_high <- summed_top5 %>%
  filter(Fire_Int_Groups == "High")

dominanceplot_high_EU <- ggplot(subset(summed_top5_high, Continent.x == "Eurasia"), aes(x = YSF_interval, y = meancov, 
                                                                                            color = genus)) +
  geom_point(size = 3) +
  geom_line(aes(color = genus), size = 1) +
  theme_bw() +
  scale_color_discrete() +
  labs(color = "Genus") +
  ggtitle("Europe", subtitle = "High intensity fire") +
  xlab("Years since fire") +
  ylab("Mean cover") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 

dominanceplot_high_EU

dominanceplot_high_NA <- ggplot(subset(summed_top5_high, Continent.x == "North_America"), aes(x = YSF_interval, y = meancov, 
                                                                                                  color = genus)) +
  geom_point(size = 3) +
  geom_line(aes(color = genus), size = 1) +
  theme_bw() +
  scale_color_discrete() +
  labs(color = "Genus") +
  ggtitle("North America", subtitle = "High intensity fire") +
  xlab("Years since fire") +
  ylab("Mean cover") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 

dominanceplot_high_NA

CompleteDominancePlot <- (dominanceplot_Low_EU|dominanceplot_Low_NA)/
  (dominanceplot_medium_EU|dominanceplot_medium_NA)/
  (dominanceplot_high_EU|dominanceplot_high_NA)

CompleteDominancePlot

ggsave(plot = CompleteDominancePlot, filename = "dominancePlots.png",
       dpi = 300, width = 13, height = 15.78)
