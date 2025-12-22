#Step 1. Load packages
library(tidyverse)


#Step 2. Load raw data and divide into metadata and species matrix
df <- read.csv ("Raw_Data.csv", sep = ";")

metadata <- df %>%
  select(-contains("postfire"))

species_raw <- df %>%
  select(StudyID, contains("postfire")) %>%
  select(
    StudyID,
    matches("_postfire$|_postfire_cover$")
  )

species_raw <- species_raw %>%
  mutate(SampleID = row_number())

species_names <- species_raw %>%
  select(SampleID, StudyID, matches("_postfire$")) %>%
  pivot_longer(
    cols = matches("_postfire$"),
    names_to = "base",
    names_pattern = "(.*)_postfire$",
    values_to = "species"
  )

species_cover <- species_raw %>%
  select(SampleID, StudyID, matches("_postfire_cover$")) %>%
  pivot_longer(
    cols = matches("_postfire_cover$"),
    names_to = "base",
    names_pattern = "(.*)_postfire_cover$",
    values_to = "cover"
  )

species_long <- left_join(
  species_names,
  species_cover,
  by = c("SampleID", "StudyID", "base")
)

species_long <- species_long %>%
  mutate(cover = as.numeric(unlist(cover)))


##Here things go wrong. My matrix collapses multiple combinations of 
#species = 0 and cover = 0 into one column. But that loses information.

species_matrix_subplot <- species_long %>%
  pivot_wider(
    id_cols = c(SampleID, StudyID),
    names_from = base,      # base = original Dominant_X_Y column
    values_from = cover,
    values_fill = 0
  )

species_mat <- species_matrix %>%
  select(-StudyID) %>%
  as.matrix()

#Step 3. Load climate and trait data
TRY_Traits <- read.csv ("TRY_Traits_20251219.csv", sep = ";")

Perc_Data <- read.csv("site_climatological_percipitation_WorldClim_v2.csv",
                      sep = ";")
Temp_Data <- read.csv("site_climatological_temperature_WorldClim_v2.csv",
                      sep = ";")