#Step 1. Load packages
library(tidyverse)

#Step 2. Load raw data
df <- read.csv ("Raw_Data.csv", sep = ";")


# Summarize number of studies per year per fire intensity
summary_df <- df %>%
  group_by(Years_since_fire, Fire_Int_Groups) %>%
  summarise(num_studies = n(),
            sqrm = sum(Sample_size * Plot_size),
            .groups = 'drop')

#Step 3. Plot number of studies per year
studies_plot <- ggplot(summary_df, aes(x = Years_since_fire, 
                                       y = num_studies,
                                       color = Fire_Int_Groups)) +
  geom_point(aes(fill = Fire_Int_Groups)) +
  geom_line(aes(color = Fire_Int_Groups)) +
  theme_bw()

studies_plot

#Step 4. Plot sqr meters per year

sqrm_plot <- ggplot(summary_df, aes(x = Years_since_fire, 
                                       y = log10(sqrm),
                                       color = Fire_Int_Groups)) +
  geom_point(aes(fill = Fire_Int_Groups)) +
  geom_line(aes(color = Fire_Int_Groups)) +
  theme_bw()

sqrm_plot
