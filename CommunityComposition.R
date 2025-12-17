#Step 1. Load packages
library(tidyverse)


#Step 2. Load matrix
matrix <- read.csv ("Species_Matrix.csv", sep = ";")

#Step 3. Load climate and trait data
TRY_Traits <- read.csv ("TRY_Data.csv", sep = ";")

Climate_Data <- read.csv ("Climate_Data.csv", sep = ";")