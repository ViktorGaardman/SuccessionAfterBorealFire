library(tidyverse)
library(maps)
library(dbscan)

set.seed(123)  # reproducibility

dataset <- read.csv ("Coordinate_Data_20251217.csv", sep = ";")

worldmap <- map_data("world")

#Make a circle of points for overlapping or
#nearly overlapping study points

coords <- as.matrix(dataset[, c("Longitude", "Latitude")])

# tolerance in degrees (≈ 10 km at mid-latitudes)
eps <- 0.5

clusters <- dbscan(coords, eps = eps, minPts = 1)

dataset_clustered <- dataset %>%
  mutate(cluster = clusters$cluster)

sunflower <- dataset_clustered %>%
  group_by(cluster) %>%
  mutate(
    n = n(),
    id = row_number(),
    angle = ifelse(n > 1, 2 * pi * (id - 1) / n, 0),
    radius = ifelse(n > 1, 0.7 + 0.02 * sqrt(n), 0),
    Long_sf = Longitude + radius * cos(angle),
    Lat_sf  = Latitude  + radius * sin(angle)
  ) %>%
  ungroup()

map <- ggplot() +
  geom_map(
    data = worldmap, map = worldmap,
    aes(long, lat, map_id = region),
    color = "black", fill = "lightgray", linewidth = 0.1
  ) + coord_fixed(xlim = c(-170, 110), ylim = c(30, 80), expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
  )

Studypoint_map <- map + geom_point(data=sunflower, 
                 aes(x=Long_sf, y=Lat_sf),color = "firebrick",
                 fill="gold", pch=21, size = 0.7, stroke = 0.2)

ggsave(filename = "Studypoint_map.TIFF", plot = Studypoint_map,
       height = 2.63, width = 6.5, dpi = 300)
