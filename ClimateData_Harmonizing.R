library(terra)

# ---- A. stack monthly rasters ------------------------------------------------
tif_dir <- "C:/Users/vikto/OneDrive/Skrivbord/Vetenskapliga studier/Fire Ecology/SuccessionAfterBorealFire/WorldClim2_Data/AverageTemp"            # adjust to your folder

tifs  <- list.files(pattern = "\\.tif$", 
                    path = tif_dir,
                    recursive = TRUE, full.names = TRUE)

tifs  <- sort(tifs)                          # guarantees 1970‑01 … 2000‑12 order

rast_stack <- rast(tifs)                     # 372‑layer SpatRaster

# annotate temporal dimension (optional but handy)
time(rast_stack) <- seq(as.Date("1970-01-15"),
                        by   = "month",
                        length.out = nlyr(rast_stack))

nlyr(rast_stack)      # should be 372
unique(format(time(rast_stack), "%Y"))


# ---- B. read your sampling sites --------------------------------------------
sites <- read.csv("Coordinate_Data_20251217.csv", sep=";")               # site_id, lat, lon (WGS84)

names(sites)

pts   <- vect(sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# ---- C. extract cell values --------------------------------------------------
# choose method = "bilinear" if you want a smooth estimate instead of the
# center‑of‑cell value
vals <- extract(rast_stack, pts, method = "bilinear")
# 'vals' columns: ID (row number) + 372 temperature values

plot(rast_stack[[1]], main = "Raster with correctly placed points")
points(pts, col = "red", pch = 20)


# keep Title instead of numeric ID
sites <- sites[!duplicated(sites$Title), ]
vals$ID <- sites$Title
rownames(vals) <- vals$ID
vals$ID <- NULL

MAT <- rowMeans(vals, na.rm = TRUE)
SD_monthly <- apply(vals, 1, sd, na.rm = TRUE)
CV_monthly <- SD_monthly / MAT

out <- data.frame(
  site_id = sites$Title,
  MAT     = MAT,
  SD      = SD_monthly,
  CV      = CV_monthly
)

write.csv(out, "site_climatological_temperature_WorldClim_v2.csv",
          row.names = FALSE)


####Percipitation#####
# stack precip
tif_dir_p <- "C:/Users/vikto/OneDrive/Skrivbord/Vetenskapliga studier/Fire Ecology/SuccessionAfterBorealFire/WorldClim2_Data/Percipitation"            # adjust to your folder


tif_p <- list.files(pattern="\\.tif$", 
                    path = tif_dir_p,
                    full.names=TRUE)

rast_prec <- rast(tif_p)

# extract
pvals <- extract(rast_prec, pts, method = "bilinear")
pvals$ID <- sites$Title

# monthly to annual climatology
pvals_num <- pvals[, sapply(pvals, is.numeric)]
AvgPer <- rowSums(pvals_num[ , -1], na.rm = TRUE)
SD_monthly <- apply(pvals_num, 1, sd, na.rm = TRUE)
CV_monthly <- SD_monthly / AvgPer

out_per <- data.frame(
  site_id = sites$Title,
  AvgPer     = AvgPer,
  SD      = SD_monthly,
  CV      = CV_monthly
)

write.csv(out_per, "site_climatological_percipitation_WorldClim_v2.csv",
          row.names = FALSE)
