library(terra)
getwd()
# ---- A. stack monthly rasters ------------------------------------------------
tif_dir <- "/WorldClim2_Data/AverageTemp"            # adjust to your folder

tifs  <- list.files(pattern = "\\.tif$", 
                    recursive = TRUE, full.names = TRUE)
tifs  <- sort(tifs)                          # guarantees 1970‑01 … 2000‑12 order

rast_stack <- rast(tifs)                     # 372‑layer SpatRaster

# annotate temporal dimension (optional but handy)
time(rast_stack) <- seq(as.Date("1970-01-15"),
                        by   = "month",
                        length.out = nlyr(rast_stack))

# ---- B. read your sampling sites --------------------------------------------
sites <- read.csv("Coordinate_Data.csv", sep=";")               # site_id, lat, lon (WGS84)

names(sites)

pts   <- vect(sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# ---- C. extract cell values --------------------------------------------------
# choose method = "bilinear" if you want a smooth estimate instead of the
# center‑of‑cell value
vals <- extract(rast_stack, pts, method = "bilinear")
# 'vals' columns: ID (row number) + 372 temperature values

plot(rast_stack[[1]], main = "Raster with correctly placed points")
points(pts, col = "red", pch = 20)


# keep site_id instead of numeric ID
vals$ID <- sites$site_id
rownames(vals) <- vals$ID
vals$ID <- NULL

# ---- D. aggregate monthly → annual -------------------------------------------
dates <- time(rast_stack)
yrs   <- format(dates, "%Y")                 # vector length 372, e.g. "1970"


# split the 372 columns by year, then rowMeans
yearly <- sapply(split.default(vals, yrs), 
                 FUN = rowMeans, na.rm = TRUE)

# 'yearly' dims: sites × 31 years (1970…2000)
yearly <- as.data.frame(yearly)
yearly$site_id <- rownames(yearly)
yearly <- yearly[, c("site_id", sort(names(yearly)[-ncol(yearly)]))]

write.csv(yearly, "site_mean_temps_1970_2000.csv", row.names = FALSE)



####Percipitation#####
# stack precip
tif_p <- list.files(pattern="\\.tif$", full.names=TRUE)
rast_prec <- rast(tif_p)

# extract
pvals <- extract(rast_prec, pts, method = "bilinear")
pvals$ID <- sites$Title

# monthly to annual climatology
pvals$annual_mm <- rowSums(pvals[ , -1], na.rm = TRUE)
write.csv(pvals, "site_precip_1970_2000.csv", row.names = FALSE)
