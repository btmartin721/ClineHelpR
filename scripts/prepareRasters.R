
# Example data can be found in a Dryad Digital Repository (doi: XXXX)

library("raster")
library("rgdal")

#### Get mean among 12 months of raster layers. ####
####
####

# Mean annual solar radiation; worldclim.org
setwd("bioclim_R/layers/wc2.1_30s_srad/")

# Load into raster stack
solar.files <- list.files(pattern = "*.tif", full.names = TRUE)
solar.rad <- raster::stack(solar.files)

setwd("bioclim_R/layers/wc2.1_30s_wind/") # found at worldclim.org; wind speed layer

wind.files <- list.files(pattern = "*.tif", full.names = TRUE)
wind <- raster::stack(wind.files)

# Calculate annual means
solar.mean <- raster::calc(solar.rad, fun = mean, na.rm = TRUE)
wind.mean <- raster::calc(wind, fun = mean, na.rm = TRUE)

raster::plot(solar.mean)
raster::plot(wind.mean)

raster::writeRaster(x = solar.mean,
                    filename = "maxent_rasters/solarRAD_annualMean.tif",
                    driver = "GeoTiff")

raster::writeRaster(x = wind.mean,
                    filename = "maxent_rasters/windSpeed_annualMean.tif",
                    driver = "GeoTiff")

# List all rasters in maxent_rasters to load all of them.
files <-
  list.files(path = "maxent_rasters/",
             pattern = "wc*",
             full.names = TRUE)

# Create a raster stack from the list of files.
wc <- raster::stack(files)

samples <-
  read.csv(
    "sample_localities_maxent_southeast.csv",
    header = TRUE,
    stringsAsFactors = FALSE
  )

# Get CRS (coordinate reference system)
mycrs <- raster::crs(wc[[1]])

# Make sure values are numeric.
samples[,3] <- as.numeric(as.character(samples[,3]))
samples[,4] <- as.numeric(as.character(samples[,4]))

# Get dataframe of longitude, latitude.
coords <- data.frame(samples[,4], samples[,3])

# Create spatialpoints object.
p <- sp::SpatialPoints(coords = coords, proj4string=mycrs)



# Get the bounding box of the points
bb <- raster::bbox(p)

# Add bb.buffer * resolution (in arc-seconds) to bb for cropping raster layer.
envs.buf <- raster::extent(bb[1]-1,
                           bb[3]+1,
                           bb[2]-1,
                           bb[4]+1)


# Crop raster extent to sample bounding box + bb.buffer
envs.cropped <- raster::crop(wc, envs.buf)

raster::writeRaster(
  envs.cropped,
  filename = file.path(
    "maxent_rasters",
    paste0("crop_", names(envs.cropped))
  ),
  bylayer = TRUE,
  format = "GTiff"
)

