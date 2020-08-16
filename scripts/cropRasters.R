
# Example data can be found in a Dryad Digital Repository (doi: XXXX)

# This script will prepare specific rasters downloaded from https://worldclim.org
# and the National Land Cover Database, https://www.mrlc.gov/
# Essentially, it calculates the annual mean among 12 months of raster data for
# Wind speed and solar radiation from worldclim.org,
# and crops them all to the same extent.
# The extent it uses is defined by the sampling extent plus a buffer
# The cropped rasters are then saved for input into run_enmevalPipeline.R
# See ClinePlotR/scripts/run_enmevalPipeline.R
# If you aren't using the specific layers listed above, this script doesn't
# need to be run.

##########*******Not all rasters were included here because of their large size (>150GB).
##########*******Only fully cropped rasters were included in DRYAD exampleData

library("raster")
library("rgdal")

###############################################################################
#### Get mean among 12 months of raster layers. ####
###############################################################################

# Set working directory here
# Solar radiation for each of 12 months; https://worldclim.org
setwd("exampleData/ENMeval_bioclim/rasterLayers/original/wc2.1_30s_srad/")

# Directory to save output raster files.
rasterDIR <- "./output_rasters"

# Load into raster stack
solar.files <- list.files(pattern = "*.tif", full.names = TRUE)
solar.rad <- raster::stack(solar.files)

# Wind speed for each of 12 months; found at https://worldclim.org;
setwd("../wc2.1_30s_wind/")

wind.files <- list.files(pattern = "*.tif", full.names = TRUE)
wind <- raster::stack(wind.files)

# Calculate annual means
solar.mean <- raster::calc(solar.rad, fun = mean, na.rm = TRUE)
wind.mean <- raster::calc(wind, fun = mean, na.rm = TRUE)

# Plot the rasters
raster::plot(solar.mean)
raster::plot(wind.mean)

# Write the rasters to file in rasterDIR
raster::writeRaster(x = solar.mean,
                    filename = file.path(rasterDIR, "solarRAD_annualMean.tif"),
                    driver = "GeoTiff")

raster::writeRaster(x = wind.mean,
                    filename = file.path(rasterDIR, "windSpeed_annualMean.tif"),
                    driver = "GeoTiff")

# National Land Cover Database.
# This one was obtained from https://www.mrlc.gov/
# Not included in DRYAD, so if you want to use this script you'll have to
# download this file.
#nlcd <-
#  raster::raster(
#    "nlcd/a_nlcd2011.tif")
#  )

# List all rasters in rasterDIR to load all of them.
# Lists all files that start with wc*
# If your filenames are different, change the pattern argument
files <-
  list.files(path = rasterDIR,
             pattern = "wc*",
             full.names = TRUE)

# Create a raster stack from the list of files.
wc <- raster::stack(files)

setwd("../../../../../")

samples <-
  read.csv(
    file.path("exampleData",
              "ENMeval_bioclim",
              "localityInfo",
              "sample_localities_maxent_southeast.csv"),
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



# Get the bounding box of the sampling points
bb <- raster::bbox(p)

# Add bb.buffer * resolution (in arc-seconds) to bb for cropping raster layer.
envs.buf <- raster::extent(bb[1]-1,
                           bb[3]+1,
                           bb[2]-1,
                           bb[4]+1)


# Crop raster extent to sample bounding box + bb.buffer
envs.cropped <- raster::crop(wc, envs.buf)

# Write the cropped rasters to new files.
raster::writeRaster(
  envs.cropped,
  filename = file.path(rasterDIR,
                       paste0("crop_", names(envs.cropped))),
  bylayer = TRUE,
  format = "GTiff"
)

## I prepared NLCD raster in QGIS because it was taking a long time in R.
## I ran the raster warp function to resample at same extent and resolution as
## the worldclim layers. The pixel size was:
## 0.008333333333333333218,-0.008333333333333373116

# The next step is to run the enmeval pipeline
# Can use the run_enmevalPipeline.R script in ClinePlotR/scripts directory
# I made the categorical raster load first by making the filename
# have an earlier letter.
# That way I can easily specify which one is categorical when running ENMeval.

