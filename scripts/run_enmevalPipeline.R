#!/usr/bin/env Rscript

# I first prepared the rasters in prepareRasters.R
# Now I run the ENMeval pipeline

plotDIR <- "../../../../Dissertation/BOX/gis/bioclim_R/plots"

envList <-
  prepare_rasters(
    raster.dir = "../../../../Dissertation/BOX/gis/bioclim_R/layers/maxent_rasters/crop",
    sample.file = "../../../../Dissertation/BOX/gis/bioclim_R/sample_localities_maxent_southeast_noNA.csv",
    header = TRUE,
    bb.buffer = 0.5,
   plotDIR = plotDIR
   )

bg <- partition_raster_bg(env.list = envList, plotDIR = "../../../../Dissertation/BOX/gis/bioclim_R/plots/")

# Saved to file and reloaded because rJava ran out of memory
saveRDS(bg, file = "../../../../Dissertation/BOX/gis/bioclim_R/gis_Routput/bg.rds")
saveRDS(envList, file = "../../../../Dissertation/BOX/gis/bioclim_R/gis_Routput/envList.rds")

envList <-
  readRDS(
    file = "../../../../Dissertation/BOX/gis/bioclim_R/gis_Routput/envList.rds"
    )

envs.fg <- envList[[1]]
coords <- envList[[3]]
rm(envList)
gc()

bg <- readRDS("../../../../Dissertation/BOX/gis/bioclim_R/gis_Routput/bg.rds")

# Give rJava more memory.
options(java.parameters = "-Xmx24g")

eval.par <- runENMeval(envs.fg = envs.fg,
                       bg = bg, parallel = FALSE,
                       categoricals = 1,
                       partition.method = "checkerboard1",
                       coords = coords)

saveRDS(eval.par, "../../../../Dissertation/BOX/gis/bioclim_R/gis_Routput/enm_eval.rds")

summarize_ENMeval(
  eval.par = eval.par,
  plotDIR = plotDIR,
  minLat = 25,
  maxLat = 45,
  minLon = -100,
  maxLon = -45,
  imp.margins = c(15.1, 4.1, 3.1, 2.1))

# Use delta.AICc to find best model and plot response curves for it
# For this dataset it was features classes LQH and RM = 2.5
pdf(file = file.path("../../../../Dissertation/BOX/gis/bioclim_R/plots/",
                     "ENMeval_responseCurves.pdf"),
    width=7,
    height=7,
    onefile=TRUE)
# Plot response curves.
dismo::response(eval.par@models[[28]])
dev.off()

# I want to save the best model so I can open it in QGIS
# Which was LQH feature classes with RM = 2.5
raster::writeRaster(eval.par@predictions@layers[[28]], filename = "../../../../Dissertation/BOX/gis/bioclim_R/plots/LQH_2.5.tif", format = "GTiff", device = "GTiff")
