
# Example data can be found in a Dryad Digital Repository (doi: XXXX)

# This script will run Introgress, then use the output from
# run_enmevalPipeline.R (see ClinePlotR/scripts directory)
# to correlate genomic clines and hybrid indices with
# latitude, longitude, and the environmental variables.

library("ClinePlotR")

# Obtained from https://worldclim.org; finest resolution available.
rasterDIR <- "../clinehelpr_analysis/data/exampleData/ENMeval_bioclim/"

# Introrgress input data directory
dataDIR <- "../clinehelpr_analysis/data/exampleData/introgress"

# Object saved in run_enmevalPipeline.R script
# See ClinePlotR/scripts directory
envList <- readRDS("exampleData/ENMeval_bioclim/Robjects/envList.rds")

# Extract raster values at each sample locality.
rasterPoint.list <- extractPointValues(envList)

# Make directory to store extracted raster values per sample locality
dir.create(file.path(dataDIR, "extractedRasterValues"), showWarnings = FALSE)

# Write them to file for each raster layer
for (i in 1:length(rasterPoint.list)) {
  write.table(
    rasterPoint.list[[i]],
    file.path(
      dataDIR,
      "extractedRasterValues",
      paste0("rasterPoints_tutorial",
             colnames(rasterPoint.list[[i]][4]),
             ".csv")
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = ","
  )
}

plotDIR <- "../clinehelpr_analysis/notebooks/plots"

genind <- adegenet::read.structure(
  file = "../clinehelpr_analysis/data/exampleData/bgc/bgc_structureFiles/eatt.str",
  n.ind = 101,
  n.loc = 233,
  onerowperind = FALSE,
  col.lab = 1,
  col.pop = 2,
  row.marknames = 1,
  NA.char = -9,
  ask = FALSE)

genind2introgress(gen = genind,
                  p1 = "2",
                  p2 = "3",
                  admix = "1",
                  missingPerLoc = 0.5,
                  missingPerInd = 0.9,
                  prefix=file.path(plotDIR, "eatt"))

# Run introgress.
# There are many parameters that can be used to fine-tune the analysis.
# E.g. minDelt determines the allele frequency threshold (delta)
# required for a locus to be tested in Introgress.
# can be lowered to e.g., 0.7 or 0.6, which will include
# more SNPs.

eatt <- runIntrogress(
  p1.file = file.path(plotDIR, "eatt_p1data.txt"),
  p2.file = file.path(plotDIR, "eatt_p2data.txt"),
  admix.file = file.path(plotDIR, "eatt_admix.txt"),
  loci.file = file.path(plotDIR, "eatt_loci.txt"),
  clineLabels = c("EA", "Het", "TT"),
  minDelt = 0.8,
  prefix = "EATT",
  outputDIR = plotDIR,
  sep = "\t",
  fixed = FALSE,
  pop.id = FALSE,
  ind.id = FALSE
)

# Save Introgress output to R object.
saveRDS(eatt, file = file.path(dataDIR, "outputRobject", "eatt_introgress.rds"))

# Subset individuals for only the populations I want
# This was done because the locality data has more individuals
# than were used for introgress
rasterPoint.list.subset <-
  lapply(rasterPoint.list,
         subsetIndividuals,
         file.path(dataDIR, "inputFiles", "eatt_inds.txt"))

# Correlate genomic clines/hybrid index with environment/lat/lon
# Can use different correlation methods. E.g. pearson or kendall
clinesXenvironment(
  clineList = eatt,
  rasterPointValues = rasterPoint.list.subset,
  clineLabels = c("EA", "Het", "TT"),
  outputDIR = file.path(dataDIR, "outputFiles", "clines"),
  clineMethod = "permutation",
  prefix = "eatt",
  cor.method = "spearman"
)


