
# Example data can be found in a Dryad Digital Repository (doi: XXXX)

# This script will run Introgress, then use the output from
# run_enmevalPipeline.R (see ClinePlotR/scripts directory)
# to correlate genomic clines and hybrid indices with
# latitude, longitude, and the environmental variables.

library("ClinePlotR")

# Obtained from https://worldclim.org; finest resolution available.
rasterDIR <- "exampleData/ENMeval_bioclim/"

# Introrgress input data directory
dataDIR <- "exampleData/introgress"

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
      paste0("rasterPoints_",
             colnames(rasterPoint.list[[i]][4]),
             ".csv")
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = ","
  )
}

# Run introgress.
# There are many parameters that can be used to fine-tune the analysis.
# E.g. minDelt determines the allele frequency threshold (delta)
# required for a locus to be tested in Introgress.
# can be lowered to e.g., 0.7 or 0.6, which will include
# more SNPs.
eatt <- runIntrogress(
  p1.file = file.path(dataDIR, "inputFiles", "EATT_p1data.txt"),
  p2.file = file.path(dataDIR, "inputFiles", "EATT_p2data.txt"),
  admix.file = file.path(dataDIR, "inputFiles", "EATT_admix.txt"),
  loci.file = file.path(dataDIR, "inputFiles", "EATT_loci.txt"),
  clineLabels = c("EA", "Het", "TT"),
  minDelt = 0.8,
  prefix = "EATT",
  outputDIR = file.path(dataDIR, "outputFiles", "introgress_plots"),
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


