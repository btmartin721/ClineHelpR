
# Example data can be found in a Dryad Digital Repository (doi: XXXX)

library("ClinePlotR")

rasterDIR <- "bioclim_R" # Obtained from worldclim.org
dataDIR <- "introgress/"

envList <- readRDS(file.path(rasterDIR, "gis_Routput", "envList.rds"))

rasterPoint.list <- extractPointValues(envList)

dir.create(file.path(dataDIR, "rasterPoints"), showWarnings = FALSE)

for (i in 1:length(rasterPoint.list)) {
  write.table(
    rasterPoint.list[[i]],
    file.path(
      rasterDIR,
      "rasterPoints",
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

gutt <- runIntrogress(
  p1.file = file.path(dataDIR, "GUTT_p1data.txt"),
  p2.file = file.path(dataDIR, "GUTT_p2data.txt"),
  admix.file = file.path(dataDIR, "GUTT_admix.txt"),
  loci.file = file.path(dataDIR, "GUTT_loci.txt"),
  clineLabels = c("GU", "Het", "TT"),
  minDelt = 0.8,
  prefix = "GUTT",
  outputDIR = file.path(dataDIR, "outputFiles"),
  sep = "\t",
  fixed = FALSE,
  pop.id = FALSE,
  ind.id = FALSE
)

saveRDS(gutt, file = file.path(dataDIR, "rawRoutput", "gutt_introgress.rds"))

# Subset individuals for only the populations I want
rasterPoint.list.subset <-
  lapply(rasterPoint.list,
         subsetIndividuals,
         file.path(dataDIR, "gutt_inds.txt"))

# Correlate genomic clines/hybrid index with environment/lat/lon
clinesXenvironment(
  clineList = gutt,
  rasterPointValues = rasterPoint.list.subset,
  clineLabels = c("GU", "Het", "TT"),
  outputDIR = file.path(dataDIR, "outputFiles", "GUTT"),
  clineMethod = "permutation",
  prefix = "GUTT",
  cor.method = "spearman"
)


