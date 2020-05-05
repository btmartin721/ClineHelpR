#!/usr/bin/env Rscript

rasterDIR <- "../../../../Dissertation/BOX/gis/bioclim_R"
dataDIR <- "../../../../Dissertation/BOX/introgress/"

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

# EA X GU
eagu <- runIntrogress(
  p1.file = file.path(dataDIR, "EAGU_p1data.txt"),
  p2.file = file.path(dataDIR, "EAGU_p2data.txt"),
  admix.file = file.path(dataDIR, "EAGU_admix_nocoordscut.txt"),
  loci.file = file.path(dataDIR, "EAGU_loci.txt"),
  clineLabels = c("EA", "Het", "GU"),
  minDelt = 0.7,
  prefix = "EAGU",
  outputDIR = file.path(dataDIR, "outputFiles"),
  sep = "\t",
  fixed = FALSE,
  pop.id = FALSE,
  ind.id = FALSE
)

# Subset individuals for only the populations I want
rasterPoint.list.subset <-
  lapply(rasterPoint.list,
         subsetIndividuals,
         file.path(dataDIR, "eagu_inds.txt"))

# Correlate genomic clines/hybrid index with environment/lat/lon
clinesXenvironment(
  clineList = eagu,
  rasterPointValues = rasterPoint.list.subset,
  clineLabels = c("EA", "Het", "GU"),
  outputDIR = file.path(dataDIR, "outputFiles", "EAGU"),
  clineMethod = "permutation",
  prefix = "EAGU",
)

dir.create(file.path(dataDIR, "rawRoutput"), showWarnings = FALSE)
saveRDS(eagu, file = file.path(dataDIR, "rawRoutput", "eagu_introgress.rds"))

#########################################################
## EA X TT ##
#########################################################
eatt <- runIntrogress(
  p1.file = file.path(dataDIR, "EATT_p1data_nocoordscut.txt"),
  p2.file = file.path(dataDIR, "EATT_p2data.txt"),
  admix.file = file.path(dataDIR, "EATT_admix_nocoordscut.txt"),
  loci.file = file.path(dataDIR, "EATT_loci.txt"),
  clineLabels = c("EA", "Het", "TT"),
  minDelt = 0.8,
  prefix = "EATT",
  outputDIR = file.path(dataDIR, "outputFiles"),
  sep = "\t",
  fixed = FALSE,
  pop.id = FALSE,
  ind.id = FALSE
)

saveRDS(eatt, file = file.path(dataDIR, "rawRoutput", "eatt_introgress.rds"))

# Subset individuals for only the populations I want
rasterPoint.list.subset <-
  lapply(rasterPoint.list,
         subsetIndividuals,
         file.path(dataDIR, "eatt_inds.txt"))

# Correlate genomic clines/hybrid index with environment/lat/lon
clinesXenvironment(
  clineList = eatt,
  rasterPointValues = rasterPoint.list.subset,
  clineLabels = c("EA", "Het", "TT"),
  outputDIR = file.path(dataDIR, "outputFiles", "EATT"),
  clineMethod = "permutation",
  prefix = "EATT"
)

#########################################################
## GU X TT ##
#########################################################
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
  cor.method = "auto"
)


