
# Example data can be obtained from a Dryad Digital Repository (doi: XXXX)

library("ClinePlotR")

plotDIR <- "./plots"
genesDIR <- "genes"
fullDIR <- "fulldataset"

pafInfo  <- "refmap_asm20.scaffolds.tdt"

prefix <- "gutt"
admixPop <- "GUTT"
popmap <- "popmaps/gutt.bgc.popmap.txt"
genes.loci.file <- "genes/gutt_bgc_loci.txt"
full.loci.file <- "fulldataset/gutt_full_bgc_loci.txt"
scafInfoDIR <- "./scafInfo"
gene.size <- 4e6
other.size <- 1e6

# GFF file is the only one that can't be found in vignettes/exampledata.
# It can be found at: https://www.ncbi.nlm.nih.gov/genome/?term=Terrapene
gff <- "gff/genes.gff"



bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = prefix)

plot_traces(df.list = bgc.genes,
            prefix = prefix,
            plotDIR = plotDIR)

gene.outliers <-
  get_bgc_outliers(
    df.list = bgc.genes,
    admix.pop = admixPop,
    popmap = popmap,
    loci.file = genes.loci.file
  )

bgc.full <-
  combine_bgc_output(results.dir = fullDIR,
                     prefix = prefix)

plot_traces(
  df.list = bgc.full,
  prefix = paste0(prefix, "_full"),
  plotDIR = plotDIR
)

full.outliers <-
  get_bgc_outliers(
    df.list = bgc.full,
    admix.pop = admixPop,
    popmap = popmap,
    loci.file = full.loci.file
  )

ab.range <-
  data.frame(
    "full.alpha" = c(
      "min" = min(full.outliers[[1]]$alpha),
      "max" = max(full.outliers[[1]]$alpha)
    ),
    "full.beta" = c(
      "min" = min(full.outliers[[1]]$beta),
      "max" = max(full.outliers[[1]]$beta)
    ),
    "genes.alpha" = c(
      "min" = min(gene.outliers[[1]]$alpha),
      "max" = max(gene.outliers[[1]]$alpha)
    ),
    "genes.beta" = c(min(gene.outliers[[1]]$beta), max(gene.outliers[[1]]$beta))
  )

write.table(
  data.frame("Header" = rownames(ab.range), ab.range),
  file = file.path(plotDIR, paste0(prefix, "_ab.ranges.csv")),
  sep = ",",
  row.names = F,
  col.names = T,
  quote = F
)

phiPlot(
  outlier.list = gene.outliers,
  popname = paste0(admixPop, " Genes"),
  line.size = 0.35,
  saveToFile = paste0(prefix, "_genes"),
  plotDIR = plotDIR,
  both.outlier.tests = FALSE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)

phiPlot(
  outlier.list = full.outliers,
  popname = paste0(admixPop, " All"),
  line.size = 0.2,
  saveToFile = paste0(prefix, "_genome"),
  plotDIR = plotDIR,
  both.outlier.tests = FALSE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)

alphaBetaPlot(
  gene.outliers,
  alpha.color = "cornflowerblue",
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = prefix,
  plotDIR = "./plots", padding = 0.2,
)

gff <- parseGFF(gff.filepath = gff)

genes.annotated <-
  join_bgc_gff(
    prefix = prefix,
    outlier.list = gene.outliers,
    gff.data = gff,
    scafInfoDIR = scafInfoDIR
  )

ref.trans.info <-
  plot_outlier_ideogram(
    prefix = prefix,
    outliers.genes = genes.annotated,
    outliers.full.scaffolds = full.outliers,
    pafInfo = pafInfo,
    plotDIR = plotDIR,
    gene.size = gene.size,
    other.size = other.size,
    both.outlier.tests = FALSE
  )

write.table(
  genes.annotated,
  file = file.path(plotDIR, paste0(prefix, "_genesAnnotated.csv")),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = ","
)

write.table(
  ref.trans.info,
  file = file.path(plotDIR, paste0(prefix, "_refTransInfo.csv")),
  quote = TRUE,
  row.names = FALSE,
  col.names = TRUE,
  sep = ","
)



