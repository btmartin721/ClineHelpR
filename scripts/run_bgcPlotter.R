#!/usr/bin/env Rscript

library("ClinePlotR")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  writeLines("Usage: ./run_bgcPlotter.R [file_prefix] [pop name for plots]")
  stop("Error: Exactly two argument must be supplied.n",
       call.=FALSE)
}

prefix <- args[1]
admixPop <- args[2]
plotDIR <- "../../plots_maf"
genesDIR <- "../../bgc_results_genes"
fullDIR <-
  "D:/Dissertation/BOX/bgc_annotations/full_dataset/bgc_results_full/"

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
    popmap = file.path("../../popmaps/", paste0(prefix, ".bgc.popmap_final.txt")),
    loci.file = file.path("../../data_vcf_maf/", paste0(prefix, "_bgc_loci.txt"))
  )

rm(bgc.genes)
gc()

bgc.full <-
  combine_bgc_output(results.dir = fullDIR,
                     prefix = prefix, discard = 2500)

plot_traces(df.list = bgc.full,
         prefix = paste0(prefix, "_full"),
         plotDIR = plotDIR)

full.outliers <-
  get_bgc_outliers(
    df.list = bgc.full,
    admix.pop = admixPop,
    popmap = file.path("../../popmaps/", paste0(prefix, ".bgc.popmap_final.txt")),
    loci.file = file.path(fullDIR,
                          paste0(prefix, "_bgc_loci.txt"))
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

rm(bgc.full)
gc()

phiPlot(outlier.list = gene.outliers,
        popname = paste0(admixPop, " Genes"),
        line.size = 0.25,
        saveToFile = paste0(prefix, "_genes"),
        plotDIR = plotDIR,
        both.outlier.tests = FALSE,
        hist.y.origin = 1.2,
        hist.height = 1.8,
        margins = c(160.0, 5.5, 5.5, 5.5),
        hist.binwidth = 0.05)

gc()

phiPlot(outlier.list = full.outliers,
        popname = paste0(admixPop, " All"),
        line.size = 0.25,
        saveToFile = paste0(prefix, "_genome"),
        plotDIR = plotDIR,
        both.outlier.tests = FALSE,
        hist.y.origin = 1.2,
        hist.height = 1.8,
        margins = c(160.0, 5.5, 5.5, 5.5),
        hist.binwidth = 0.05)

gc()

gff <- parseGFF(gff.filepath = "../../gff/genes.gff")


genes.annotated <-
  join_bgc_gff(prefix = prefix,
               outlier.list = gene.outliers,
               gff.data = gff,
               scafInfoDIR = "../../scaffold_info")

plot_outlier_ideogram(
  prefix = prefix,
  outliers.genes = genes.annotated,
  outliers.full.scaffolds = full.outliers,
  pafInfo = "../../pafscaff/refmap_asm20_scafTrans_pafscaff.scaffolds.tdt",
  plotDIR = plotDIR,
  missing.chrs = c("chr11", "chr21", "chr25"),
  miss.chr.length = c(4997863, 1374423, 1060959),
  gene.size = 1e6,
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
