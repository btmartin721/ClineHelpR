#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  writeLines("Usage: ./run_bgcPlotter.R [file_prefix] [pop name for plots]")
  stop("Error: Exactly two argument must be supplied.n",
       call.=FALSE)
}

prefix <- args[1]
admixPop <- args[2]

bgc.genes <-
  combine_bgc_output(results.dir = "../../bgc_results_genes/",
                                     prefix = prefix,
                                     discard = 3000)

plot_traces(df.list = bgc.genes,
         prefix = prefix,
         plotDIR = "../../plots")

gene.outliers <-
  get_bgc_outliers(
    df.list = bgc.genes,
    admix.pop = admixPop,
    popmap = file.path("../../popmaps/", paste0(prefix, ".bgc.popmap_final.txt")),
    loci.file = file.path("../../loci_files/", paste0(prefix, "_bgc_loci.txt"))
  )

rm(bgc.genes)
gc()

bgc.full <-
  combine_bgc_output(results.dir = "../../bgc_results_full/",
                     prefix = prefix, discard = 3000)

plot_traces(df.list = bgc.full,
         prefix = paste0(prefix, "_full"),
         plotDIR = "../../plots")

full.outliers <-
  get_bgc_outliers(
    df.list = bgc.full,
    admix.pop = admixPop,
    popmap = file.path("../../popmaps/", paste0(prefix, ".bgc.popmap_final.txt")),
    loci.file = file.path("../../bgc_results_full/",
                          paste0(prefix, "_bgc_loci.txt"))
  )

rm(bgc.full)
gc()

phiPlot(outlier.list = gene.outliers,
        popname = paste0(admixPop, " Genes"),
        line.size = 0.25,
        saveToFile = paste0(prefix, "_genes"),
        plotDIR = "../../plots/",
        both.outlier.tests = TRUE,
        hist.y.origin = 1.2,
        hist.height = 1.8,
        margins = c(160.0, 5.5, 5.5, 5.5),
        hist.binwidth = 0.05)

gc()

phiPlot(outlier.list = full.outliers,
        popname = paste0(admixPop, " All"),
        line.size = 0.25,
        saveToFile = paste0(prefix, "_genome"),
        plotDIR = "../../plots/",
        both.outlier.tests = TRUE,
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
  plotDIR = "../../plots",
  missing.chrs = c("chr11", "chr21", "chr25"),
  miss.chr.length = c(4997863, 1374423, 1060959)
)
