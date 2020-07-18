#!/usr/bin/env Rscript

library("ClinePlotR")

plotBGCresults <-
  function(prefix,
           admixPop,
           plotDIR,
           genesDIR,
           fullDIR,
           pafInfo,
           gene.size,
           other.size) {
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
        popmap = file.path(
          "../bgcPlotter/popmaps/",
          paste0(prefix, ".bgc.popmap_final.txt")
        ),
        loci.file = file.path(
          "../bgcPlotter/data_vcf_maf_genes/",
          paste0(prefix, "_bgc_loci.txt")
        )
      )

    rm(bgc.genes)
    gc()

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
        popmap = file.path(
          "../bgcPlotter/popmaps/",
          paste0(prefix, ".bgc.popmap_final.txt")
        ),
        loci.file = file.path(
          "../bgcPlotter/data_vcf_maf_full/",
          paste0(prefix, "_bgc_loci.txt")
        )
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

    gc()

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

    gc()

    gff <- parseGFF(gff.filepath = "../bgcPlotter/gff/genes.gff")


    genes.annotated <-
      join_bgc_gff(
        prefix = prefix,
        outlier.list = gene.outliers,
        gff.data = gff,
        scafInfoDIR = "../../scaffold_info_trachemys"
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

    return(ref.trans.info)
  }

plotDIR <- "../plots_maf_trachemys"
genesDIR <- "../bgcPlotter/bgc_results_genes"
fullDIR <-
  "../bgcPlotter/bgc_results_full/"

pafInfo  <-
  "../bgcPlotter/pafscaff_trachemys/pafscaff_asm20_scafTrans_tscripta.scaffolds.tdt"

eagu.prefix <- "eagu"
eagu.admixPop <- "EAGU"

eatt.prefix <- "eatt"
eatt.admixPop <- "EATT"

gutt.prefix <- "gutt"
gutt.admixPop <- "GUTT"


eagu.ref.trans.info <-
  plotBGCresults(eagu.prefix,
                 eagu.admixPop,
                 plotDIR,
                 genesDIR,
                 fullDIR,
                 pafInfo,
                 4e6,
                 1e6)

eatt.ref.trans.info <-
  plotBGCresults(eatt.prefix,
                 eatt.admixPop,
                 plotDIR,
                 genesDIR,
                 fullDIR,
                 pafInfo,
                 4e6,
                 1e6)

gutt.ref.trans.info <-
  plotBGCresults(gutt.prefix,
                 gutt.admixPop,
                 plotDIR,
                 genesDIR,
                 fullDIR,
                 pafInfo,
                 4e6,
                 1e6)

eatt.bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = eatt.prefix)

eatt.gene.outliers <-
  get_bgc_outliers(
    df.list = eatt.bgc.genes,
    admix.pop = eatt.admixPop,
    popmap = file.path(
      "../bgcPlotter/popmaps/",
      paste0(eatt.prefix, ".bgc.popmap_final.txt")
    ),
    loci.file = file.path(
      "../bgcPlotter/data_vcf_maf_genes/",
      paste0(eatt.prefix, "_bgc_loci.txt")
    )
  )

eagu.bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = eagu.prefix)

eagu.gene.outliers <-
  get_bgc_outliers(
    df.list = eagu.bgc.genes,
    admix.pop = eagu.admixPop,
    popmap = file.path(
      "../bgcPlotter/popmaps/",
      paste0(eagu.prefix, ".bgc.popmap_final.txt")
    ),
    loci.file = file.path(
      "../bgcPlotter/data_vcf_maf_genes/",
      paste0(eagu.prefix, "_bgc_loci.txt")
    )
  )

gutt.bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = gutt.prefix)

gutt.gene.outliers <-
  get_bgc_outliers(
    df.list = gutt.bgc.genes,
    admix.pop = gutt.admixPop,
    popmap = file.path(
      "../bgcPlotter/popmaps/",
      paste0(gutt.prefix, ".bgc.popmap_final.txt")
    ),
    loci.file = file.path(
      "../bgcPlotter/data_vcf_maf_genes/",
      paste0(gutt.prefix, "_bgc_loci.txt")
    )
  )


alphaBetaPlot(
  eagu.gene.outliers,
  alpha.color = "cornflowerblue",
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = eagu.prefix,
  plotDIR = "../plots_maf_trachemys/", padding = 0.2
)

alphaBetaPlot(
  eatt.gene.outliers,
  alpha.color = "cornflowerblue",
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = eatt.prefix,
  plotDIR = "../plots_maf_trachemys/", padding = 0.2,
)

alphaBetaPlot(
  gutt.gene.outliers,
  alpha.color = "cornflowerblue",
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = gutt.prefix,
  plotDIR = "../plots_maf_trachemys/", padding = 0.2,
)
