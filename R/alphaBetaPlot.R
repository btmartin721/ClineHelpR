#' Makes Alpha ~ Beta Plots for BGC results
#'
#' This function generates an Alpha x Beta 2D density plot
#' with polygon hulls encapsulating alpha and beta outliers.
#' The plot can be printed to the plot
#' window in Rstudio (default) or saved to file by specifying a file prefix
#' with the saveToFile option.
#'
#' Generates Alpha x Beta Density Plots, Highlighting Outliers.
#' @param outlier.list List of data.frames returned from get_bgc_outliers().
#'                     see ?get_bgc_outliers for more info
#' @param overlap.zero Boolean. If TRUE, outliers include when the credible
#'                     doesn't contain 0
#' @param qn.interval Boolean. If TRUE, outliers = outside quantile interval
#' @param both.outlier.tests If TRUE, outliers = both overlap.zero and
#'                           qn.interval
#' @param line.size Size of regression lines in plot. Default = 0.25
#' @param saveToFile If specified, saves plots to file. The value
#'                   for saveToFile should be the prefix for the filename
#'                   you want to save to
#' @param plotDIR Directory path. If saveToFile is specified, plots are
#'                saved in this directory
#' @param device File format to save plot. Supports ggplot2::ggsave devices
#' @param neutral.color Color for non-outlier loci. Default = "gray"
#' @param alpha.color Color for alpha outlier loci. Default = "blue"
#' @param beta.color Color for beta outlier loci. Default = "red"
#' @param both.color Color for loci that are alpha and beta outliers. Default = "purple"
#' @param margins Vector of margins for phi plot: c(Top, Right, Bottom, Left)
#'                Top is extended to include the Hybrid Index histogram
#' @param margin.units Units for margins parameter. Default = "points"
#' @param height Height for plot. Default = 7
#' @param width Width for plot. Default = 7
#' @param dim.units Units for height and width. Default = "in"
#' @param text.size Size for plot text. Default = 18
#' @param dpi DPI for saving plot. Default = 300
#' @param padding Padding to add to X and Y axes. Default=0.1
#' @param alpha_limits Limits (e.g., c(1,1) for plotting alpha. Default = NULL
#' @param beta_limits Limits (e.g., c(1,1) for plotting beta. Default = NULL
#' @export
#' @examples
#' alphaBetaPlot(outlier.list = outliers, saveToFile = "pop1_phiPlot",
#' plotDIR = "./bgc_plots", device = "png")
#'
#' phiPlot(outlier.list = gene.outliers, ,
#' line.size = 0.1, qn.interval = FALSE, overlap.zero = FALSE,
#' both.outlier.tests = TRUE, saveToFile = "mus_phiPlot")
#'
#' phiPlot(outlier.list = full.outliers,  overlap.zero = FALSE)
alphaBetaPlot <- function(outlier.list,
                    overlap.zero = TRUE,
                    qn.interval = TRUE,
                    both.outlier.tests = FALSE,
                    saveToFile = NULL,
                    plotDIR = "./plots",
                    device = "pdf",
                    neutral.color = "gray",
                    alpha.color = "blue",
                    beta.color = "red",
                    both.color = "purple",
                    margins = c(150.0, 5.5, 5.5, 5.5),
                    margin.units = "points",
                    padding=0.1,
                    height = 7,
                    width = 7,
                    dim.units = "in",
                    text.size = 18,
                    dpi = 300,
                    alpha_limits=NULL,
                    beta_limits=NULL){

  # To use the dplyr pipe.
  `%>%` <- dplyr::`%>%`

  snps <- outlier.list[[1]] # SNPs data.frame
  hi <- outlier.list[[3]] # hybrid index data.frame

  rm(outlier.list)
  gc()

  # If only want crazy.a and crazy.b outliers
  if (both.outlier.tests){
    overlap.zero = FALSE
    qn.interval = FALSE
    writeLines("\n\nboth.outlier.tests is set.")
    writeLines("Ignoring overlap.zero qn.interval settings\n")
    snps$alpha.signif <-
      snps$crazy.a == TRUE

    snps$beta.signif <-
      snps$crazy.b == TRUE
  }

  # TRUE if alpha/beta excess/outliers.
  if (overlap.zero & qn.interval){
    snps$alpha.signif <-
      (!is.na(snps$alpha.excess) | !is.na(snps$alpha.outlier))
    snps$beta.signif <-
      (!is.na(snps$beta.excess) | !is.na(snps$beta.outlier))

  } else if (overlap.zero & !qn.interval){
    snps$alpha.signif <- (!is.na(snps$alpha.excess))
    snps$beta.signif <- (!is.na(snps$beta.excess))

  } else if (!overlap.zero & qn.interval){
    snps$alpha.signif <- (!is.na(snps$alpha.outlier))
    snps$beta.signif <- (!is.na(snps$beta.outlier))

  } else if (!overlap.zero & !qn.interval & !both.outlier.tests){
    mywarning <-
      paste(
        "\n\noverlap.zero, qn.interval, and both.outlier.tests",
        "were all FALSE. Setting both.outlier.tests to TRUE\n\n"
      )
    warning(paste(strwrap(mywarning), collapse = "\n"))

    snps$alpha.signif <-
      snps$crazy.a == TRUE

    snps$beta.signif <-
      snps$crazy.b == TRUE
  }

  isAlphaOutliers <- TRUE
  isBetaOutliers <- TRUE

  # If there aren't any alpha or beta outliers.
  if (all(snps$alpha.signif == FALSE)){
    isAlphaOutliers <- FALSE
    warning("\n\nWarning: No alpha outliers were identified.\n")
    if (both.outlier.tests){
      writeLines("\n\nTry setting both.outlier.tests to FALSE")
    }
  }

  if (all(snps$beta.signif == FALSE)){
    isBetaOutliers <- FALSE
    warning("\n\nWarning: No beta outliers were identified.\n\n")
    if (both.outlier.tests){
      writeLines("\n\nTry setting both.outlier.tests to FALSE")
    }
  }

  snps[snps$alpha.signif == TRUE & snps$alpha > 0.0, "alpha.signif"] <- "pos"
  snps[snps$alpha.signif == TRUE & snps$alpha < 0.0, "alpha.signif"] <- "neg"
  snps[snps$beta.signif == TRUE & snps$beta > 0.0, "beta.signif"] <- "pos"
  snps[snps$beta.signif == TRUE & snps$beta < 0.0, "beta.signif"] <- "neg"

  snps["fill_color"] <- neutral.color
  snps[snps$alpha.signif != FALSE, "fill_color"] <- alpha.color
  snps[snps$beta.signif != FALSE, "fill_color"] <- beta.color
  snps[snps$alpha.signif != FALSE & snps$beta.signif != FALSE, "fill_color"] <- both.color

  alpha_hulls <- snps[snps$alpha.signif != FALSE,]
  beta_hulls <- snps[snps$beta.signif != FALSE,]


  ab.plot <- ggplot2::ggplot(snps) +
    ggplot2::stat_density_2d(ggplot2::aes(x=alpha, y=beta, fill=..level..), geom="polygon") +
    ggplot2::scale_fill_distiller(palette= "Spectral", direction=1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.line = ggplot2::element_line(colour = "black"),
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = text.size,
                                          colour = "black"),
        text = ggplot2::element_text(size = text.size,
                                     colour = "black"),
        axis.ticks = ggplot2::element_line(size = ggplot2::rel(1.5)),
        plot.margin = ggplot2::unit(margins,
                                    margin.units)
    ) +
    ggplot2::geom_point(ggplot2::aes(x=alpha, y=beta, col=fill_color), size=2) +
    ggplot2::scale_colour_identity()


  if (nrow(beta_hulls) > 0){
    ab.plot <- ab.plot +
      ggforce::geom_mark_hull(data=beta_hulls, mapping=ggplot2::aes(x=alpha, y=beta,
                                                   group=beta.signif,
                                                   color=beta.color),
                     concavity=5, expand=0, radius=0, alpha=0.7)
  }

  if (nrow(alpha_hulls) > 0){
    ab.plot <- ab.plot +
      ggforce::geom_mark_hull(data=alpha_hulls, mapping=ggplot2::aes(x=alpha, y=beta,
                                                  group=alpha.signif,
                                                  color=alpha.color),
                     concavity=5, expand=0, radius=0, alpha=0.7)
  }

  if (length(alpha_limits) == 2){
    ab.plot <- ab.plot + ggplot2::xlim(alpha_limits)
  }else{
    ab.plot <- ab.plot + ggplot2::xlim(min(snps$alpha)-padding, max(snps$alpha)+padding)
  }

  if (length(beta_limits) == 2){
    ab.plot <- ab.plot + ggplot2::ylim(beta_limits)
  }else{
    ab.plot <- ab.plot + ggplot2::ylim(min(snps$beta)-padding, max(snps$beta)+padding)
  }

  if (!is.null(saveToFile)){
    # If saveToFile is specified, save plot as PDF.
    # Create plotDIR if doesn't already exist.
    dir.create(plotDIR, showWarnings = FALSE)

    ggplot2::ggsave(
      filename = paste0(saveToFile, "_alphaXbeta.pdf"),
      plot = ab.plot,
      device = device,
      plotDIR,
      width = width,
      height = height,
      dpi = dpi,
      units = dim.units
    )

  } else{
    # If saveToFile is not specified by user.
    print(ab.plot)
  }

}
