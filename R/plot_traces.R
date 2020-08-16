#' Plot LnL Traces to Assess Convergence for BGC Log-likelihoods
#'
#' Wrapper function to make BGC parameter trace plots
#' (LnL, alpha, beta, qa, qb, hi). It assumes that the BGC runs have been
#' aggregated. See ?combine_bgc_output. The plot should be inspected for
#' convergence before proceeding further in the pipeline.
#'
#' @param df.list List of data.frames containing aggregated BGC results.
#'                See ?combine_bgc_output
#' @param prefix Prefix used for input files
#' @param plotDIR Directory to save the plots
#' @export
#' @examples
#' plot_traces(df.list = aggregated.results, prefix = "eatt",
#'          plotDIR = "./bgc_plots")
plot_traces <- function(df.list, prefix, plotDIR = "./plots"){

  # Create plotDIR if doesn't already exist.
  dir.create(plotDIR, showWarnings = FALSE)

  plotLnL(df.list[[1]], plotDIR, prefix)
  plotBGCparams(df = df.list[[2]],
                ylab = "Alpha",
                bgc_param = "alpha",
                prefix = prefix,
                plotDIR = plotDIR)
  plotBGCparams(df = df.list[[3]],
                ylab = "Beta",
                bgc_param = "beta",
                prefix = prefix,
                plotDIR = plotDIR)
  plotBGCparams(df = df.list[[4]],
                ylab = "Gamma-Quantile",
                bgc_param = "qa",
                prefix = prefix,
                plotDIR = plotDIR)
  plotBGCparams(df = df.list[[5]],
                ylab = "Zeta-Quantile",
                bgc_param = "qb",
                prefix = prefix,
                plotDIR = plotDIR)
  plotBGCparams(df = df.list[[6]],
                ylab = "Hybrid Index",
                bgc_param = "hi",
                prefix = prefix,
                plotDIR = plotDIR)

  writeLines(paste("Saved trace plots to", plotDIR))

}

#' Function to plot BGC parameter traces
#' @param df Data.frame for BGC log-likelihoods
#' @param plotDIR Directory path to save plots in
#' @param prefix Character string; Prefix for output files
#' @noRd
plotLnL <- function(df, plotDIR, prefix){

  # Load in and transpose LnL data.
  lnl.df <- as.data.frame(t(df))

  # Set column name.
  colnames(lnl.df) <- c("lnl")

  # Round to nearest integer.
  lnl.df$lnl <- round(lnl.df$lnl)

  # Make scatterplot of thinned lnl samples.
  lnl <- data.frame(samp=1:length(lnl.df$lnl),
                    lnl=lnl.df$lnl)
  p1 <- ggplot2::ggplot(lnl,
                        ggplot2::aes(samp, lnl)) +
    ggplot2::geom_point(pch = 21) +
    ggplot2::theme_classic() +
    ggplot2::xlab('MCMC Sample') +
    ggplot2::ylab('Log Likelihood')

  ggplot2::ggsave(
    filename = file.path(plotDIR,
                         paste0(prefix,
                                '_bgc_LnL_convergence.pdf')),
    plot = p1,
    width = 6,
    height = 4,
    useDingbats = FALSE
  )
}

#' Function to plot BGC parameter traces
#' @param df Data.frame for one of the BGC parameters
#' @param ylab Character string indicating Y-axis label for trace plot
#' @param bgc_param Character string or BGC parameter to plot:
#'                  LnL, alpha, beta, qa, qb, or hi
#' @param prefix Character string; prefix for output files
#' @param plotDIR Directory path to save plots in
#' @noRd
plotBGCparams <- function(df, ylab, bgc_param, prefix, plotDIR){

  # Get mean of parameter for each MCMC sample
  df.mean <- lapply(df, mean)
  param.df <- data.frame(samp = 1:ncol(df),
                         param = unlist(df.mean),
                         stringsAsFactors = FALSE)

  p1 <- ggplot2::ggplot(param.df,
                        ggplot2::aes(samp, param)) +
    ggplot2::geom_line() +
    ggplot2::theme_classic() +
    ggplot2::xlab('MCMC Sample') +
    ggplot2::ylab(ylab)

  ggplot2::ggsave(
    filename = file.path(plotDIR,
                         paste0(prefix,
                                "_",
                                bgc_param,
                                "_convergence",
                                ".pdf")),
    plot = p1,
    width = 6,
    height = 4,
    useDingbats = FALSE
  )
}
