#' Plot LnL Traces to Assess Convergence for BGC Log-likelihoods
#'
#' This function makes a log-likelihood trace plot. The thin parameter should
#' equal the thinning parameter used for BGC. This function assumes that the
#' BGC runs have been aggregated. See ?combine_bgc_output for help. The plot
#' should be inspected for convergence before proceeding.
#' Requires ggplot2 package.
#'
#' @param df.list List of data.frames containing aggregated BGC results.
#'                See ?combine_bgc_output
#' @param prefix Prefix used for input files
#' @param thin Thinning parameter used in BGC
#' @param plotDIR Directory to save the plots
#' @export
#' @examples
#' plot_lnl(df.list = aggregated.results, prefix = "population1",
#'          thin = 40, plotDIR = "./bgc_plots")
plot_lnl <- function(df.list, prefix, thin, plotDIR = "./plots"){

  # Create plotDIR if doesn't already exist.
  dir.create(plotDIR, showWarnings = FALSE)

  # Load in and transpose LnL data.
  lnl.df <- as.data.frame(t(df.list[[1]]))

  # Set column name.
  colnames(lnl.df) <- c("lnl")

  # Round to nearest integer.
  lnl.df$lnl <- round(lnl.df$lnl)

  # Make scatterplot of thinned lnl samples.
  lnl <- data.frame(samp=1:length(lnl.df$lnl),
                    lnl=lnl.df$lnl)
  p1 <- ggplot2::ggplot(lnl,
                        ggplot2::aes(samp,lnl)) +
    ggplot2::geom_point(pch=21) +
    ggplot2::theme_classic() +
    ggplot2::xlab('MCMC Sample') +
    ggplot2::ylab('Log Likelihood')
  ggplot2::ggsave(filename = file.path(plotDIR,
                                       paste0(prefix,
                                              '_bgc_LnL_convergence.pdf')),
                  plot=p1,
                  width=6,
                  height=4,
                  useDingbats=FALSE)

  writeLines(paste0("Saved LnL trace plots to ", plotDIR))

}
