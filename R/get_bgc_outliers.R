#' Identify Outlier BGC SNPs Using Two Criteria
#'
#' These functions identify alpha and beta BGC outliers using two criteria.
#' First, if a SNP's 95% credible interval for alpha or beta doesn't overlap
#' with 0, it is considered to have excess ancestry (positive or negative).
#' Second, if a SNP's alpha or beta falls outside the quantile interval
#' qn/2 and (qn-1)/2, it is considered an outlier (positive or negative).
#' If both criteriea are met for alpha or beta, crazy.a or crazy.b are TRUE.
#' These criteria are discussed in the Gompert BGC papers.
#'
#' Get credible intervals of a BGC parameter
#' @param df data.frame of one BGC parameter (e.g. alpha)
#' @return data.frame with credible intervals for the parameter
#' @noRd
get_cis <- function(df){
  # Get credible interval for parameters.
  cis <- apply(df, 1, function(x) bayestestR::ci(x = x, ci=0.95))
  cis <- do.call("rbind", cis)
  param.mean <- rowMeans(as.matrix(df))
  param.median <- apply(df, 1, function(x) median(x))

  df2 <- data.frame(mean=param.mean,
                    median=param.median,
                    lb=cis$CI_low,
                    ub=cis$CI_high,
                    stringsAsFactors = FALSE)

  # Clear memory.
  rm(param.mean, param.median, cis)
  gc()
  return(df2)
}

#' Get alpha and beta outliers
#' @param a.out data.frame of alpha BGC results
#' @param b.out data.frame of beta BGC results
#' @param qa.out data.frame of qa (gamma-quantile) BGC results
#' @param qb.out data.frame of qb (zeta-quantile) BGC results
#' @param loci Path to two-column file containing loci names (in BGC order)
#' @param qn Upper quantile interval value
#' @return data.frame containing SNP outlier info
#' @noRd
get_ab_outliers <- function(a.out, b.out, qa.out, qb.out, loci, qn){
  # Get excess ancestry and outliers for alpha and beta parameters.
  names(a.out)[3:4] <- c('lb','ub')
  names(b.out)[3:4] <- c('lb','ub')
  names(qa.out)[3:4] <- c('lb','ub')
  names(qb.out)[3:4] <- c('lb','ub')

  # Excess ancestry compared to genome-wide average introgression
  # If lb > 0 or ub < 0, it has excess ancestry for P0 or P1, respectively.
  # Do alpha.
  a.out$excess <- NA
  a.out$excess[a.out$lb > 0] <- 'pos'
  a.out$excess[a.out$ub < 0] <- 'neg'

  # Beta.
  b.out$excess <- NA
  b.out$excess[b.out$lb > 0] <- 'pos'
  b.out$excess[b.out$ub < 0] <- 'neg'

  # Add in quantiles
  a.out$q <- qa.out$mean
  b.out$q <- qb.out$mean

  # qnorm takes SD, so take the square root of the quantile estimate.
  # This basically generates the interval to check if alpha and beta are
  # outliers.
  # The interval is defined by qn = n/2 and qn = 1-n/2.
  # If the median for alpha or beta fall outside that interval, it's an
  # outlier.
  # By default, qn is set to 0.025 and 0.975
  qn_lower <- 1 - qn

  a.out$qlb <- qnorm(qn_lower,0,sqrt(a.out$q))
  a.out$qub <- qnorm(qn,0,sqrt(a.out$q))

  b.out$qlb <- qnorm(qn_lower,0,sqrt(b.out$q))
  b.out$qub <- qnorm(qn,0,sqrt(b.out$q))

  # Check if median falls outside qn interval.
  a.out$outlier <- NA
  a.out$outlier[a.out$median > a.out$qub] <- 'pos'
  a.out$outlier[a.out$median < a.out$qlb] <- 'neg'

  b.out$outlier <- NA
  b.out$outlier[b.out$median > b.out$qub] <- 'pos'
  b.out$outlier[b.out$median < b.out$qlb] <- 'neg'

  # Read loci file generated using my vcf2bgc.py script.
  # Contains actual scaffold names.
  snps <- read.delim(loci, sep = " ", stringsAsFactors = FALSE)

  ### Copy the data to snps data.frame.
  # Mean a and b
  snps$alpha <- a.out$mean
  snps$beta <- b.out$mean

  # Excess ancestry and outliers.
  snps$alpha.excess <- a.out$excess
  snps$beta.excess <- b.out$excess
  snps$alpha.outlier <- a.out$outlier
  snps$beta.outlier <- b.out$outlier

  # If has both excess ancestry AND is an outlier.
  snps$crazy.a <- !is.na(snps$alpha.excess) & !is.na(snps$alpha.outlier)
  snps$crazy.b <- !is.na(snps$beta.excess) & !is.na(snps$beta.outlier)

  # Clear memory.
  rm(a.out, b.out, qa.out, qb.out, loci)
  gc()

  return(snps)
}

#' Gets hybrid index mean across MCMC iterations for each individual
#' @param df data.frame of BGC hybrid index results
#' @return Matrix containing MCMC hybrid index means
#' @noRd
get_hi <- function(df){
  hi.mean <- rowMeans(as.matrix(df))
  return(hi.mean)
}

#' Parses BGC output and identifies outlier SNPs
#' @param df.list List of data.frame objects from combine_bgc_output()
#' @param admix.pop String identifying admixed population in population map
#'                  file
#' @param popmap Population map file consisting of two tab-separated columns:
#'               indID   popID (no header line)
#' @param qn Upper quantile interval value for determining outlier SNPs.
#'           Default = 0.975
#' @param loci.file Path to two-column tab-separated file containing loci names
#'                  (in same order as BGC input). Column 1 should be locus
#'                  name, column2 should be the position of the SNP.
#'                  header line should be: #CHROM POS
#' @param save.obj File path. If supplied, saves outlier info as .RDS object.
#'                 The value for save.obj should be the filename for the
#'                 RDS object
#' @return List containing 3 data.frame objects. 1) outlier info,
#'         2)popmap info, 3) hybrid index info
#' @export
#' @examples
#' outliers <- get_bgc_outliers(df.list = aggregated.results,
#'                              admix.pop = "population1",
#'                              popmap = "./popmap.txt",
#'                              qn = 0.975,
#'                              loci.file = "./bgc_loci.txt",
#'                              save.obj = "./outlierInfo.rds")
get_bgc_outliers <- function(df.list,
                             admix.pop,
                             popmap,
                             qn = 0.975,
                             loci.file,
                             save.obj = NULL){

  # Get 95% credible intervals for alpha, beta, gamma-quantile (qa),
  # zeta-quantile (qb).
  a <- get_cis(df.list[[2]])
  gc()
  b <- get_cis(df.list[[3]])
  gc()
  qa <- get_cis(df.list[[4]])
  gc()
  qb <- get_cis(df.list[[5]])
  gc()

  if (!file.exists(loci.file)){
    stop(paste0("Error: The loci file ", loci.file, " does not exist!"))
  }

  if (!file.exists(popmap)){
    stop(paste0("Error: The popmap file ", popmap, " does not exist!"))
  }

  # Find SNPs with excess ancestry and outliers SNPs.
  snps <- get_ab_outliers(a, b, qa, qb, loci.file, qn)
  gc()

  # Get hybrid indexes.
  hi <- get_hi(df.list[[6]])

  # Read in popmap.
  # Tab-separated, two-columns: indID\tpopID
  popmap <- read.table(file = popmap, stringsAsFactors = FALSE)

  # Get data.frame with indIDs of admixed individuals and hybrid index.
  admix <- popmap[ which(popmap$V2 == admix.pop),]
  hi.df <- data.frame(id=admix$V1, hi=hi, stringsAsFactors = FALSE)

  res <- list(snps, popmap, hi.df)

  rm(a, b, qa, qb, snps, hi, popmap, admix, hi.df)
  gc()
  # If user specified filepath: save res as rds object.
  if(!is.null(save.obj)){
    saveRDS(res, save.obj)
  }

  gc()
  return(res)

}

