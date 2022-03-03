#' Aggregate BGC Output
#'
#' This function aggregates multiple BGC runs into one file. Runs are joined
#' by columns and must all contain the same number of post-burnin MCMC
#' iterations. It assumes that input files are contained in a single directory
#' and have specific file suffixes. See the README.md file for details on
#' file suffixes. Files can have population prefixes.
#'
#' @param results.dir Path to a directory containing all the BGC results
#' @param prefix Prefix to the input files; e.g., "population1"
#' @param thin Thin to every n MCMC samples. E.g. thin=2 cuts samples in half
#' @param discard Discard first N samples from each BGC run
#' @return A list of data.frames (each data.frame is a BGC parameter)
#' @export
#' @examples
#' aggregate.results <- combine_bgc_output(results.dir = "./results",
#'                                         prefix = "population1")
#'
#' aggregate.results <- combine_bgc_output(results.dir = "./results",
#'                                         prefix = "population2",
#'                                         thin = 2)
#'
#' aggregate.results <- combine_bgc_output(results.dir = "./results",
#'                                         prefix = "pop3",
#'                                         discard = 2000)
combine_bgc_output <- function(results.dir,
                               prefix,
                               thin = NULL,
                               discard = NULL){
  
  writeLines(paste0("\n\nLoading input files with prefix ", prefix, "...\n"))
  
  lnl <- list.files(path = results.dir,
                    pattern = paste0(prefix, ".*LnL_\\d+$"),
                    full.names = TRUE)
  check_if_files(lnl)
  
  a0 <- list.files(path = results.dir,
                   pattern = paste0(prefix, ".*a0_\\d+$"),
                   full.names = TRUE)
  check_if_files(a0)
  
  b0 <- list.files(path = results.dir,
                   pattern = paste0(prefix, ".*b0_\\d+$"),
                   full.names = TRUE)
  check_if_files(b0)
  
  qa <- list.files(path = results.dir,
                   pattern = paste0(prefix, ".*qa_\\d+$"),
                   full.names = TRUE)
  check_if_files(qa)
  
  qb <- list.files(path = results.dir,
                   pattern = paste0(prefix, ".*qb_\\d+$"),
                   full.names = TRUE)
  check_if_files(qb)
  
  hi <- list.files(path = results.dir,
                   pattern = paste0(prefix, ".*hi_\\d+$"),
                   full.names = TRUE)
  check_if_files(hi)
  
  gc()
  
  if (!countBGCruns(list(lnl, a0, b0, qa, qb, hi))){
    stop("Error: All parameters must have the same number of BGC runs!")
  }
  
  if (!is.null(thin)){
    writeLines(paste0("\n\nThinning MCMC samples with frequency: ", thin))
  }
  if (!is.null(discard)){
    writeLines(paste0("\n\nDiscarding samples as burnin: ", discard))
  }
  
  # Aggregate the BGC runs into one data.frame
  lnl.df <- bind_bgc_runs(lnl, thin = thin, discard = discard)
  gc()
  a0.df <- bind_bgc_runs(a0, thin = thin, discard = discard)
  gc()
  b0.df <- bind_bgc_runs(b0, thin = thin, discard = discard)
  gc()
  qa.df <- bind_bgc_runs(qa, thin = thin, discard = discard)
  gc()
  qb.df <- bind_bgc_runs(qb, thin = thin, discard = discard)
  gc()
  hi.df <- bind_bgc_runs(hi, thin = thin, discard = discard)
  gc()
  
  writeLines(paste0("\n\n # MCMC samples = ", ncol(lnl.df)))
  

writeLines(paste0("\n\nDone! Found ", length(lnl), " BGC runs.\n\n"))

df.list <-
  list(lnl=lnl.df, a0=a0.df, b0=b0.df, qa=qa.df, qb=qb.df, hi=hi.df)

return(df.list)
}

#' Validates that input files were found
#' @param files A vector of file paths
#' @noRd
check_if_files <- function(files){
  # Make sure list.files found the correct input files.
  if (length(files) == 0){
    stop(paste0("Error: The files ",
                files,
                " were not found in the input directory."))
  }
}

#' Combine bgc runs together by binding columns (uses dplyr and data.table)
#' @param files A vector of file paths
#' @param thin Boolean. If TRUE, thins MCMC to every nth sample
#' @noRd
bind_bgc_runs <- function(files, thin = NULL, discard = NULL){
  
  data_list <- lapply(files, data.table::fread, stringsAsFactors = FALSE)
  
  if (!is.null(thin)){
    data_list <- lapply(data_list, function(x) thin_MCMC(x, thin))
  }
  if (!is.null(discard)){
    data_list <- lapply(data_list, function(x) discardNsamples(x, discard))
  }
  
  df <- dplyr::bind_cols(data_list)
  
  rm(data_list)
  gc()
  return(df)
}

#' Validate that all parameters have same number of BGC runs
#' @param l List of filepath vectors
#' @noRd
countBGCruns <- function(l){
  return(length(unique(sapply(l, length))) == 1)
}

#' Thin MCMC iterations every Nth iteration.
#' @param df data.frame with each MCMC sample as a column
#' @param n Only keep every nth iteration
#' @noRd
thin_MCMC <- function(df, n){
  df.new <- df[ ,seq(from = 1, to = ncol(df), by = n)]
  rm(df)
  gc()
  return(df.new)
}

#' Discard first N samples of each run
#' @param df data.frame containing
#' @noRd
discardNsamples <- function(df, n){
  df.new <- df[,(n+1):ncol(df)]
  rm(df)
  gc()
  return(df.new)
}
