#' Reads in BGC Outlier Data and Joins with GFF Data
#'
#' This function joins BGC outlier data with info from a GFF file.
#' It outputs a joined data.frame of BGC loci that match to GFF genes.
#' The BGC outlier loci names must be GenBank transcript IDs that can be
#' related to the transcript IDs in the GFF file The BGC outlier object must
#' also be generated using upstream methods. See ?get_bgc_outliers for
#' more info.


#' Read In BGC Data and Join With GFF Annotations
#' @param prefix Prefix for output file
#' @param outlier.list List generated from get_bgc_outliers().
#'                     See ?get_bgc_outliers.
#' @param gff.data Object output from read_parse_GFF(). See ?read_parse_GFF
#' @param scafInfoDIR Directory to save output from this function.
#' @return Data.frame containing outliers genes joined with GFF data
#' @export
#' @examples
#' outlier.genes <- join_bgc_gff(prefix = "population1",
#'                               outlier.list = outliers,
#'                               gff.data = gff.df,
#'                               scafInfoDIR = "./scaffold_info")
join_bgc_gff <- function(prefix,
                         outlier.list,
                         gff.data,
                         scafInfoDIR = "./scaffold_info"){

  # Function to join outlier BGC data and gene scaffold info.

  # Get data.frame of outlier snps.
  snps <- outlier.list[[1]]

  # Get df with only outlier loci (alpha or beta)
  outliers <-
    snps[ which(!is.na(snps$alpha.excess) |
                  !is.na(snps$alpha.outlier) |
                  !is.na(snps$beta.excess) |
                  !is.na(snps$beta.outlier)), ]

  # Remove version number from transcript ids and save to new column.
  outliers$TranscriptId <- gsub("(.*)\\..*", "\\1", outliers$X.CHROM)

  # Inner Join dataframes with matching TranscriptId columns.
  outliers.genes <- merge(outliers, gff.data, by = "TranscriptId", all = FALSE)

  # Get only rows with unique transcript ids.
  unique_rows <- !duplicated(outliers.genes["TranscriptId"])
  outliers.genes <- outliers.genes[unique_rows,]

  # Print out gene names for matched outlier loci.
  writeLines("\n\nOutlier genes matching to GFF genes:\n\n")
  print(outliers.genes$gene)

  # Create output directory if doesn't already exist.
  dir.create(path = scafInfoDIR, showWarnings = FALSE)

  # Save the output to csv file.
  write.csv(x = outliers.genes,
            file = file.path(scafInfoDIR,
                          paste0(prefix, "_scaffold_info.csv")))

  return(outliers.genes)
}
