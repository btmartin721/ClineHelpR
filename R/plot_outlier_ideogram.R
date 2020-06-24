#' Plot an Ideogram for BGC Outlier Genes
#'
#' For this function to work, you must have a closely related reference genome
#' and must run minimap2 externally (tested with v2.17). minimap2 maps assembly
#' scaffolds to a reference genome.
#' This includes all the unplaced assembly
#' scaffolds as well as the transcriptome combined into one assembly fasta
#' file. minimap2 should be run with the --cs and -N 100 options with default
#' PAF file output format. You will need to adjust the asm value to suit how
#' closely related your reference genome is (i.e. sequence divergence).
#' See the minimap2 documentation for more info:
#' https://github.com/lh3/minimap2. Then, you need to run PAFScaff on the PAF
#' file output from minimap2 (https://github.com/slimsuite/pafscaff). PAFScaff
#' parses the minimap2 output and improves the mapping. The scaffolds.tdt.file
#' value is the *.scaffolds.tdt output file from PAFScaff, and you need it
#' for this function. Ref chromosomes should be named integers in PAFScaff
#' (1, 2, etc.). This funtion also depends on previous output from
#' get_bgc_outliers() and join_bgc_gff(). See ?get_bgc_outliers and
#' ?join_bgc_gff. You should run BGC and get_bgc_outliers() separately
#' for both transcriptome-aligned data and all other scaffolded loci. Then
#' run join_bgc_gff() on only the transcriptome-aligned data.
#' Both are required as input to make the ideograms. You can also see
#' https://cran.r-project.org/web/packages/RIdeogram/index.html for more
#' info on plotting the RIdeograms.
#'
#'
#' Function to plot outlier BGC loci as heatmaps on chromosome ideograms.
#' @param prefix Prefix for output files
#' @param outliers.genes.annotated List containing gene outlier data from
#'                                 get_bgc_outliers(). See ?get_bgc_outliers.
#'                                 This must be outliers from a transcriptome
#'                                 alignment
#' @param outliers.full.scaffolds List containing outlier data from
#'                                get_bgc_outliers(). See ?get_bgc_outliers.
#'                                This must be loci aligned to full scaffolds
#' @param pafInfo Path to *.scaffolds.tdt file output from PAFScaff
#' @param plotDIR Directory to save output plots
#' @param both.outlier.tests Boolean; If TRUE, scaffold outliers must meet both the
#'                           overlap.zero and qn.interval criteria
#' @param both.outlier.tests.genes Boolean; If TRUE, gene outliers must meet
#'                                 both the overlap.zero.genes and
#'                                 qn.interval.genes criteria
#' @param overlap.zero Boolean; If TRUE, scaffold outliers are SNPs whose
#'                     credible interval does not contain zero
#' @param overlap.zero.genes Boolean; If TRUE, gene outliers are SNPs whose
#'                           credible interval does not contain zero
#' @param qn.interval Boolean; If TRUE, scaffold outliers fall outside the
#'                    quantile interval qn/2 and 1-qn/2
#' @param qn.interval.genes Boolean; If TRUE, gene outliers fall outside the
#'                          quantile interval qn/2 and 1-qn/2
#' @param missing.chrs If specified, must be character vector of missing
#'                     chromosome names. Chromosome numbers should be prefixed
#'                     with "chr". I.e., c("chr3", "chr6"). If some chromosomes
#'                     don't get plotted, use this option
#' @param missing.chr.length Vector of integer lengths (in bp) of
#'                           missing chromosomes. Must also be specified
#'                           if missing.chrs is used. Vector must also be
#'                           the same length as missing.chrs
#' @param gene.size Adjust the size for each outlier transcriptome gene
#'                  on the ideogram. If the loci appear too small or large
#'                  on the ideogram, adjust just this parameter
#'
#' @param other.size Adjust the size for each outlier scaffold gene on the
#'                   ideogram
#' @param convert_svg Device to convert SVG output plot. Default is pdf, but
#'                    you can use png or other commonly used devices
#' @param colorset1 Vector of colors for RIdeogram alpha heatmap.
#'                  Default is the same as the RIdeogram defaults
#' @param colorset2 Vector of colors for RIdeogram beta heatmap.
#'                  Default is the same as the RIdeogram defaults
#' @param chrnum.prefix Prefix for chromosome numbers on ideaogram plot
#' @param genes.only Boolean; If TRUE, only include known genes on ideogram
#' @param linked.only Boolean; If TRUE, only include non-genes on ideogram
#' @return Data.frame containing reference info for gene outliers.
#' @export
#' @examples
#' plot_outlier_ideogram(prefix = "population1",
#'                       outliers.genes.annotated = outliers.genes.annotated,
#'                       outliers.full.scaffolds = outliers.full.scaffolds,
#'                       pafInfo = "./population1.scaffolds.tdt",
#'                       plotDIR = "./plots",
#'                       both.outlier.tests = TRUE,
#'                       missing.chrs = c("chr11", "chr21", "chr25"),
#'                       miss.chr.length = c(4997863, 1374423, 1060959),
#'                       gene.size = 1e6,
#'                       other.size = 5e5,
#'                       convert_svg = "png",
#'                       colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
#'                       colorset2 = c("#fc8d59", "#ffffbf", "#91bfdb")
plot_outlier_ideogram <- function(prefix,
                                  outliers.genes,
                                  outliers.full.scaffolds,
                                  pafInfo,
                                  plotDIR = "./plots",
                                  both.outlier.tests = FALSE,
                                  both.outlier.tests.genes = FALSE,
                                  overlap.zero = TRUE,
                                  overlap.zero.genes = TRUE,
                                  qn.interval = TRUE,
                                  qn.interval.genes = TRUE,
                                  missing.chrs = NULL,
                                  miss.chr.length = NULL,
                                  gene.size = 5e5,
                                  other.size = 1e5,
                                  convert_svg = "pdf",
                                  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
                                  colorset2 = c("#4575b4", "#ffffbf", "#d73027"),
                                  chrnum.prefix = NULL,
                                  genes.only = FALSE,
                                  linked.only = FALSE){

  ###############################################################
  ### Read input objects and files.
  ###############################################################

  # Make output directory if it doesn't already exist.
  dir.create(plotDIR, showWarnings = FALSE)

  # Transcript outliers that were present in the GFF file.
  scaffold.trans.info <- outliers.genes

  # Copy full dataset info. Obtained using PAFScaff and minimap
  scaffold.ref.info <- clean_sort_pafscaff(pafInfo)

  # Copy outlier data.
  full.outliers <- outliers.full.scaffolds
  full.scaf <- full.outliers[[1]]

  #############################################################
  ### Parse and trim data for joins below.
  #############################################################

  ## Remove GenBank version numbers from scaffold ids.
  ## It can result in df joins missing some scaffolds if versions differ.
  full.scaf$ScafId <- gsub("\\..*", "", full.scaf$X.CHROM)
  scaffold.ref.info$ScafId <- gsub("\\..*", "", scaffold.ref.info$Scaffold)
  scaffold.trans.info$ScafId <- gsub("\\..*", "", scaffold.trans.info$seqid)

  # Subset transcript scaffold info if starts with "X" to get just transcripts.
  transcript.ref.info <-
    scaffold.ref.info[gdata::startsWith(scaffold.ref.info$ScafId, "X"),]

  #################################################################
  ### Join data frames to associate assembly and ref scaffolds.
  #################################################################

  # Join dataframes to where all info is contained in one df.
  # Preserve all ref data.
  ref.merged <-
    merge(scaffold.ref.info, full.scaf, by = "ScafId", all = FALSE)

  ref.trans.merged <- merge(scaffold.trans.info,
                            transcript.ref.info,
                            by.x = "TranscriptId",
                            by.y = "ScafId",
                            all = FALSE)

  # Get SNP position on scaffold:
  ref.trans.merged$snpPOS <-
    ref.trans.merged$RefStart + ref.trans.merged$POS

  if (both.outlier.tests.genes){

    overlap.zero.genes <- FALSE
    qn.interval.genes <- FALSE

    writeLines("\n\nboth.outlier.tests is TRUE for genes.")
    writeLines("Ignoring overlap.zero qn.interval settings\n")

    # Subset extreme outliers
    # (had both alpha excess AND was alpha outlier. Likewise for beta)
    ref.trans.merged.alpha <- ref.trans.merged[ref.trans.merged$crazy.a == TRUE,]
    ref.trans.merged.beta <- ref.trans.merged[ref.trans.merged$crazy.b == TRUE,]
    ref.trns.merged.either <-
      ref.trans.merged[(ref.trans.merged$crazy.a == TRUE |
                          ref.trans.merged$crazy.b == TRUE),]

    if (nrow(ref.trans.merged.alpha) == 0){
      stop("No alpha outliers in full dataset. Set both.outlier.tests = FALSE")
    }

    if (nrow(ref.trans.merged.beta) == 0){
      stop("No beta outliers in full dataset. Set both.outlier.tests = FALSE.")
    }
  }

  # TRUE if alpha/beta excess/outliers.
  if (overlap.zero.genes & qn.interval.genes){
    ref.trans.merged.alpha <-
      ref.trans.merged[(!is.na(ref.trans.merged$alpha.excess) |
                    !is.na(ref.trans.merged$alpha.outlier)),]
    ref.trans.merged.beta <-
      ref.trans.merged[(!is.na(ref.trans.merged$beta.excess) |
                    !is.na(ref.trans.merged$beta.outlier)),]

    ref.trans.merged.either <-
      ref.trans.merged[(
        !is.na(ref.trans.merged$alpha.excess) |
          !is.na(ref.trans.merged$alpha.outlier) |
          !is.na(ref.trans.merged$beta.excess) |
          !is.na(ref.trans.merged$beta.outlier)
      ), ]

  } else if (overlap.zero.genes & !qn.interval.genes){
    ref.trans.merged.alpha <- ref.trans.merged[!is.na(ref.trans.merged$alpha.excess),]
    ref.trans.merged.beta <- ref.trans.merged[!is.na(ref.trans.merged$beta.excess),]
    ref.trans.merged.either <-
      ref.trans.merged[(!is.na(ref.trans.merged$alpha.excess) |
                          !is.na(ref.trans.merged$beta.excess)), ]

  } else if (!overlap.zero.genes & qn.interval.genes){
    ref.trans.merged.alpha <- ref.trans.merged[!is.na(ref.trans.merged$alpha.outlier),]
    ref.trans.merged.beta <- ref.trans.merged[!is.na(ref.trans.merged$beta.outlier),]

    ref.trans.merged.either <-
      ref.trans.merged[(!is.na(ref.trans.merged$alpha.outlier) |
                          !is.na(ref.trans.merged$beta.outlier)), ]

  } else if (!overlap.zero.genes & !qn.interval.genes & !both.outlier.tests.genes){
    mywarning <-
      paste(
        "\n\noverlap.zero, qn.interval, and both.outlier.tests",
        "were all FALSE. Setting both.outlier.tests to TRUE\n\n"
      )
    warning(paste(strwrap(mywarning), collapse = "\n"))

    # Subset extreme outliers
    # (had both alpha excess AND was alpha outlier. Likewise for beta)
    ref.trans.merged.alpha <- ref.trans.merged[ref.trans.merged$crazy.a == TRUE,]
    ref.trans.merged.beta <- ref.trans.merged[ref.trans.merged$crazy.b == TRUE,]
    ref.trans.merged.either <-
      ref.trans.merged[(ref.trans.merged$crazy.a == TRUE |
                          ref.trans.merged$crazy.b == TRUE),]
  }

  # Remove rows with NA in snpPOS. NAs caused issues with RIdeogram.
  ref.trans.merged.alpha <-
    ref.trans.merged.alpha[!is.na(ref.trans.merged.alpha$snpPOS), ]
  ref.trans.merged.beta <-
    ref.trans.merged.beta[!is.na(ref.trans.merged.beta$snpPOS), ]
  ref.trans.merged.either <-
    ref.trans.merged.either[!is.na(ref.trans.merged.either$snpPOS), ]

  ############################################################
  ### Run RIDEOGRAM
  ############################################################

  ############################################################
  ### MAKE KARYOTYPE FILE
  ############################################################

  # Columns for chromosomes are supposed to be as follows:
  # Chr: a character representing the chromosome/contig/region name like
  #      'chr1' or '1' or 'ch1'
  # Start: a numeric value to specify chromosome (or chromosome region)
  #        start position. If you are considering entire chromosome this
  #        value is typically 0.
  # End: a numeric value specifying chromosome/contig/region end position.
  # Again, if you are considering entire chromosome,
  #       then this value is the length of chromosome.

  # Columns for heatmap annotations are supposed to be as follows:
  # Chr: a character specifying the chromosome name.
  #        [NOTE: the chromosome names should be consistent in chromosome
  #        and data files.]
  # Start: A numeric specifying element start position.
  # End: A numeric specifying element end position.
  # Value: A numeric specifying the data value.

  if (!is.null(chrnum.prefix)){
    # Create a data.frame with the chromosome information.
    chrom.df <- data.frame(Chr=paste0(chrnum.prefix, ref.merged$Ref),
                           Start=as.integer(0),
                           End=ref.merged$RefLen,
                           stringsAsFactors = FALSE)
  } else if (is.null(chrnum.prefix)){
    # Create a data.frame with the chromosome information.
    chrom.df <- data.frame(Chr=ref.merged$Ref,
                           Start=as.integer(0),
                           End=ref.merged$RefLen,
                           stringsAsFactors = FALSE)
  }

  chrom.df <- chrom.df[!duplicated(chrom.df$Chr),]
  chrom.df <- chrom.df[order(chrom.df$Chr),]

  # If missing.chrs and miss.chr.length are defined.
  if (!is.null(missing.chrs) & !is.null(miss.chr.length)){

      if (length(missing.chrs) != length(miss.chr.length)){
        stop("Error: missing.chrs and miss.chr.length must both
             be the same length")
      }

      # Make new df of missing chromosomes.
      missing.chr.df <- data.frame(Chr=missing.chrs,
                                   Start=as.integer(0),
                                   End=as.integer(miss.chr.length),
                                   stringsAsFactors = FALSE)

      chrom.df <- rbind(chrom.df, missing.chr.df)

  } else if (!is.null(missing.chrs) & is.null(miss.chr.length)){
      stop("Error: missing.chrs and miss.chr.length must have
           both defined or neither defined.")
  } else if (!is.null(miss.chr.length) & is.null(missing.chrs)){
      stop("Error: missing.chrs and miss.chr.length must have
           both defined or neither defined.")
  }

  # Sort numerically instead of alphabetically.
  # R sorts strings alphabetically by default,
  # and it was sorting the chroms out of order.
  # e.g. chr1 chr10 etc.
  chrom.df$Chr <- as.character(chrom.df$Chr)
  chrom.df <- chrom.df[gtools::mixedorder(as.character(chrom.df$Chr)),]


  ########################################################
  ### MAKE HEATMAP FILES - GENES ONLY
  ########################################################


  # Make input objects for RIdeogram.
  # This one is for genes only.
  # Made the genes bigger than the genomic loci for clarity.
  # Added gene.size length to each SNP.
  # Do both alpha and beta outliers.

  if (!is.null(chrnum.prefix)){

  heatmap.alpha <-
      data.frame(
        Chr = as.character(paste0(chrnum.prefix, ref.trans.merged.alpha$Ref)),
        Start = ref.trans.merged.alpha$RefStart,
        End = ref.trans.merged.alpha$RefStart+gene.size,
        Value = ref.trans.merged.alpha$alpha,
        stringsAsFactors = FALSE
      )

  # Genes only.
  heatmap.beta <-
    data.frame(
      Chr = as.character(paste0(chrnum.prefix, ref.trans.merged.beta$Ref)),
      Start = ref.trans.merged.beta$RefStart,
      End = ref.trans.merged.beta$RefStart+gene.size,
      Value = ref.trans.merged.beta$beta,
      stringsAsFactors = FALSE
    )

  } else if (is.null(chrnum.prefix)){
    heatmap.alpha <-
      data.frame(
        Chr = as.character(ref.trans.merged.alpha$Ref),
        Start = ref.trans.merged.alpha$RefStart,
        End = ref.trans.merged.alpha$RefStart+gene.size,
        Value = ref.trans.merged.alpha$alpha,
        stringsAsFactors = FALSE
      )

    heatmap.beta <-
      data.frame(
        Chr = as.character(ref.trans.merged.beta$Ref),
        Start = ref.trans.merged.beta$RefStart,
        End = ref.trans.merged.beta$RefStart+gene.size,
        Value = ref.trans.merged.beta$beta,
        stringsAsFactors = FALSE
      )
  }




  # Filter out outliers that went beyond chromosome length.
  # Probably because many of the chromosomes are incomplete.
  heatmap.alpha <- filter_overhang(heatmap.alpha, chrom.df)
  heatmap.beta <- filter_overhang(heatmap.beta, chrom.df)

  ########################################################
  ### MAKE HEATMAP FILES - GENOMIC LOCI
  ########################################################

  # Get positions of actual SNPs.
  ref.merged$snpPOS <- ref.merged$RefStart + ref.merged$POS

  if (both.outlier.tests){

    overlap.zero = FALSE
    qn.interval = FALSE

    writeLines("\n\nboth.outlier.tests is TRUE.")
    writeLines("Ignoring overlap.zero qn.interval settings\n")

    # Subset extreme outliers
    # (had both alpha excess AND was alpha outlier. Likewise for beta)
    ref.merged.alpha <- ref.merged[ref.merged$crazy.a == TRUE,]
    ref.merged.beta <- ref.merged[ref.merged$crazy.b == TRUE,]

    if (nrow(ref.merged.alpha) == 0){
      stop("No alpha outliers in full dataset. Set both.outlier.tests = FALSE")
    }

    if (nrow(ref.merged.beta) == 0){
      stop("No beta outliers in full dataset. Set both.outlier.tests = FALSE.")
    }
  }

  # TRUE if alpha/beta excess/outliers.
  if (overlap.zero & qn.interval){
    ref.merged.alpha <-
      ref.merged[(!is.na(ref.merged$alpha.excess) |
                   !is.na(ref.merged$alpha.outlier)),]
    ref.merged.beta <-
      ref.merged[(!is.na(ref.merged$beta.excess) |
                    !is.na(ref.merged$beta.outlier)),]

  } else if (overlap.zero & !qn.interval){
      ref.merged.alpha <- ref.merged[!is.na(ref.merged$alpha.excess),]
      ref.merged.beta <- ref.merged[!is.na(ref.merged$beta.excess),]

  } else if (!overlap.zero & qn.interval){
      ref.merged.alpha <- ref.merged[!is.na(ref.merged$alpha.outlier),]
      ref.merged.beta <- ref.merged[!is.na(ref.merged$beta.outlier),]

  } else if (!overlap.zero & !qn.interval & !both.outlier.tests){
    mywarning <-
      paste(
        "\n\noverlap.zero, qn.interval, and both.outlier.tests",
        "were all FALSE. Setting both.outlier.tests to TRUE\n\n"
      )
    warning(paste(strwrap(mywarning), collapse = "\n"))

    # Subset extreme outliers
    # (had both alpha excess AND was alpha outlier. Likewise for beta)
    ref.merged.alpha <- ref.merged[ref.merged$crazy.a == TRUE,]
    ref.merged.beta <- ref.merged[ref.merged$crazy.b == TRUE,]
  }

  # Remove rows with NA in snpPOS. NAs caused issues with RIdeogram.
  ref.merged.alpha <- ref.merged.alpha[!is.na(ref.merged.alpha$snpPOS),]
  ref.merged.beta <- ref.merged.beta[!is.na(ref.merged.beta$snpPOS),]

  # Now make the input objects for all the other outliers.
  # Made length longer for clarity.
  # Made length snpPOS + 1e5. Smaller than genes though.

  if (!is.null(chrnum.prefix)){
    linked.alpha <-
      data.frame(
        Chr = as.character(paste0(chrnum.prefix, ref.merged.alpha$Ref)),
        Start = as.integer(ref.merged.alpha$snpPOS),
        End = as.integer(ref.merged.alpha$snpPOS + other.size),
        Value = as.numeric(as.character(ref.merged.alpha$alpha)),
        stringsAsFactors = FALSE
      )

    linked.beta <-
      data.frame(
        Chr = as.character(paste0(chrnum.prefix, ref.merged.beta$Ref)),
        Start = as.integer(ref.merged.beta$snpPOS),
        End = as.integer(ref.merged.beta$snpPOS + other.size),
        Value = as.numeric(as.character(ref.merged.beta$beta)),
        stringsAsFactors = FALSE
      )

  } else if (is.null(chrnum.prefix)){

    linked.alpha <-
      data.frame(
        Chr = as.character(ref.merged.alpha$Ref),
        Start = as.integer(ref.merged.alpha$snpPOS),
        End = as.integer(ref.merged.alpha$snpPOS+other.size),
        Value = as.numeric(as.character(ref.merged.alpha$alpha)),
        stringsAsFactors = FALSE
      )

    linked.beta <-
      data.frame(
        Chr = as.character(ref.merged.beta$Ref),
        Start = as.integer(ref.merged.beta$snpPOS),
        End = as.integer(ref.merged.beta$snpPOS+other.size),
        Value = as.numeric(as.character(ref.merged.beta$beta)),
        stringsAsFactors = FALSE
      )

  }

  # Remove loci that exteded past chromosome length.
  linked.alpha <- filter_overhang(linked.alpha, chrom.df)
  linked.beta <- filter_overhang(linked.beta, chrom.df)


  ########## ALL LOCI PLOTTED TOGETHER #################

  if (!genes.only & !linked.only){

    all.alpha.loci <- rbind(heatmap.alpha, linked.alpha)
    all.beta.loci <- rbind(heatmap.beta, linked.beta)

  } else if (genes.only & !linked.only){

    all.alpha.loci <- heatmap.alpha
    all.beta.loci <- heatmap.beta

  } else if (!genes.only & linked.only){
    all.alpha.loci <- linked.alpha
    all.beta.loci <- linked.beta

  } else if (genes.only & linked.only){
    warning("genes.only and linked.only were both set to TRUE. Including both")
    all.alpha.loci <- rbind(heatmap.alpha, linked.alpha)
    all.beta.loci <- rbind(heatmap.beta, linked.beta)
  }


  plotdf <-
    data.frame(
      alpha.min = min(all.alpha.loci$Value),
      alpha.max = max(all.alpha.loci$Value),
      beta.min = min(all.beta.loci$Value),
      beta.max = max(all.beta.loci$Value),
      stringsAsFactors = FALSE
    )

  # Run RIdeogram and make SVD and PDF plots.
  RIdeogram::ideogram(
    karyotype = chrom.df,
    overlaid = all.alpha.loci,
    label = all.beta.loci,
    label_type = "heatmap",
    colorset1 = colorset1,
    colorset2 = colorset2,
    output = file.path(plotDIR, paste0(prefix, "_ideogram.svg"))
  )

  # Converts SVG to specified format.
  RIdeogram::convertSVG(file.path(plotDIR,
                                  paste0(prefix, "_ideogram.svg")),
                        device = convert_svg)

  # Return data.frame with ref info for gene outliers.
  return(ref.trans.merged.either)

}

#' Function to Clean and Parse PAFScaff *.scaffolds.tdt file
#' @param scaffolds.tdt.file filepath to PAFScaff output file
#' @return data.frame with cleaned and sorted PAFScaff data
#' @noRd
clean_sort_pafscaff <- function(scaffolds.tdt.file){

  pafInfo <- read.table(file = scaffolds.tdt.file,
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        blank.lines.skip = TRUE,
                        quote="",
                        fill = FALSE)
  # Remove RevComp prefix from Description column.
  pafInfo$Description <-
    gsub(pattern = "^RevComp ", replacement = "", x = pafInfo$Description)

  # Remove unplaced scaffolds.
  pafInfo <- pafInfo[!grepl("^unplaced", pafInfo$Qry),]

  # Sort by Ref (Chromosome), RefStart, then RefEnd columns.
  pafInfo <- pafInfo[with(pafInfo, order(Ref, RefStart, RefEnd)), ]

  return(pafInfo)
}

#' Function to filter out loci that fall outside the chromosome length.
#' E.g., if the reference chromosomes are incomplete.
#' @param df data.frame with BGC outlier info formatted for RIdeogram
#' @param chrom.df data.frame with karyotype info formatted for RIdeogram
#' @return data.frame with loci exceeding the chromosome length removed
#' @noRd
filter_overhang <- function(df, chrom.df){
  df2 <- merge(df, chrom.df, by="Chr", all.y=TRUE)

  `%>%` <- dplyr::`%>%`

  tmp <- df2 %>%
    dplyr::group_by(Chr) %>%
    dplyr::filter(End.x < End.y)

  df3 <- data.frame(Chr=tmp$Chr, Start=tmp$Start.x,
                    End=tmp$End.x, Value=tmp$Value,
                    stringsAsFactors = FALSE)
  return(df3)
}
