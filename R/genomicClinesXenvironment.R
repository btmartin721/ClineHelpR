#' Functions to run INTROGRESS and correlate genomic clines X environment
#'
#' These functions run INTROGRESS and generate genomic clines. The genomic
#' clines are then correlated with the raster values at each sampling locality
#' for the raster layers determined to be the best predictors for species
#' distributions in MAXENT. See ?prepare_rasters, ?partition_raster_bg,
#' ?runENMeval, ?summarize_ENMeval, and ?extractPointValues.
#'
#' Function to run INTROGRESS
#' @param p1.input Character string or genind2introgress object;
#'                 File path for INTROGRESS p1 input file or
#'                 p1 object from genind2introgress
#' @param p2.file Character string or genind2introgress object;
#'                 File path for INTROGRESS p2 input file or p2 object from
#'                 genind2introgress
#' @param admix.file Character string or genind2introgress object;
#'                   File path for INTROGRESS admixed file or admix object
#'                   from genind2intrgoress
#' @param loci.file Character string or genind2introgress object;
#'                  File path for INTROGRESS loci file or
#'                  loci object from genind2introgress
#' @param clineLabels Character vector of length == 3 for c(P1, Het, P2)
#'                    populations
#' @param minDelt Numeric; Minimum allele frequency delta to retain loci for
#'                INTROGRESS
#' @param prefix Character string; Population prefix for input and output files
#' @param outputDIR Character string; File path for INTROGRESS output
#' @param showPLOTS Boolean; Whether to print the plots to the screen
#' @param sep Character string; Column delimiter for INTROGRESS input files
#' @param pop.id Boolean; Value for INTROGRESS parameter pop.id.
#'               See ?prepare.data in introgress R pacakge
#' @param ind.id Boolean; Value for INTROGRESS parameter ind.id
#'               See ?prepare.data in introgress R package
#' @param fixed Boolean; TRUE if input data are fixed SNPs
#' @return List object containing genomic cline output
#' @export
#' @examples
#' eatt <- runIntrogress(p1.file = "EATT_p1data.txt",
#'                       p2.file = "EATT_p2data.txt",
#'                       admix.file = "EATT_admix.txt",
#'                       loci.file = "EATT_loci.txt",
#'                       clineLabels = c("EA", "Het", "TT")
#'                       minDelt = 0.8,
#'                       prefix = "eatt",
#'                       outputDIR = "./introgress_plots",
#'                       sep = "\t",
#'                       pop.id = FALSE,
#'                       ind.id = FALSE,
#'                       fixed = FALSE)
runIntrogress <- function(p1.input,
                          p2.input,
                          admix.input,
                          loci.input,
                          clineLabels = c("P1", "Het", "P2"),
                          minDelt = 0.8,
                          prefix = "population1",
                          outputDIR = "./plots",
                          showPLOTS = FALSE,
                          sep = "\t",
                          pop.id = FALSE,
                          ind.id = FALSE,
                          fixed = FALSE){


  if (!requireNamespace("introgress", quietly = TRUE)){
    warning("introgress package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("rgdal", quietly = TRUE)){
    warning("rgdal package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("raster", quietly = TRUE)){
    warning("raster package must be installed to use this functionality")
    return(NULL)
  }

  if (length(clineLabels) != 3){
    stop("Error: clineLabels parameter must be character vector of length 3")
  }

  # Create outputDIR if doesn't exist.
  dir.create(outputDIR, showWarnings = FALSE)

  if (file.exists(admix.input) &&
      file.exists(loci.input) &&
      file.exists(p1.input) &&
      file.exists(p2.input)){

    # Read in the data from files
    gen.data <- read.table(file = admix.input, header = FALSE, sep = sep)
    loci.data <- read.table(file = loci.input, header = TRUE, sep = sep)
    p1 <- read.table(file = p1.input, header = FALSE, sep = sep)
    p2 <- read.table(file = p2.input, header = FALSE, sep = sep)

  } else {
    # Read in the objects
    gen.data <- admix.input
    loci.data <- loci.input
    p1 <- p1.input
    p2 <- p2.input
  }

  # Set the filenames of the PDF outfiles generated in script
  # This first block is for the permutation genomic cline method
  hIndex_histPDF <- file.path(outputDIR,
                              paste0(prefix, "_hIndexHistogram.pdf"))

  compositePDF <- file.path(outputDIR, paste0(prefix, "_composite.pdf"))
  singleClinesPDF <- file.path(outputDIR, paste0(prefix, "_indClines.pdf"))
  imagePDF <- file.path(outputDIR, paste0(prefix, "_image.pdf"))

  # Contains summary genomic cline data,
  summaryFile <- file.path(outputDIR, paste0(prefix, "_summaryData.txt"))


  # This block is for the parametric genomic cline method
  parametric_compositePDF <- file.path(outputDIR,
                                       paste0(prefix,
                                              "_parametric_composite.pdf"))
  parametric_singleClinesPDF <- file.path(outputDIR,
                                          paste0(prefix,
                                                 "_parametric_indClines.pdf"))
  parametric_summaryFile <- file.path(outputDIR,
                                      paste0(prefix,
                                             "_parametric_summaryData.txt"))

  # If fixed == TRUE this will be plotted.
  triPlotPDF <- file.path(outputDIR, paste0(prefix, "_trianglePlot.pdf"))

  # Set the delta value for keeping only loci with allele frequency
  # differentials > minDelt
  minDelt <- minDelt

  # Prep data for analysis
  count.matrix <- introgress::prepare.data(
    admix.gen = gen.data,
    loci.data = loci.data,
    parental1 = p1,
    parental2 = p2,
    pop.id = FALSE,
    ind.id = FALSE
  )


  # Calulate hybrid index for each individual using maximum-likelihood
  hi.index <-
    introgress::est.h(introgress.data = count.matrix,
                      loci.data = loci.data,
                      fixed = fixed)

  # Plot distribution of hybrid indices
  pdf(hIndex_histPDF)
    hist(hi.index$h,breaks = 20)
  dev.off()

  if (isTRUE(showPLOTS)){
    print(hist(hi.index$h,breaks = 20))
  }

  if (fixed){
    # Plot H.index x interspecifit heterozygosity
    int.h <- introgress::calc.intersp.het(introgress.dat = count.matrix)
    par(mfrow = c(1, 1))

    if (isTRUE(showPLOTS)){
      introgress::triangle.plot(hi.index = hi.index,
                                int.het = int.h,
                                pdf = FALSE)
    }

    introgress::triangle.plot(hi.index = hi.index,
                  int.het = int.h,
                  pdf = TRUE,
                  out.file = triPlotPDF)
  }

  # Calculate a measure of allele frequency different between parental pops
  # this is why you need to make sure no hybrids are in the parental pop
  # references
  deltas <- vector()
  for (i in 1:length(loci.data$locus)) {
    deltas[[i]] <-
      introgress::delta(count.matrix$Parental1.allele.freq[i, 1:2],
                         count.matrix$Parental2.allele.freq[i, 1:2])
    }

   # Estimate genomic clines ONLY for loci with high delta distance between
  # Parent1 and Parent2
  perm_clines <- introgress::genomic.clines(
    introgress.data = count.matrix,
    hi.index = hi.index,
    loci.data = loci.data,
    sig.test = TRUE,
    classification = TRUE,
    loci.touse = which(deltas > minDelt)
  )

  if (isTRUE(showPLOTS)) {
    #Plot genomic clines
    introgress::composite.clines(perm_clines,
                                 pdf = FALSE,
                                 labels = clineLabels)

    # There was a bug with the mk.image feature, so I removed it.
    # Not sure what was wrong with it.
  }

  # Plot genomic clines
  introgress::composite.clines(perm_clines,
                   pdf = TRUE,
                   out.file = compositePDF,
                   labels = clineLabels)


  # There was a bug with the mk.image feature, so I removed it.
  # Not sure what was wrong with it.

  # Pvalues for genomic clines.
  pValues <- perm_clines$Summary.data[,4]

  write.table(pValues,
              file.path(outputDIR, paste0(prefix, "_pvalues.csv")),
              quote = FALSE,
              row.names = FALSE)

  if (isTRUE(showPLOTS)) {
    # Plot the clines.
    introgress::clines.plot(
      cline.data = perm_clines,
      rplots = 3,
      cplots = 3,
      pdf = FALSE,
      quantiles = TRUE,
      cd = clineLabels
    )
  }

  # Plot the clines.
  introgress::clines.plot(
    cline.data = perm_clines,
    rplots = 3,
    cplots = 3,
    pdf = TRUE,
    out.file = singleClinesPDF,
    quantiles = TRUE,
    cd = clineLabels
  )

  # Esimate genomic clines ONLY for loci with high delta distance between
  # Parent1 and Parent2
  parametric_clines <- introgress::genomic.clines(
    introgress.data = count.matrix,
    hi.index = hi.index,
    loci.data = loci.data,
    sig.test = TRUE,
    classification = TRUE,
    loci.touse = which(deltas > minDelt),
    method = "parametric"
  )

  if (isTRUE(showPLOTS)) {
    # Plot parametric composite clines
    introgress::composite.clines(
      parametric_clines,
      pdf = FALSE,
      labels = clineLabels
    )

    # Plot parametric clines.
    introgress::clines.plot(
      cline.data = parametric_clines,
      rplots = 3,
      cplots = 3,
      pdf = FALSE,
      quantiles = TRUE,
      cd = clineLabels
    )
  }

  # Plot parametric composite clines
  introgress::composite.clines(
    parametric_clines,
    pdf = TRUE,
    out.file = parametric_compositePDF,
    labels = clineLabels
  )

  # Plot parametric clines.
  introgress::clines.plot(
    cline.data = parametric_clines,
    rplots = 3,
    cplots = 3,
    pdf = TRUE,
    out.file = parametric_singleClinesPDF,
    quantiles = TRUE,
    cd = clineLabels
  )

  # Write summary permutation clines
  write.table(
    perm_clines$Summary.data,
    file = summaryFile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  # Write summary parametric clines
  write.table(
    parametric_clines$Summary.data,
    file = parametric_summaryFile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  clineList <-
    list(
      permutation = perm_clines,
      parametric = parametric_clines,
      hi = hi.index,
      deltas = deltas,
      minDelt = minDelt
    )

  return(clineList)
}


#' Function to correlate INTROGRESS genomic clines with raster data
#' @param clineList List object returned from runIntrogress
#' @param rasterPointValues List of data.frames with raster values at each
#'                     sample locality. See ?extractPointValues
#' @param rastersToUse Integer vector indicating which raster layers to use for
#'                     genomic cline X environment plots. If NULL, uses all
#'                     rasters
#' @param sample.info File path to sample info file with 4 columns:
#'                    indID,popID,latitude,longitude
#' @param header Boolean specifying if sample.info has a header line
#' @param coords data.frame with coordinates. Can be taken from env.list,
#'               an object returned from the prepare_rasters() function.
#'               See ?prepare_rasters(). If used, value is env.list[[3]]
#' @param clineLabels Character vector of length == 3 designating population
#'                    names for c(P1, Het, P2)
#' @param outputDIR File path to directory for outputting plots
#' @param showPLOTS Boolean; Whether to print the plots to the screen
#' @param clineMethod Character string indicating desired method for generating
#'                    genomic clines. Must be either "permutation" or
#'                    "parametric". See ?introgress
#' @param prefix Character prefix for output files
#' @param cor.method Method for correlations c("pearson", "kendall", "spearman)
#' @export
#' @examples
#' clinesXenvironment(clineList = eatt,
#'                    rasterPointValues = rasterPoints,
#'                    rasterToUse = NULL,
#'                    clineLabels = c("EA", "Het", "TT"),
#'                    outputDIR = "./cline_plots",
#'                    clineMethod = "permutation",
#'                    prefix = "eatt",
#'                    cor.method = "auto")
clinesXenvironment <- function(clineList,
                               rasterPointValues,
                               rastersToUse = NULL,
                               clineLabels = c("P1, Het", "P2"),
                               outputDIR = "./plots",
                               showPLOTS = FALSE,
                               clineMethod = "permutation",
                               prefix = "population1",
                               cor.method = "auto"){

    clines <- clineList[[clineMethod]]
    hi <- clineList[["hi"]]

    # Load coordinate data
    coords <- data.frame(rasterPointValues[[1]]$lng,
                         rasterPointValues[[1]]$lat,
                         stringsAsFactors = FALSE)

    colnames(coords) <- c("lng", "lat")

  # Make sure coords are numeric.
  coords$lng <- as.numeric(as.character(coords$lng))
  coords$lat <- as.numeric(as.character(coords$lat))

  lat_norm <- range01(coords$lat)
  lon_norm <- range01(coords$lng)

  # Get +- 10% lower/upper bounds for lat_norm.
  l <- data.frame(col1 = ((lat_norm) - (lat_norm * 0.1)),
                  col2 = lat_norm,
                  col3 = ((lat_norm) + (lat_norm * 0.1)))

  # Change column names.
  colnames(l) <- c("lower", "h", "upper")

  #### Now do it for longitude. ####

  # Get +- 10% lower/upper bounds for lon_norm.
  lon <- data.frame(col1 = ((lon_norm) - (lon_norm * 0.1)),
                    col2 = lon_norm,
                    col3 = ((lon_norm) + (lon_norm * 0.1)))

  colnames(lon)<-c("lower", "h", "upper")

  # Cline X Environment Plots

  ## latitude ##

  # Create output directory if it doesn't already exist.
  dir.create(outputDIR, showWarnings = FALSE)

  pdf(
    file = file.path(outputDIR, paste0(prefix, "_clinesXLatLon.pdf")),
    width = 11,
    height = 7,
    onefile = TRUE
  )
  plot_lm_clines(
    l$h,
    clines$Fitted.AA,
    "Latitude",
    paste0("Genomic Clines (", clineLabels[1], " X ", clineLabels[3], ")")
  )

  plot_lm_hindex(
    l$h,
    hi$h,
    "Latitude",
    paste0("Hybrid Index (", clineLabels[1], " X ", clineLabels[3], ")")
  )

  plot_lm_clines(
    lon$h,
    clines$Fitted.AA,
    "Longitude",
    paste0("Genomic Clines (", clineLabels[1], " X ", clineLabels[3], ")")
  )

  plot_lm_hindex(
    lon$h,
    hi$h,
    "Longitude",
    paste0("Hybrid Index (", clineLabels[1], " X ", clineLabels[3], ")")
  )
  dev.off()

  if (isTRUE(showPLOTS)) {
    print(plot_lm_clines(
      l$h,
      clines$Fitted.AA,
      "Latitude",
      paste0("Genomic Clines (", clineLabels[1], " X ", clineLabels[3], ")")
    ))

    print(plot_lm_hindex(
      l$h,
      hi$h,
      "Latitude",
      paste0("Hybrid Index (", clineLabels[1], " X ", clineLabels[3], ")")
    ))

    print(plot_lm_clines(
      lon$h,
      clines$Fitted.AA,
      "Longitude",
      paste0("Genomic Clines (", clineLabels[1], " X ", clineLabels[3], ")")
    ))

    print(plot_lm_hindex(
      lon$h,
      hi$h,
      "Longitude",
      paste0("Hybrid Index (", clineLabels[1], " X ", clineLabels[3], ")")
    ))
  }

  #### Now do environments (e.g. temp, precip) on X-axis ####

  # If subset not specified, use all raster layers
  if (is.null(rastersToUse)){
    importantLayers <- rasterPointValues
    rastersToUse <- c(1:length(rasterPointValues))
  } else{
    # Load in most important environment data.frames from list of data.frames
    importantLayers <- rasterPointValues[rastersToUse]
  }

  # Do Hybrid Index ~ Lat/Lon correlations
  latCorrSummary <- hiXenvCorr(l$h, hi$h, cor.method, "Latitude")

  colnames(latCorrSummary) <- c("CorrelationStats", "Values", "EnvVar")

  lonCorrSummary <- hiXenvCorr(lon$h, hi$h, cor.method, "Longitude")

  # Write correlations to file.
  # latitude
  write.table(
    x = latCorrSummary,
    file = file.path(outputDIR, paste0(
      prefix,
      "_corrSummary.csv"
    )),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    append = FALSE
  )

  # longitude
  write.table(
    x = lonCorrSummary,
    file = file.path(outputDIR, paste0(
      prefix,
      "_corrSummary.csv")),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
  pdf(
    file = file.path(outputDIR, paste0(prefix, "_clinesXenv.pdf")),
    width = 11,
    height = 7,
    onefile = TRUE
  )
  for (i in 1:length(importantLayers)) {
    env <- importantLayers[[i]][,4]

    # Genomic Clines ~ Environment
    plot_lm_clines(
      env,
      clines$Fitted.AA,
      paste("Raster Layer", colnames(importantLayers[[i]][4]), sep = " "),
      paste0("Genomic Clines (",
             clineLabels[1],
             " X ",
             clineLabels[3], ")")
    )

    # HybridIndex ~ Environment
    plot_lm_hindex(
      env,
      hi$h,
      paste("Raster Layer", colnames(importantLayers[[i]][4]), sep = " "),
      paste0("Hybrid Index (", clineLabels[1], " X ", clineLabels[3], ")")
    )
  }
  dev.off()

  if (isTRUE(showPLOTS)) {
    for (i in 1:length(importantLayers)) {
      env <- importantLayers[[i]][,4]

      # Genomic Clines ~ Environment
      print(plot_lm_clines(
        env,
        clines$Fitted.AA,
        paste("Raster Layer", colnames(importantLayers[[i]][4]), sep = " "),
        paste0("Genomic Clines (",
               clineLabels[1],
               " X ",
               clineLabels[3], ")")
      ))

      # HybridIndex ~ Environment
      print(plot_lm_hindex(
        env,
        hi$h,
        paste("Raster Layer", colnames(importantLayers[[i]][4]), sep = " "),
        paste0("Hybrid Index (", clineLabels[1], " X ", clineLabels[3], ")")
      ))
    }
  }

  for (i in 1:length(importantLayers)) {
    env <- importantLayers[[i]][,4]

    # Do hybrid index ~ env correlations.
    corrSummary <- hiXenvCorr(env, hi$h, cor.method, colnames(importantLayers[[i]][4]))

    # Write correlations summaries to file.
    write.table(
      x = corrSummary,
      file = file.path(outputDIR, paste0(
        prefix,
        "_corrSummary.csv"
      )),
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

#' Function to get range of values normalized from 0 to 1
#' @param x Vector to normalize
#' @noRd
range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}

#' Function to plot smoothed linear models for each
#'  genomic cline ~ environment/lat/lon.
#' @param x Numeric X-axis environmental variable
#' @param y Numeric Y-axis genomic cline variable
#' @param xlab Character string for X-axis label
#' @param ylab Character string for Y-axis label
#' @noRd
plot_lm_clines <- function(x, y, xlab, ylab) {

  # First, get the genomic clines into list of data.frames.
  dflist <- list()
  for (i in 1:nrow(y)) {
    dflist[[i]] <- data.frame(y[i, ], x)
  }

  # Get the locus names.
  locinames <- data.frame(attributes(y)[[2]][[1]], stringsAsFactors = FALSE)
  colnames(locinames) <- c("loci")

  #return(locinames.vec)

  # Set column names. y=genomic clines, x = environment.
  for (i in 1:nrow(y)){
    colnames(dflist[[i]]) <- c("y", "x")
  }

  # Reorder by "x" column (environment).
  dflist <- lapply(dflist, function(i) i[order(i[[2]]),])


  # Plot regressions using linear model and 3rd order polynomial.
  lm <- ggplot2::ggplot(dplyr::bind_rows(dflist, .id = "df"),
                    ggplot2::aes(x, y, colour = df)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
    ggplot2::xlab(label = xlab) + ggplot2::ylab(label = ylab) +
    ggplot2::theme_bw() + ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text = ggplot2::element_text(size = 18),
      text = ggplot2::element_text(size = 18),
      axis.ticks = ggplot2::element_line(size = ggplot2::rel(1.5))
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
    ggplot2::scale_color_hue(name = "Loci", labels = locinames$loci)
  print(lm)
}

#' Function to do linear model for and plot Hybrid Index ~ environment.
#' Only slightly different from plot_lm_clines.
#' @param x Numeric X-axis variable
#' @param y Numeric Y-axis variable
#' @param xlab Character string for X-axis plot label
#' @param ylab Character string for Y-axis plot label
#' @noRd
plot_lm_hindex <- function(x, y, xlab, ylab){

  df <- data.frame(x = x, y = y)

  # Plot regressions using linear model and 3rd order polynomial.
  lm <-
    ggplot2::ggplot(data = df, ggplot2::aes(x, y)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
    ggplot2::xlab(label = xlab) + ggplot2::ylab(label = ylab) +
    ggplot2::theme_bw() + ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text = ggplot2::element_text(size = 18),
      text = ggplot2::element_text(size = 18),
      axis.ticks = ggplot2::element_line(size = ggplot2::rel(1.5))
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), oob = scales::squish)
  print(lm)
}

#' Function to subset individuals from the extractRasterPoints data.frame
#' @param rasterPoint.df data.frame from extractRasterPoints list
#' @param inds File path to one-column input file containing subset of
#'             individuals you want to use for clinesXenvironment()
#' @return data.frame with subset of individuals
#' @export
subsetIndividuals <- function(rasterPoint.df,
                              inds){

  popInds <- read.table(inds, header = FALSE, stringsAsFactors = FALSE)
  rasterPoint.df <- subset(rasterPoint.df,
                           subset = inds %in% popInds$V1)
  return(rasterPoint.df)
}

#' Function to perform hybrid index ~ environment correlations
#' @param env Numeric vector containing environmental values on X-axis
#' @param hi Numeric vector containing hybrid indices
#' @param cor.method Character string indicating correlation method
#'                   c("pearson", "kendall", "spearman")
#' @param envName Character string indicating name of environmental variable
#' @return data.frame with correlation output
#' @noRd
hiXenvCorr <- function(env, hi, cor.method, envName){

  # Checks for normality. If P < 0.05, variable is not normally distributed
  # So it uses kendall instead of pearson
  if (cor.method == "auto"){
    sw.env <- shapiro.test(env)
    sw.hi <- shapiro.test(hi)

    if (sw.env$p.value <= 0.05 | sw.hi$p.value <= 0.05){
      cor.method <- "kendall"
      est <- "tau"
    } else{
      cor.method <- "pearson"
      est <- "R2"
    }
  } else if (cor.method == "pearson"){
    est <- "R2"
  } else if (cor.method == "spearman"){
    est <- "rho"
  } else if (cor.method == "kendall"){
    est <- "tau"
  } else{
    stop("cor.method must be one of c('auto','pearson','kendall','spearman')")
  }

  # Correlations.
  corrCoeff <- cor.test(env, hi, method = cor.method)

  # Make info into data.frame.
  corrSummary <-
    data.frame(
      "CorrelationStats" = c("Method",
        "TestStatistic",
        "Pvalue",
        est,
        "Hypothesis"),
     "Values" = c(
        corrCoeff$method,
        corrCoeff$statistic,
        corrCoeff$p.value,
        corrCoeff$estimate,
        corrCoeff$alternative),
      "Environment" = c(rep(envName, 5))
    )

  return(corrSummary)
}

