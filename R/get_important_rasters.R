#' Script to extract raster data for all sampling localities
#' and find the most important features.
#'
#' This script takes as input a directory of rasters, crops them to the
#' sampling extent, finds the raster values at each sample locality,
#' and uses MAXENT and ENMeval to determine the most important raster
#' layers (i.e. features).
#'
#'
#' Function to load in and crop the raster layers to the sampling extent
#'
#' This function takes as input a directory of rasters. Only the desired
#' rasters should be included in raster.dir. You also will need a
#' comma-delimited sample.file that has four columns in a specific order:
#' (sampleIDs,populationIDs,latitude,longitude). You can choose if a header
#' line is present. The rasters will all be put into a stack and cropped to
#' the extent of the sample localities + bb.buffer.
#' @param raster.dir Directory of rasters to load and crop
#' @param sample.file CSV file with sample information (sampleID,popID,lat,lon)
#' @param header Boolean; Does sample.file have a header line?
#' @param bb.buffer Integer; Buffer around sample bounding box.
#'                  bb.buffer = bb.buffer * resolution (in arc-seconds)
#' @param plotDIR Directory to save plots to
#' @param showPLOTS Boolean; Whether to print plots to screen
#' @return List with cropped rasters and other info
#' @export
#' @examples
#' envList <- prepare_rasters(raster.dir = "uncroppedLayers",
#'                            sample.file = file.path("exampleData",
#'                                          "ENMeval_bioclim",
#'                                          "localityInfo",
#'                            "sample_localities_maxent_southeast_noNA.csv"),
#'                            header = TRUE,
#'                            bb.buffer = 10,
#'                            plotDIR = "./plots",
#'                            showPLOTS = TRUE)
prepare_rasters <- function(raster.dir,
                            sample.file,
                            header = TRUE,
                            bb.buffer = 10,
                            plotDIR = "./plots",
                            showPLOTS = FALSE){

  if (!requireNamespace("raster", quietly = TRUE)){
    warning("The raster package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("sp", quietly = TRUE)){
    warning("The sp package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("dismo", quietly = TRUE)){
    warning("The dismo package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("rJava", quietly = TRUE)){
    warning("The rJava package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("ENMeval", quietly = TRUE)){
    warning("The ENMeval package must be installed to use this functionality")
    return(NULL)
  }

  dir.create(plotDIR, showWarnings = FALSE)

  samples <-
    read.csv(file = sample.file,
             header = header,
             stringsAsFactors = FALSE)

  # Make sure values are numeric.
  samples[,3] <- as.numeric(as.character(samples[,3]))
  samples[,4] <- as.numeric(as.character(samples[,4]))

  # Get dataframe of longitude, latitude.
  coords <- data.frame(samples[,4], samples[,3])

  # Make into factor for bg color assignment later.
  pops <- base::as.factor(samples[,2])

  # Change column names.
  colnames(coords) <- c("lng", "lat")

  # Get raster filenames (all in one directory)
  files <- list.files(path = raster.dir, full.names = T)

  writeLines("\n\nOrder of raster files: \n")

  # Order alphanumerically.
  files <- gtools::mixedsort(files)

  for (i in 1:length(files)){
    writeLines(paste0(i, ": ", files[i]))
  }

  # Put the rasters into a RasterStack.
  envs <- raster::stack(files)

  writeLines(paste0("\n\nLoaded ", raster::nlayers(envs),
                    " raster layers.."))

  # Plot first raster in the stack, bio1.
  #raster::plot(envs[[1]], main=names(envs)[1])

  # Get CRS (coordinate reference system)
  mycrs <- raster::crs(envs[[1]])

  # Create spatialpoints object.
  p <- sp::SpatialPoints(coords = coords, proj4string=mycrs)

  # Get the bounding box of the points
  bb <- raster::bbox(p)

  # Add (bb.buffer * resolution / 2) to sample bounds for background extent.
  bb.buf <- raster::extent(bb[1]-bb.buffer,
                           bb[3]+bb.buffer,
                           bb[2]-bb.buffer,
                           bb[4]+bb.buffer)

  # Add bb.buffer * resolution (in arc-seconds) to bb for cropping raster layer.
  envs.buf <- raster::extent(bb[1]-(bb.buffer/2),
                             bb[3]+(bb.buffer/2),
                             bb[2]-(bb.buffer/2),
                             bb[4]+(bb.buffer/2))

  writeLines(paste0("\n\nCropping to samples with buffer of ",
                    (bb.buffer/2),
                    " degrees\n"))

  # Crop raster extent to sample bounding box + bb.buffer
  envs.cropped <- raster::crop(envs, envs.buf)

  writeLines(paste0("\nCropping background layers with buffer of ",
                    bb.buffer,
                    " degrees\n"))

  # Crop environmental layers to tmatch the study extent.
  # Used for background data later.
  envs.backg <- raster::crop(envs, bb.buf)

  counter <- 1
  # Save raster plots with sample points layered on top
  pdf(file = file.path(plotDIR, "croppedRasterPlots.pdf"),
      width = 7,
      height = 7,
      onefile = T)
    for (i in 1:raster::nlayers(envs.backg)){
      writeLines(paste0("Saving raster plot ", i, " to disk..."))
      raster::plot(envs.backg[[i]])
      dismo::points(coords, pch=21, bg=pops)
      counter <- counter+1
    }
  dev.off()

  if (isTRUE(showPLOTS)) {
    for (i in 1:raster::nlayers(envs.backg)){
      writeLines(paste0("Saving raster plot ", i, " to disk..."))
      raster::plot(envs.backg[[i]])
      dismo::points(coords, pch=21, bg=pops)
      counter <- counter+1
    }
  }

  envList <- list(envs.cropped, envs.backg, coords, p, pops, samples[,1])

  rm(envs, p, bb, bb.buf, envs.buf, envs.cropped, envs.backg)
  gc(verbose = FALSE)

  return(envList)
}

#' Function to make partitions from background and foreground rasters.
#'
#' This function uses the output from prepare_raster() and makes
#' background partitions using several methods.
#' @param env.list Object output from prepare_rasters() function
#' @param categoricals Character vector of categorical layer names
#' @param number.bg.points Number of background points to generate
#' @param bg.raster Raster to use for background points
#' @param agg.factor A vector of 1 or 2 numbers for the checkerboard
#'                           methods
#' @param plotDIR Directory to save the plots to
#' @param showPLOTS Boolean; Whether to print plots to screen
#' @return Object with background points
#' @export
#' @examples
#' bg <- partition_raster_bg(env.list = envList,
#'                           number.bg.points = 10000,
#'                           bg.raster = 1,
#'                           plotDIR = "./plots",
#'                           showPLOTS = TRUE,
#'                           agg.factor = 2)
partition_raster_bg <- function(env.list,
                                categoricals = NULL,
                                number.bg.points = 10000,
                                bg.raster = 1,
                                plotDIR = "./plots",
                                showPLOTS = FALSE,
                                agg.factor = 2){

  if (!requireNamespace("raster", quietly = TRUE)){
    warning("The raster package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("sp", quietly = TRUE)){
    warning("The sp package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("dismo", quietly = TRUE)){
    warning("The dismo package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("rJava", quietly = TRUE)){
    warning("The rJava package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("ENMeval", quietly = TRUE)){
    warning("The ENMeval package must be installed to use this functionality")
    return(NULL)
  }

  dir.create(plotDIR, showWarnings = FALSE)

  bg.raster <- as.integer(as.character(bg.raster))

  envs.cropped <- env.list[[1]]
  envs.backg <- env.list[[2]]
  coords <- env.list[[3]]
  p <- env.list[[4]]
  pops <- env.list[[5]]
  inds <- env.list[[6]]

  rm(env.list)
  invisible(gc(verbose = FALSE))

  writeLines(paste0("\n\nGenerating ", number.bg.points, " random points..."))

  # Randomly sample number.bg.points background points from
  # one background extent raster
  # (only one per cell without replacement).
  # Note: If the raster has number.bg.points pixels,
  # you'll get a warning and all pixels will be used for background.
  bg <- dismo::randomPoints(envs.backg[[bg.raster]], n=number.bg.points)
  bg <- raster::as.data.frame(bg)

  # Change column names.
  colnames(coords) <- c("lng", "lat")
  colnames(bg) <- colnames(coords)

  rm(pops)
  invisible(gc())

  writeLines("\nSaving background partition plots...")
  writeLines("\nPerforming block method...")
  ### Block Method
  blocks <- ENMeval::get.block(coords, bg)
  plot_partitions(envs.fg = envs.cropped,
                  envs.bg = envs.backg,
                  bg = bg,
                  coords = coords,
                  partition = blocks,
                  partition.type = "Block",
                  categoricals = categoricals,
                  plotDIR = plotDIR,
                  showPLOTS = showPLOTS)
  # pdf(file = file.path(plotDIR, "bgPartitions_blocks.pdf"),
  #     width = 7,
  #     height = 7,
  #     onefile = TRUE)

    # p1 <- ENMeval::evalplot.grps(pts = coords,
    #                        pts.grp = blocks$occs.grp,
    #                        envs = envs.backg) +
    #   ggplot2::ggtitle("Spatial Block Partitions: Occurrences")


  rm(blocks)
  invisible(gc(verbose = FALSE))

  writeLines("Performing checkerboard1 method...")
  ### Checkerboard1 Method
  check1 <- ENMeval::get.checkerboard1(coords,
                                       envs.backg,
                                       bg,
                                       aggregation.factor = agg.factor)

  plot_partitions(envs.fg = envs.cropped,
                  envs.bg = envs.backg,
                  bg = bg,
                  coords = coords,
                  partition = check1,
                  partition.type = "Checkerboard1",
                  categoricals = categoricals,
                  plotDIR = plotDIR,
                  showPLOTS = showPLOTS)
    # pdf(file = file.path(plotDIR, "bgPartitions_checkerboard1.pdf"),
    #     width = 7,
    #     height = 7,
    #     onefile = TRUE)
      # raster::plot(
      #   envs.backg[[bg.raster]],
      #   col="gray",
      #   main = "Partitions into Training and Test Data (Aggregation Factor=5)",
      #   sub="Checkerboard1 Method",
      #   xlab = "Longitude",
      #   ylab = "Latitude"
      #   )
    #   dismo::points(bg, pch=21, bg=check1$occ.grp)
    # dev.off()

    # if (isTRUE(showPLOTS)) {
    #   raster::plot(
    #     envs.backg[[bg.raster]],
    #     col="gray",
    #     main = "Partitions into Training and Test Data (Aggregation Factor=5)",
    #     sub="Checkerboard1 Method",
    #     xlab = "Longitude",
    #     ylab = "Latitude"
    #   )
    #   dismo::points(bg, pch=21, bg=check1$occ.grp)
    # }

    rm(check1)
    invisible(gc(verbose = FALSE))

    writeLines("Performing checkerboard2 method...")
    ### Using the Checkerboard2 method.
    check2 <-
      ENMeval::get.checkerboard2(coords,
                                 envs.backg,
                                 bg,
                                 aggregation.factor = c(agg.factor,
                                                      agg.factor))

    plot_partitions(envs.fg = envs.cropped,
                    envs.bg = envs.backg,
                    bg = bg,
                    coords = coords,
                    partition = check2,
                    partition.type = "Checkerboard2",
                    categoricals = categoricals,
                    plotDIR = plotDIR,
                    showPLOTS = showPLOTS)
  #   pdf(file = file.path(plotDIR, "bgPartitions_checkerboard2.pdf"),
  #       width = 7,
  #       height = 7,
  #       onefile = TRUE)
  #   raster::plot(envs.backg[[bg.raster]], col="gray",
  #        main = "Partitions into Training and Test Data",
  #        sub="Checkerboard2 Method",
  #        xlab = "Longitude",
  #        ylab = "Latitude")
  #   dismo::points(bg, pch=21, bg=check2$bg.grp)
  #   dismo::points(coords, pch=21,
  #                 bg=check2$occ.grp,
  #                 col="white",
  #                 cex=1.5)
  # dev.off()

  # if (isTRUE(showPLOTS)) {
  #   raster::plot(envs.backg[[bg.raster]], col="gray",
  #                main = "Partitions into Training and Test Data",
  #                sub="Checkerboard2 Method",
  #                xlab = "Longitude",
  #                ylab = "Latitude")
  #   dismo::points(bg, pch=21, bg=check2$bg.grp)
  #   dismo::points(coords, pch=21,
  #                 bg=check2$occ.grp,
  #                 col="white",
  #                 cex=1.5)
  # }

  rm(check2)
  invisible(gc())

  writeLines("\nDone!")

  return(bg)
}

#' Function to plot partitions and environmental differences and similarities
#' @param envs.fg Raster stack object (envList[[1]])
#' @param envs.bg Background points on raster object (envList[[2]])
#' @param bg Background object
#' @param coords DataFrame of coordinates (occs) (envList[[3]])
#' @param partitions Partition object
#' @param partition.type string; Type of partition being used
#' @param categoricals Character vector; names of categorical layers
#' @param plotDIR Directory to save plots
#' @param showPLOTS Whether to print plots to the screen.
#' @noRd
plot_partitions <- function(envs.fg,
                            envs.bg,
                            bg,
                            coords,
                            partitions,
                            partition.type,
                            categoricals = NULL,
                            plotDIR = "./",
                            showPLOTS = FALSE){

  # Plot partitioning of points.
  p1 <- ENMeval::evalplot.grps(pts = coords,
                               pts.grp = partitions$occs.grp,
                               envs = envs.bg) +
    ggplot2::ggtitle(paste0("Spatial ", partition.type, " Partitions: Occurrences"))

  #bg <- bg[[which(!is.na(envs.bg[[1]]@data@values))]]

  p2 <- ENMeval::evalplot.grps(pts = bg,
                               pts.grp = partitions$bg.grp,
                               envs = envs.bg) +
    ggplot2::ggtitle(paste0("Spatial ", partition.type, " Partitions: Background"))

  occs.z <- cbind(coords, raster::extract(envs.fg, coords))
  bg.z <- cbind(bg, raster::extract(envs.fg, bg))

  p3 <- ENMeval::evalplot.envSim.hist(sim.type = "mess",
                                      ref.data = "occs",
                                      occs.z = occs.z,
                                      bg.z = bg.z,
                                      occs.grp = partitions$occs.grp,
                                      bg.grp = partitions$bg.grp,
                                      categoricals = categoricals)

  p4 <- ENMeval::evalplot.envSim.hist(sim.type = "most_diff",
                                      ref.data = "occs",
                                      occs.z = occs.z,
                                      bg.z = bg.z,
                                      occs.grp = partitions$occs.grp,
                                      bg.grp = partitions$bg.grp,
                                      categoricals = categoricals)

  p5 <- ENMeval::evalplot.envSim.hist(sim.type = "most_sim",
                                      ref.data = "occs",
                                      occs.z = occs.z,
                                      bg.z = bg.z,
                                      occs.grp = partitions$occs.grp,
                                      bg.grp = partitions$bg.grp,
                                      categoricals = categoricals)


  ggplot2::ggsave(filename = file.path(plotDIR, paste0("bgPartitions_",
                                    partition.type,
                                    "_occurrences.pdf")),
                  plot = p1,
                  device = "pdf",
                  width = 7,
                  height = 7,
                  units = "in")

  ggplot2::ggsave(filename = file.path(plotDIR, paste0("bgPartitions_",
                                                       partition.type,
                                                       "_background.pdf")),
                  plot = p2,
                  device = "pdf",
                  width = 7,
                  height = 7,
                  units = "in")

  ggplot2::ggsave(filename = file.path(plotDIR, paste0("bgPartitions_",
                                    partition.type,
                                    "_multivar_envSimilarity.pdf")),
                  plot = p3,
                  device = "pdf",
                  width = 7,
                  height = 10,
                  units = "in")

  ggplot2::ggsave(filename = file.path(plotDIR, paste0("bgPartitions_",
                                    partition.type,
                                    "_most_different.pdf")),
                  plot = p4,
                  device = "pdf",
                  width = 7,
                  height = 10,
                  units = "in")

  ggplot2::ggsave(filename = file.path(plotDIR, paste0("bgPartitions_",
                                    partition.type,
                                    "_most_similar.pdf")),
                  plot = p5,
                  device = "pdf",
                  width = 7,
                  height = 10,
                  units = "in")

  # raster::plot(envs.backg[[bg.raster]], col="gray",
  #      main = "Partitions into Training and Test Data",
  #      sub="Block Method",
  #      xlab = "Longitude",
  #      ylab = "Latitude")
  # dismo::points(coords, pch=21, bg=blocks$occ.grp)
  # dev.off()

  if (isTRUE(showPLOTS)) {
    # Plot partitioning of points.
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
  }
}

#' Function to run ENMeval and MAXENT on the raster layers.
#'
#' This function takes as input the object output from preparte_rasters() and
#' the background partitioning method you would like to use
#' (see ENMeval documentation). You can visualize how the partitioning
#' methods will look by viewing the PDF output from partition_raster_bg.
#' See ?partition_raster_bg for more info.
#' @param env.list Object output from prepare_rasters() function
#' @param bg Object with background points; generated from partition_raster_bg.
#'           If NULL, the background layer will be calculated for you.
#' @param partition.method Method used for background point partitioning
#' @param parallel If TRUE, ENMeval is run parallelized with np CPU cores
#' @param np Number of parallel cores to use if parallel = TRUE
#' @param RMvalues Vector of non-negative regularization multiplier values.
#'                 Higher values impose a stronger penalty on model complexity
#' @param feature.classes Character vector of feature classes to be used
#' @param categoricals Vector indicating which (if any) of the input
#'                     environmental layers are categorical.
#' @param algorithm Character vector. Defaults to dismo's maxent.jar
#' @return ENMeval object
#' @export
#' @examples
#' eval.par <- runENMeval(env.list = envList,
#'                        bg = bg,
#'                        partition.method = "checkerboard1",
#'                        parallel = FALSE,
#'                        RMvalues = seq(0.5, 4, 0.5),
#'                        feature.classes = c("L",
#'                                            "LQ",
#'                                            "H",
#'                                            "LQH",
#'                                            "LQHP",
#'                                            "LQHPT")
#'                        categoricals = NULL,
#'                        algorithm = "maxent.jar")
runENMeval <- function(env.list,
                       bg,
                       partition.method,
                       parallel = FALSE,
                       np = 2,
                       RMvalues = seq(0.5, 4, 0.5),
                       feature.classes = c("L",
                                           "LQ",
                                           "H",
                                           "LQH",
                                           "LQHP",
                                           "LQHPT"),
                       categoricals = NULL,
                       algorithm = "maxent.jar"
                       ){

  if (!requireNamespace("raster", quietly = TRUE)){
    warning("The raster package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("sp", quietly = TRUE)){
    warning("The sp package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("dismo", quietly = TRUE)){
    warning("The dismo package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("rJava", quietly = TRUE)){
    warning("The rJava package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("ENMeval", quietly = TRUE)){
    warning("The ENMeval package must be installed to use this functionality")
    return(NULL)
  }

  if (isTRUE(parallel)){
    numCores <- np
  } else if (isFALSE(parallel)) {
    numCores = 1
  } else {
    stop("parallel argument must be either TRUE or FALSE")
  }

  envs.fg <- env.list[[1]]
  coords <- env.list[[3]]

  rm(env.list)
  invisible(gc())

  # NOTE: Had to make some changes when ENVeval switched to version 2.0.0.
  # Some arguments were deprecated, others were changed.
  # occ changed to occs
  # env changed to envs
  # bg.coords changed to bg
  # NOTE: If bg is NULL, it calculates the random points for you.
  # feature.classes and RMvalues got combined into a named list.
  # methods became partitions
  # categoricals no longer supports integer to choose categorical layers;
  # Now have to specify a character vector of the layer name.
  # aggregation.factor was deprecated.
  eval.par <-
    ENMeval::ENMevaluate(
      occs = coords,
      envs = envs.fg,
      bg = bg,
      tune.args = list(fc = feature.classes, rm = RMvalues),
      partitions = partition.method,
      parallel = parallel,
      algorithm = algorithm,
      categoricals = categoricals,
      numCores = np
    )

  return(eval.par)
}

#' Function to summarize ENMeval output
#' @param eval.par Object returned from runENMeval
#' @param minLat Minimum latitude to plot for predictions
#' @param maxLat Maximum latitude to plot for predictions
#' @param minLon Minimum longitude to plot for predictions
#' @param maxLon Maximum longitude to plot for predictions
#' @param examine.predictions Character vector of feature classes to examine
#'                            how complexity affects predictions
#' @param examine.stats Character vector of stats to evaluate
#' @param RMvalues Vector of non-negative RM values to examine how
#'                 complexity affects predictions
#' @param nullmodel.iter Integer; Number of iterations to perform with null
#'                       model.
#' @param plotDIR Directory to save plots to
#' @param showPLOTS Boolean; Whether to print plots to screen
#' @param niche.overlap Boolean. If TRUE, calculates pairwise niche overlap
#'                      matrix
#' @param plot.width Integer. Specify plot widths
#' @param plot.height Integer. Specify plot heights
#' @param imp.margins Integer vector. Margins of variable importance barplot.
#'                    c(bottom, left, top, right)
#' @export
#' @examples
#' summarize_ENMeval(eval.par = eval.par,
#'                  minLat = 20,
#'                  maxLat = 50,
#'                  minLon = -110,
#'                  maxLon = -40,
#'                  examine.predictions = c("L",
#'                                          "LQ",
#'                                          "H",
#'                                          "LQH",
#'                                          "LQHP",
#'                                          "LQHPT"),
#'                  RMvalues = seq(0.5, 4, 0.5),
#'                  plotDIR = "./plots",
#'                  showPLOTS = TRUE,
#'                  niche.overlap = FALSE,
#'                  plot.width = 7,
#'                  plot.height = 7,
#'                  imp.margins = c(10.0, 4.1, 4.1, 2.1))
summarize_ENMeval <- function(eval.par,
                              minLat = 20,
                              maxLat = 50,
                              minLon = -110,
                              maxLon = -40,
                              examine.predictions = c("L",
                                                      "LQ",
                                                      "H",
                                                      "LQH",
                                                      "LQHP",
                                                      "LQHPT"),
                              examine.stats = c("auc.val",
                                                "auc.diff",
                                                "cbi.val",
                                                "or.mtp",
                                                "or.10p"),
                              RMvalues = seq(0.5, 4, 0.5),
                              nullmodel.iter = 100,
                              plotDIR = "./plots",
                              showPLOTS = FALSE,
                              niche.overlap = FALSE,
                              plot.width = 7,
                              plot.height = 7,
                              imp.margins = c(10.0, 4.1, 4.1, 2.1)){

  dir.create(plotDIR, showWarnings = FALSE)

  `%>%` <- dplyr::`%>%`

  if (!requireNamespace("raster", quietly = TRUE)){
    warning("The raster package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("sp", quietly = TRUE)){
    warning("The sp package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("dismo", quietly = TRUE)){
    warning("The dismo package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("rJava", quietly = TRUE)){
    warning("The rJava package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("ENMeval", quietly = TRUE)){
    warning("The ENMeval package must be installed to use this functionality")
    return(NULL)
  }

  if (niche.overlap){
    # Calculate niche overlap.
    overlap <- ENMeval::calc.niche.overlap(eval.par@predictions, stat="D")

    # Write niche overlap to a file.
    write.csv(x = overlap, file = file.path(plotDIR,
                                            "ENMeval_nicheOverlap.csv"))
  }

  # Write ENMeval results to file.
  write.csv(x = eval.par@results, file = file.path(plotDIR,
                                                   "ENMeval_results.csv"))

  # Get best model.
  bestAICc <- eval.par@results[which(eval.par@results$delta.AICc==0),]

  # Write to file.
  write.csv(bestAICc, file = file.path(plotDIR, "ENMeval_bestAICc.csv"))

  # Save it to file.
  pdf(file = file.path(plotDIR,
                       "ENMeval_relativeOccurenceRate_bestModel.pdf"),
      width = plot.width,
      height = plot.height,
      onefile = TRUE)

    raster::plot(eval.par@predictions[[which(eval.par@results$delta.AICc==0)]],
               main="Relative Occurrence Rate")
  dev.off()

  if (isTRUE(showPLOTS)) {
    raster::plot(eval.par@predictions[[which(eval.par@results$delta.AICc==0)]],
                 main="Relative Occurrence Rate")
  }

  # Look at the model object for our "AICc optimal" model:
  aic.opt <- eval.par@models[[which(eval.par@results$delta.AICc==0)]]

  # Write it to file.
  write.csv(x = aic.opt@results,
            file = file.path(plotDIR,
                             "ENMeval_maxentModel_aicOptimal.csv"))

  # Get a data.frame of two variable importance metrics,
  # percent contribution and permutation importance:
  # For the model with AICc == 0.
  # NOTE: var.importance is a deprecated function
  # It is now saved as an object inside the results list.
  # This was a change with ENMeval version 2.0.0
  varImportance <-
    eval.par@variable.importance[[which(eval.par@results$delta.AICc==0)]]

  # Write them to file.
  write.csv(varImportance, file = file.path(plotDIR,
                                            "ENMeval_varImportance.csv"))

  # The "lambdas" slot shows which variables were included in the model.
  # If the coefficient is 0, that variable was not included in the model.
  lambdas <- aic.opt@lambdas

  write.csv(lambdas, file = file.path(plotDIR, "ENMeval_lambdas.csv"))

  pdf(file = file.path(plotDIR, "ENMeval_AUC_importancePlots.pdf"),
      width = plot.width,
      height = plot.height,
      onefile = TRUE)
    p1 <- ENMeval::evalplot.stats(e = eval.par,
                            stats = examine.stats,
                            color = "fc", x.var = "rm")

    p2 <- ENMeval::evalplot.stats(e = eval.par,
                                  stats = examine.stats,
                                  color = "rm", x.var = "fc")

    print(p1)
    print(p2)

    # Plot permutation importance.
    df <- varImportance
    par(mar=imp.margins)
    p3 <- raster::barplot(df$permutation.importance,
                    names.arg = df$variable,
                    las = 2,
                    ylab = "Permutation Importance")

    print(p3)

  dev.off()

  if (isTRUE(showPLOTS)) {
    print(p1)
    print(p2)
    print(p3)
  }

  # Let's see how model complexity changes the predictions in our example
  pdf(file = file.path(plotDIR, "model_Eval_and_Predictions.pdf"),
      width = plot.width,
      height = plot.height,
      onefile = TRUE)

    res <- ENMeval::eval.results(eval.par)

    opt.aicc <- res %>% dplyr::filter(delta.AICc == 0)

    opt.seq <- res %>%
      dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
      dplyr::filter(auc.val.avg == max(auc.val.avg))

    mod.seq <- ENMeval::eval.models(eval.par)[[opt.seq$tune.args]]
    raster::plot(mod.seq)

    # View the response curves
    #mod.seq <- dismo::response(ENMeval::eval.models(eval.par)[[opt.seq$tune.args]])
    #raster::plot(mod.seq, main="MAXENT Response Curves")
    #dev.off()

    # Assess predictions for optimal model
    pred.seq <- ENMeval::eval.predictions(eval.par)[[opt.seq$tune.args]]

    # Plot predictions for optimal model
    raster::plot(pred.seq, main=paste0("Optimal Model Prediction (",
                                       opt.seq$tune.args,
                                       ")"))

    # Plot binned background points
    points(ENMeval::eval.bg(eval.par),
           pch = 21,
           bg = ENMeval::eval.bg.grp(eval.par))

    # Plot occurence points on top
    points(ENMeval::eval.occs(eval.par),
           pch = 21,
           bg = ENMeval::eval.occs.grp(eval.par))

  # Plot predictions for all models
    for (i in 1:length(examine.predictions)){
      for (j in 1:length(RMvalues)){
        raster::plot(ENMeval::eval.predictions(eval.par)[[paste0(
          "fc.",
          examine.predictions[i],
          "_rm.",
          RMvalues[j])]],
                     ylim=c(minLat,maxLat),
                     xlim=c(minLon, maxLon),
                     main=paste0("fc.",
                                  examine.predictions[i],
                                  "_rm.",
                                  RMvalues[j],
                                 " Prediction"),
        )
      }
    }

    ### Compare against null models for significance
    mod.null <- ENMeval::ENMnulls(e = eval.par,
                                  mod.settings = list(fc = "LQ",
                                                      rm = tail(RMvalues, n=1)),
                                  no.iter = nullmodel.iter,
                                  eval.stats = examine.stats)

    pnull.hist <- ENMeval::evalplot.nulls(mod.null,
                                     stats = examine.stats,
                                     plot.type = "histogram")

    # Violin plot
    pnull.violin <- ENMeval::evalplot.nulls(mod.null,
                                            stats = examine.stats,
                                            plot.type = "violin")

    print(pnull.hist)
    print(pnull.violin)

  dev.off()

  if (isTRUE(showPLOTS)){
    raster::plot(mod.seq, main="MAXENT Response Curves")
    #dismo::response(ENMeval::eval.models(eval.par)[[opt.seq$tune.args]])

    # Plot predictions
    # Plot predictions for optimal model
    raster::plot(pred.seq, main=paste0("Optimal Model Prediction (",
                                       opt.seq$tune.args,
                                       ")"))

    # Plot binned background points
    points(ENMeval::eval.bg(eval.par),
           pch = 21,
           bg = ENMeval::eval.bg.grp(eval.par))

    # Plot occurence points on top
    points(ENMeval::eval.occs(eval.par),
           pch = 21,
           bg = ENMeval::eval.occs.grp(eval.par))

    # Plot predictions for all models
    for (i in 1:length(examine.predictions)){
      for (j in 1:length(RMvalues)){
        raster::plot(ENMeval::eval.predictions(eval.par)[[paste0(
          "fc.",
          examine.predictions[i],
          "_rm.",
          RMvalues[j])]],
          ylim=c(minLat,maxLat),
          xlim=c(minLon, maxLon),
          main=paste0("fc.",
                      examine.predictions[i],
                      "_rm.",
                      RMvalues[j],
                      " Prediction"),
        )
      }
    }

    print(pnull.hist)
    print(pnull.violin)

  }
}

#' Function to extract raster values at each sample point for raster stack
#' @param env.list List object generated from prepare_rasters()
#' @return List of data.frames with raster values at each sample locality
#'         for each raster layer, and list of raster names
#' @export
#' @examples
#' rasterPoints <- extractPointValues(envList)
extractPointValues <- function(env.list){

  if (!requireNamespace("raster", quietly = TRUE)){
    warning("The raster package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("sp", quietly = TRUE)){
    warning("The sp package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("dismo", quietly = TRUE)){
    warning("The dismo package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("rJava", quietly = TRUE)){
    warning("The rJava package must be installed to use this functionality")
    return(NULL)
  }

  if(!requireNamespace("ENMeval", quietly = TRUE)){
    warning("The ENMeval package must be installed to use this functionality")
    return(NULL)
  }


  envs.cropped <- env.list[[1]]
  p <- env.list[[4]]
  inds <- env.list[[6]]

  raster.points <- list()

  writeLines("\nExtracting raster values at sample points")

  # Extract raster values at each sample point.
  for (i in 1:raster::nlayers(envs.cropped)){
    raster.points[[i]] <- raster::extract(envs.cropped[[i]], p)
  }

  writeLines("\nAdding sampleIDs to raster dataframe")

  # Add sample IDs to point/raster value dataframe.
  rp.df <- lapply(raster.points,
                  function(x) {
                    cbind.data.frame(
                      inds,
                      raster::coordinates(p),
                      x
                    )
                  }
  )

  rasterNames <- env.list[[1]]@data@names

  # Use raster name as x colname
  for (i in 1:length(rasterNames)){
    colnames(rp.df[[i]]) <- c("inds", "lng", "lat", rasterNames[i])
  }

  rm(raster.points)
  gc(verbose = FALSE)

  # Get which samples are on NA raster values.
  test4na <- lapply(rp.df, function(x) x[raster::rowSums(is.na(x)) > 0,])

  # If some samples are on NA raster values: Print them and stop with error.
  allEmpty <- lapply(test4na, function(x) nrow(x) == 0)
  allEmpty <- unlist(allEmpty)

  if(!all(allEmpty)){
    writeLines("\n\nThe following samples are located on NA raster values:\n")
    print(test4na)
    stop("The above samples have NA raster values. Please remove the them.")
  }

  rm(test4na)
  gc(verbose = FALSE)

  return(rp.df)
}
