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
#' @return List with cropped rasters and other info
#' @export
prepare_rasters <- function(raster.dir,
                            sample.file,
                            header = TRUE,
                            bb.buffer = 10,
                            plotDIR = "./plots"){

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
#' @param number.bg.points Number of background points to generate
#' @param bg.raster Raster to use for background points
#' @param agg.factor A vector of 1 or 2 numbers for the checkerboard
#'                           methods
#' @param plotDIR Directory to save the plots to
#' @return Object with background points
#' @export
partition_raster_bg <- function(env.list,
                                number.bg.points = 10000,
                                bg.raster = 1,
                                plotDIR = "./plots",
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
  gc(verbose = FALSE)

  writeLines(paste0("\n\nGenerating ", number.bg.points, " random points..."))

  # Randomly sample 10,000 background points from one background extent raster
  # (only one per cell without replacement).
  # Note: If the raster has <10,000 pixels,
  # you'll get a warning and all pixels will be used for background.
  bg <- dismo::randomPoints(envs.backg[[bg.raster]], n=number.bg.points)
  bg <- raster::as.data.frame(bg)

  # Change column names.
  colnames(coords) <- c("lng", "lat")

  rm(pops)
  gc(verbose = FALSE)

  writeLines("\nSaving background partition plots...")
  writeLines("\nPerforming block method...")
  ### Block Method
  blocks <- ENMeval::get.block(coords, bg)
  pdf(file = file.path(plotDIR, "bgPartitions_blocks.pdf"),
      width = 7,
      height = 7,
      onefile = TRUE)
    raster::plot(envs.backg[[bg.raster]], col="gray",
         main = "Partitions into Training and Test Data",
         sub="Block Method",
         xlab = "Longitude",
         ylab = "Latitude")
    dismo::points(coords, pch=21, bg=blocks$occ.grp)
  dev.off()

  rm(blocks)
  gc(verbose = FALSE)

  writeLines("Performing checkerboard1 method...")
  ### Checkerboard1 Method
  check1 <- ENMeval::get.checkerboard1(coords,
                                       envs.cropped,
                                       bg,
                                       aggregation.factor = agg.factor)
    pdf(file = file.path(plotDIR, "bgPartitions_checkerboard1.pdf"),
        width = 7,
        height = 7,
        onefile = TRUE)
      raster::plot(
        envs.backg[[bg.raster]],
        col="gray",
        main = "Partitions into Training and Test Data (Aggregation Factor=5)",
        sub="Checkerboard1 Method",
        xlab = "Longitude",
        ylab = "Latitude"
        )
      dismo::points(bg, pch=21, bg=check1$occ.grp)
    dev.off()

    rm(check1)
    gc(verbose = FALSE)

    writeLines("Performing checkerboard2 method...")
    ### Using the Checkerboard2 method.
    check2 <-
      ENMeval::get.checkerboard2(coords,
                                 envs.cropped,
                                 bg,
                                 aggregation.factor = c(agg.factor,
                                                      agg.factor))
    pdf(file = file.path(plotDIR, "bgPartitions_checkerboard2.pdf"),
        width = 7,
        height = 7,
        onefile = TRUE)
    raster::plot(envs.backg[[bg.raster]], col="gray",
         main = "Partitions into Training and Test Data",
         sub="Checkerboard2 Method",
         xlab = "Longitude",
         ylab = "Latitude")
    dismo::points(bg, pch=21, bg=check2$bg.grp)
    dismo::points(coords, pch=21,
                  bg=check2$occ.grp,
                  col="white",
                  cex=1.5)
  dev.off()

  rm(check2)
  gc(verbose = FALSE)

  writeLines("\nDone!")

  return(bg)
}

#' Function to run ENMeval and MAXENT on the raster layers.
#'
#' This function takes as input the object output from preparte_rasters() and
#' the background partitioning method you would like to use
#' (see ENMeval documentation). You can visualize how the partitioning
#' methods will look by viewing the PDF output from partition_raster_bg.
#' See ?partition_raster_bg for more info.
#' @param envs.fg First element of envs.list returned from prepare_raster_bg()
#' @param bg Object with background points; generated from partition_raster_bg
#' @param coords Data.Frame with (lon,lat). 3rd element of envs.list
#' @param partition.method Method used for background point partitioning
#' @param parallel If TRUE, ENMeval is run parallelized with np CPU cores
#' @param np Number of parallel cores to use if parallel = TRUE
#' @param RMvalues Vector of non-negative regularization multiplier values.
#'                 Higher values impose a stronger penalty on model complexity
#' @param feature.classes Character vector of feature classes to be used
#' @param categoricals Vector indicating which (if any) of the input
#'                     environmental layers are categorical.
#' @param agg.factor Aggregation factor(s) for checkerboard
#'                           partitioning method.
#' @param algorithm Character vector. Defaults to dismo's maxent.jar
#' @return ENMeval object
#' @export
runENMeval <- function(envs.fg,
                       bg,
                       coords,
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
                       agg.factor = 2,
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

  if (partition.method == "checkerboard2"){
    agg.factor <- c(agg.factor, agg.factor)
  }

  if (parallel == TRUE){
    eval.par <-
      ENMeval::ENMevaluate(
        occ = coords,
        env = envs.fg,
        bg.coords = bg,
        method = partition.method,
        RMvalues = RMvalues,
        fc = feature.classes,
        parallel = parallel,
        algorithm = algorithm,
        categoricals = categoricals,
        aggregation.factor = agg.factor,
        numCores = np
      )
  } else if (parallel == FALSE){
    eval.par <-
      ENMeval::ENMevaluate(
        occ = coords,
        env = envs.fg,
        bg.coords = bg,
        method = partition.method,
        RMvalues = RMvalues,
        fc = feature.classes,
        parallel = parallel,
        algorithm = algorithm,
        categoricals = categoricals,
        aggregation.factor = agg.factor
      )
  } else {
    stop("Parallel must be either TRUE or FALSE")
  }

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
#' @param RMvalues Vector of non-negative RM values to examine how
#'                 complexity affects predictions
#' @param plotDIR Directory to save plots to
#' @param niche.overlap Boolean. If TRUE, calculates pairwise niche overlap
#'                      matrix
#' @param plot.width Integer. Specify plot widths
#' @param plot.height Integer. Specify plot heights
#' @param imp.margins Integer vector. Margins of variable importance barplot.
#'                    c(bottom, left, top, right)
#' @export
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
                              RMvalues = seq(0.5, 4, 0.5),
                              plotDIR = "./plots",
                              niche.overlap = FALSE,
                              plot.width = 7,
                              plot.height = 7,
                              imp.margins = c(10.0, 4.1, 4.1, 2.1)){

  dir.create(plotDIR, showWarnings = FALSE)


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
  eval.par@results[which(eval.par@results$delta.AICc==0),]
  eval.par@results[which(eval.par@results$delta.AICc==0),]
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

  # Look at the model object for our "AICc optimal" model:
  aic.opt <- eval.par@models[[which(eval.par@results$delta.AICc==0)]]

  # Write it to file.
  write.csv(x = aic.opt@results,
            file = file.path(plotDIR,
                             "ENMeval_maxentModel_aicOptimal.csv"))

  # Get a data.frame of two variable importance metrics:
  # percent contribution and permutation importance.
  varImportance <- ENMeval::var.importance(aic.opt)

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
    ENMeval::eval.plot(eval.par@results)
    ENMeval::eval.plot(eval.par@results, 'avg.test.AUC',
                       variance ='var.test.AUC')
    ENMeval::eval.plot(eval.par@results, "avg.diff.AUC",
                       variance="var.diff.AUC")

    # Plot permutation importance.
    df <- ENMeval::var.importance(aic.opt)
    par(mar=imp.margins)
    raster::barplot(df$permutation.importance,
                    names.arg = df$variable,
                    las = 2,
                    ylab = "Permutation Importance")
  dev.off()

  # Let's see how model complexity changes the predictions in our example
  pdf(file = file.path(plotDIR, "modelPredictions.pdf"),
      width = plot.width,
      height = plot.height,
      onefile = TRUE)
    for (i in 1:length(examine.predictions)){
      for (j in 1:length(RMvalues)){
        raster::plot(eval.par@predictions[[paste(examine.predictions[i],
                                                    RMvalues[j],
                                                    sep = "_")]],
                       ylim=c(minLat,maxLat),
                       xlim=c(minLon, maxLon),
                       main=paste0(examine.predictions[i],
                                   "_",
                                   RMvalues[j],
                                   " Prediction"),
        )
      }
    }

  dev.off()
}

#' Function to extract raster values at each sample point for raster stack
#' @param env.list List object generated from prepare_rasters()
#' @return List of data.frames with raster values at each sample locality
#'         for each raster layer, and list of raster names
#' @export
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
