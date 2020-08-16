
# Example data can be obtained from a Dryad Digital Repository (doi: XXXX)
# To run this script, please download the example data.
# Change the paths to the input files as needed below.

library("ClinePlotR")

# Set working directory here.
# setwd()

###############################################################################
# Set input and output paths.
# Also set various parameters you'll need to run the script.
# Adjust these as necessary if you use diferent paths or filenames.
###############################################################################

# output directory to save plots in; will create if doesn't exist
plotDIR <- "./plots"

# Directory with BGC output files; transcriptomic alignment
genesDIR <- "exampleData/bgc/bgc_outputFiles/genes"

# Directory with BGC output files; scaffold alignment
fullDIR <- "exampleData/bgc/bgc_outputFiles/fulldataset"

# PAFScaff output file; needed for ideogram
pafInfo  <-
  "exampleData/PAFScaff/pafscaff_asm20_scafTrans_tscripta.scaffolds.tdt"

# prefix for BGC output files
prefix <- "eatt"

# admixed population from popmap file
admixPop <- "EATT"

# Path to population map (popmap) file
popmap <- "exampleData/popmaps/bgc/eatt.bgc.popmap_final.txt"

# File with locus names (transcriptomic alignment)
genes.loci.file <- "exampleData/bgc/bgc_lociFiles/genes/eatt_bgc_loci.txt"

# File with locus names (scaffold alignment)
full.loci.file <- "exampleData/bgc/bgc_lociFiles/fulldataset/eatt_bgc_loci.txt"

# Directory to save scaffold info files; Will create if doesnt't exist
scafInfoDIR <- "./scafInfo"

# Set size of heatmap bands on ideogram. Measured in nucleotide bases
# Strictly for better visualization. Can be adjusted as needed
gene.size <- 4e6 # known genes
other.size <- 1e6 # scaffold alignment

# GFF file is from the Terrapene mexicana triunguis scaffold-level genome
# It can be found at: https://www.ncbi.nlm.nih.gov/genome/?term=Terrapene
gff <- "exampleData/gff/genes_Terrapene.gff"

###############################################################################
# Run BGC plotting functions
###############################################################################

##### Transcripomic alignment ######

# Combine multiple independent BGC runs together
bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = prefix)

# Plot the likelihood and parameter traces
plot_traces(df.list = bgc.genes,
            prefix = prefix,
            plotDIR = plotDIR)

# Detect BGC outliers for known genes
gene.outliers <-
  get_bgc_outliers(
    df.list = bgc.genes,
    admix.pop = admixPop,
    popmap = popmap,
    loci.file = genes.loci.file
  )


##### Scaffold alignment #####

# Aggregate runs together
bgc.full <-
  combine_bgc_output(results.dir = fullDIR,
                     prefix = prefix)

# Plot likelihood and parameter traces
plot_traces(
  df.list = bgc.full,
  prefix = paste0(prefix, "_full"),
  plotDIR = plotDIR
)

# Detect outliers in scaffold alignment
full.outliers <-
  get_bgc_outliers(
    df.list = bgc.full,
    admix.pop = admixPop,
    popmap = popmap,
    loci.file = full.loci.file
  )

# Not ClinePlotR functions.
# This just saves the minimum and maximum alpha and beta outlier values
# Can be used to add to the legend
ab.range <-
  data.frame(
    "full.alpha" = c(
      "min" = min(full.outliers[[1]]$alpha),
      "max" = max(full.outliers[[1]]$alpha)
    ),
    "full.beta" = c(
      "min" = min(full.outliers[[1]]$beta),
      "max" = max(full.outliers[[1]]$beta)
    ),
    "genes.alpha" = c(
      "min" = min(gene.outliers[[1]]$alpha),
      "max" = max(gene.outliers[[1]]$alpha)
    ),
    "genes.beta" = c(min(gene.outliers[[1]]$beta), max(gene.outliers[[1]]$beta))
  )

# Write alpha beta ranges to file.
write.table(
  data.frame("Header" = rownames(ab.range), ab.range),
  file = file.path(plotDIR, paste0(prefix, "_ab.ranges.csv")),
  sep = ",",
  row.names = F,
  col.names = T,
  quote = F
)

# Make the phi plot for the transcriptomic alignment
# Any of these parameters can be adjusted as needed.
# Here, both.outlier.tests is FALSE
# This means that outliers are flagged if they are significant in either method
phiPlot(
  outlier.list = gene.outliers,
  popname = paste0(admixPop, " Genes"),
  line.size = 0.35,
  saveToFile = paste0(prefix, "_genes"),
  plotDIR = plotDIR,
  both.outlier.tests = FALSE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)

# Here both.outlier.tests is TRUE
# This means that both outlier tests have to be significant to flag the outlier
phiPlot(
  outlier.list = gene.outliers,
  popname = paste0(admixPop, " Genes"),
  line.size = 0.35,
  saveToFile = paste0(prefix, "_genes_bothOutlierTests"),
  plotDIR = plotDIR,
  both.outlier.tests = TRUE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)

# Phi plot for scaffold alignment
phiPlot(
  outlier.list = full.outliers,
  popname = paste0(admixPop, " All"),
  line.size = 0.2,
  saveToFile = paste0(prefix, "_scaffold"),
  plotDIR = plotDIR,
  both.outlier.tests = FALSE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)

# alphabetaplot
# 2-D contour plot with hulls for outliers
alphaBetaPlot(
  gene.outliers,
  alpha.color = "cornflowerblue",
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = prefix,
  plotDIR = "./plots",
  padding = 0.2,
)

# Read and parse GFF file
gff <- parseGFF(gff.filepath = gff)

# Get annotation info for BGC outliers (transcriptomic alignment)
genes.annotated <-
  join_bgc_gff(
    prefix = prefix,
    outlier.list = gene.outliers,
    gff.data = gff,
    scafInfoDIR = scafInfoDIR
  )

# Plot ideogram
# Here, both.outlier.tests is FALSE
# This means that the SNP is an outlier if flagged by either outlier test
ref.trans.info <-
  plot_outlier_ideogram(
    prefix = paste0(prefix, "_eitherOutlierTests"),
    outliers.genes = genes.annotated,
    outliers.full.scaffolds = full.outliers,
    pafInfo = pafInfo,
    plotDIR = plotDIR,
    gene.size = gene.size,
    other.size = other.size,
    both.outlier.tests = FALSE
  )

# Write annotation information to file so it can be used later
write.table(
  genes.annotated,
  file = file.path(plotDIR, paste0(prefix, "_genesAnnotated.csv")),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = ","
)

# Write annotation information to file so it can be used later
# Has some extra columns compared to the genesAnnotated.csv file
write.table(
  ref.trans.info,
  file = file.path(plotDIR, paste0(prefix, "_refTransInfo.csv")),
  quote = TRUE,
  row.names = FALSE,
  col.names = TRUE,
  sep = ","
)



