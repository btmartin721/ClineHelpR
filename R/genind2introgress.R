#' Generate inputs for Introgress from genind object
#'
#' This function parses a genind object to create input data frames for
#' running genomic cline analysis in the Introgress R package
#'
#' @param gen Genind object containing data for parental and admixed populations
#' @param p1 Character or vector containing population names for the p1 group
#' @param p2 Character or vector containing population names for the p2 group
#' @param admix Character or vector containing population names for the admixed individuals
#' @param missingPerInd Proportion of missing genotypes allowed to retain an individual
#' @param missingPerLoc Proportion of missing genotypes allowed to retain a locus
#' @param subset An optional parameter defining a number of loci to randomly retain
#' @param drop Boolean value defining whether or not to remove monomorphic loci
#' @param prefix Prefix for output Introgress files (default="introgress")
#' @export
#' @examples
#' introgress_input <- genind2introgress(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"))
#'
#' introgress_input <- genind2introgress(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"),
#'                                         missingPerInd=0.75, missingPerLoc=0.25,
#'                                         drop=TRUE)
#'
#' introgress_input <- genind2introgress(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"),
#'                                         subset=500)
genind2introgress <- function(
	gen,
	p1,
	p2,
	admix,
	missingPerLoc=0.5,
	missingPerInd=0.5,
	subset=NULL,
	drop=TRUE,
	prefix="introgress"
){

	ret<-list()

	all_pops<-c(p1, p2, admix)

	missingPerLoc = 1.0 - missingPerLoc
	missingPerInd = 1.0 - missingPerInd
	# drop individuals not in selected pops
	sub<-gen[adegenet::pop(gen) %in% all_pops,]

	# drop loci having more than missingPerLoc missing data
	sub<-filter_missingByLoc(sub, prop=missingPerLoc)

	# drop individuals having more than missingPerInd missing data
	# NOTE: Also dropping loci that become monomorphic
	sub<-filter_missingByInd(sub, prop=missingPerInd, drop=drop)
  #print(sub)
	# if necessary, randomly subset loci to the number requested
	if (!is.null(subset)){
		sub<-filter_randomSubsetLoci(sub, sample=subset)
	}
  #print(sub)
	# begin parsing remaining data to create inputs for introgress
	# set up loci information table ()
	loci.data<-data.frame(cbind(c(names(sub$loc.n.all)), sub$type), stringsAsFactors = FALSE)
	colnames(loci.data)<-c("locus", "type")
	loci.data[loci.data["type"]=="codom", "type"]<-"C"
	loci.data[loci.data["type"]=="dom", "type"]<-"D"

	# genotype table for p1
	p1_genind<-sub[adegenet::pop(sub) %in% p1,]
	p1.data<-adegenet::genind2df(p1_genind, sep="/")
	p1.data<-subset(p1.data, select=-c(pop))
	p1.data[is.na(p1.data)] <- "NA/NA"

	# genotype table for p2
	p2_genind<-sub[adegenet::pop(sub) %in% p2,]
	p2.data<-adegenet::genind2df(p2_genind, sep="/")
	p2.data<-subset(p2.data, select=-c(pop))
	p2.data[is.na(p2.data)] <- "NA/NA"

	# genotype table for admix population
	admix_genind<-sub[adegenet::pop(sub) %in% admix,]
	admix.data<-adegenet::genind2df(admix_genind, sep="/")
	pop<-admix.data$pop
	admix.data<-subset(admix.data, select=-c(pop))
	admix.data[is.na(admix.data)] <- "NA/NA"
	admix.data <- rbind(pop, row.names(admix.data), admix.data)

	# write output files
	print(paste0("Writing outputs with prefix: ", prefix))
	if (file.exists(paste0(prefix, "_p1data.txt"))) {
	  # Delete file if it exists
	  file.remove(paste0(prefix, "_p1data.txt"))
	}
	write.table(data.table::transpose(p1.data),
	            file=paste0(prefix, "_p1data.csv"),
	            sep="\t", col.names=FALSE, row.names=FALSE,
	            quote = FALSE)

	if (file.exists(paste0(prefix, "_p2data.txt"))) {
	  # Delete file if it exists
	  file.remove(paste0(prefix, "_p2data.txt"))
	}
	write.table(data.table::transpose(p2.data),
	            file=paste0(prefix, "_p2data.csv"),
	            sep="\t", col.names=FALSE, row.names=FALSE,
	            quote = FALSE)

	if (file.exists(paste0(prefix, "_loci.txt"))) {
	  # Delete file if it exists
	  file.remove(paste0(prefix, "_loci.txt"))
	}
	#l<-as.data.frame(t(loci.data))
  write.table(loci.data, file=paste0(prefix, "_loci.txt"),
            sep="\t", col.names=TRUE, row.names=FALSE,
            quote = FALSE)

  if (file.exists(paste0(prefix, "_admix.txt"))) {
    # Delete file if it exists
    file.remove(paste0(prefix, "_admix.txt"))
  }
  write.table(data.table::transpose(admix.data),
              file=paste0(prefix, "_admix.txt"),
              sep="\t", col.names=FALSE, row.names=FALSE,
              quote = FALSE)
}

###############################################
#' Utility function to filter missing data by locus.
#' @param gen DataFrame Genotypes to filter.
#' @param prop Float Missing data threshold for which to filter loci.
#' @return DataFrame A filtered DataFrame where loci with missing data > prop are removed.
#' @noRd
filter_missingByLoc<-function(gen, prop=0.5){
  propTyped<-adegenet::propTyped(gen, by="loc")
	#print(missing)
	#print(missing>prop)
	#print(genind[loc=c(missing>prop)])
	return(gen[loc=c(propTyped>prop)])
}

#' Utility function to filter missing data by individual.
#' @param gen DataFrame Genotypes to filter.
#' @param prop Float Proportion of missing data threshold.
#' @param drop Boolean Whether to drop monomorphic loci.
#' @return DataFrame Individuals with missing data > prop filtered out.
#' @noRd
filter_missingByInd<-function(gen, prop=0.5, drop=TRUE){
	propTyped<-adegenet::propTyped(gen, by="ind")
	#print(propTyped)
	#print(prop)
	#print(propTyped)
	#print(propTyped>prop)
	return(gen[c(propTyped>prop), drop=drop])
}

#' Utility function to randomly subset loci.
#' @param gen DataFrame Genotypes to randomly subset.
#' @param sample Integer Number of loci to subset.
#' @return DataFrame Random subset of genotypes.
#' @noRd
filter_randomSubsetLoci<-function(gen, sample){
	loc<-adegenet::nLoc(gen)
	if (sample < loc){
		locs <- sort(sample(1:loc, sample, replace=FALSE))
		return(gen[loc=c(locs)])
	}
}

