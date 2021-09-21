#!/usr/local/bin/R
#'Generate inputs for BGC without genotype uncertainty
#'
#' This function parses a genind object to create input data frames for 
#' running genomic cline analysis in BGC
#'
#' @param gen Genind object containing data for parental and admixed populations
#' @param p1 Character or vector containing population names for the p1 group
#' @param p2 Character or vector containing population names for the p2 group
#' @param admix Character or vector containing population names for the admixed individuals
#' @param missingPerInd Proportion of missing genotypes allowed to retain an individual 
#' @param missingPerLoc Proportion of missing genotypes allowed to retain a locus
#' @param subset An optional parameter defining a number of loci to randomly retain
#' @param prefix A string giving the prefix for bgc input file names
#' @examples
#' bgc_input <- genind2bgc(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"))
#'
#' bgc_input <- genind2bgc(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"),
#'                                         missingPerInd=0.75, missingPerLoc=0.25)
#'
#' bgc_input <- genind2bgc(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"),
#'                                         subset=500)
genind2bgc <- function(
	gen,
	p1,
	p2,
	admix,
	missingPerLoc=0.5,
	missingPerInd=0.5,
	subset=NULL,
	drop=TRUE,
	prefix='bgc'
){
	
	ret<-list()
	
	all_pops<-c(p1, p2, admix)
	
	#drop individuals not in selected pops
	sub<-gen[pop(gen) %in% all_pops,]
	
	#drop loci having more than missingPerLoc missing data
	sub<-filter_missingByLoc(sub, prop=missingPerLoc)
	
	#drop individuals having more than missingPerInd missing data
	#NOTE: Also dropping loci that become monomorphic
	sub<-filter_missingByInd(sub, prop=missingPerInd, drop=drop)
	
	#if necessary, randomly subset loci to the number requested
	if (!is.null(subset)){
		sub<-filter_randomSubsetLoci(sub, sample=subset)
	}
	
	#get locus names 
	loci.data<-data.frame(cbind(c(names(sub$loc.n.all))), stringsAsFactors = FALSE)
	colnames(loci.data)<-c("locus")
	
	#genotype table for p1
	p1_genind<-sub[pop(sub) %in% p1,]
	p1.data<-adegenet::genind2df(p1_genind, sep="/")
	p1.data<-subset(p1.data, select=-c(pop))
	#print(p1.data)
	
	#genotype table for p2
	p2_genind<-sub[pop(sub) %in% p2,]
	p2.data<-adegenet::genind2df(p2_genind, sep="/")
	p2.data<-subset(p2.data, select=-c(pop))

	# #genotype table for admix population
	admix_genind<-sub[pop(sub) %in% admix,]
	admix.data<-adegenet::genind2df(admix_genind, sep="/")
	pops<-unique(admix.data$pop)

	#loop through and gather allele counts 
	p1.out <- list()
	p2.out <- list()
	admix.out <- list()
	
	for (l in 1:length(loci.data$locus)){
	  #print(loci.data$locus[l])
	  alleles <- getRefAlt(c(p1.data[,l], p2.data[,l], admix.data[,l+1]))
	  if (length(alleles) != 2){
	    print(paste0("Locus ", loci.data$locus[l], " not biallelic. Skipping."))
	    next
	  }
	  #append locus names
	  locnum <- paste0("locus_", gsub("L", "", loci.data$locus[l]))
	  p1.out <- c(p1.out, locnum)
	  p2.out <- c(p2.out, locnum)
	  admix.out <- c(admix.out, locnum)
	  
	  #get allele counts for parental pops
	  p1.g <- unlist(lapply(p1.data[,l], function(x) strsplit(x, "/")))
	  p1.out <- c(p1.out, paste(getCounts(p1.g, alleles), collapse=" "))
	  p2.g <- unlist(lapply(p2.data[,l], function(x) strsplit(x, "/")))
	  p2.out <- c(p2.out, paste(getCounts(p2.g, alleles), collapse=" "))
	  
	  #get genotypes for admix pops
	  for (p in pops){
	    admix.out <- c(admix.out, paste0("pop_", p))
	    pop_dat <- admix.data[admix.data$pop==p,l+1]
	    for (i in pop_dat){
	      g <- unlist(lapply(i, function(x) strsplit(x, "/")))
	      admix.out <- c(admix.out, paste(getCounts(g, alleles), collapse=" "))
	    }
	  }

	}
  #write output files
	print(paste0("Writing outputs with prefix: ", prefix))
	if (file.exists(paste0(prefix, "_p0in.txt"))) {
	  #Delete file if it exists
	  file.remove(paste0(prefix, "_p0in.txt"))
	}
	t<- lapply(p1.out, write, paste0(prefix, "_p0in.txt"), append=TRUE)
	if (file.exists(paste0(prefix, "_p1in.txt"))) {
	  #Delete file if it exists
	  file.remove(paste0(prefix, "_p1in.txt"))
	}
	t<- lapply(p2.out, write, paste0(prefix, "_p1in.txt"), append=TRUE)
	if (file.exists(paste0(prefix, "_admixedin.txt"))) {
	  #Delete file if it exists
	  file.remove(paste0(prefix, "_admixedin.txt"))
	}
	t<-lapply(admix.out, write, paste0(prefix, "_admixedin.txt"), append=TRUE)
  #return(NULL)
}

###############################################
#Utility functions
getCounts<-function(genotypes, alleles){
  counts<-summary(as.factor(genotypes))
  ret <-c(as.numeric(counts[alleles[1]]), as.numeric(counts[alleles[2]]))
  ret[is.na(ret)] <- 0
  return(ret)
}

getRefAlt<-function(genotypes){
  l <- unlist(lapply(genotypes[!is.na(genotypes)], function(x) strsplit(x, "/")))
  n <- names(sort(summary(as.factor(l)), decreasing=T))
  return(n)
}

filter_missingByLoc<-function(gen, prop=0.5){
	missing<-adegenet::propTyped(gen, by="loc")
	#print(genind[loc=c(missing>prop)])
	return(gen[loc=c(missing>prop)])
}

filter_missingByInd<-function(gen, prop=0.5, drop=TRUE){
	missing<-adegenet::propTyped(gen, by="ind")
	#print(missing>prop)
	return(gen[c(missing>prop), drop=drop])
}

filter_randomSubsetLoci<-function(gen, sample){
	loc<-adegenet::nLoc(gen)
	if (sample < loc){
		locs <- sort(sample(1:loc, sample, replace=FALSE))
		return(gen[loc=c(locs)])
	}
}

